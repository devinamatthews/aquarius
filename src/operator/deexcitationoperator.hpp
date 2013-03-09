/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_OPERATOR_DEEXCITATIONOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_DEEXCITATIONOPERATOR_HPP_

#include <exception>

#include "tensor/dist_tensor.hpp"
#include "tensor/spinorbital.hpp"
#include "tensor/util.h"
#include "util/util.h"
#include "util/distributed.hpp"
#include "scf/scf.hpp"

#include "mooperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T, int np, int nh> class ScaledDeexcitationOperator;

template <typename T, int np, int nh=np>
class DeexcitationOperator : public MOOperator<T>
{
    protected:
        const int spin;
        std::vector< tensor::SpinorbitalTensor< tensor::DistTensor<T> >* > r;

    public:
        DeexcitationOperator(const scf::UHF<T>& uhf, const int spin=0)
        : MOOperator<T>(uhf), spin(spin),
          r(std::min(np,nh)+1, (tensor::SpinorbitalTensor< tensor::DistTensor<T> >*)NULL)
        {
            if (abs(spin%2) != (np+nh)%2 || abs(spin) > np+nh) throw std::logic_error("incompatible spin");

            int N = uhf.getMolecule().getNumOrbitals();
            int nI = uhf.getMolecule().getNumAlphaElectrons();
            int ni = uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;

            for (int ex = 0;ex <= std::min(np,nh);ex++)
            {
                int npex = (np > nh ? ex+np-nh : ex);
                int nhex = (nh > np ? ex+nh-np : ex);

                r[ex] = new tensor::SpinorbitalTensor< tensor::DistTensor<T> >(0, autocc::Manifold(npex,nhex), -spin);

                for (int pspin = std::max(-npex,spin-nhex);pspin <= std::min(npex,spin+nhex);pspin++)
                {
                    if (abs(pspin%2) != npex%2) continue;

                    int hspin = pspin-spin;
                    int pa = (npex+pspin)/2;
                    int pb = (npex-pspin)/2;
                    int ha = (nhex+hspin)/2;
                    int hb = (nhex-hspin)/2;

                    int ndim = npex+nhex;
                    std::vector<int> len(ndim);
                    std::vector<int> sym(ndim, AS);

                    int i = 0;
                    for (int j = 0;j < ha;j++,i++) len[i] = nI;
                    for (int j = 0;j < hb;j++,i++) len[i] = ni;
                    for (int j = 0;j < pa;j++,i++) len[i] = nA;
                    for (int j = 0;j < pb;j++,i++) len[i] = na;

                    if (pa > 0) sym[pa-1] = NS;
                    if (pa+pb > 0) sym[pa+pb-1] = NS;
                    if (pa+pb+ha > 0) sym[pa+pb+ha-1] = NS;
                    if (pa+pb+ha+hb > 0) sym[pa+pb+ha+hb-1] = NS;

                    r[ex]->addSpinCase(new tensor::DistTensor<T>(uhf.ctf, ndim, len.data(), sym.data(), true),
                                       0, autocc::Manifold(pa, ha));
                }
            }
        }

        T operator*(const DeexcitationOperator<T,np,nh>& other) const
        {
            T res = (T)0;

            for (int ex = 0;ex <= std::min(np,nh);ex++)
            {
                res += tensor::scalar((*r[ex])*(*other.r[ex]))/factorial(ex)/factorial(ex+abs(np-nh));
            }

            return res;
        }

        DeexcitationOperator<T,np,nh>& operator=(const T other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] = other;
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator+=(const T other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] += other;
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator-=(const T other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] -= other;
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator*=(const T other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] *= other;
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator/=(const T other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] /= other;
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator=(const DeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] = *other.r[ex];
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator+=(const DeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] += *other.r[ex];
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator-=(const DeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] -= *other.r[ex];
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator*=(const DeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] *= *other.r[ex];
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator/=(const DeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] /= *other.r[ex];
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator=(const ScaledDeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] = *other.op.r[ex]*other.factor;
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator+=(const ScaledDeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] += *other.op.r[ex]*other.factor;
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator-=(const ScaledDeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] -= *other.op.r[ex]*other.factor;
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator*=(const ScaledDeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] *= *other.op.r[ex]*other.factor;
            return *this;
        }

        DeexcitationOperator<T,np,nh>& operator/=(const ScaledDeexcitationOperator<T,np,nh>& other)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++) *r[ex] /= *other.op.r[ex]*other.factor;
            return *this;
        }

        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& operator[](int i)
        {
            assert(i >= abs(np-nh) && i <= std::max(np,nh));
            return *r[i-abs(np-nh)];
        }

        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& operator[](int i) const
        {
            assert(i >= abs(np-nh) && i <= std::max(np,nh));
            return *r[i-abs(np-nh)];
        }

        const ScaledDeexcitationOperator<T,np,nh> operator*(const T factor) const
        {
            return ScaledDeexcitationOperator<T,np,nh>(*this, factor);
        }

        friend ScaledDeexcitationOperator<T,np,nh> operator*(const T factor, const DeexcitationOperator<T,np,nh>& other)
        {
            return ScaledDeexcitationOperator<T,np,nh>(other, factor);
        }
};

template <typename T, int np, int nh>
class ScaledDeexcitationOperator
{
    friend class DeexcitationOperator<T,np,nh>;

    protected:
        const DeexcitationOperator<T,np,nh>& op;
        T factor;

        ScaledDeexcitationOperator(const DeexcitationOperator<T,np,nh>& op, T factor)
        : op(op), factor(factor) {}
};

}
}

#endif
