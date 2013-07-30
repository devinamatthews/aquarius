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

#include "tensor/composite_tensor.hpp"
#include "tensor/spinorbital_tensor.hpp"

#include "mooperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T, int np, int nh=np>
class DeexcitationOperator
: public MOOperator<T>,
  public tensor::CompositeTensor< DeexcitationOperator<T,np,nh>,
                                  tensor::SpinorbitalTensor<T>, T >
{
    INHERIT_FROM_COMPOSITE_TENSOR(CONCAT(DeexcitationOperator<T,np,nh>),
                                  tensor::SpinorbitalTensor<T>, T)

    protected:
        const int spin;

    public:
        DeexcitationOperator(const DeexcitationOperator<T,np,nh>& other)
        : MOOperator<T>(other.uhf),
          tensor::CompositeTensor< DeexcitationOperator<T,np,nh>,
           tensor::SpinorbitalTensor<T>, T >(other),
          spin(other.spin) {}

        DeexcitationOperator(const scf::UHF<T>& uhf, const int spin=0)
        : MOOperator<T>(uhf),
          tensor::CompositeTensor< DeexcitationOperator<T,np,nh>,
           tensor::SpinorbitalTensor<T>, T >(std::max(np,nh)+1),
          spin(spin)
        {
            if (abs(spin%2) != (np+nh)%2 || abs(spin) > np+nh) throw std::logic_error("incompatible spin");

            int N = uhf.getMolecule().getNumOrbitals();
            int nI = uhf.getMolecule().getNumAlphaElectrons();
            int ni = uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;

            for (int ex = 0;ex <= std::min(np,nh);ex++)
            {
                int idx = ex+std::abs(np-nh);

                int npex = (np > nh ? ex+np-nh : ex);
                int nhex = (nh > np ? ex+nh-np : ex);

                tensors_[idx].tensor_ = new tensor::SpinorbitalTensor<T>(0, autocc::Manifold(npex,nhex), -spin);

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

                    if (ha > 0) sym[ha-1] = NS;
                    if (ha+hb > 0) sym[ha+hb-1] = NS;
                    if (ha+hb+pa > 0) sym[ha+hb+pa-1] = NS;
                    if (ha+hb+pa+pb > 0) sym[ha+hb+pa+pb-1] = NS;

                    tensors_[idx].tensor_->addSpinCase(new tensor::DistTensor<T>(uhf.ctf, ndim, len.data(), sym.data(), true),
                                                       0, autocc::Manifold(pa, ha));
                }
            }
        }

        T dot(bool conja, const op::DeexcitationOperator<T,np,nh>& A, bool conjb) const
        {
            T s = (T)0;

            for (int i = 0;i <= std::min(np,nh);i++)
            {
                s += (*this)(i).dot(conja, A(i), conjb)/(T)factorial(i)/(T)factorial(i+std::abs(np-nh));
            }

            return s;
        }
};

}
}

#endif
