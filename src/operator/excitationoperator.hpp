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

#ifndef _AQUARIUS_OPERATOR_EXCITATIONOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_EXCITATIONOPERATOR_HPP_

#include "util/stl_ext.hpp"
#include "tensor/composite_tensor.hpp"
#include "tensor/spinorbital_tensor.hpp"

#include "mooperator.hpp"
#include "denominator.hpp"

namespace aquarius
{
namespace op
{

template <typename T, int np, int nh=np>
class ExcitationOperator
: public MOOperator,
  public tensor::CompositeTensor< ExcitationOperator<T,np,nh>,
                                  tensor::SpinorbitalTensor<T>, T >
{
    INHERIT_FROM_COMPOSITE_TENSOR(CONCAT(ExcitationOperator<T,np,nh>),
                                  tensor::SpinorbitalTensor<T>, T)

    protected:
        const int spin;

    public:
        ExcitationOperator(const Arena& arena, const Space& occ, const Space& vrt, int spin=0)
        : MOOperator(arena, occ, vrt),
          tensor::CompositeTensor< ExcitationOperator<T,np,nh>,
           tensor::SpinorbitalTensor<T>, T >(std::max(np,nh)+1),
          spin(spin)
        {
            for (int ex = 0;ex <= std::min(np,nh);ex++)
            {
                int nv = ex+(np > nh ? np-nh : 0);
                int no = ex+(nh > np ? nh-np : 0);

                tensors[ex+std::abs(np-nh)].tensor =
                    new tensor::SpinorbitalTensor<T>(arena, std::vec(vrt,occ), std::vec(nv,0), std::vec(0,no), spin);
            }
        }

        void weight(const Denominator<T>& d)
        {
            std::vector<const std::vector<T>*> da = vec(&d.getDA(), &d.getDI());
            std::vector<const std::vector<T>*> db = vec(&d.getDa(), &d.getDi());

            for (int ex = 0;ex <= std::min(np,nh);ex++)
            {
                if (ex== 0 && np == nh) continue;
                tensors[ex+std::abs(np-nh)].tensor->weight(da, db);
            }
        }

        T dot(bool conja, const op::ExcitationOperator<T,np,nh>& A, bool conjb) const
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
