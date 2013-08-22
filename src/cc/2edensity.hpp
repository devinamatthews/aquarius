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

#ifndef _AQUARIUS_CC_DENSITY_HPP_
#define _AQUARIUS_CC_DENSITY_HPP_

#include "operator/2eoperator.hpp"

#include "1edensity.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class TwoElectronDensity : public op::TwoElectronOperator<U>
{
    public:
        /*
         * Form the SCF density
         *
         * P1_ab = D_ab (2D?)
         *
         * P2_abcd = P1_ab P1_cd + 1/2 P1_ac P1_bd
         */
        TwoElectronDensity(const op::MOSpace<U>& occ, const op::MOSpace<U>& vrt,
                           const tensor::DistTensor<U>& Da, const tensor::DistTensor<U>& Db);

        /*
         * Form the unrelaxed CCSD Density
         */
        TwoElectronDensity(const op::ExcitationOperator<U,2>& T);

        /*
         * Form the mixed perturbed CCSD Density
         */
        TwoElectronDensity(const op::DeexcitationOperator<U,2>& L,
                           const op::ExcitationOperator<U,2>& T,
                           const op::ExcitationOperator<U,2>& TA);

        /*
         * Form the relaxed CCSD Density
         */
        TwoElectronDensity(const op::DeexcitationOperator<U,2>& L,
                           const op::ExcitationOperator<U,2>& T);

        /*
         * Form the relaxed perturbed CCSD Density
         */
        TwoElectronDensity(const op::DeexcitationOperator<U,2>& L,
                           const op::DeexcitationOperator<U,2>& LA,
                           const op::ExcitationOperator<U,2>& T,
                           const op::ExcitationOperator<U,2>& TA);

        double getS2() const;

        double getProjectedMultiplicity() const;
};

}
}

#endif
