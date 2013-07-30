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

#ifndef _AQUARIUS_CC_1EDENSITY_HPP_
#define _AQUARIUS_CC_1EDENSITY_HPP_

#include "operator/exponentialoperator.hpp"
#include "operator/deexcitationoperator.hpp"
#include "operator/1eoperator.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class OneElectronDensity : public op::OneElectronOperator<U>
{
    public:
        /*
         * Form the SCF density
         */
        OneElectronDensity(const scf::UHF<U>& uhf);

        /*
         * Form the unrelaxed CCSD density
         */
        OneElectronDensity(const op::ExponentialOperator<U,2>& T);

        /*
         * Form the partial perturbed CCSD Density
         */
        OneElectronDensity(const op::DeexcitationOperator<U,2>& L,
                           const op::ExponentialOperator<U,2>& T,
                           const op::ExcitationOperator<U,2>& TA);

        /*
         * Form the relaxed CCSD density
         */
        OneElectronDensity(const op::DeexcitationOperator<U,2>& L,
                           const op::ExponentialOperator<U,2>& T);

        /*
         * Form the relaxed perturbed CCSD Density
         */
        OneElectronDensity(const op::DeexcitationOperator<U,2>& L,
                           const op::DeexcitationOperator<U,2>& LA,
                           const op::ExponentialOperator<U,2>& T,
                           const op::ExcitationOperator<U,2>& TA);
};

}
}

#endif
