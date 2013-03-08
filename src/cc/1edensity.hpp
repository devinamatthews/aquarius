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

#include "util/distributed.hpp"
#include "tensor/spinorbital.hpp"
#include "tensor/dist_tensor.hpp"
#include "scf/scf.hpp"
#include "operator/excitationoperator.hpp"
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
        OneElectronDensity(const scf::UHF<U>& uhf)
        : op::OneElectronOperator<U>(uhf)
        {
            //TODO: SCF density
            assert(0);
        }

        OneElectronDensity(const DeexcitationOperator<U,2>& L, const ExponentialOperator<U,2>& T)
        : op::OneElectronOperator<U>(L.getSCF())
        {
            this->ia["ijab"] = -L[1]["ia"];

            this->ab["ab"] = -0.5*T[2]["aemn"]*L[2]["mnbe"];

            this->ij["ij"] = T[1]["ei"]*L[1]["je"];
            this->ij["ij"] += 0.5*T[2]["efim"]*L[2]["jmef"];

            this->ai["ai"] = T[1]["ai"];
            this->ai["ai"] -= T[2]["aeim"]*L[1]["me"];
            this->ai["ai"] += this->ij["mi"]*T[1]["am"];
            this->ai["ai"] -= this->ab["ae"]*T[1]["ei"];

            this->ab["ab"] -= T[1]["am"]*L[1]["mb"];
        }
};

}
}

#endif
