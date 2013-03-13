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

#ifndef _AQUARIUS_CC_HBAR_HPP_
#define _AQUARIUS_CC_HBAR_HPP_

#include "operator/2eoperator.hpp"
#include "operator/exponentialoperator.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class Hbar : public op::TwoElectronOperator<U>
{
    public:
        Hbar(const op::TwoElectronOperator<U>& moints, const op::ExponentialOperator<U,2>& T)
        : op::TwoElectronOperator<U>(const_cast<op::TwoElectronOperator<U>&>(moints), op::TwoElectronOperator<U>::ALL)
        {
            tensor::SpinorbitalTensor< tensor::DistTensor<U> > Tau(T(2));
            Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

            this->ai = 0;
            this->abij = 0;

            this->ia["me"] = this->ijab["mnef"]*T(1)["fn"];

            this->ij["mi"] += 0.5*this->ijab["nmef"]*T(2)["efni"];
            this->ij["mi"] += this->ia["me"]*T(1)["ei"];
            this->ij["mi"] += this->ijka["mnif"]*T(1)["fn"];

            this->ab["ae"] -= 0.5*this->ijab["mnfe"]*T(2)["famn"];
            this->ab["ae"] -= this->ia["me"]*T(1)["am"];
            this->ab["ae"] += this->aibc["anef"]*T(1)["fn"];

            this->ijkl["mnij"] += 0.5*this->ijab["mnef"]*Tau["efij"];
            this->ijkl["mnij"] += this->ijka["mnie"]*T(1)["ej"];

            this->ijka["mnie"] += this->ijab["mnfe"]*T(1)["fi"];

            this->iajk["mbij"] += 0.5*this->aibc["bmfe"]*Tau["efij"];
            this->iajk["mbij"] -= this->aibj["bmej"]*T(1)["ei"];
            this->iajk["mbij"] += this->ijka["mnie"]*T(2)["bejn"];
            this->iajk["mbij"] -= this->ijkl["mnij"]*T(1)["bn"];
            this->iajk["mbij"] += this->ia["me"]*T(2)["ebij"];

            this->aibj["amei"] -= this->ijab["mnef"]*T(2)["afin"];
            this->aibj["amei"] -= this->aibc["amfe"]*T(1)["fi"];
            this->aibj["amei"] += 0.5*this->ijka["nmie"]*T(1)["an"];

            this->abci["abej"] += 0.5*this->ijka["nmje"]*T(2)["abmn"];
            this->abci["abej"] -= this->aibj["amej"]*T(1)["bm"];
            this->abci["abej"] += this->aibc["amef"]*T(2)["fbmj"];
            this->abci["abej"] += this->abcd["abef"]*T(1)["fj"];
            this->abci["abej"] -= this->ia["me"]*T(2)["abmj"];

            this->aibj["amei"] += 0.5*this->ijka["nmie"]*T(1)["an"];

            this->abcd["abef"] += 0.5*this->ijab["mnef"]*Tau["abmn"];
            this->abcd["abef"] -= this->aibc["amef"]*T(1)["bm"];

            this->aibc["amef"] -= this->ijab["nmef"]*T(1)["an"];
        }
};

}
}

#endif
