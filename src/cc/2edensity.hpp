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

#include "util/distributed.hpp"
#include "tensor/spinorbital.hpp"
#include "tensor/dist_tensor.hpp"
#include "scf/scf.hpp"
#include "operator/2eoperator.hpp"

#include "lambdaccsd.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class TwoElectronDensity : public op::TwoElectronOperator<U>
{
    protected:
        OneElectronDensity<U> D1;

    public:
        TwoElectronDensity(const scf::UHF<U>& uhf)
        : D1(uhf), op::TwoElectronOperator(D1, op::TwoElectronOperator::NONE)
        {
            //TODO: SCF density
            assert(0);
        }

        TwoElectronDensity(const DeexcitationOperator& L, const ExponentialOperator& T)
        : D1(L, T), op::TwoElectronOperator<U>(D1, op::TwoElectronOperator::NONE)
        {
            tensor::SpinorbitalTensor< tensor::DistTensor<U> > aitmp(this->ai);

            tensor::SpinorbitalTensor< tensor::DistTensor<U> > Tau(T[2]);
            Tau["abij"] += 0.5*T[1]["ai"]*T[1]["bj"];

            this->ijab["ijab"] = -L[2]["ijab"];

            this->ijka["ijka"] = L[2]["ijea"]*T[1]["ek"];

            this->aibc["aibc"] = -L[2]["mibc"]*T[1]["am"];

            this->ijkl["ijkl"] = -0.5*L[2]["ijef"]*Tau["efkl"];

            this->abcd["abcd"] = -0.5*L[2]["mncd"]*Tau["abmn"];

            this->aibj["aibj"] = L[1]["ib"]*T[1]["aj"];
            this->aibj["aibj"] += L[2]["imbe"]*T[2]["aejm"];
            this->aibj["aibj"] -= this->ijka["mijb"]*T[1]["am"];

            this->iajk["iajk"] = this->ij["ij"]*T[1]["ak"];

            this->abci["abci"] = -L[1]["mc"]*T[2]["abmi"];
            this->abci["abci"] += this->ab["ac"]*T[1]["bi"];
            this->abci["abci"] += this->aibc["amce"]*T[2]["beim"];
            this->abci["abci"] -= 0.5*this->ijka["nmic"]*Tau["abmn"];

            this->abij["abij"] = Tau["abij"];
            this->abij["abij"] += this->ijkl["mnij"]*Tau["abmn"];
            this->abij["abij"] += this->ij["mi"]*T[2]["abmj"];
            this->abij["abij"] -= this->ab["ae"]*T[2]["ebij"];
            this->abij["abij"] -= this->aibj["amej"]*T[2]["ebmj"];
            this->abij["abij"] += this->iajk["mbij"]*T[1]["am"];
            aitmp["ai"] = -this->aibc["amef"]*T[2]["efim"];
            this->abij["abij"] += aitmp["ai"]*T[1]["bj"];

            this->iajk["iajk"] += L[1]["ie"]*T[2]["eajk"];
            this->iajk["iajk"] += this->ijka["imje"]*T[2]["aekm"];
            this->iajk["iajk"] += 0.5*this->aibc["aife"]*Tau["efjk"];
        }
};

}
}

#endif
