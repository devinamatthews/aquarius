/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must reRin the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documenRtion and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANRBILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENRL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DAR, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_OPERATOR_PERTURBEDST2EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_PERTURBEDST2EOPERATOR_HPP_

#include "2eoperator.hpp"
#include "st2eoperator.hpp"
#include "exponentialoperator.hpp"

namespace aquarius
{
namespace op
{

/*
 *      _A    -T  A  T    _   A      A  T     _  A
 * Form X  = e   X  e  + [X, T ] = (X  e ) + (X T ) , up to two-electron terms
 *                                        c        c
 */
template <typename U, int nex> class PerturbedSTTwoElectronOperator;

template <typename U>
class PerturbedSTTwoElectronOperator<U,2> : public STTwoElectronOperator<U,2>
{
    protected:
        const STTwoElectronOperator<U,2>& X;
        const ExcitationOperator<U,2>& TA;

        void initialize(const STTwoElectronOperator<U,2>& X,
                        const ExponentialOperator<U,2>& T,
                        const ExcitationOperator<U,2>& TA)
        {
            OneElectronOperator<U> I(this->uhf);

            tensor::SpinorbitalTensor<U>& IMI = I.getIJ();
            tensor::SpinorbitalTensor<U>& IAE = I.getAB();
            tensor::SpinorbitalTensor<U>& IME = I.getIA();

            IME["me"]  = X.ijab["mnef"]*TA(1)["fn"];

            IMI["mi"]  = X.ijka["mnie"]*TA(1)["en"];
            IMI["mi"] += 0.5*X.ijab["mnef"]*TA(2)["efin"];

            IAE["ae"]  = X.aibc["amef"]*TA(1)["fm"];
            IAE["ae"] -= 0.5*X.ijab["mnef"]*TA(2)["afmn"];

            this->ia["ia"] += IME["ia"];

            this->ij["ij"] += X.ia["ie"]*TA(1)["ej"];
            this->ij["ij"] += IMI["ij"];

            this->ab["ab"] -= X.ia["mb"]*TA(1)["am"];
            this->ab["ab"] += IAE["ab"];

            this->ai["ai"] += X.ab["ae"]*TA(1)["ei"];
            this->ai["ai"] -= X.ij["mi"]*TA(1)["am"];
            this->ai["ai"] += X.ia["me"]*TA(2)["aeim"];
            this->ai["ai"] -= X.aibj["amei"]*TA(1)["em"];
            this->ai["ai"] -= 0.5*X.ijka["mnie"]*TA(2)["aemn"];
            this->ai["ai"] += 0.5*X.aibc["amef"]*TA(2)["efim"];

            this->ijka["ijka"] += X.ijab["ijea"]*TA(1)["ek"];

            this->aibc["aibc"] -= X.ijab["mibc"]*TA(1)["am"];

            this->ijkl["ijkl"] += X.ijka["ijke"]*TA(1)["el"];
            this->ijkl["ijkl"] += 0.5*X.ijab["ijef"]*TA(2)["efkl"];

            this->abcd["abcd"] -= X.aibc["amcd"]*TA(1)["bm"];
            this->abcd["abcd"] += 0.5*X.ijab["mncd"]*TA(2)["abmn"];

            this->aibj["aibj"] += X.aibc["aibe"]*TA(1)["ej"];
            this->aibj["aibj"] -= X.ijka["imjb"]*TA(1)["am"];
            this->aibj["aibj"] -= X.ijab["mibe"]*TA(2)["aemj"];

            this->iajk["iajk"] += IME["ie"]*T(2)["eajk"];
            this->iajk["iajk"] += X.aibj["aiej"]*TA(1)["ek"];
            this->iajk["iajk"] -= X.ijkl["imjk"]*TA(1)["am"];
            this->iajk["iajk"] += X.ijka["iejm"]*TA(2)["aekm"];
            this->iajk["iajk"] += 0.5*X.aibc["aife"]*TA(2)["efjk"];

            this->abci["abci"] -= IME["mc"]*T(2)["abmi"];
            this->abci["abci"] -= X.aibj["amci"]*TA(1)["bm"];
            this->abci["abci"] += X.abcd["abce"]*TA(1)["ei"];
            this->abci["abci"] += X.aibc["amce"]*TA(2)["beim"];
            this->abci["abci"] += 0.5*X.ijka["nmic"]*TA(2)["abmn"];

            this->abij["abij"] += IAE["ae"]*T(2)["ebij"];
            this->abij["abij"] -= IMI["mi"]*T(2)["abmj"];
            this->abij["abij"] += X.abci["abej"]*TA(1)["ei"];
            this->abij["abij"] -= X.iajk["mbij"]*TA(1)["am"];
            this->abij["abij"] += 0.5*X.abcd["abef"]*TA(2)["efij"];
            this->abij["abij"] += 0.5*X.ijkl["mnij"]*TA(2)["abmn"];
            this->abij["abij"] -= X.aibj["amei"]*TA(2)["ebmj"];
        }

    public:
        PerturbedSTTwoElectronOperator(const STTwoElectronOperator<U,2>& X, const OneElectronOperator<U>& XA,
                                       const ExponentialOperator<U,2>& T, const ExcitationOperator<U,2>& TA)
        : STTwoElectronOperator<U,2>(XA, T), X(X), TA(TA)
        {
            initialize(X, T, TA);
        }

        PerturbedSTTwoElectronOperator(const STTwoElectronOperator<U,2>& X, const TwoElectronOperator<U>& XA,
                const ExponentialOperator<U,2>& T, const ExcitationOperator<U,2>& TA)
        : STTwoElectronOperator<U,2>(XA, T), X(X), TA(TA)
        {
            initialize(X, T, TA);
        }

        void contract(const ExcitationOperator<U,2>& R, ExcitationOperator<U,2>& Z, bool connected=true) const
        {
            STTwoElectronOperator<U,2>::contract(R, Z, connected);

            op::OneElectronOperator<U> I(this->uhf);

            tensor::SpinorbitalTensor<U>& IMI = I.getIJ();
            tensor::SpinorbitalTensor<U>& IAE = I.getAB();

            IMI["mi"] = X.getIJKA()["mnie"]*R(1)["en"];
            IAE["ae"] = X.getAIBC()["amef"]*R(1)["fm"];

            Z(2)["abij"] += IAE["ae"]*TA(2)["ebij"];
            Z(2)["abij"] -= IMI["mi"]*TA(2)["abmj"];
        }

        void contract(const DeexcitationOperator<U,2>& L, DeexcitationOperator<U,2>& Z, bool connected=false) const
        {
            STTwoElectronOperator<U,2>::contract(L, Z, connected);

            op::OneElectronOperator<U> I(this->uhf);

            tensor::SpinorbitalTensor<U>& IMN = I.getIJ();
            tensor::SpinorbitalTensor<U>& IEF = I.getAB();

            IMN["mn"] =  0.5*L(2)["moef"]*TA(2)["efno"];
            IEF["ef"] = -0.5*L(2)["mnfg"]*TA(2)["egmn"];

            Z(1)["ia"] -= IMN["mn"]*X.getIJKA()["nima"];
            Z(1)["ia"] -= IEF["ef"]*X.getAIBC()["fiea"];
        }
};

}
}

#endif
