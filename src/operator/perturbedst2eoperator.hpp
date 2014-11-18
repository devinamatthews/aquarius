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

namespace aquarius
{
namespace op
{

/*
 *      _A    -T  A  T    _   A      A  T     _  A
 * Form X  = e   X  e  + [X, T ] = (X  e ) + (X T ) , up to two-electron terms
 *                                        c        c
 */
template <typename U>
class PerturbedSTTwoElectronOperator : public STTwoElectronOperator<U>
{
    protected:
        const STTwoElectronOperator<U>& X;
        const ExcitationOperator<U,2>& TA;

    public:
        template <int N>
        PerturbedSTTwoElectronOperator(const std::string& name, const STTwoElectronOperator<U>& X, const OneElectronOperator<U>& XA,
                                       const ExcitationOperator<U,N>& T, const ExcitationOperator<U,N>& TA)
        : STTwoElectronOperator<U>(name, XA, T), X(X), TA(TA)
        {
            OneElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

            tensor::SpinorbitalTensor<U>& IMI = I.getIJ();
            tensor::SpinorbitalTensor<U>& IAE = I.getAB();
            tensor::SpinorbitalTensor<U>& IME = I.getIA();

            IME["me"]  = X.getIJAB()["mnef"]*TA(1)["fn"];

            IMI["mi"]  = X.getIJAK()["nmei"]*TA(1)["en"];
            IMI["mi"] += 0.5*X.getIJAB()["mnef"]*TA(2)["efin"];

            IAE["ae"]  = X.getAIBC()["amef"]*TA(1)["fm"];
            IAE["ae"] -= 0.5*X.getIJAB()["mnef"]*TA(2)["afmn"];

            this->ia["ia"] += IME["ia"];

            this->ij["ij"] += X.getIA()["ie"]*TA(1)["ej"];
            this->ij["ij"] += IMI["ij"];

            this->ab["ab"] -= X.getIA()["mb"]*TA(1)["am"];
            this->ab["ab"] += IAE["ab"];

            this->ai["ai"] += X.getAB()["ae"]*TA(1)["ei"];
            this->ai["ai"] -= X.getIJ()["mi"]*TA(1)["am"];
            this->ai["ai"] += X.getIA()["me"]*TA(2)["aeim"];
            this->ai["ai"] -= X.getAIBJ()["amei"]*TA(1)["em"];
            this->ai["ai"] -= 0.5*X.getIJAK()["nmei"]*TA(2)["aemn"];
            this->ai["ai"] += 0.5*X.getAIBC()["amef"]*TA(2)["efim"];

            this->getIJAK()["ijak"] += X.getIJAB()["ijae"]*TA(1)["ek"];

            this->getAIBC()["aibc"] -= X.getIJAB()["mibc"]*TA(1)["am"];

            this->getIJKL()["ijkl"] += X.getIJAK()["jiek"]*TA(1)["el"];
            this->getIJKL()["ijkl"] += 0.5*X.getIJAB()["ijef"]*TA(2)["efkl"];

            this->getABCD()["abcd"] -= X.getAIBC()["amcd"]*TA(1)["bm"];
            this->getABCD()["abcd"] += 0.5*X.getIJAB()["mncd"]*TA(2)["abmn"];

            this->getAIBJ()["aibj"] += X.getAIBC()["aibe"]*TA(1)["ej"];
            this->getAIBJ()["aibj"] -= X.getIJAK()["mibj"]*TA(1)["am"];
            this->getAIBJ()["aibj"] -= X.getIJAB()["mibe"]*TA(2)["aemj"];

            this->getAIJK()["aijk"] += IME["ie"]*T(2)["aejk"];
            this->getAIJK()["aijk"] += X.getAIBJ()["aiek"]*TA(1)["ej"];
            this->getAIJK()["aijk"] -= X.getIJKL()["mijk"]*TA(1)["am"];
            this->getAIJK()["aijk"] += X.getIJAK()["miek"]*TA(2)["aejm"];
            this->getAIJK()["aijk"] += 0.5*X.getAIBC()["aief"]*TA(2)["efjk"];

            this->getABCI()["abci"] -= IME["mc"]*T(2)["abmi"];
            this->getABCI()["abci"] -= X.getAIBJ()["amci"]*TA(1)["bm"];
            this->getABCI()["abci"] += X.getABCD()["abce"]*TA(1)["ei"];
            this->getABCI()["abci"] += X.getAIBC()["amce"]*TA(2)["beim"];
            this->getABCI()["abci"] += 0.5*X.getIJAK()["mnci"]*TA(2)["abmn"];

            this->abij["abij"] += IAE["ae"]*T(2)["ebij"];
            this->abij["abij"] -= IMI["mi"]*T(2)["abmj"];
            this->abij["abij"] += X.getABCI()["abej"]*TA(1)["ei"];
            this->abij["abij"] -= X.getAIJK()["bmji"]*TA(1)["am"];
            this->abij["abij"] += 0.5*X.getABCD()["abef"]*TA(2)["efij"];
            this->abij["abij"] += 0.5*X.getIJKL()["mnij"]*TA(2)["abmn"];
            this->abij["abij"] -= X.getAIBJ()["amei"]*TA(2)["ebmj"];
        }

        template <int N>
        PerturbedSTTwoElectronOperator(const std::string& name, const STTwoElectronOperator<U>& X, const TwoElectronOperator<U>& XA,
                                       const ExcitationOperator<U,N>& T, const ExcitationOperator<U,N>& TA)
        : STTwoElectronOperator<U>(name, XA, T), X(X), TA(TA)
        {
            OneElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

            tensor::SpinorbitalTensor<U>& IMI = I.getIJ();
            tensor::SpinorbitalTensor<U>& IAE = I.getAB();
            tensor::SpinorbitalTensor<U>& IME = I.getIA();

            IME["me"]  = X.getIJAB()["mnef"]*TA(1)["fn"];

            IMI["mi"]  = X.getIJAK()["nmei"]*TA(1)["en"];
            IMI["mi"] += 0.5*X.getIJAB()["mnef"]*TA(2)["efin"];

            IAE["ae"]  = X.getAIBC()["amef"]*TA(1)["fm"];
            IAE["ae"] -= 0.5*X.getIJAB()["mnef"]*TA(2)["afmn"];

            this->ia["ia"] += IME["ia"];

            this->ij["ij"] += X.getIA()["ie"]*TA(1)["ej"];
            this->ij["ij"] += IMI["ij"];

            this->ab["ab"] -= X.getIA()["mb"]*TA(1)["am"];
            this->ab["ab"] += IAE["ab"];

            this->ai["ai"] += X.getAB()["ae"]*TA(1)["ei"];
            this->ai["ai"] -= X.getIJ()["mi"]*TA(1)["am"];
            this->ai["ai"] += X.getIA()["me"]*TA(2)["aeim"];
            this->ai["ai"] -= X.getAIBJ()["amei"]*TA(1)["em"];
            this->ai["ai"] -= 0.5*X.getIJAK()["nmei"]*TA(2)["aemn"];
            this->ai["ai"] += 0.5*X.getAIBC()["amef"]*TA(2)["efim"];

            this->getIJAK()["ijak"] += X.getIJAB()["ijae"]*TA(1)["ek"];

            this->getAIBC()["aibc"] -= X.getIJAB()["mibc"]*TA(1)["am"];

            this->getIJKL()["ijkl"] += X.getIJAK()["jiek"]*TA(1)["el"];
            this->getIJKL()["ijkl"] += 0.5*X.getIJAB()["ijef"]*TA(2)["efkl"];

            this->getABCD()["abcd"] -= X.getAIBC()["amcd"]*TA(1)["bm"];
            this->getABCD()["abcd"] += 0.5*X.getIJAB()["mncd"]*TA(2)["abmn"];

            this->getAIBJ()["aibj"] += X.getAIBC()["aibe"]*TA(1)["ej"];
            this->getAIBJ()["aibj"] -= X.getIJAK()["mibj"]*TA(1)["am"];
            this->getAIBJ()["aibj"] -= X.getIJAB()["mibe"]*TA(2)["aemj"];

            this->getAIJK()["aijk"] += IME["ie"]*T(2)["aejk"];
            this->getAIJK()["aijk"] += X.getAIBJ()["aiek"]*TA(1)["ej"];
            this->getAIJK()["aijk"] -= X.getIJKL()["mijk"]*TA(1)["am"];
            this->getAIJK()["aijk"] += X.getIJAK()["miek"]*TA(2)["aejm"];
            this->getAIJK()["aijk"] += 0.5*X.getAIBC()["aief"]*TA(2)["efjk"];

            this->getABCI()["abci"] -= IME["mc"]*T(2)["abmi"];
            this->getABCI()["abci"] -= X.getAIBJ()["amci"]*TA(1)["bm"];
            this->getABCI()["abci"] += X.getABCD()["abce"]*TA(1)["ei"];
            this->getABCI()["abci"] += X.getAIBC()["amce"]*TA(2)["beim"];
            this->getABCI()["abci"] += 0.5*X.getIJAK()["mnci"]*TA(2)["abmn"];

            this->abij["abij"] += IAE["ae"]*T(2)["ebij"];
            this->abij["abij"] -= IMI["mi"]*T(2)["abmj"];
            this->abij["abij"] += X.getABCI()["abej"]*TA(1)["ei"];
            this->abij["abij"] -= X.getAIJK()["bmji"]*TA(1)["am"];
            this->abij["abij"] += 0.5*X.getABCD()["abef"]*TA(2)["efij"];
            this->abij["abij"] += 0.5*X.getIJKL()["mnij"]*TA(2)["abmn"];
            this->abij["abij"] -= X.getAIBJ()["amei"]*TA(2)["ebmj"];
        }
};

}
}

#endif
