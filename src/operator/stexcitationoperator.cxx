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

#include "stexcitationoperator.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::tensor;

template <typename U>
STExcitationOperator<U,2>::STExcitationOperator(const OneElectronOperator<U>& X, const ExcitationOperator<U,2>& T)
: ExcitationOperator<U,2>(X.arena, X.occ, X.vrt)
{
    transform(X, T, *this);
}

template <typename U>
STExcitationOperator<U,2>::STExcitationOperator(const TwoElectronOperator<U>& X, const ExcitationOperator<U,2>& T)
: ExcitationOperator<U,2>(X.arena, X.occ, X.vrt)
{
    transform(X, T, *this);
}

template <typename U>
void STExcitationOperator<U,2>::transform(const OneElectronOperator<U>& X,
                                          const ExcitationOperator<U,2>& T,
                                                ExcitationOperator<U,2>& Z)
{
    OneElectronOperator<U> W(const_cast<OneElectronOperator<U>&>(X),
                             OneElectronOperator<U>::AB|OneElectronOperator<U>::IJ);

    SpinorbitalTensor<U>& FAI = W.getAI();
    SpinorbitalTensor<U>& FME = W.getIA();
    SpinorbitalTensor<U>& FAE = W.getAB();
    SpinorbitalTensor<U>& FMI = W.getIJ();

    Z(0) = FME["me"]*T(1)["em"];

    FMI["mi"] += FME["me"]*T(1)["ei"];

    Z(1)["ai"]  = FAI["ai"];
    Z(1)["ai"] += T(2)["aeim"]*FME["me"];
    Z(1)["ai"] += T(1)["ei"]*FAE["ae"];
    Z(1)["ai"] -= T(1)["am"]*FMI["mi"];

    FAE["ae"] -= FME["me"]*T(1)["am"];

    Z(2)["abij"]  = FAE["af"]*T(2)["fbij"];
    Z(2)["abij"] -= FMI["ni"]*T(2)["abnj"];
}

template <typename U>
void STExcitationOperator<U,2>::transform(const TwoElectronOperator<U>& X,
                                          const ExcitationOperator<U,2>& T,
                                                ExcitationOperator<U,2>& Z)
{
    TwoElectronOperator<U> W(const_cast<TwoElectronOperator<U>&>(X), TwoElectronOperator<U>::AB|
                                                                     TwoElectronOperator<U>::IJ|
                                                                     TwoElectronOperator<U>::IA|
                                                                     TwoElectronOperator<U>::IJKL|
                                                                     TwoElectronOperator<U>::IJAK|
                                                                     TwoElectronOperator<U>::AIJK|
                                                                     TwoElectronOperator<U>::AIBJ);

    SpinorbitalTensor<U>& FAI = W.getAI();
    SpinorbitalTensor<U>& FME = W.getIA();
    SpinorbitalTensor<U>& FAE = W.getAB();
    SpinorbitalTensor<U>& FMI = W.getIJ();
    SpinorbitalTensor<U>& WABIJ = W.getABIJ();
    SpinorbitalTensor<U>& WMNEF = W.getIJAB();
    SpinorbitalTensor<U>& WAMEF = W.getAIBC();
    SpinorbitalTensor<U>& WABEJ = W.getABCI();
    SpinorbitalTensor<U>& WABEF = W.getABCD();
    SpinorbitalTensor<U>& WMNIJ = W.getIJKL();
    SpinorbitalTensor<U>& WMNEJ = W.getIJAK();
    SpinorbitalTensor<U>& WAMIJ = W.getAIJK();
    SpinorbitalTensor<U>& WAMEI = W.getAIBJ();

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    Z(0) = (U)0.0;
    //Z(0)  = FME["me"]*T(1)["em"];
    //Z(0) += 0.25*WMNEF["mnef"]*Tau["efmn"];

    FME["me"] += WMNEF["mnef"]*T(1)["fn"];

    FMI["mi"] += 0.5*WMNEF["mnef"]*T(2)["efin"];
    FMI["mi"] += FME["me"]*T(1)["ei"];
    FMI["mi"] += WMNEJ["nmfi"]*T(1)["fn"];

    WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*Tau["efij"];
    WMNIJ["mnij"] += WMNEJ["mnej"]*T(1)["ei"];

    WMNEJ["mnej"] += WMNEF["mnef"]*T(1)["fj"];

    Z(1)["ai"]  = FAI["ai"];
    Z(1)["ai"] -= T(1)["em"]*WAMEI["amei"];
    Z(1)["ai"] += 0.5*WAMEF["amef"]*Tau["efim"];
    Z(1)["ai"] -= 0.5*WMNEJ["mnei"]*T(2)["eamn"];
    Z(1)["ai"] += T(2)["aeim"]*FME["me"];
    Z(1)["ai"] += T(1)["ei"]*FAE["ae"];
    Z(1)["ai"] -= T(1)["am"]*FMI["mi"];

    FAE["ae"] -= 0.5*WMNEF["mnef"]*T(2)["afmn"];
    FAE["ae"] -= FME["me"]*T(1)["am"];
    FAE["ae"] += WAMEF["amef"]*T(1)["fm"];

    WAMIJ["amij"] += 0.5*WAMEF["amef"]*Tau["efij"];
    WAMIJ["amij"] += WAMEI["amej"]*T(1)["ei"];

    WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
    WAMEI["amei"] += WAMEF["amef"]*T(1)["fi"];
    WAMEI["amei"] -= WMNEJ["nmei"]*T(1)["an"];

    Z(2)["abij"]  = WABIJ["abij"];
    Z(2)["abij"] += FAE["af"]*T(2)["fbij"];
    Z(2)["abij"] -= FMI["ni"]*T(2)["abnj"];
    Z(2)["abij"] += WABEJ["abej"]*T(1)["ei"];
    Z(2)["abij"] -= WAMIJ["amij"]*T(1)["bm"];
    Z(2)["abij"] += 0.5*WABEF["abef"]*Tau["efij"];
    Z(2)["abij"] += 0.5*WMNIJ["mnij"]*Tau["abmn"];
    Z(2)["abij"] -= WAMEI["amei"]*T(2)["ebmj"];
}

INSTANTIATE_SPECIALIZATIONS_2(STExcitationOperator,2);
