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
STExcitationOperator<U,2>::STExcitationOperator(const std::string& name, const OneElectronOperator<U>& X, const ExcitationOperator<U,2>& T)
: ExcitationOperator<U,2>(name, X.arena, X.occ, X.vrt)
{
    transform(X, T, *this);
}

template <typename U>
STExcitationOperator<U,2>::STExcitationOperator(const std::string& name, const TwoElectronOperator<U>& X, const ExcitationOperator<U,2>& T)
: ExcitationOperator<U,2>(name, X.arena, X.occ, X.vrt)
{
    transform(X, T, *this);
}

template <typename U>
void STExcitationOperator<U,2>::transform(const OneElectronOperator<U>& X,
                                          const ExcitationOperator<U,2>& T,
                                                ExcitationOperator<U,2>& Z)
{
    OneElectronOperator<U> W("W", const_cast<OneElectronOperator<U>&>(X),
                             OneElectronOperator<U>::AB|OneElectronOperator<U>::IJ);
    transform(X, T, Z, W);
}

template <typename U>
void STExcitationOperator<U,2>::transform(const OneElectronOperator<U>& X,
                                          const ExcitationOperator<U,2>& T,
                                                ExcitationOperator<U,2>& Z,
                                                OneElectronOperator<U>& W)
{
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
    TwoElectronOperator<U> W("W", const_cast<TwoElectronOperator<U>&>(X), TwoElectronOperator<U>::AB|
                                                                     TwoElectronOperator<U>::IJ|
                                                                     TwoElectronOperator<U>::IA|
                                                                     TwoElectronOperator<U>::IJKL|
                                                                     TwoElectronOperator<U>::IJAK|
                                                                     TwoElectronOperator<U>::AIJK|
                                                                     TwoElectronOperator<U>::AIBJ);
    SpinorbitalTensor<U> Tau("Tau", T(2));
    transform(X, T, Tau, Z, W);
}

template <typename U>
void STExcitationOperator<U,2>::transform(const TwoElectronOperator<U>& X,
                                          const ExcitationOperator<U,2>& T,
                                                SpinorbitalTensor<U>& Tau,
                                                ExcitationOperator<U,2>& Z,
                                                TwoElectronOperator<U>& W)
{
    W.getAB() = X.getAB();
    W.getIJ() = X.getIJ();
    W.getIA() = X.getIA();
    W.getIJKL() = X.getIJKL();
    W.getIJAK() = X.getIJAK();
    W.getAIJK() = X.getAIJK();
    W.getAIBJ() = X.getAIBJ();

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

    Tau = T(2);
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

template <typename U>
void STExcitationOperator<U,2>::mytransform(const TwoElectronOperator<U>& X,
                                          const ExcitationOperator<U,2>& T,
                                                SpinorbitalTensor<U>& Tau,
                                                ExcitationOperator<U,2>& Z,
                                                TwoElectronOperator<U>& W)
{
    W.getAB() = X.getAB();
    W.getIJ() = X.getIJ();
    W.getIA() = X.getIA();
    W.getIJKL() = X.getIJKL();
    W.getIJAK() = X.getIJAK();
    W.getAIJK() = X.getAIJK();
    W.getAIBJ() = X.getAIBJ();

    SpinorbitalTensor<U>& FBC = W.getAB();
    SpinorbitalTensor<U>& FKJ = W.getIJ();
    SpinorbitalTensor<U>& WIJAB = W.getIJAB();
    SpinorbitalTensor<U>& WABCD = W.getABCD();
    SpinorbitalTensor<U>& WKLIJ = W.getIJKL();
    SpinorbitalTensor<U>& WBKCJ = W.getAIBJ();
    SpinorbitalTensor<U>& FAI = W.getAI();
    SpinorbitalTensor<U>& FKC = W.getIA();    
    SpinorbitalTensor<U>& WABIJ = W.getABIJ();
    SpinorbitalTensor<U>& WAKCD = W.getAIBC();
    SpinorbitalTensor<U>& WABCJ = W.getABCI();
    SpinorbitalTensor<U>& WKLCI = W.getIJAK();
    SpinorbitalTensor<U>& WAKIJ = W.getAIJK();
    

    Tau = T(2);
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    Z(0) = (U)0.0;

    FKC["kc"] += WIJAB["klcd"]*T(1)["dl"]; // CCS 14, through CCS 11 and CCS 6. CCS 10 through CCS 2. CCSD 17, through CCS 11 and CCD 3. CCSD 18 through CCSD 4 and CCD 2

    FKJ["kj"] += 0.5*WIJAB["klcd"]*T(2)["dclj"]; // CCD 9, through CCD 3 and CCS 9 through CCS 6
    FKJ["kj"] += FKC["kc"]*T(1)["cj"]; // CCS 11, through CCS 6. CCSD 3, through CCD 3
    FKJ["kj"] -= WKLCI["klcj"]*T(1)["cl"]; // CCS 13, through CCS 6. CCSD 10, through CCD 3. Sign change because <kl||ic> = -<kl||ci>

    WKLIJ["klij"] += 0.5*WIJAB["klcd"]*Tau["cdij"]; // CCD 7, CCSD 14, 15, 21 through CCD 5. 
    WKLIJ["klij"] += WKLCI["klcj"]*T(1)["ci"]; // CCSD 8, 20 through CCD 5

    WKLCI["klci"] += WIJAB["klcd"]*T(1)["di"]; // CCS 8, through CCS 4 and CCSD 16 through CCSD 6 and CCD 6.

    Z(1)["ai"]  = FAI["ai"]; // CCS 1
    Z(1)["ai"] -= T(1)["ck"]*WBKCJ["akci"]; // CCS 7
    Z(1)["ai"] += 0.5*WAKCD["akcd"]*Tau["cdik"]; // CCS 3, 12. Factor of 0.5 doesn't make sense for 12. Makes even less sense considering the 0.5 factor in the T(1)*T(1) part of Tau.
    Z(1)["ai"] += 0.5*WKLCI["klci"]*T(2)["ackl"]; // CCS 4. Sign change because <kl||ic> = -<kl||ci>
    Z(1)["ai"] += T(2)["acik"]*FKC["kc"]; // CCS 2
    Z(1)["ai"] += T(1)["ci"]*FBC["ac"]; // CCS 5
    Z(1)["ai"] -= T(1)["ak"]*FKJ["ki"]; // CCS 6

    FBC["bc"] -= 0.5*WIJAB["klcd"]*T(2)["dblk"]; // CCD 10, through CCD 2
    FBC["bc"] -= FKC["kc"]*T(1)["bk"]; // CCSD 4 through CCD 2
    FBC["bc"] += WAKCD["bkcd"]*T(1)["dk"]; // CCSD 9 through CCD 2

    WAKIJ["akij"] += 0.5*WAKCD["akcd"]*Tau["cdij"]; // CCSD 7, 19, through CCCD 2. Sign on my 19 seems wrong. 
    WAKIJ["akij"] -= WBKCJ["akci"]*T(1)["cj"]; // CCSD 13 through CCSD 2

    WBKCJ["bkcj"] -= 0.5*WIJAB["klcd"]*T(2)["bdjl"]; // CCD 8, through CCD 6
    WBKCJ["bkcj"] += WAKCD["bkcd"]*T(1)["dj"]; // CCSD 5, through CCD 6
    WBKCJ["bkcj"] += WKLCI["klcj"]*T(1)["bl"]; // CCSD 6, through CCD 6

    Z(2)["abij"]  = WABIJ["abij"]; // CCD 1
    Z(2)["abij"] += FBC["bc"]*T(2)["acij"]; // CCD 2
    Z(2)["abij"] -= FKJ["kj"]*T(2)["abik"]; // CCD 3
    Z(2)["abij"] += WABCJ["abcj"]*T(1)["ci"]; // CCSD 1
    Z(2)["abij"] -= WAKIJ["akij"]*T(1)["bk"]; // CCSD 2
    Z(2)["abij"] += 0.5*WABCD["abcd"]*Tau["cdij"]; // CCD 4, CCSD 11
    Z(2)["abij"] += 0.5*WKLIJ["klij"]*Tau["abkl"]; // CCD 5, CCSD 12
    Z(2)["abij"] -= WBKCJ["bkcj"]*T(2)["acik"]; // CCD 6
}

INSTANTIATE_SPECIALIZATIONS_2(STExcitationOperator,2);
