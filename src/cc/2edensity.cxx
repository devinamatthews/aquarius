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

#include "2edensity.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::op;
using namespace aquarius::scf;
using namespace aquarius::tensor;

template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const UHF<U>& uhf)
: TwoElectronOperator<U>(OneElectronDensity<U>(uhf))
{
    this->abcd["abcd"]  =     this->ab["ab"]*this->ab["cd"];
    this->abcd["abcd"] += 0.5*this->ab["ac"]*this->ab["bd"];

    this->abci["abci"]  =     this->ab["ab"]*this->ai["ci"];
    this->abci["abci"] += 0.5*this->ab["ac"]*this->ai["bi"];

    this->abij["abij"]  =     this->ab["ab"]*this->ij["ij"];
    this->abij["abij"] += 0.5*this->ai["ai"]*this->ai["bj"];

    this->aibj["aibj"]  =     this->ai["ai"]*this->ai["bj"];
    this->aibj["aibj"] += 0.5*this->ab["ab"]*this->ij["ij"];

    this->ijka["ijka"]  =     this->ij["ij"]*this->ia["ka"];
    this->ijka["ijka"] += 0.5*this->ij["ik"]*this->ia["ja"];

    this->ijkl["ijkl"]  =     this->ij["ij"]*this->ij["kl"];
    this->ijkl["ijkl"] += 0.5*this->ij["ik"]*this->ij["jl"];
}

/*
 * Form the unrelaxed CCSD Density
 */
template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const ExponentialOperator<U,2>& T)
: TwoElectronOperator<U>(OneElectronDensity<U>(T))
{
    this->abij["abij"]  = T(2)["abij"];
    this->abij["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];
}

/*
 * Form the mixed perturbed CCSD Density
 */
template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const DeexcitationOperator<U,2>& L,
                                          const ExponentialOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: TwoElectronOperator<U>(OneElectronDensity<U>(L, T, TA))
{
    TwoElectronOperator<U> I(this->uhf);

    SpinorbitalTensor<U>& IIJ = I.getIJ();
    SpinorbitalTensor<U>& IAB = I.getAB();
    SpinorbitalTensor<U>& IAI = I.getAI();

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    SpinorbitalTensor<U> TauA(TA(2));
    TauA["abij"] += T(1)["ai"]*TA(1)["bj"];

    this->ijka["ijka"] = -L(2)["ijea"]*TA(1)["ek"];

    this->aibc["aibc"] = L(2)["mibc"]*TA(1)["am"];

    this->ijkl["ijkl"] = 0.5*L(2)["ijef"]*TauA["efkl"];

    this->abcd["abcd"] = 0.5*L(2)["mncd"]*TauA["abmn"];

    this->aibj["aibj"] -= L(1)["ib"]*TA(1)["aj"];
    this->aibj["aibj"] -= L(2)["imbe"]*TA(2)["aejm"];

    IAB["ab"]  = 0.5*T(2)["aemn"]*L(2)["mnbe"];
    IIJ["ij"]  = 0.5*T(2)["efim"]*L(2)["jmef"];
    IAI["ai"]  = IAB["ae"]*TA(1)["ei"];
    IAI["ai"] += IIJ["mi"]*TA(1)["am"];
    IAI["ai"] += L(1)["me"]*TA(2)["aeim"];

    this->abci["abci"] += L(1)["mc"]*TA(2)["abmi"];
    this->abci["abci"] += this->ab["ac"]*T(1)["bi"];
    this->abci["abci"] += IAB["ac"]*TA(1)["bi"];
    this->abci["abci"] += this->aibc["amce"]*T(2)["beim"];
    this->abci["abci"] -= this->aibj["amci"]*T(1)["bm"];

    this->iajk["iajk"] -= L(1)["ie"]*TA(2)["eajk"];
    this->iajk["iajk"] += this->ij["ij"]*T(1)["ak"];
    this->iajk["iajk"] -= IIJ["ij"]*TA(1)["ak"];
    this->iajk["iajk"] += this->ijka["imje"]*T(2)["aekm"];

    this->abij["abij"] += TauA["abij"];
    this->abij["abij"] += this->ijkl["mnij"]*Tau["abmn"];
    this->abij["abij"] += this->abcd["abef"]*Tau["efij"];
    this->abij["abij"] += this->ij["mi"]*T(2)["abmj"];
    this->abij["abij"] -= this->ab["ae"]*T(2)["ebij"];
    this->abij["abij"] -= IIJ["mi"]*TA(2)["abmj"];
    this->abij["abij"] -= IAB["ae"]*TA(2)["ebij"];
    this->abij["abij"] -= this->aibj["amej"]*T(2)["ebmj"];
    this->abij["abij"] += this->iajk["mbij"]*T(1)["am"];
    this->abij["abij"] -= this->abci["abej"]*T(1)["ei"];
    this->abij["abij"] += IAI["ai"]*T(1)["bj"];

    this->abci["abci"] -= 0.5*this->ijka["nmic"]*Tau["abmn"];
    this->abci["abci"] -= this->abcd["abce"]*T(1)["ei"];

    this->iajk["iajk"] += 0.5*this->aibc["aife"]*Tau["efjk"];
    this->iajk["iajk"] += this->aibj["aiej"]*T(1)["ek"];
    this->iajk["iajk"] += this->ijkl["imjk"]*T(1)["am"];

    this->aibj["aibj"] += this->ijka["imjb"]*T(1)["am"];
    this->aibj["aibj"] -= this->aibc["aibe"]*T(1)["ej"];
}

/*
 * Form the relaxed CCSD Density
 */
template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const DeexcitationOperator<U,2>& L,
                                          const ExponentialOperator<U,2>& T)
: TwoElectronOperator<U>(OneElectronDensity<U>(L, T))
{
    SpinorbitalTensor<U> aitmp(this->ai);

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    this->ijab["ijab"] =  L(2)["ijab"];

    this->ijka["ijka"] = -L(2)["ijea"]*T(1)["ek"];

    this->aibc["aibc"] =  L(2)["mibc"]*T(1)["am"];

    this->ijkl["ijkl"] =  0.5*L(2)["ijef"]*Tau["efkl"];

    this->abcd["abcd"] =  0.5*L(2)["mncd"]*Tau["abmn"];

    this->aibj["aibj"]  = -L(1)["ib"]*T(1)["aj"];
    this->aibj["aibj"] -=  L(2)["imbe"]*T(2)["aejm"];
    this->aibj["aibj"] +=  this->ijka["imjb"]*T(1)["am"];

    this->iajk["iajk"]  = this->ij["ij"]*T(1)["ak"];

    this->abci["abci"]  = L(1)["mc"]*T(2)["abmi"];
    this->abci["abci"] += this->ab["ac"]*T(1)["bi"];
    this->abci["abci"] += this->aibc["amce"]*T(2)["beim"];
    this->abci["abci"] -= 0.5*this->ijka["nmic"]*Tau["abmn"];

    this->abij["abij"]  = Tau["abij"];
    this->abij["abij"] += this->ijkl["mnij"]*Tau["abmn"];
    this->abij["abij"] += this->ij["mi"]*T(2)["abmj"];
    this->abij["abij"] -= this->ab["ae"]*T(2)["ebij"];
    this->abij["abij"] -= this->aibj["amej"]*T(2)["ebmj"];
    this->abij["abij"] += this->iajk["mbij"]*T(1)["am"];
    aitmp["ai"] = -this->aibc["amef"]*T(2)["efim"];
    this->abij["abij"] += aitmp["ai"]*T(1)["bj"];

    this->iajk["iajk"] -= L(1)["ie"]*T(2)["eajk"];
    this->iajk["iajk"] += this->ijka["imje"]*T(2)["aekm"];
    this->iajk["iajk"] += 0.5*this->aibc["aife"]*Tau["efjk"];
}

/*
 * Form the relaxed perturbed CCSD Density
 */
template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const DeexcitationOperator<U,2>& L,
                                          const DeexcitationOperator<U,2>& LA,
                                          const ExponentialOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: TwoElectronOperator<U>(OneElectronDensity<U>(L, LA, T, TA))
{
    //TODO
    assert(0);
}

template <typename U>
double TwoElectronDensity<U>::getS2() const
{
    //TODO
    assert(0);
}

template <typename U>
double TwoElectronDensity<U>::getProjectedMultiplicity() const
{
    return sqrt(1+4*getS2());
}

INSTANTIATE_SPECIALIZATIONS(TwoElectronDensity);
