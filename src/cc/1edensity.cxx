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

#include "1edensity.hpp"

using namespace std;
using namespace aquarius::cc;
using namespace aquarius::op;
using namespace aquarius::scf;
using namespace aquarius::tensor;

template <typename U>
OneElectronDensity<U>::OneElectronDensity(const UHF<U>& uhf)
: OneElectronOperator<U>(uhf)
{
    int N = uhf.getMolecule().getNumOrbitals();
    int nI = uhf.getMolecule().getNumAlphaElectrons();
    int ni = uhf.getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;

    vector<int> NN = vec(NS,NS);
    vector<int> sizeII = vec(nI,nI);
    vector<int> sizeii = vec(ni,ni);
    vector<int> sizeAI = vec(nA,nI);
    vector<int> sizeai = vec(na,ni);

    const DistTensor<U>& CA = uhf.getCA();
    const DistTensor<U>& Ca = uhf.getCa();
    const DistTensor<U>& CI = uhf.getCI();
    const DistTensor<U>& Ci = uhf.getCi();

    DistTensor<U> DIJ(this->arena, 2, sizeII, NN);
    DistTensor<U> Dij(this->arena, 2, sizeii, NN);
    DistTensor<U> DAI(this->arena, 2, sizeAI, NN);
    DistTensor<U> Dai(this->arena, 2, sizeai, NN);

    DIJ["IJ"] = CI["pI"]*CI["pJ"];
    Dij["ij"] = Ci["pi"]*Ci["pj"];
    DAI["AI"] = CA["pA"]*CI["pI"];
    Dai["ai"] = Ca["pa"]*Ci["pi"];

    this->ab(1,0,1,0)["AB"] = DAI["AI"]*DAI["BI"];
    this->ab(0,0,0,0)["ab"] = Dai["ai"]*Dai["bi"];

    this->ij(0,1,0,1)["IJ"] = DIJ["IK"]*DIJ["JK"];
    this->ij(0,0,0,0)["ij"] = Dij["ik"]*Dij["jk"];

    this->ai(1,0,0,1)["AI"] = DAI["AJ"]*DIJ["IJ"];
    this->ai(0,0,0,0)["ai"] = Dai["aj"]*Dij["ij"];
}

/*
 * Form the unrelaxed CCSD density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const ExponentialOperator<U,2>& T)
: OneElectronOperator<U>(T.getSCF(), false)
{
    this->ai["ai"] = T(1)["ai"];
}

/*
 * Form the partial perturbed CCSD Density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const DeexcitationOperator<U,2>& L,
                                          const ExponentialOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: OneElectronOperator<U>(L.getSCF(), false)
{
    OneElectronOperator<U> I(this->uhf);

    SpinorbitalTensor<U>& IIJ = I.getIJ();
    SpinorbitalTensor<U>& IAB = I.getAB();

    IAB["ab"] = 0.5*T(2)["aemn"]*L(2)["mnbe"];
    IIJ["ij"] = 0.5*T(2)["efim"]*L(2)["jmef"];

    this->ab["ab"] += TA(1)["am"]*L(1)["mb"];
    this->ab["ab"] += 0.5*TA(2)["aemn"]*L(2)["mnbe"];

    this->ij["ij"] -= TA(1)["ei"]*L(1)["je"];
    this->ij["ij"] -= 0.5*TA(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] += TA(1)["ai"];
    this->ai["ai"] += TA(2)["aeim"]*L(1)["me"];
    this->ai["ai"] += this->ij["mi"]*T(1)["am"];
    this->ai["ai"] -= this->ab["ae"]*T(1)["ei"];
    this->ai["ai"] -= IIJ["mi"]*TA(1)["am"];
    this->ai["ai"] -= IAB["ae"]*TA(1)["ei"];
}

/*
 * Form the relaxed CCSD density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const DeexcitationOperator<U,2>& L,
                                          const ExponentialOperator<U,2>& T)
: OneElectronOperator<U>(L.getSCF(), false)
{
    this->ia["ia"] += L(1)["ia"];

    this->ab["ab"] += 0.5*T(2)["aemn"]*L(2)["mnbe"];

    this->ij["ij"] -= T(1)["ei"]*L(1)["je"];
    this->ij["ij"] -= 0.5*T(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] += T(1)["ai"];
    this->ai["ai"] += T(2)["aeim"]*L(1)["me"];
    this->ai["ai"] += this->ij["mi"]*T(1)["am"];
    this->ai["ai"] -= this->ab["ae"]*T(1)["ei"];

    this->ab["ab"] += T(1)["am"]*L(1)["mb"];
}

/*
 * Form the relaxed perturbed CCSD Density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const DeexcitationOperator<U,2>& L,
                                          const DeexcitationOperator<U,2>& LA,
                                          const ExponentialOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: OneElectronOperator<U>(L.getSCF(), false)
{
    OneElectronOperator<U> I(this->uhf);

    SpinorbitalTensor<U>& IIJ = I.getIJ();
    SpinorbitalTensor<U>& IAB = I.getAB();

    this->ia["ia"] += LA(1)["ia"];

    this->ab["ab"] += 0.5*T(2)["aemn"]*LA(2)["mnbe"];

    this->ij["ij"] -= T(1)["ei"]*LA(1)["je"];
    this->ij["ij"] -= 0.5*T(2)["efim"]*LA(2)["jmef"];

    this->ai["ai"] += T(1)["ai"];
    this->ai["ai"] += T(2)["aeim"]*LA(1)["me"];
    this->ai["ai"] += this->ij["mi"]*T(1)["am"];
    this->ai["ai"] -= this->ab["ae"]*T(1)["ei"];

    this->ab["ab"] += T(1)["am"]*LA(1)["mb"];

    IAB["ab"]  = TA(1)["am"]*L(1)["mb"];
    IAB["ab"] += 0.5*TA(2)["aemn"]*L(2)["mnbe"];

    this->ab["ab"] += IAB["ab"];
    this->ij["ij"] -= IIJ["ij"];

    IIJ["ij"]  = TA(1)["ei"]*L(1)["je"];
    IIJ["ij"] += 0.5*TA(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] += TA(1)["ai"];
    this->ai["ai"] += TA(2)["aeim"]*L(1)["me"];
    this->ai["ai"] -= IIJ["mi"]*T(1)["am"];
    this->ai["ai"] -= IAB["ae"]*T(1)["ei"];

    IAB["ab"] = 0.5*T(2)["aemn"]*L(2)["mnbe"];
    IIJ["ij"] = 0.5*T(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] -= IIJ["mi"]*TA(1)["am"];
    this->ai["ai"] -= IAB["ae"]*TA(1)["ei"];
}

INSTANTIATE_SPECIALIZATIONS(OneElectronDensity);
