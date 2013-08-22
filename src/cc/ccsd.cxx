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

#include "ccsd.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::scf;
using namespace aquarius::input;
using namespace aquarius::tensor;

template <typename U>
CCSD<U>::CCSD(const Config& config, TwoElectronOperator<U>& H)
: Iterative(config), ExponentialOperator<U,2>(H.getSCF()),
  T(*this), D(H.getSCF()), Z(H.getSCF()),
  H(H), diis(config.get("diis"))
{
    D(0) = (U)1.0;
    D(1)["ai"]  = H.getIJ()["ii"];
    D(1)["ai"] -= H.getAB()["aa"];
    D(2)["abij"]  = H.getIJ()["ii"];
    D(2)["abij"] += H.getIJ()["jj"];
    D(2)["abij"] -= H.getAB()["aa"];
    D(2)["abij"] -= H.getAB()["bb"];

    D = 1/D;

    Z(0) = (U)0.0;
    T(0) = (U)0.0;
    T(1) = H.getAI()*D(1);
    T(2) = H.getABIJ()*D(2);

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    energy = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));

    conv =          T(1)(0).reduce(CTF_OP_MAXABS);
    conv = max(conv,T(1)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,T(2)(0).reduce(CTF_OP_MAXABS));
    conv = max(conv,T(2)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,T(2)(2).reduce(CTF_OP_MAXABS));
}

template <typename U>
void CCSD<U>::_iterate()
{
    STExcitationOperator<U,2>::transform(H, T, Z);
    //Z(0) = (U)0.0;
    T.set_name("T");
    Z.set_name("Z");
    H.set_name("H");

    Z *= D;
    T += Z;

    SpinorbitalTensor<U> Tau(T(2));
    Tau.set_name("Tau");
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    energy = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));

    conv =          Z(1)(0).reduce(CTF_OP_MAXABS);
    conv = max(conv,Z(1)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(0).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(2).reduce(CTF_OP_MAXABS));

    diis.extrapolate(T, Z);
}

template <typename U>
double CCSD<U>::getProjectedS2(const UHF<U>& uhf,
                               const SpinorbitalTensor<U>& T1,
                               const SpinorbitalTensor<U>& T2)
{
    tCTF_World<U>& ctf = T1(0).ctf;

    int N = uhf.getMolecule().getNumOrbitals();
    int nI = uhf.getMolecule().getNumAlphaElectrons();
    int ni = uhf.getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;

    int shapeNN[] = {NS,NS};
    int shapeNNNN[] = {NS,NS,NS,NS};
    int sizeAI[] = {nA,nI};
    int sizeAi[] = {nA,ni};
    int sizeaI[] = {na,nI};
    int sizeai[] = {na,ni};
    int sizeIi[] = {nI,ni};
    int sizeIn[] = {nI,N};
    int sizein[] = {ni,N};
    int sizeAaIi[] = {nA,na,nI,ni};

    const DistTensor<U>& CA = uhf.getCA();
    const DistTensor<U>& Ca = uhf.getCa();
    const DistTensor<U>& CI = uhf.getCI();
    const DistTensor<U>& Ci = uhf.getCi();
    const DistTensor<U>& S = uhf.getOverlap();
    DistTensor<U> DAI(ctf, 2, sizeAI, shapeNN, false);
    DistTensor<U> DAi(ctf, 2, sizeAi, shapeNN, false);
    DistTensor<U> DaI(ctf, 2, sizeaI, shapeNN, false);
    DistTensor<U> Dai(ctf, 2, sizeai, shapeNN, false);
    DistTensor<U> DIj(ctf, 2, sizeIi, shapeNN, false);
    DistTensor<U> DAbIj(ctf, 4, sizeAaIi, shapeNNNN, false);
    DistTensor<U> tmp1(ctf, 2, sizeIn, shapeNN, false);
    DistTensor<U> tmp2(ctf, 2, sizein, shapeNN, false);

    tmp1["Iq"] = CI["pI"]*S["pq"];
    DIj["Ij"] = tmp1["Iq"]*Ci["qj"];
    DaI["aI"] = tmp1["Iq"]*Ca["qa"];

    tmp2["iq"] = Ci["pi"]*S["pq"];
    DAi["Ai"] = tmp2["iq"]*CA["qA"];

    DAI["AI"] = DAi["Aj"]*DIj["Ij"];
    Dai["ai"] = DaI["aJ"]*DIj["Ji"];
    DAbIj["AbIj"] = DAi["Aj"]*DaI["bI"];

    const DistTensor<U>& T1A = T1(1,0,0,1);
    const DistTensor<U>& T1B = T1(0,0,0,0);
    DistTensor<U> TauAB(T2(1,0,0,1));

    TauAB["AbIj"] += T1A["AI"]*T1B["bj"];

    U S2 = uhf.getS2();

    U S2T11 = -scalar(DAI*T1A);
    U S2T12 = -scalar(Dai*T1B);
    U S2T2 = -scalar(DAbIj*TauAB);

    return abs(S2+S2T11+S2T12+S2T2);
}

template <typename U>
double CCSD<U>::getProjectedS2() const
{
    return getProjectedS2(this->uhf, T(1), T(2));
}

template <typename U>
double CCSD<U>::getProjectedMultiplicity(const UHF<U>& uhf,
                                         const SpinorbitalTensor<U>& T1,
                                         const SpinorbitalTensor<U>& T2)
{
    return sqrt(1+4*getProjectedS2(uhf, T1, T2));
}

template <typename U>
double CCSD<U>::getProjectedMultiplicity() const
{
    return getProjectedMultiplicity(this->uhf, T(1), T(2));
}

INSTANTIATE_SPECIALIZATIONS(CCSD);
