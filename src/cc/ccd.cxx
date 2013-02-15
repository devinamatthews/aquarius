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

#include "cc.hpp"

#include <cfloat>

using namespace std;
using namespace libtensor;
using namespace aquarius::autocc;
using namespace aquarius::input;
using namespace aquarius::time;
using namespace aquarius::scf;

namespace aquarius
{
namespace cc
{

CCD::CCD(const Config& config, MOIntegrals& moints)
: Iterative(config), moints(moints),
  T2("ab,ij"), E2("ab,ij"), D2("ab,ij"), Z2("ab,ij")
{
    int N = moints.getSCF().getMolecule().getNumOrbitals();
    int nI = moints.getSCF().getMolecule().getNumAlphaElectrons();
    int ni = moints.getSCF().getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;

    int sizeAAII[] = {nA, nA, nI, nI};
    int sizeAaIi[] = {nA, na, nI, ni};
    int sizeaaii[] = {na, na, ni, ni};

    int shapeNNNN[] = {NS, NS, NS, NS};
    int shapeANAN[] = {AS, NS, AS, NS};

    DistWorld *dw = moints.getFAB().getSpinCase(0).dw;

    T2.addSpinCase(new DistTensor(4, sizeAAII, shapeANAN, dw, false), "AB,IJ", "ABIJ");
    T2.addSpinCase(new DistTensor(4, sizeAaIi, shapeNNNN, dw, false), "Ab,Ij", "AbIj");
    T2.addSpinCase(new DistTensor(4, sizeaaii, shapeANAN, dw, false), "ab,ij", "abij");

    E2.addSpinCase(new DistTensor(4, sizeAAII, shapeANAN, dw, false), "AB,IJ", "ABIJ");
    E2.addSpinCase(new DistTensor(4, sizeAaIi, shapeNNNN, dw, false), "Ab,Ij", "AbIj");
    E2.addSpinCase(new DistTensor(4, sizeaaii, shapeANAN, dw, false), "ab,ij", "abij");

    D2.addSpinCase(new DistTensor(4, sizeAAII, shapeANAN, dw, false), "AB,IJ", "ABIJ");
    D2.addSpinCase(new DistTensor(4, sizeAaIi, shapeNNNN, dw, false), "Ab,Ij", "AbIj");
    D2.addSpinCase(new DistTensor(4, sizeaaii, shapeANAN, dw, false), "ab,ij", "abij");

    Z2.addSpinCase(new DistTensor(4, sizeAAII, shapeANAN, dw, false), "AB,IJ", "ABIJ");
    Z2.addSpinCase(new DistTensor(4, sizeAaIi, shapeNNNN, dw, false), "Ab,Ij", "AbIj");
    Z2.addSpinCase(new DistTensor(4, sizeaaii, shapeANAN, dw, false), "ab,ij", "abij");

    D2["abij"]  = moints.getFIJ()["ii"];
    D2["abij"] += moints.getFIJ()["jj"];
    D2["abij"] -= moints.getFAB()["aa"];
    D2["abij"] -= moints.getFAB()["bb"];

    int size;
    double * data;
    data = D2.getSpinCase(0).getRawData(&size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D2.getSpinCase(1).getRawData(&size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D2.getSpinCase(2).getRawData(&size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }

    T2["abij"] = moints.getVABIJ()["abij"]*D2["abij"];

    energy = 0.25*scalar(moints.getVABIJ()["efmn"]*T2["efmn"]);

    conv =          T2.getSpinCase(0).reduce(CTF_OP_MAXABS);
    conv = max(conv,T2.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,T2.getSpinCase(2).reduce(CTF_OP_MAXABS));
}

void CCD::_iterate()
{
    Hamiltonian H(moints, Hamiltonian::FAE|
                          Hamiltonian::FMI|
                          Hamiltonian::WMNIJ|
                          Hamiltonian::WMBEJ);

    SpinorbitalTensor<DistTensor>& FAE = H.getFAE();
    SpinorbitalTensor<DistTensor>& FMI = H.getFMI();
    SpinorbitalTensor<DistTensor>& WMNEF = H.getWMNEF();
    SpinorbitalTensor<DistTensor>& WABEF = H.getWABEF();
    SpinorbitalTensor<DistTensor>& WMNIJ = H.getWMNIJ();
    SpinorbitalTensor<DistTensor>& WMBEJ = H.getWMBEJ();

    //FAE["aa"] = 0.0;
    //FMI["ii"] = 0.0;

    FAE = 0.0;
    FMI = 0.0;

    PROFILE_SECTION(calc_WMNEF)
    Z2["abij"] = WMNEF["ijab"];
    PROFILE_STOP

    PROFILE_SECTION(calc_FAE)
    FAE["ae"] = -0.5*WMNEF["mnef"]*T2["afmn"];
    Z2["abij"] += FAE["af"]*T2["fbij"];
    PROFILE_STOP

    PROFILE_SECTION(calc_FMI)
    FMI["mi"] = 0.5*WMNEF["mnef"]*T2["efin"];
    Z2["abij"] -= FMI["ni"]*T2["abnj"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WABEF)
    Z2["abij"] += 0.5*WABEF["abef"]*T2["efij"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WMNIJ)
    WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*T2["efij"];
    Z2["abij"] += 0.5*WMNIJ["mnij"]*T2["abmn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WMBEJ)
    WMBEJ["maei"] += 0.5*WMNEF["mnef"]*T2["afin"];
    Z2["abij"] += WMBEJ["maei"]*T2["ebmj"];
    PROFILE_STOP

    PROFILE_SECTION(calc_EN)
    E2["abij"]  = Z2["abij"]*D2["abij"];
    E2["abij"] -= T2["abij"];
    T2["abij"] += E2["abij"];
    energy      = 0.25*scalar(WMNEF["mnef"]*T2["efmn"]);
    PROFILE_STOP

    conv =           E2.getSpinCase(0).reduce(CTF_OP_MAXABS);
    conv = max(conv, E2.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv, E2.getSpinCase(2).reduce(CTF_OP_MAXABS));
}

}
}
