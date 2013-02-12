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
 * ARE DISCLAIMED. IN NO EVENT SHALL EDGAR SOLOMONIK BE LIABLE FOR ANY
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
using namespace aquarius::time;
using namespace aquarius::input;
using namespace aquarius::scf;

namespace aquarius
{
namespace cc
{

CCSD::CCSD(const Config& config, MOIntegrals& moints)
: Iterative(config), moints(moints),
  T1(moints.getFAI()), E1("a,i"), D1("a,i"), Z1("a,i"),
  T2("ac,ij"), E2("ad,ij"), D2("af,ij"), Z2("ae,ij")
{
    int N = moints.getSCF().getMolecule().getNumOrbitals();
    int nI = moints.getSCF().getMolecule().getNumAlphaElectrons();
    int ni = moints.getSCF().getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;

    int sizeAI[] = {nA, nI};
    int sizeai[] = {na, ni};
    int sizeAAII[] = {nA, nA, nI, nI};
    int sizeAaIi[] = {nA, na, nI, ni};
    int sizeaaii[] = {na, na, ni, ni};

    int shapeNN[] = {NS, NS};
    int shapeNNNN[] = {NS, NS, NS, NS};
    int shapeANAN[] = {AS, NS, AS, NS};

    DistWorld *dw = moints.getFAB().getSpinCase(0).dw;

    E1.addSpinCase(new DistTensor(2, sizeAI, shapeNN, dw, false), "A,I", "AI");
    E1.addSpinCase(new DistTensor(2, sizeai, shapeNN, dw, false), "a,i", "ai");

    D1.addSpinCase(new DistTensor(2, sizeAI, shapeNN, dw, false), "A,I", "AI");
    D1.addSpinCase(new DistTensor(2, sizeai, shapeNN, dw, false), "a,i", "ai");

    Z1.addSpinCase(new DistTensor(2, sizeAI, shapeNN, dw, false), "A,I", "AI");
    Z1.addSpinCase(new DistTensor(2, sizeai, shapeNN, dw, false), "a,i", "ai");

    T2.addSpinCase(new DistTensor(4, sizeAAII, shapeANAN, dw, false), "AC,IJ", "ACIJ");
    T2.addSpinCase(new DistTensor(4, sizeAaIi, shapeNNNN, dw, false), "Ac,Ij", "AcIj");
    T2.addSpinCase(new DistTensor(4, sizeaaii, shapeANAN, dw, false), "ac,ij", "acij");

    T2["abij"] = moints.getVABIJ()["abij"];

    E2.addSpinCase(new DistTensor(4, sizeAAII, shapeANAN, dw, false), "AD,IJ", "ADIJ");
    E2.addSpinCase(new DistTensor(4, sizeAaIi, shapeNNNN, dw, false), "Ad,Ij", "AdIj");
    E2.addSpinCase(new DistTensor(4, sizeaaii, shapeANAN, dw, false), "ad,ij", "adij");

    D2.addSpinCase(new DistTensor(4, sizeAAII, shapeANAN, dw, false), "AF,IJ", "AFIJ");
    D2.addSpinCase(new DistTensor(4, sizeAaIi, shapeNNNN, dw, false), "Af,Ij", "AfIj");
    D2.addSpinCase(new DistTensor(4, sizeaaii, shapeANAN, dw, false), "af,ij", "afij");

    Z2.addSpinCase(new DistTensor(4, sizeAAII, shapeANAN, dw, false), "AE,IJ", "AEIJ");
    Z2.addSpinCase(new DistTensor(4, sizeAaIi, shapeNNNN, dw, false), "Ae,Ij", "AeIj");
    Z2.addSpinCase(new DistTensor(4, sizeaaii, shapeANAN, dw, false), "ae,ij", "aeij");

    D1["ai"]  = moints.getFIJ()["ii"];
    D1["ai"] -= moints.getFAB()["aa"];

    D2["abij"]  = moints.getFIJ()["ii"];
    D2["abij"] += moints.getFIJ()["jj"];
    D2["abij"] -= moints.getFAB()["aa"];
    D2["abij"] -= moints.getFAB()["bb"];

    int size;
    double * data;
    data = D1.getSpinCase(0).getRawData(&size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D1.getSpinCase(1).getRawData(&size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
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
}

void CCSD::_iterate()
{
    Hamiltonian H(moints, Hamiltonian::FAE|
                          Hamiltonian::FMI|
                          Hamiltonian::FME|
                          Hamiltonian::WMNIJ|
                          Hamiltonian::WMNIE|
                          Hamiltonian::WMBIJ|
                          Hamiltonian::WMBEJ);

    SpinorbitalTensor<DistTensor>& FME = H.getFME();
    SpinorbitalTensor<DistTensor>& FAE = H.getFAE();
    SpinorbitalTensor<DistTensor>& FMI = H.getFMI();
    SpinorbitalTensor<DistTensor>& WMNEF = H.getWMNEF();
    SpinorbitalTensor<DistTensor>& WAMEF = H.getWAMEF();
    SpinorbitalTensor<DistTensor>& WABEJ = H.getWABEJ();
    SpinorbitalTensor<DistTensor>& WABEF = H.getWABEF();
    SpinorbitalTensor<DistTensor>& WMNIJ = H.getWMNIJ();
    SpinorbitalTensor<DistTensor>& WMNIE = H.getWMNIE();
    SpinorbitalTensor<DistTensor>& WMBIJ = H.getWMBIJ();
    SpinorbitalTensor<DistTensor>& WMBEJ = H.getWMBEJ();

    //FAE["aa"] = 0.0;
    //FMI["ii"] = 0.0;

    FAE = 0.0;
    FMI = 0.0;

    SpinorbitalTensor<DistTensor> Tau(T2);
    Tau["abij"] += 0.5*T1["ai"]*T1["bj"];

    /**************************************************************************
     *
     * Intermediates for T1->T1 and T2->T1
     */
    PROFILE_SECTION(calc_FEM)
    FME["me"] = WMNEF["mnef"]*T1["fn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_FMI)
    FMI["mi"] = 0.5*WMNEF["mnef"]*T2["efin"];
    FMI["mi"] += FME["me"]*T1["ei"];
    FMI["mi"] += WMNIE["mnif"]*T1["fn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WMNIE)
    WMNIE["mnie"] += WMNEF["mnfe"]*T1["fi"];
    PROFILE_STOP
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * T1->T1 and T2->T1
     */
    PROFILE_SECTION(calc_FAI)
    Z1["ai"] = moints.getFAI()["ai"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T1_IN_T1_RING)
    //Z1["ai"] -= T1["em"]*WMBEJ["amei"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T1_ABCI)
    //Z1["ai"] += 0.5*WAMEF["amef"]*Tau["efim"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T1_IJKA)
    //Z1["ai"] -= 0.5*WMNIE["mnie"]*T2["aemn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T1_FME)
    //Z1["ai"] += T2["aeim"]*FME["me"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T1_IN_T1_FAE)
    //Z1["ai"] += T1["ei"]*FAE["ae"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T1_IN_T1_FMI)
    //Z1["ai"] -= T1["am"]*FMI["mi"];
    PROFILE_STOP
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * Intermediates for T1->T2 and T2->T2
     */
    PROFILE_SECTION(calc_FAE)
    FAE["ae"] = -0.5*WMNEF["mnef"]*T2["afmn"];
    FAE["ae"] -= FME["me"]*T1["am"];
    FAE["ae"] += WABEJ["efan"]*T1["fn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WIJKL)
    WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*Tau["efij"];
    WMNIJ["mnij"] += WMNIE["mnie"]*T1["ej"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WAIJK)
    WMBIJ["mbij"] += 0.5*WABEJ["febm"]*Tau["efij"];
    WMBIJ["mbij"] += WMBEJ["amej"]*T1["ei"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WMBEJ)
    WMBEJ["maei"] += 0.5*WMNEF["mnef"]*T2["afin"];
    WMBEJ["maei"] += WABEJ["feam"]*T1["fi"];
    WMBEJ["maei"] -= WMNIE["nmie"]*T1["an"];
    PROFILE_STOP
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * T1->T2 and T2->T2
     */
    PROFILE_SECTION(calc_WMNEF)
    Z2["abij"] = WMNEF["ijab"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T2_FAE)
    Z2["abij"] += FAE["af"]*T2["fbij"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T2_FMI)
    Z2["abij"] -= FMI["ni"]*T2["abnj"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T1_IN_T2_ABCI)
    Z2["abij"] += WABEJ["abej"]*T1["ei"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T1_IN_T2_IJKA)
    Z2["abij"] -= WMBIJ["mbij"]*T1["am"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T2_ABCD)
    Z2["abij"] += 0.5*WABEF["abef"]*Tau["efij"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T2_IJKL)
    Z2["abij"] += 0.5*WMNIJ["mnij"]*Tau["abmn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T2_RING)
    Z2["abij"] += WMBEJ["maei"]*T2["ebmj"];
    PROFILE_STOP
    /*
     *************************************************************************/

    PROFILE_SECTION(calc_EN)
    E1["ai"]     = Z1["ai"]*D1["ai"];
    E1["ai"]    -= T1["ai"];
    T1["ai"]    += E1["ai"];
    E2["abij"]   = Z2["abij"]*D2["abij"];
    E2["abij"]  -= T2["abij"];
    T2["abij"]  += E2["abij"];
    Tau["abij"]  = T2["abij"];
    Tau["abij"] += 0.5*T1["ai"]*T1["bj"];
    energy = 0.25*scalar(WMNEF["mnef"]*Tau["efmn"]);
    PROFILE_STOP

    conv =          E1.getSpinCase(0).reduce(CTF_OP_MAXABS);
    conv = max(conv,E1.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,E2.getSpinCase(0).reduce(CTF_OP_MAXABS));
    conv = max(conv,E2.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,E2.getSpinCase(2).reduce(CTF_OP_MAXABS));
}

}
}
