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
using namespace aquarius::time;
using namespace aquarius::scf;
using namespace aquarius::input;

namespace aquarius
{
namespace cc
{

CCSDT::CCSDT(const Config& config, MOIntegrals& moints)
: Iterative(config), moints(moints),
  T1("a,i"), E1("a,i"), D1("a,i"), Z1("a,i"),
  T2("ab,ij"), E2("ab,ij"), D2("ab,ij"), Z2("ab,ij"),
  T3("abc,ijk"), E3("abc,ijk"), D3("abc,ijk"), Z3("abc,ijk")
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
    int sizeAAAIII[] = {nA, nA, nA, nI, nI, nI};
    int sizeAAaIIi[] = {nA, nA, na, nI, nI, ni};
    int sizeAaaIii[] = {nA, na, na, nI, ni, ni};
    int sizeaaaiii[] = {na, na, na, ni, ni, ni};

    int shapeNN[] = {NS, NS};
    int shapeNNNN[] = {NS, NS, NS, NS};
    int shapeANAN[] = {AS, NS, AS, NS};
    int shapeANNANN[] = {AS, NS, NS, AS, NS, NS};
    int shapeNANNAN[] = {NS, AS, NS, NS, AS, NS};
    int shapeAANAAN[] = {AS, AS, NS, AS, AS, NS};

    DistWorld *dw = moints.getFAB().getSpinCase(0).dw;

    T1.addSpinCase(new DistTensor(2, sizeAI, shapeNN, dw, false), "A,I", "AI");
    T1.addSpinCase(new DistTensor(2, sizeai, shapeNN, dw, false), "a,i", "ai");

    E1.addSpinCase(new DistTensor(2, sizeAI, shapeNN, dw, false), "A,I", "AI");
    E1.addSpinCase(new DistTensor(2, sizeai, shapeNN, dw, false), "a,i", "ai");

    D1.addSpinCase(new DistTensor(2, sizeAI, shapeNN, dw, false), "A,I", "AI");
    D1.addSpinCase(new DistTensor(2, sizeai, shapeNN, dw, false), "a,i", "ai");

    Z1.addSpinCase(new DistTensor(2, sizeAI, shapeNN, dw, false), "A,I", "AI");
    Z1.addSpinCase(new DistTensor(2, sizeai, shapeNN, dw, false), "a,i", "ai");

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

    T3.addSpinCase(new DistTensor(6, sizeAAAIII, shapeAANAAN, dw, true), "ABC,IJK", "ABCIJK");
    T3.addSpinCase(new DistTensor(6, sizeAAaIIi, shapeANNANN, dw, true), "ABc,IJk", "ABcIJk");
    T3.addSpinCase(new DistTensor(6, sizeAaaIii, shapeNANNAN, dw, true), "Abc,Ijk", "AbcIjk");
    T3.addSpinCase(new DistTensor(6, sizeaaaiii, shapeAANAAN, dw, true), "abc,ijk", "abcijk");

    E3.addSpinCase(new DistTensor(6, sizeAAAIII, shapeAANAAN, dw, false), "ABC,IJK", "ABCIJK");
    E3.addSpinCase(new DistTensor(6, sizeAAaIIi, shapeANNANN, dw, false), "ABc,IJk", "ABcIJk");
    E3.addSpinCase(new DistTensor(6, sizeAaaIii, shapeNANNAN, dw, false), "Abc,Ijk", "AbcIjk");
    E3.addSpinCase(new DistTensor(6, sizeaaaiii, shapeAANAAN, dw, false), "abc,ijk", "abcijk");

    D3.addSpinCase(new DistTensor(6, sizeAAAIII, shapeAANAAN, dw, false), "ABC,IJK", "ABCIJK");
    D3.addSpinCase(new DistTensor(6, sizeAAaIIi, shapeANNANN, dw, false), "ABc,IJk", "ABcIJk");
    D3.addSpinCase(new DistTensor(6, sizeAaaIii, shapeNANNAN, dw, false), "Abc,Ijk", "AbcIjk");
    D3.addSpinCase(new DistTensor(6, sizeaaaiii, shapeAANAAN, dw, false), "abc,ijk", "abcijk");

    Z3.addSpinCase(new DistTensor(6, sizeAAAIII, shapeAANAAN, dw, false), "ABC,IJK", "ABCIJK");
    Z3.addSpinCase(new DistTensor(6, sizeAAaIIi, shapeANNANN, dw, false), "ABc,IJk", "ABcIJk");
    Z3.addSpinCase(new DistTensor(6, sizeAaaIii, shapeNANNAN, dw, false), "Abc,Ijk", "AbcIjk");
    Z3.addSpinCase(new DistTensor(6, sizeaaaiii, shapeAANAAN, dw, false), "abc,ijk", "abcijk");

    D1["ai"]  = moints.getFIJ()["ii"];
    D1["ai"] -= moints.getFAB()["aa"];

    D2["abij"]  = moints.getFIJ()["ii"];
    D2["abij"] += moints.getFIJ()["jj"];
    D2["abij"] -= moints.getFAB()["aa"];
    D2["abij"] -= moints.getFAB()["bb"];

    D3["abcijk"]  = moints.getFIJ()["ii"];
    D3["abcijk"] += moints.getFIJ()["jj"];
    D3["abcijk"] += moints.getFIJ()["kk"];
    D3["abcijk"] -= moints.getFAB()["aa"];
    D3["abcijk"] -= moints.getFAB()["bb"];
    D3["abcijk"] -= moints.getFAB()["cc"];

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
    data = D3.getSpinCase(0).getRawData(&size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D3.getSpinCase(1).getRawData(&size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D3.getSpinCase(2).getRawData(&size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D3.getSpinCase(3).getRawData(&size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }

    T1["ai"] = moints.getFAI()["ai"]*D1["ai"];
    T2["abij"] = moints.getVABIJ()["abij"]*D2["abij"];

    SpinorbitalTensor<DistTensor> Tau(T2);
    Tau["abij"] += 0.5*T1["ai"]*T1["bj"];

    energy = 0.25*scalar(moints.getVABIJ()["efmn"]*Tau["efmn"]);

    conv =          T1.getSpinCase(0).reduce(CTF_OP_MAXABS);
    conv = max(conv,T1.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,T2.getSpinCase(0).reduce(CTF_OP_MAXABS));
    conv = max(conv,T2.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,T2.getSpinCase(2).reduce(CTF_OP_MAXABS));
}

void CCSDT::_iterate()
{
    Hamiltonian H(moints, Hamiltonian::FAE|
                          Hamiltonian::FMI|
                          Hamiltonian::FME|
                          Hamiltonian::WAMEF|
                          Hamiltonian::WABEJ|
                          Hamiltonian::WABEF|
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
    FMI["mi"] = 0.5*WMNEF["nmef"]*T2["efni"];
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
    Z1["ai"] -= T1["em"]*WMBEJ["amei"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T1_ABCI)
    Z1["ai"] += 0.5*WAMEF["amef"]*Tau["efim"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T1_IJKA)
    Z1["ai"] -= 0.5*WMNIE["mnie"]*T2["aemn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T1_FME)
    Z1["ai"] += T2["aeim"]*FME["me"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T1_IN_T1_FAE)
    Z1["ai"] += T1["ei"]*FAE["ae"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T1_IN_T1_FMI)
    Z1["ai"] -= T1["am"]*FMI["mi"];
    PROFILE_STOP
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * Intermediates for T1->T2 and T2->T2
     */
    PROFILE_SECTION(calc_FAE)
    FAE["ae"] = -0.5*WMNEF["mnfe"]*T2["famn"];
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

    /**************************************************************************
     *
     * Intermediates for CCSDT
     */
    PROFILE_SECTION(calc_WAIJK2)
    WMBIJ["mbij"] += WMNIE["mnie"]*T2["bejn"];
    WMBIJ["mbij"] -= WMNIJ["mnij"]*T1["bn"];
    WMBIJ["mbij"] += FME["me"]*T2["ebij"];
    //WMBIJ["mbij"] += 0.5*WMNEF["mnef"]*T3["efbinj"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WMBEJ2)
    WMBEJ["maei"] *= 2;
    WMBEJ["maei"] += moints.getVAIBJ()["amei"];
    WMBEJ["maei"] -= WAMEF["amfe"]*T1["fi"];
    WMBEJ["maei"] += 1.5*WMNIE["nmie"]*T1["an"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WABCI)
    WABEJ["abej"] += WAMEF["amef"]*T2["fbmj"];
    WABEJ["abej"] += 0.5*WMNIE["nmje"]*T2["abmn"];
    WABEJ["abej"] += WABEF["abef"]*T1["fj"];
    WABEJ["abej"] -= WMBEJ["mbej"]*T1["am"];
    //WABEJ["abej"] -= 0.5*WMNEF["mnef"]*T3["afbmnj"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WMBEJ3)
    WMBEJ["maei"] -= 0.5*WMNIE["nmie"]*T1["an"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WABCD)
    WABEF["abef"] -= WAMEF["amef"]*T1["bm"];
    WABEF["abef"] += 0.5*WMNEF["mnef"]*Tau["abmn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WAIBC)
    WAMEF["amef"] -= WMNEF["nmef"]*T1["an"];
    PROFILE_STOP
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * CCSDT Iteration
     */
    //PROFILE_SECTION(calc_T2_IN_T3_ABCI)
    //Z3["abcijk"] = WABEJ["bcek"]*T2["aeij"];
    //PROFILE_STOP

    PROFILE_SECTION(calc_T2_IN_T3_IJKA)
    //Z3["abcijk"] -= WAMIJ["cmkj"]*T2["abim"];
    Z3["abcijk"] = moints.getVIJKA()["jkmc"]*T2["abim"];
    PROFILE_STOP

    //PROFILE_SECTION(calc_T3_IN_T2_ABCI)
    //Z2["abij"] += 0.5*WAMEF["bmef"]*T3["aefijm"];
    //PROFILE_STOP

    PROFILE_SECTION(calc_T3_IN_T2_IJKA)
    //Z2["abij"] -= 0.5*WMNIE["mnje"]*T3["abeimn"];
    Z2["abij"] -= 0.5*moints.getVIJKA()["mnje"]*T3["abeimn"];
    PROFILE_STOP

    /*
    PROFILE_SECTION(calc_T3_IN_T2_FME)
    Z2["abij"] += FME["me"]*T3["abeijm"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T3_IN_T1)
    Z1["ai"] += 0.25*WMNEF["mnef"]*T3["aefimn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T3_IN_T3_FAE)
    Z3["abcijk"] += FAE["ce"]*T3["abeijk"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T3_IN_T3_FMI)
    Z3["abcijk"] -= FMI["mk"]*T3["abcijm"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T3_IN_T3_ABCD)
    Z3["abcijk"] += 0.5*WABEF["abef"]*T3["efcijk"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T3_IN_T3_IJKL)
    Z3["abcijk"] += 0.5*WMNIJ["mnij"]*T3["abcmnk"];
    PROFILE_STOP

    PROFILE_SECTION(calc_T3_IN_T3_RING)
    Z3["abcijk"] += WMBEJ["maei"]*T3["ebcmjk"];
    PROFILE_STOP
    */
    /*
     **************************************************************************/

    PROFILE_SECTION(calc_EN)
    E1["ai"]      = Z1["ai"]*D1["ai"];
    E1["ai"]     -= T1["ai"];
    T1["ai"]     += E1["ai"];
    E2["abij"]    = Z2["abij"]*D2["abij"];
    E2["abij"]   -= T2["abij"];
    T2["abij"]   += E2["abij"];
    E3["abcijk"]  = Z3["abcijk"]*D3["abcijk"];
    E3["abcijk"] -= T3["abcijk"];
    T3["abcijk"] += E3["abcijk"];
    Tau["abij"]   = T2["abij"];
    Tau["abij"]  += 0.5*T1["ai"]*T1["bj"];
    energy        = 0.25*scalar(WMNEF["mnef"]*Tau["efmn"]);
    PROFILE_STOP

    conv =           E1.getSpinCase(0).reduce(CTF_OP_MAXABS);
    conv = max(conv, E1.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv, E2.getSpinCase(0).reduce(CTF_OP_MAXABS));
    conv = max(conv, E2.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv, E2.getSpinCase(2).reduce(CTF_OP_MAXABS));
    conv = max(conv, E3.getSpinCase(0).reduce(CTF_OP_MAXABS));
    conv = max(conv, E3.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv, E3.getSpinCase(2).reduce(CTF_OP_MAXABS));
    conv = max(conv, E3.getSpinCase(3).reduce(CTF_OP_MAXABS));
}

}
}
