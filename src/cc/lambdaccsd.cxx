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
#include <iostream>

using namespace std;
using namespace aquarius::tensor;
using namespace aquarius::time;
using namespace aquarius::input;
using namespace aquarius::scf;

namespace aquarius
{
namespace cc
{

LambdaCCSD::LambdaCCSD(const Config& config, CCSD& ccsd)
: Distributed<double>(ccsd.ctf), Iterative(config), moints(ccsd.moints), ccsd(ccsd),
  L1("i,a"), E1("i,a"), D1("i,a"), Z1("i,a"),
  L2("ij,ab"), E2("ij,ab"), D2("ij,ab"), Z2("ij,ab"),
  H(moints, Hamiltonian::FAE|
            Hamiltonian::FMI|
            Hamiltonian::FME|
            Hamiltonian::WMNIJ|
            Hamiltonian::WMNIE|
            Hamiltonian::WMBIJ|
            Hamiltonian::WMBEJ),
  diis(config.get("diis"), 2, 2)
{
    int N = moints.getSCF().getMolecule().getNumOrbitals();
    int nI = moints.getSCF().getMolecule().getNumAlphaElectrons();
    int ni = moints.getSCF().getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;

    int sizeIA[] = {nI, nA};
    int sizeia[] = {ni, na};
    int sizeIIAA[] = {nI, nI, nA, nA};
    int sizeIiAa[] = {nI, ni, nA, na};
    int sizeiiaa[] = {ni, ni, na, na};

    int shapeNN[] = {NS, NS};
    int shapeNNNN[] = {NS, NS, NS, NS};
    int shapeANAN[] = {AS, NS, AS, NS};

    L1.addSpinCase(new DistTensor<double>(ctf, 2, sizeIA, shapeNN, false), "I,A", "IA");
    L1.addSpinCase(new DistTensor<double>(ctf, 2, sizeia, shapeNN, false), "i,a", "ia");

    E1.addSpinCase(new DistTensor<double>(ctf, 2, sizeIA, shapeNN, false), "I,A", "IA");
    E1.addSpinCase(new DistTensor<double>(ctf, 2, sizeia, shapeNN, false), "i,a", "ia");

    D1.addSpinCase(new DistTensor<double>(ctf, 2, sizeIA, shapeNN, false), "I,A", "IA");
    D1.addSpinCase(new DistTensor<double>(ctf, 2, sizeia, shapeNN, false), "i,a", "ia");

    Z1.addSpinCase(new DistTensor<double>(ctf, 2, sizeIA, shapeNN, false), "I,A", "IA");
    Z1.addSpinCase(new DistTensor<double>(ctf, 2, sizeia, shapeNN, false), "i,a", "ia");

    L2.addSpinCase(new DistTensor<double>(ctf, 4, sizeIIAA, shapeANAN, false), "IJ,AB", "IJAB");
    L2.addSpinCase(new DistTensor<double>(ctf, 4, sizeIiAa, shapeNNNN, false), "Ij,Ab", "IjAb");
    L2.addSpinCase(new DistTensor<double>(ctf, 4, sizeiiaa, shapeANAN, false), "ij,ab", "ijab");

    E2.addSpinCase(new DistTensor<double>(ctf, 4, sizeIIAA, shapeANAN, false), "IJ,AB", "IJAB");
    E2.addSpinCase(new DistTensor<double>(ctf, 4, sizeIiAa, shapeNNNN, false), "Ij,Ab", "IjAb");
    E2.addSpinCase(new DistTensor<double>(ctf, 4, sizeiiaa, shapeANAN, false), "ij,ab", "ijab");

    D2.addSpinCase(new DistTensor<double>(ctf, 4, sizeIIAA, shapeANAN, false), "IJ,AB", "IJAB");
    D2.addSpinCase(new DistTensor<double>(ctf, 4, sizeIiAa, shapeNNNN, false), "Ij,Ab", "IjAb");
    D2.addSpinCase(new DistTensor<double>(ctf, 4, sizeiiaa, shapeANAN, false), "ij,ab", "ijab");

    Z2.addSpinCase(new DistTensor<double>(ctf, 4, sizeIIAA, shapeANAN, false), "IJ,AB", "IJAB");
    Z2.addSpinCase(new DistTensor<double>(ctf, 4, sizeIiAa, shapeNNNN, false), "Ij,Ab", "IjAb");
    Z2.addSpinCase(new DistTensor<double>(ctf, 4, sizeiiaa, shapeANAN, false), "ij,ab", "ijab");

    D1["ia"]  = moints.getFIJ()["ii"];
    D1["ia"] -= moints.getFAB()["aa"];

    D2["ijab"]  = moints.getFIJ()["ii"];
    D2["ijab"] += moints.getFIJ()["jj"];
    D2["ijab"] -= moints.getFAB()["aa"];
    D2["ijab"] -= moints.getFAB()["bb"];

    int64_t size;
    double * data;
    data = D1.getSpinCase(0).getRawData(size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D1.getSpinCase(1).getRawData(size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D2.getSpinCase(0).getRawData(size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D2.getSpinCase(1).getRawData(size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }
    data = D2.getSpinCase(2).getRawData(size);
    for (int i=0; i<size; i++){
      if (fabs(data[i]) > DBL_MIN)
        data[i] = 1./data[i];
    }

    SpinorbitalTensor< DistTensor<double> >& T1 = ccsd.T1;
    SpinorbitalTensor< DistTensor<double> >& T2 = ccsd.T2;
    SpinorbitalTensor< DistTensor<double> >& FME = H.getFME();
    SpinorbitalTensor< DistTensor<double> >& FAE = H.getFAE();
    SpinorbitalTensor< DistTensor<double> >& FMI = H.getFMI();
    SpinorbitalTensor< DistTensor<double> >& WMNEF = H.getWMNEF();
    SpinorbitalTensor< DistTensor<double> >& WAMEF = H.getWAMEF();
    SpinorbitalTensor< DistTensor<double> >& WABEJ = H.getWABEJ();
    SpinorbitalTensor< DistTensor<double> >& WABEF = H.getWABEF();
    SpinorbitalTensor< DistTensor<double> >& WMNIJ = H.getWMNIJ();
    SpinorbitalTensor< DistTensor<double> >& WMNIE = H.getWMNIE();
    SpinorbitalTensor< DistTensor<double> >& WMBIJ = H.getWMBIJ();
    SpinorbitalTensor< DistTensor<double> >& WMBEJ = H.getWMBEJ();

    SpinorbitalTensor< DistTensor<double> > Tau(T2);
    Tau["abij"] += 0.5*T1["ai"]*T1["bj"];

    PROFILE_SECTION(calc_FEM)
    FME["me"] = WMNEF["mnef"]*T1["fn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_FMI)
    FMI["mi"] = 0.5*WMNEF["nmef"]*T2["efni"];
    FMI["mi"] += FME["me"]*T1["ei"];
    FMI["mi"] += WMNIE["mnif"]*T1["fn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_FAE)
    FAE["ae"] = -0.5*WMNEF["mnfe"]*T2["famn"];
    FAE["ae"] -= FME["me"]*T1["am"];
    FAE["ae"] += WABEJ["efan"]*T1["fn"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WMNIE)
    WMNIE["mnie"] += WMNEF["mnfe"]*T1["fi"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WIJKL)
    WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*Tau["efij"];
    WMNIJ["mnij"] += WMNIE["mnie"]*T1["ej"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WAIJK)
    WMBIJ["mbij"] += 0.5*WABEJ["febm"]*Tau["efij"];
    WMBIJ["mbij"] += WMBEJ["mbej"]*T1["ei"];
    WMBIJ["mbij"] += WMNIE["mnie"]*T2["bejn"];
    WMBIJ["mbij"] -= WMNIJ["mnij"]*T1["bn"];
    WMBIJ["mbij"] += FME["me"]*T2["ebij"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WMBEJ)
    WMBEJ["maei"] += WMNEF["mnef"]*T2["afin"];
    WMBEJ["maei"] += WABEJ["feam"]*T1["fi"];
    WMBEJ["maei"] -= 0.5*WMNIE["nmie"]*T1["an"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WABCI)
    WABEJ["abej"] += 0.5*WMNIE["nmje"]*T2["abmn"];
    WABEJ["abej"] -= WMBEJ["mbej"]*T1["am"];
    WABEJ["abej"] += WAMEF["amef"]*T2["fbmj"];
    WABEJ["abej"] += WABEF["abef"]*T1["fj"];
    WABEJ["abej"] -= FME["me"]*T2["abmj"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WMBEJ3)
    WMBEJ["maei"] -= 0.5*WMNIE["nmie"]*T1["an"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WABCD)
    WABEF["abef"] += 0.5*WMNEF["mnef"]*Tau["abmn"];
    WABEF["abef"] -= WAMEF["amef"]*T1["bm"];
    PROFILE_STOP

    PROFILE_SECTION(calc_WAIBC)
    WAMEF["amef"] -= WMNEF["nmef"]*T1["an"];
    PROFILE_STOP

    L1["ia"] = -FME["ia"]*D1["ia"];
    L2["ijab"] = -WMNEF["ijab"]*D2["ijab"];
}

void LambdaCCSD::_iterate()
{
    Hamiltonian D(moints, Hamiltonian::FAE|
                          Hamiltonian::FMI|
                          Hamiltonian::FME|
                          Hamiltonian::WMNIJ|
                          Hamiltonian::WMNIE|
                          Hamiltonian::WMBIJ|
                          Hamiltonian::WMBEJ);

    SpinorbitalTensor< DistTensor<double> >& T1 = ccsd.T1;
    SpinorbitalTensor< DistTensor<double> >& T2 = ccsd.T2;
    SpinorbitalTensor< DistTensor<double> >& FME = H.getFME();
    SpinorbitalTensor< DistTensor<double> >& FAE = H.getFAE();
    SpinorbitalTensor< DistTensor<double> >& FMI = H.getFMI();
    SpinorbitalTensor< DistTensor<double> >& WMNEF = H.getWMNEF();
    SpinorbitalTensor< DistTensor<double> >& WAMEF = H.getWAMEF();
    SpinorbitalTensor< DistTensor<double> >& WABEJ = H.getWABEJ();
    SpinorbitalTensor< DistTensor<double> >& WABEF = H.getWABEF();
    SpinorbitalTensor< DistTensor<double> >& WMNIJ = H.getWMNIJ();
    SpinorbitalTensor< DistTensor<double> >& WMNIE = H.getWMNIE();
    SpinorbitalTensor< DistTensor<double> >& WMBIJ = H.getWMBIJ();
    SpinorbitalTensor< DistTensor<double> >& WMBEJ = H.getWMBEJ();
    SpinorbitalTensor< DistTensor<double> >& DMN = D.getFMI();
    SpinorbitalTensor< DistTensor<double> >& DEF = D.getFAE();

    SpinorbitalTensor< DistTensor<double> > Tau(T2);
    Tau["abij"] += 0.5*T1["ai"]*T1["bj"];

    /**************************************************************************
     *
     * "Density-like" intermediates
     */
    PROFILE_SECTION(calc_DMN)
    DMN["mn"] = 0.5*T2["efno"]*L2["moef"];
    PROFILE_STOP

    PROFILE_SECTION(calc_DEF)
    DEF["ef"] = -0.5*L2["mnfg"]*T2["egmn"];
    PROFILE_STOP
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * L1->L1 and L2->L1
     */
    PROFILE_SECTION(calc_MFEM)
    Z1["ia"] = -FME["ia"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L1_IN_L1_RING)
    Z1["ia"] += L1["me"]*WMBEJ["ieam"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L1_ABCI)
    Z1["ia"] += 0.5*WABEJ["efam"]*L2["imef"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L1_IJKA)
    Z1["ia"] -= 0.5*WMBIJ["iemn"]*L2["mnae"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L1_IN_L1_FAE)
    Z1["ia"] += L1["ie"]*FAE["ea"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L1_IN_L1_FMI)
    Z1["ia"] -= L1["ma"]*FMI["im"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L1_DMN)
    Z1["ia"] -= DMN["mn"]*WMNIE["nima"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L1_DEF)
    Z1["ia"] -= DEF["ef"]*WAMEF["fiea"];
    PROFILE_STOP
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * L1->L2 and L2->L2
     */
    PROFILE_SECTION(calc_WMNEF)
    Z2["ijab"] = -WMNEF["ijab"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L2_FAE)
    Z2["ijab"] += FAE["fa"]*L2["ijfb"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L2_FMI)
    Z2["ijab"] -= FMI["in"]*L2["njab"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L1_IN_L2_FME)
    Z2["ijab"] += FME["ia"]*L1["jb"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L1_IN_L2_ABCI)
    Z2["ijab"] += WAMEF["ejab"]*L1["ie"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L1_IN_L2_IJKA)
    Z2["ijab"] -= WMNIE["ijmb"]*L1["ma"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L2_ABCD)
    Z2["ijab"] += 0.5*WABEF["efab"]*L2["ijef"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L2_IJKL)
    Z2["ijab"] += 0.5*WMNIJ["ijmn"]*L2["mnab"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L2_RING)
    Z2["ijab"] += WMBEJ["ieam"]*L2["mjeb"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L2_DMN)
    Z2["ijab"] -= DMN["im"]*WMNEF["mjab"];
    PROFILE_STOP

    PROFILE_SECTION(calc_L2_IN_L2_DEF)
    Z2["ijab"] += DEF["ae"]*WMNEF["ijeb"];
    PROFILE_STOP
    /*
     *************************************************************************/

    PROFILE_SECTION(calc_EN)

    E1["ia"]     = Z1["ia"]*D1["ia"];
    E2["ijab"]   = Z2["ijab"]*D2["ijab"];

    SpinorbitalTensor< DistTensor<double> > L1H(L1);
    SpinorbitalTensor< DistTensor<double> > L2H(L2);

    L1H["ia"] -= L1["ia"]*D1["ia"];
    L1H["ia"] += L1["ia"]*ccsd.getEnergy();
    L2H["ijab"] -= L2["ijab"]*D2["ijab"];
    L2H["ijab"] += L2["ijab"]*ccsd.getEnergy();

    double LL = 1 + scalar(L1*L1) + 0.25*scalar(L2*L2);
    double LHL = ccsd.getEnergy() + scalar(L1H*L1) + 0.25*scalar(L2H*L2);

    energy = LHL/LL;

    E1["ia"]    -= L1["ia"];
    L1["ia"]    += E1["ia"];
    E2["ijab"]  -= L2["ijab"];
    L2["ijab"]  += E2["ijab"];

    //energy = ccsd.energy;

    conv =          E1.getSpinCase(0).reduce(CTF_OP_MAXABS);
    conv = max(conv,E1.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,E2.getSpinCase(0).reduce(CTF_OP_MAXABS));
    conv = max(conv,E2.getSpinCase(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,E2.getSpinCase(2).reduce(CTF_OP_MAXABS));

    E2 *= 0.5;

    vector<SpinorbitalTensor< DistTensor<double> >*> L(2);
    L[0] = &L1;
    L[1] = &L2;
    vector<SpinorbitalTensor< DistTensor<double> >*> E(2);
    E[0] = &E1;
    E[1] = &E2;
    //diis.extrapolate(L, E);

    PROFILE_STOP
}

}
}
