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

#include "choleskymoints.hpp"

using namespace aquarius::tensor;

namespace aquarius
{
namespace scf
{

CholeskyMOIntegrals::CholeskyMOIntegrals(const CholeskyUHF& uhf)
: MOIntegrals(uhf)
{
    doTransformation(uhf.chol);
}

void CholeskyMOIntegrals::doTransformation(const CholeskyIntegrals& chol)
{
    const DistTensor<double>& cA = uhf.getCA();
    const DistTensor<double>& ca = uhf.getCa();
    const DistTensor<double>& cI = uhf.getCI();
    const DistTensor<double>& ci = uhf.getCi();
    const DistTensor<double>& Lpq = chol.getL();
    const DistTensor<double>& D = chol.getD();

    int N = uhf.getMolecule().getNumOrbitals();
    int nI = uhf.getMolecule().getNumAlphaElectrons();
    int ni = uhf.getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;
    int R = chol.getRank();

    int sizeIIR[] = {nI, nI, R};
    int sizeiiR[] = {ni, ni, R};
    int sizeAAR[] = {nA, nA, R};
    int sizeaaR[] = {na, na, R};
    int sizeAIR[] = {nA, nI, R};
    int sizeaiR[] = {na, ni, R};

    int shapeNNN[] = {NS, NS, NS};

    DistTensor<double> LIJ(ctf, 3, sizeIIR, shapeNNN, false);
    DistTensor<double> Lij(ctf, 3, sizeiiR, shapeNNN, false);
    DistTensor<double> LAB(ctf, 3, sizeAAR, shapeNNN, false);
    DistTensor<double> Lab(ctf, 3, sizeaaR, shapeNNN, false);
    DistTensor<double> LAI(ctf, 3, sizeAIR, shapeNNN, false);
    DistTensor<double> Lai(ctf, 3, sizeaiR, shapeNNN, false);

    {
        int sizeNIR[] = { N, nI, R};
        int sizeNiR[] = { N, ni, R};
        int sizeNAR[] = { N, nA, R};
        int sizeNaR[] = { N, na, R};

        DistTensor<double> LpI(ctf, 3, sizeNIR, shapeNNN, false);
        DistTensor<double> Lpi(ctf, 3, sizeNiR, shapeNNN, false);
        DistTensor<double> LpA(ctf, 3, sizeNAR, shapeNNN, false);
        DistTensor<double> Lpa(ctf, 3, sizeNaR, shapeNNN, false);

        LpI["pIR"] = Lpq["pqR"]*cI["qI"];
        Lpi["piR"] = Lpq["pqR"]*ci["qi"];
        LpA["pAR"] = Lpq["pqR"]*cA["qA"];
        Lpa["paR"] = Lpq["pqR"]*ca["qa"];

        LIJ["IJR"] = LpI["pJR"]*cI["pI"];
        Lij["ijR"] = Lpi["pjR"]*ci["pi"];
        LAI["AIR"] = LpI["pIR"]*cA["pA"];
        Lai["aiR"] = Lpi["piR"]*ca["pa"];
        LAB["ABR"] = LpA["pBR"]*cA["pA"];
        Lab["abR"] = Lpa["pbR"]*ca["pa"];
    }

    DistTensor<double> LDIJ(ctf, 3, sizeIIR, shapeNNN, false);
    DistTensor<double> LDij(ctf, 3, sizeiiR, shapeNNN, false);
    DistTensor<double> LDAB(ctf, 3, sizeAAR, shapeNNN, false);
    DistTensor<double> LDab(ctf, 3, sizeaaR, shapeNNN, false);
    DistTensor<double> LDAI(ctf, 3, sizeAIR, shapeNNN, false);
    DistTensor<double> LDai(ctf, 3, sizeaiR, shapeNNN, false);

    LDIJ["IJR"] = D["R"]*LIJ["IJR"];
    LDij["ijR"] = D["R"]*Lij["ijR"];
    LDAI["AIR"] = D["R"]*LAI["AIR"];
    LDai["aiR"] = D["R"]*Lai["aiR"];
    LDAB["ABR"] = D["R"]*LAB["ABR"];
    LDab["abR"] = D["R"]*Lab["abR"];

    (*IJKL_)["IJKL"] = 0.5*LDIJ["IKR"]*LIJ["JLR"];
    (*IjKl_)["IjKl"] =     LDIJ["IKR"]*Lij["jlR"];
    (*ijkl_)["ijkl"] = 0.5*LDij["ikR"]*Lij["jlR"];

    (*IJKA_)["IJKA"] = LDIJ["IKR"]*LAI["AJR"];
    (*IjKa_)["IjKa"] = LDIJ["IKR"]*Lai["ajR"];
    (*iJkA_)["iJkA"] = LDij["ikR"]*LAI["AJR"];
    (*ijka_)["ijka"] = LDij["ikR"]*Lai["ajR"];

    (*ABIJ_)["ABIJ"] = 0.5*LDAI["AIR"]*LAI["BJR"];
    (*AbIj_)["AbIj"] =     LDAI["AIR"]*Lai["bjR"];
    (*abij_)["abij"] = 0.5*LDai["aiR"]*Lai["bjR"];

    (*AIBJ_)["AIBJ"]  = LDAB["ABR"]*LIJ["IJR"];
    (*AIBJ_)["AIBJ"] -= LDAI["AJR"]*LAI["BIR"];
    (*AiBj_)["AiBj"]  = LDAB["ABR"]*Lij["ijR"];
    (*aIbJ_)["aIbJ"]  = LDab["abR"]*LIJ["IJR"];
    (*aibj_)["aibj"]  = LDab["abR"]*Lij["ijR"];
    (*aibj_)["aibj"] -= LDai["ajR"]*Lai["biR"];

    (*ABCI_)["ABCI"] = LDAB["ACR"]*LAI["BIR"];
    (*AbCi_)["AbCi"] = LDAB["ACR"]*Lai["biR"];
    (*aBcI_)["aBcI"] = LDab["acR"]*LAI["BIR"];
    (*abci_)["abci"] = LDab["acR"]*Lai["biR"];

    (*ABCD_)["ABCD"] = 0.5*LDAB["ACR"]*LAB["BDR"];
    (*AbCd_)["AbCd"] =     LDAB["ACR"]*Lab["bdR"];
    (*abcd_)["abcd"] = 0.5*LDab["acR"]*Lab["bdR"];
}

}
}
