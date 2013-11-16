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

using namespace std;
using namespace aquarius;
using namespace aquarius::scf;
using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::op;
using namespace aquarius::task;
using namespace aquarius::symmetry;

template <typename T>
CholeskyMOIntegrals<T>::CholeskyMOIntegrals(const string& name, const Config& config)
: MOIntegrals<T>("choleskymoints", name, config)
{
    this->getProduct("H").addRequirement(Requirement("cholesky","cholesky"));
}

template <typename T>
void CholeskyMOIntegrals<T>::run(TaskDAG& dag, const Arena& arena)
{
    const MOSpace<T>& occ = this->template get<MOSpace<T> >("occ");
    const MOSpace<T>& vrt = this->template get<MOSpace<T> >("vrt");

    const SymmetryBlockedTensor<T>& Fa = this->template get<SymmetryBlockedTensor<T> >("Fa");
    const SymmetryBlockedTensor<T>& Fb = this->template get<SymmetryBlockedTensor<T> >("Fb");

    this->put("H", new TwoElectronOperator<T>("V", OneElectronOperator<T>("f", occ, vrt, Fa, Fb)));
    TwoElectronOperator<T>& H = this->template get<TwoElectronOperator<T> >("H");

    const CholeskyIntegrals<T>& chol = this->template get<CholeskyIntegrals<T> >("cholesky");

    const SymmetryBlockedTensor<T>& cA = vrt.Calpha;
    const SymmetryBlockedTensor<T>& ca = vrt.Cbeta;
    const SymmetryBlockedTensor<T>& cI = occ.Calpha;
    const SymmetryBlockedTensor<T>& ci = occ.Cbeta;
    const SymmetryBlockedTensor<T>& Lpq = chol.getL();
    const SymmetryBlockedTensor<T>& D = chol.getD();

    const PointGroup& group = occ.group;

    const vector<int>& N = occ.nao;
    const vector<int>& nI = occ.nalpha;
    const vector<int>& ni = occ.nbeta;
    const vector<int>& nA = vrt.nalpha;
    const vector<int>& na = vrt.nbeta;
    int R = chol.getRank();

    vector<vector<int> > sizeIIR = vec(nI, nI, vec(R));
    vector<vector<int> > sizeiiR = vec(ni, ni, vec(R));
    vector<vector<int> > sizeAAR = vec(nA, nA, vec(R));
    vector<vector<int> > sizeaaR = vec(na, na, vec(R));
    vector<vector<int> > sizeAIR = vec(nA, nI, vec(R));
    vector<vector<int> > sizeaiR = vec(na, ni, vec(R));

    vector<int> shapeNNN = vec(NS, NS, NS);

    SymmetryBlockedTensor<T> LIJ("LIJ", arena, group, 3, sizeIIR, shapeNNN, false);
    SymmetryBlockedTensor<T> Lij("Lij", arena, group, 3, sizeiiR, shapeNNN, false);
    SymmetryBlockedTensor<T> LAB("LAB", arena, group, 3, sizeAAR, shapeNNN, false);
    SymmetryBlockedTensor<T> Lab("Lab", arena, group, 3, sizeaaR, shapeNNN, false);
    SymmetryBlockedTensor<T> LAI("LAI", arena, group, 3, sizeAIR, shapeNNN, false);
    SymmetryBlockedTensor<T> Lai("Lai", arena, group, 3, sizeaiR, shapeNNN, false);

    {
        vector<vector<int> > sizeNIR = vec(N, nI, vec(R));
        vector<vector<int> > sizeNiR = vec(N, ni, vec(R));
        vector<vector<int> > sizeNAR = vec(N, nA, vec(R));
        vector<vector<int> > sizeNaR = vec(N, na, vec(R));

        SymmetryBlockedTensor<T> LpI("LpI", arena, group, 3, sizeNIR, shapeNNN, false);
        SymmetryBlockedTensor<T> Lpi("Lpi", arena, group, 3, sizeNiR, shapeNNN, false);
        SymmetryBlockedTensor<T> LpA("LpA", arena, group, 3, sizeNAR, shapeNNN, false);
        SymmetryBlockedTensor<T> Lpa("Lpa", arena, group, 3, sizeNaR, shapeNNN, false);

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

    SymmetryBlockedTensor<T> LDIJ("LDIJ", arena, group, 3, sizeIIR, shapeNNN, false);
    SymmetryBlockedTensor<T> LDij("LDij", arena, group, 3, sizeiiR, shapeNNN, false);
    SymmetryBlockedTensor<T> LDAB("LDAB", arena, group, 3, sizeAAR, shapeNNN, false);
    SymmetryBlockedTensor<T> LDab("LDab", arena, group, 3, sizeaaR, shapeNNN, false);
    SymmetryBlockedTensor<T> LDAI("LDAI", arena, group, 3, sizeAIR, shapeNNN, false);
    SymmetryBlockedTensor<T> LDai("LDai", arena, group, 3, sizeaiR, shapeNNN, false);

    LDIJ["IJR"] = D["R"]*LIJ["IJR"];
    LDij["ijR"] = D["R"]*Lij["ijR"];
    LDAI["AIR"] = D["R"]*LAI["AIR"];
    LDai["aiR"] = D["R"]*Lai["aiR"];
    LDAB["ABR"] = D["R"]*LAB["ABR"];
    LDab["abR"] = D["R"]*Lab["abR"];

    H.getIJKL()(vec(0,2),vec(0,2))["IJKL"] = 0.5*LDIJ["IKR"]*LIJ["JLR"];
    H.getIJKL()(vec(0,1),vec(0,1))["IjKl"] =     LDIJ["IKR"]*Lij["jlR"];
    H.getIJKL()(vec(0,0),vec(0,0))["ijkl"] = 0.5*LDij["ikR"]*Lij["jlR"];

    H.getIJAK()(vec(0,2),vec(1,1))["IJAK"] =  LDIJ["JKR"]*LAI["AIR"];
    H.getIJAK()(vec(0,1),vec(1,0))["IjAk"] =  LDij["jkR"]*LAI["AIR"];
    H.getIJAK()(vec(0,1),vec(0,1))["IjaK"] = -LDIJ["IKR"]*Lai["ajR"];
    H.getIJAK()(vec(0,0),vec(0,0))["ijak"] =  LDij["jkR"]*Lai["aiR"];

    H.getAIJK()(vec(1,1),vec(0,2))["AIJK"] = H.getIJAK()(vec(0,2),vec(1,1))["JKAI"];
    H.getAIJK()(vec(1,0),vec(0,1))["AiJk"] = H.getIJAK()(vec(0,1),vec(1,0))["JkAi"];
    H.getAIJK()(vec(0,1),vec(0,1))["aIJk"] = H.getIJAK()(vec(0,1),vec(0,1))["JkaI"];
    H.getAIJK()(vec(0,0),vec(0,0))["aijk"] = H.getIJAK()(vec(0,0),vec(0,0))["jkai"];

    H.getABIJ()(vec(2,0),vec(0,2))["ABIJ"] = 0.5*LDAI["AIR"]*LAI["BJR"];
    H.getABIJ()(vec(1,0),vec(0,1))["AbIj"] =     LDAI["AIR"]*Lai["bjR"];
    H.getABIJ()(vec(0,0),vec(0,0))["abij"] = 0.5*LDai["aiR"]*Lai["bjR"];

    H.getIJAB()(vec(0,2),vec(2,0))["IJAB"] = H.getABIJ()(vec(2,0),vec(0,2))["ABIJ"];
    H.getIJAB()(vec(0,1),vec(1,0))["IjAb"] = H.getABIJ()(vec(1,0),vec(0,1))["AbIj"];
    H.getIJAB()(vec(0,0),vec(0,0))["ijab"] = H.getABIJ()(vec(0,0),vec(0,0))["abij"];

    H.getAIBJ()(vec(1,1),vec(1,1))["AIBJ"]  = LDAB["ABR"]*LIJ["IJR"];
    H.getAIBJ()(vec(1,1),vec(1,1))["AIBJ"] -= LDAI["AJR"]*LAI["BIR"];
    H.getAIBJ()(vec(1,0),vec(1,0))["AiBj"]  = LDAB["ABR"]*Lij["ijR"];
    H.getAIBJ()(vec(0,1),vec(0,1))["aIbJ"]  = LDab["abR"]*LIJ["IJR"];
    H.getAIBJ()(vec(0,0),vec(0,0))["aibj"]  = LDab["abR"]*Lij["ijR"];
    H.getAIBJ()(vec(0,0),vec(0,0))["aibj"] -= LDai["ajR"]*Lai["biR"];
    H.getAIBJ()(vec(1,0),vec(0,1))["AibJ"]  = -H.getABIJ()(vec(1,0),vec(0,1))["AbJi"];
    H.getAIBJ()(vec(0,1),vec(1,0))["aIBj"]  = -H.getABIJ()(vec(1,0),vec(0,1))["BaIj"];

    H.getABCI()(vec(2,0),vec(1,1))["ABCI"] =  LDAB["ACR"]*LAI["BIR"];
    H.getABCI()(vec(1,0),vec(1,0))["AbCi"] =  LDAB["ACR"]*Lai["biR"];
    H.getABCI()(vec(1,0),vec(0,1))["AbcI"] = -LDab["bcR"]*LAI["AIR"];
    H.getABCI()(vec(0,0),vec(0,0))["abci"] =  LDab["acR"]*Lai["biR"];

    H.getAIBC()(vec(1,1),vec(2,0))["AIBC"] = H.getABCI()(vec(2,0),vec(1,1))["BCAI"];
    H.getAIBC()(vec(1,0),vec(1,0))["AiBc"] = H.getABCI()(vec(1,0),vec(1,0))["BcAi"];
    H.getAIBC()(vec(0,1),vec(1,0))["aIBc"] = H.getABCI()(vec(1,0),vec(0,1))["BcaI"];
    H.getAIBC()(vec(0,0),vec(0,0))["aibc"] = H.getABCI()(vec(0,0),vec(0,0))["bcai"];

    H.getABCD()(vec(2,0),vec(2,0))["ABCD"] = 0.5*LDAB["ACR"]*LAB["BDR"];
    H.getABCD()(vec(1,0),vec(1,0))["AbCd"] =     LDAB["ACR"]*Lab["bdR"];
    H.getABCD()(vec(0,0),vec(0,0))["abcd"] = 0.5*LDab["acR"]*Lab["bdR"];
}

INSTANTIATE_SPECIALIZATIONS(CholeskyMOIntegrals);
REGISTER_TASK(CholeskyMOIntegrals<double>,"choleskymoints");
