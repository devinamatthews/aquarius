#include "choleskymoints.hpp"

using namespace aquarius::scf;
using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace op
{

template <typename T>
CholeskyMOIntegrals<T>::CholeskyMOIntegrals(const string& name, Config& config)
: MOIntegrals<T>(name, config)
{
    this->getProduct("H").addRequirement("cholesky", "cholesky");
}

template <typename T>
bool CholeskyMOIntegrals<T>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& occ = this->template get<MOSpace<T>>("occ");
    const auto& vrt = this->template get<MOSpace<T>>("vrt");

    const auto& Fa = this->template get<SymmetryBlockedTensor<T>>("Fa");
    const auto& Fb = this->template get<SymmetryBlockedTensor<T>>("Fb");

    auto& H = this->put("H", new TwoElectronOperator<T>("V", OneElectronOperator<T>("f", occ, vrt, Fa, Fb)));

    const auto& chol = this->template get<CholeskyIntegrals<T>>("cholesky");

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

    vector<vector<int>> sizeIIR = {nI, nI, {R}};
    vector<vector<int>> sizeiiR = {ni, ni, {R}};
    vector<vector<int>> sizeAAR = {nA, nA, {R}};
    vector<vector<int>> sizeaaR = {na, na, {R}};
    vector<vector<int>> sizeAIR = {nA, nI, {R}};
    vector<vector<int>> sizeaiR = {na, ni, {R}};

    vector<int> shapeNNN = {NS, NS, NS};

    SymmetryBlockedTensor<T> LIJ("LIJ", arena, group, 3, sizeIIR, shapeNNN, false);
    SymmetryBlockedTensor<T> Lij("Lij", arena, group, 3, sizeiiR, shapeNNN, false);
    SymmetryBlockedTensor<T> LAB("LAB", arena, group, 3, sizeAAR, shapeNNN, false);
    SymmetryBlockedTensor<T> Lab("Lab", arena, group, 3, sizeaaR, shapeNNN, false);
    SymmetryBlockedTensor<T> LAI("LAI", arena, group, 3, sizeAIR, shapeNNN, false);
    SymmetryBlockedTensor<T> Lai("Lai", arena, group, 3, sizeaiR, shapeNNN, false);

    {
        vector<vector<int>> sizeNIR = {N, nI, {R}};
        vector<vector<int>> sizeNiR = {N, ni, {R}};
        vector<vector<int>> sizeNAR = {N, nA, {R}};
        vector<vector<int>> sizeNaR = {N, na, {R}};

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

    H.getIJKL()({0,2},{0,2})["IJKL"] = 0.5*LDIJ["IKR"]*LIJ["JLR"];
    H.getIJKL()({0,1},{0,1})["IjKl"] =     LDIJ["IKR"]*Lij["jlR"];
    H.getIJKL()({0,0},{0,0})["ijkl"] = 0.5*LDij["ikR"]*Lij["jlR"];

    H.getIJAK()({0,2},{1,1})["IJAK"] =  LDIJ["JKR"]*LAI["AIR"];
    H.getIJAK()({0,1},{1,0})["IjAk"] =  LDij["jkR"]*LAI["AIR"];
    H.getIJAK()({0,1},{0,1})["IjaK"] = -LDIJ["IKR"]*Lai["ajR"];
    H.getIJAK()({0,0},{0,0})["ijak"] =  LDij["jkR"]*Lai["aiR"];

    H.getAIJK()({1,1},{0,2})["AIJK"] = H.getIJAK()({0,2},{1,1})["JKAI"];
    H.getAIJK()({1,0},{0,1})["AiJk"] = H.getIJAK()({0,1},{1,0})["JkAi"];
    H.getAIJK()({0,1},{0,1})["aIJk"] = H.getIJAK()({0,1},{0,1})["JkaI"];
    H.getAIJK()({0,0},{0,0})["aijk"] = H.getIJAK()({0,0},{0,0})["jkai"];

    H.getABIJ()({2,0},{0,2})["ABIJ"] = 0.5*LDAI["AIR"]*LAI["BJR"];
    H.getABIJ()({1,0},{0,1})["AbIj"] =     LDAI["AIR"]*Lai["bjR"];
    H.getABIJ()({0,0},{0,0})["abij"] = 0.5*LDai["aiR"]*Lai["bjR"];

    H.getIJAB()({0,2},{2,0})["IJAB"] = H.getABIJ()({2,0},{0,2})["ABIJ"];
    H.getIJAB()({0,1},{1,0})["IjAb"] = H.getABIJ()({1,0},{0,1})["AbIj"];
    H.getIJAB()({0,0},{0,0})["ijab"] = H.getABIJ()({0,0},{0,0})["abij"];

    H.getAIBJ()({1,1},{1,1})["AIBJ"]  = LDAB["ABR"]*LIJ["IJR"];
    H.getAIBJ()({1,1},{1,1})["AIBJ"] -= LDAI["AJR"]*LAI["BIR"];
    H.getAIBJ()({1,0},{1,0})["AiBj"]  = LDAB["ABR"]*Lij["ijR"];
    H.getAIBJ()({0,1},{0,1})["aIbJ"]  = LDab["abR"]*LIJ["IJR"];
    H.getAIBJ()({0,0},{0,0})["aibj"]  = LDab["abR"]*Lij["ijR"];
    H.getAIBJ()({0,0},{0,0})["aibj"] -= LDai["ajR"]*Lai["biR"];
    H.getAIBJ()({1,0},{0,1})["AibJ"]  = -H.getABIJ()({1,0},{0,1})["AbJi"];
    H.getAIBJ()({0,1},{1,0})["aIBj"]  = -H.getABIJ()({1,0},{0,1})["BaIj"];

    H.getABCI()({2,0},{1,1})["ABCI"] =  LDAB["ACR"]*LAI["BIR"];
    H.getABCI()({1,0},{1,0})["AbCi"] =  LDAB["ACR"]*Lai["biR"];
    H.getABCI()({1,0},{0,1})["AbcI"] = -LDab["bcR"]*LAI["AIR"];
    H.getABCI()({0,0},{0,0})["abci"] =  LDab["acR"]*Lai["biR"];

    H.getAIBC()({1,1},{2,0})["AIBC"] = H.getABCI()({2,0},{1,1})["BCAI"];
    H.getAIBC()({1,0},{1,0})["AiBc"] = H.getABCI()({1,0},{1,0})["BcAi"];
    H.getAIBC()({0,1},{1,0})["aIBc"] = H.getABCI()({1,0},{0,1})["BcaI"];
    H.getAIBC()({0,0},{0,0})["aibc"] = H.getABCI()({0,0},{0,0})["bcai"];

    H.getABCD()({2,0},{2,0})["ABCD"] = 0.5*LDAB["ACR"]*LAB["BDR"];
    H.getABCD()({1,0},{1,0})["AbCd"] =     LDAB["ACR"]*Lab["bdR"];
    H.getABCD()({0,0},{0,0})["abcd"] = 0.5*LDab["acR"]*Lab["bdR"];
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::op::CholeskyMOIntegrals);
REGISTER_TASK(aquarius::op::CholeskyMOIntegrals<double>,"choleskymoints");
