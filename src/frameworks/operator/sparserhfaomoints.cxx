#include "sparserhfaomoints.hpp"

#include "time/time.hpp"

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
SparseRHFAOMOIntegrals<T>::SparseRHFAOMOIntegrals(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs += Requirement("occspace", "occ");
    reqs += Requirement("vrtspace", "vrt");
    reqs += Requirement("Ea", "E");
    reqs += Requirement("Fa", "F");
    reqs += Requirement("eri", "I");
    addProduct("mofock", "f", reqs);
    addProduct("<Ab|Cd>", "VABCD", reqs);
    addProduct("<Ab|Ci>", "VABCI", reqs);
    addProduct("<Ab|Ij>", "VABIJ", reqs);
    addProduct("<Ai|Bj>", "VAIBJ", reqs);
    addProduct("<Ai|Jb>", "VAIJB", reqs);
    addProduct("<Ai|Jk>", "VAIJK", reqs);
    addProduct("<Ij|Kl>", "VIJKL", reqs);
}

template <typename T>
bool SparseRHFAOMOIntegrals<T>::run(TaskDAG& dag, const Arena& arena)
{
    CTF_Timer_epoch ep("SparseRHFAOMOIntegrals");
    ep.begin();
    const auto& occ = this->template get<MOSpace<T>>("occ");
    const auto& vrt = this->template get<MOSpace<T>>("vrt");

    const auto& ints = this->template get<ERI>("I");

    const SymmetryBlockedTensor<T>& cA_ = vrt.Calpha;
    const SymmetryBlockedTensor<T>& cI_ = occ.Calpha;

    int n = ints.group.getNumIrreps();
    const vector<int>& N = occ.nao;
    const vector<int>& nI = occ.nalpha;
    const vector<int>& nA = vrt.nalpha;

    auto& E = this->template get<vector<vector<real_type_t<T>>>>("E");
    auto& F = this->template get<SymmetryBlockedTensor<T>>("F");

    auto&     f = this->put(    "f", new OneElectronOperator  <T>(  "f(AB)", occ, vrt, F, F));
    auto& VABCD = this->put("VABCD", new SymmetryBlockedTensor<T>("<Ab|Cd>", arena, occ.group, 4, {nA,nA,nA,nA}, {NS,NS,NS,NS}, false));
    auto& VABCI = this->put("VABCI", new SymmetryBlockedTensor<T>("<Ab|Ci>", arena, occ.group, 4, {nA,nA,nA,nI}, {NS,NS,NS,NS}, false));
    auto& VABIJ = this->put("VABIJ", new SymmetryBlockedTensor<T>("<Ab|Ij>", arena, occ.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));
    auto& VAIBJ = this->put("VAIBJ", new SymmetryBlockedTensor<T>("<Ai|Bj>", arena, occ.group, 4, {nA,nI,nA,nI}, {NS,NS,NS,NS}, false));
    auto& VAIJB = this->put("VAIJB", new SymmetryBlockedTensor<T>("<Ai|Jb>", arena, occ.group, 4, {nA,nI,nI,nA}, {NS,NS,NS,NS}, false));
    auto& VAIJK = this->put("VAIJK", new SymmetryBlockedTensor<T>("<Ai|Jk>", arena, occ.group, 4, {nA,nI,nI,nI}, {NS,NS,NS,NS}, false));
    auto& VIJKL = this->put("VIJKL", new SymmetryBlockedTensor<T>("<Ij|Kl>", arena, occ.group, 4, {nI,nI,nI,nI}, {NS,NS,NS,NS}, false));

    vector<vector<T>> cA(n), cI(n);

    /*
     * Read transformation coefficients
     */
    for (int i = 0;i < n;i++)
    {
        vector<int> irreps = {i,i};
        cA_.getAllData(irreps, cA[i]);
        assert(cA[i].size() == N[i]*nA[i]);
        cI_.getAllData(irreps, cI[i]);
        assert(cI[i].size() == N[i]*nI[i]);
    }

    /*
     * Resort integrals so that each node has (pq|r_k s_l) where pq
     * are sparse blocks for each sparse rs pair
     */
    pqrs_integrals<T> pqrs(N, ints);
    pqrs.collect(true);

    /*
     * First quarter-transformation
     */
    pqrs_integrals<T> pArs = pqrs.transform(B, nA, cA, true);
    pqrs_integrals<T> pIrs = pqrs.transform(B, nI, cI, true);
    pqrs.free();

    /*
     * Second quarter-transformation
     */
    pqrs_integrals<T> ABrs = pArs.transform(A, nA, cA, false);
    pArs.free();
    pqrs_integrals<T> AIrs = pIrs.transform(A, nA, cA, false);
    pqrs_integrals<T> IJrs = pIrs.transform(A, nI, cI, false);
    pIrs.free();

    /*
     * Make <Ab|Cd>
     */
    pqrs_integrals<T>& rsAB = ABrs;
    rsAB.transpose();
    rsAB.collect(false);

    pqrs_integrals<T> rDAB = rsAB.transform(B, nA, cA, true);
    rsAB.free();

    pqrs_integrals<T> CDAB = rDAB.transform(A, nA, cA, false);
    rDAB.free();
    CDAB.transcribe(VABCD, false, false, NONE);
    CDAB.free();

    /*
     * Make <Ab|Ci> and <Ab|Ij>
     */
    pqrs_integrals<T>& rsAI = AIrs;
    rsAI.transpose();
    rsAI.collect(false);

    pqrs_integrals<T> rCAI = rsAI.transform(B, nA, cA, true);
    pqrs_integrals<T> rJAI = rsAI.transform(B, nI, cI, true);
    rsAI.free();

    pqrs_integrals<T> BCAI = rCAI.transform(A, nA, cA, false);
    rCAI.free();
    BCAI.transcribe(VABCI, false, false, NONE);
    BCAI.free();

    pqrs_integrals<T> BJAI = rJAI.transform(A, nA, cA, false);
    rJAI.free();
    BJAI.transcribe(VABIJ, false, false, NONE);
    BJAI.free();

    /*
     * Make <Ij|Kl>, <Ij|Ka>, and <Ai|Bj>
     */
    pqrs_integrals<T>& rsIJ = IJrs;
    rsIJ.transpose();
    rsIJ.collect(false);

    pqrs_integrals<T> rBIJ = rsIJ.transform(B, nA, cA, true);
    pqrs_integrals<T> rLIJ = rsIJ.transform(B, nI, cI, true);
    rsIJ.free();

    pqrs_integrals<T> ABIJ = rBIJ.transform(A, nA, cA, false);
    rBIJ.free();
    ABIJ.transcribe(VAIBJ, false, false, NONE);
    ABIJ.free();

    pqrs_integrals<T> AKIJ = rLIJ.transform(A, nA, cA, false);
    pqrs_integrals<T> KLIJ = rLIJ.transform(A, nI, cI, false);
    rLIJ.free();
    AKIJ.transcribe(VAIJK, false, false, NONE);
    AKIJ.free();
    KLIJ.transcribe(VIJKL, false, false, NONE);
    KLIJ.free();

    /*
     * Make <Ai|Jb> = <Ab|Ji>
     */
    VAIJB["AiJb"] = VABIJ["AbJi"];

    //this->log(arena) << "ABCD: " << setprecision(15) << H.getABCD()({2,0},{2,0}).norm(2) << endl;
    //this->log(arena) << "AbCd: " << setprecision(15) << H.getABCD()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "abcd: " << setprecision(15) << H.getABCD()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "ABCI: " << setprecision(15) << H.getABCI()({2,0},{1,1}).norm(2) << endl;
    //this->log(arena) << "AbCi: " << setprecision(15) << H.getABCI()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "AbcI: " << setprecision(15) << H.getABCI()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "abci: " << setprecision(15) << H.getABCI()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "AIBC: " << setprecision(15) << H.getAIBC()({1,1},{2,0}).norm(2) << endl;
    //this->log(arena) << "AiBc: " << setprecision(15) << H.getAIBC()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "aIBc: " << setprecision(15) << H.getAIBC()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "aibc: " << setprecision(15) << H.getAIBC()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "ABIJ: " << setprecision(15) << H.getABIJ()({2,0},{0,2}).norm(2) << endl;
    //this->log(arena) << "AbIj: " << setprecision(15) << H.getABIJ()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "abij: " << setprecision(15) << H.getABIJ()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "AIBJ: " << setprecision(15) << H.getAIBJ()({1,1},{1,1}).norm(2) << endl;
    //this->log(arena) << "AiBj: " << setprecision(15) << H.getAIBJ()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "aIbJ: " << setprecision(15) << H.getAIBJ()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "AibJ: " << setprecision(15) << H.getAIBJ()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "aIBj: " << setprecision(15) << H.getAIBJ()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "aibj: " << setprecision(15) << H.getAIBJ()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IJAB: " << setprecision(15) << H.getIJAB()({0,2},{2,0}).norm(2) << endl;
    //this->log(arena) << "IjAb: " << setprecision(15) << H.getIJAB()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "ijab: " << setprecision(15) << H.getIJAB()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "AIJK: " << setprecision(15) << H.getAIJK()({1,1},{0,2}).norm(2) << endl;
    //this->log(arena) << "AiJk: " << setprecision(15) << H.getAIJK()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "aIJk: " << setprecision(15) << H.getAIJK()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "aijk: " << setprecision(15) << H.getAIJK()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IJAK: " << setprecision(15) << H.getIJAK()({0,2},{1,1}).norm(2) << endl;
    //this->log(arena) << "IjAk: " << setprecision(15) << H.getIJAK()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "IjaK: " << setprecision(15) << H.getIJAK()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "ijak: " << setprecision(15) << H.getIJAK()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IJKL: " << setprecision(15) << H.getIJKL()({0,2},{0,2}).norm(2) << endl;
    //this->log(arena) << "IjKl: " << setprecision(15) << H.getIJKL()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "ijkl: " << setprecision(15) << H.getIJKL()({0,0},{0,0}).norm(2) << endl;
    ep.end();

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::op::SparseRHFAOMOIntegrals);
REGISTER_TASK(aquarius::op::SparseRHFAOMOIntegrals<double>,"sparserhfaomoints");
