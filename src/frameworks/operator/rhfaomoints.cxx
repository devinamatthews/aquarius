#include "rhfaomoints.hpp"

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
RHFAOMOIntegrals<T>::RHFAOMOIntegrals(const string& name, Config& config)
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
bool RHFAOMOIntegrals<T>::run(TaskDAG& dag, const Arena& arena)
{
    CTF_Timer_epoch ep("RHFAOMOIntegrals");
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
     * are dense blocks for each sparse rs pair
     */
    pqrs_integrals<T> pqrs(N, ints);
    pqrs.collect(true);
    abrs_integrals<T> PQrs(pqrs, true);

    /*
     * First quarter-transformation
     */
    abrs_integrals<T> PArs = PQrs.transform(B, nA, cA);
    abrs_integrals<T> PIrs = PQrs.transform(B, nI, cI);
    PQrs.free();

    /*
     * Second quarter-transformation
     */
    abrs_integrals<T> ABrs = PArs.transform(A, nA, cA);
    PArs.free();
    abrs_integrals<T> AIrs = PIrs.transform(A, nA, cA);
    abrs_integrals<T> IJrs = PIrs.transform(A, nI, cI);
    PIrs.free();

    /*
     * Make <Ab|Cd>
     */
    pqrs_integrals<T> rsAB(ABrs);
    rsAB.collect(false);

    abrs_integrals<T> RSAB(rsAB, true);
    abrs_integrals<T> RDAB = RSAB.transform(B, nA, cA);
    RSAB.free();

    abrs_integrals<T> CDAB = RDAB.transform(A, nA, cA);
    RDAB.free();
    CDAB.transcribe(VABCD, false, false, NONE);
    CDAB.free();

    /*
     * Make <Ab|Ci> and <Ab|Ij>
     */
    pqrs_integrals<T> rsAI(AIrs);
    rsAI.collect(false);

    abrs_integrals<T> RSAI(rsAI, true);
    abrs_integrals<T> RCAI = RSAI.transform(B, nA, cA);
    abrs_integrals<T> RJAI = RSAI.transform(B, nI, cI);
    RSAI.free();

    abrs_integrals<T> BCAI = RCAI.transform(A, nA, cA);
    RCAI.free();
    BCAI.transcribe(VABCI, false, false, NONE);
    BCAI.free();

    abrs_integrals<T> BJAI = RJAI.transform(A, nA, cA);
    RJAI.free();
    BJAI.transcribe(VABIJ, false, false, NONE);
    BJAI.free();

    /*
     * Make <Ij|Kl>, <Ai|Jk>, and <Ai|Bj>
     */
    pqrs_integrals<T> rsIJ(IJrs);
    rsIJ.collect(false);

    abrs_integrals<T> RSIJ(rsIJ, true);
    abrs_integrals<T> RBIJ = RSIJ.transform(B, nA, cA);
    abrs_integrals<T> RLIJ = RSIJ.transform(B, nI, cI);
    RSIJ.free();

    abrs_integrals<T> ABIJ = RBIJ.transform(A, nA, cA);
    RBIJ.free();
    ABIJ.transcribe(VAIBJ, false, false, NONE);
    ABIJ.free();

    abrs_integrals<T> AKIJ = RLIJ.transform(A, nA, cA);
    abrs_integrals<T> KLIJ = RLIJ.transform(A, nI, cI);
    RLIJ.free();
    AKIJ.transcribe(VAIJK, false, false, NONE);
    AKIJ.free();
    KLIJ.transcribe(VIJKL, false, false, NONE);
    KLIJ.free();

    /*
     * Make <Ai|Jb> = <Ab|Ji>
     */
    VAIJB["AiJb"] = VABIJ["AbJi"];

    //this->log(arena) << "AbCd: " << setprecision(15) << VABCD.norm(2) << endl;
    //this->log(arena) << "AbCi: " << setprecision(15) << VABCI.norm(2) << endl;
    //this->log(arena) << "AbIj: " << setprecision(15) << VABIJ.norm(2) << endl;
    //this->log(arena) << "AiBj: " << setprecision(15) << VAIBJ.norm(2) << endl;
    //this->log(arena) << "AiJb: " << setprecision(15) << VAIJB.norm(2) << endl;
    //this->log(arena) << "AiJk: " << setprecision(15) << VAIJK.norm(2) << endl;
    //this->log(arena) << "IjKl: " << setprecision(15) << VIJKL.norm(2) << endl;
    ep.end();

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::op::RHFAOMOIntegrals);
REGISTER_TASK(aquarius::op::RHFAOMOIntegrals<double>,"rhfaomoints");
