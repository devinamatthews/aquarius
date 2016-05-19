#include "sparseaomoints.hpp"

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
T absmax(const vector<T>& c)
{
    //return c[iamax(c.size(), c.data(), 1)];
    return nrm2(c.size(), c.data(), 1);
}

template <typename T>
SparseAOMOIntegrals<T>::SparseAOMOIntegrals(const string& name, Config& config)
: MOIntegrals<T>(name, config)
{
    this->getProduct("H").addRequirement("eri", "I");
}

template <typename T>
bool SparseAOMOIntegrals<T>::run(TaskDAG& dag, const Arena& arena)
{
    CTF_Timer_epoch ep("SparseAOMOIntegrals");
    ep.begin();
    const auto& occ = this->template get<MOSpace<T>>("occ");
    const auto& vrt = this->template get<MOSpace<T>>("vrt");

    const auto& ints = this->template get<ERI>("I");

    const SymmetryBlockedTensor<T>& cA_ = vrt.Calpha;
    const SymmetryBlockedTensor<T>& ca_ = vrt.Cbeta;
    const SymmetryBlockedTensor<T>& cI_ = occ.Calpha;
    const SymmetryBlockedTensor<T>& ci_ = occ.Cbeta;

    int n = ints.group.getNumIrreps();
    const vector<int>& N = occ.nao;
    const vector<int>& nI = occ.nalpha;
    const vector<int>& ni = occ.nbeta;
    const vector<int>& nA = vrt.nalpha;
    const vector<int>& na = vrt.nbeta;

    vector<int> nI_{sum(nI)}; nI_ += vector<int>(n-1,0);
    vector<int> ni_{sum(ni)}; ni_ += vector<int>(n-1,0);
    vector<int> nA_{sum(nA)}; nA_ += vector<int>(n-1,0);
    vector<int> na_{sum(na)}; na_ += vector<int>(n-1,0);

    auto& Ea = this->template get<vector<vector<real_type_t<T>>>>("Ea");
    auto& Eb = this->template get<vector<vector<real_type_t<T>>>>("Eb");

    auto& Fa = this->template get<SymmetryBlockedTensor<T>>("Fa");
    auto& Fb = this->template get<SymmetryBlockedTensor<T>>("Fb");

    //this->put("H", new TwoElectronOperator<T>("V", OneElectronOperator<T>("f", arena, occ, vrt)));
    auto& H = this->put("H", new TwoElectronOperator<T>("V", OneElectronOperator<T>("f", occ, vrt, Fa, Fb)));

    /*
    {
        vector<int> zero(n);
        SymmetryBlockedTensor<T> Ea_occ("Ea_occ", arena, ints.group, 1, {nI}, {NS}, false);
        SymmetryBlockedTensor<T> Eb_occ("Eb_occ", arena, ints.group, 1, {ni}, {NS}, false);
        SymmetryBlockedTensor<T> Ea_vrt("Ea_vrt", arena, ints.group, 1, {nA}, {NS}, false);
        SymmetryBlockedTensor<T> Eb_vrt("Eb_vrt", arena, ints.group, 1, {na}, {NS}, false);

        Ea_occ.slice(1.0, false, Ea, {zero}, 0.0);
        Ea_vrt.slice(1.0, false, Ea, {nI}, 0.0);
        Eb_occ.slice(1.0, false, Eb, {zero}, 0.0);
        Eb_vrt.slice(1.0, false, Eb, {ni}, 0.0);

        H.getIJ()({0,1},{0,1})["II"] = Ea_occ["I"];
        H.getIJ()({0,0},{0,0})["ii"] = Eb_occ["i"];
        H.getAB()({1,0},{1,0})["AA"] = Ea_vrt["A"];
        H.getAB()({0,0},{0,0})["aa"] = Eb_vrt["a"];
    }
    */

    SymmetryBlockedTensor<T> ABIJ__("<AB|IJ>", arena, ints.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false);
    SymmetryBlockedTensor<T> abij__("<ab|ij>", arena, ints.group, 4, {na,na,ni,ni}, {NS,NS,NS,NS}, false);

    vector<vector<T>> cA(n), ca(n), cI(n), ci(n);

    /*
     * Read transformation coefficients
     */
    for (int i = 0;i < n;i++)
    {
        vector<int> irreps = {i,i};
        cA_.getAllData(irreps, cA[i]);
        assert(cA[i].size() == N[i]*nA[i]);
        ca_.getAllData(irreps, ca[i]);
        assert(ca[i].size() == N[i]*na[i]);
        cI_.getAllData(irreps, cI[i]);
        assert(cI[i].size() == N[i]*nI[i]);
        ci_.getAllData(irreps, ci[i]);
        assert(ci[i].size() == N[i]*ni[i]);
    }

    #define SHOWIT(name) cout << #name ": " << absmax(name.ints) << endl;

    /*
     * Resort integrals so that each node has (pq|r_k s_l) where pq
     * are sparse blocks for each sparse rs pair
     */
    pqrs_integrals<T> pqrs(N, ints);
    pqrs.collect(true);

    /*
     * First quarter-transformation
     */
    //SHOWIT(PQrs);
    pqrs_integrals<T> pArs = pqrs.transform(B, nA, cA, true);
    //SHOWIT(PArs);
    pqrs_integrals<T> pars = pqrs.transform(B, na, ca, true);
    //SHOWIT(Pars);
    pqrs_integrals<T> pIrs = pqrs.transform(B, nI, cI, true);
    //SHOWIT(PIrs);
    pqrs_integrals<T> pirs = pqrs.transform(B, ni, ci, true);
    //SHOWIT(Pirs);
    pqrs.free();

    /*
     * Second quarter-transformation
     */
    pqrs_integrals<T> ABrs = pArs.transform(A, nA, cA, false);
    //SHOWIT(ABrs);
    pArs.free();
    pqrs_integrals<T> abrs = pars.transform(A, na, ca, false);
    //SHOWIT(abrs);
    pars.free();
    pqrs_integrals<T> AIrs = pIrs.transform(A, nA, cA, false);
    //SHOWIT(AIrs);
    pqrs_integrals<T> IJrs = pIrs.transform(A, nI, cI, false);
    //SHOWIT(IJrs);
    pIrs.free();
    pqrs_integrals<T> airs = pirs.transform(A, na, ca, false);
    //SHOWIT(airs);
    pqrs_integrals<T> ijrs = pirs.transform(A, ni, ci, false);
    //SHOWIT(ijrs);
    pirs.free();

    /*
     * Make <AB||CD>
     */
    pqrs_integrals<T>& rsAB = ABrs;
    rsAB.transpose();
    rsAB.collect(false);

    pqrs_integrals<T> rDAB = rsAB.transform(B, nA, cA, true);
    rsAB.free();

    pqrs_integrals<T> CDAB = rDAB.transform(A, nA, cA, false);
    //SHOWIT(CDAB);
    rDAB.free();
    CDAB.transcribe(H.getABCD()({2,0},{2,0}), true, true, NONE);
    CDAB.free();

    /*
     * Make <Ab|Cd> and <ab||cd>
     */
    pqrs_integrals<T>& rsab = abrs;
    rsab.transpose();
    rsab.collect(false);

    pqrs_integrals<T> rDab = rsab.transform(B, nA, cA, true);
    //SHOWIT(RDab);
    pqrs_integrals<T> rdab = rsab.transform(B, na, ca, true);
    //SHOWIT(Rdab);
    rsab.free();

    pqrs_integrals<T> CDab = rDab.transform(A, nA, cA, false);
    //SHOWIT(CDab);
    rDab.free();
    CDab.transcribe(H.getABCD()({1,0},{1,0}), false, false, NONE);
    CDab.free();

    pqrs_integrals<T> cdab = rdab.transform(A, na, ca, false);
    //SHOWIT(cdab);
    rdab.free();
    cdab.transcribe(H.getABCD()({0,0},{0,0}), true, true, NONE);
    cdab.free();

    /*
     * Make <AB||CI>, <Ab|cI>, and <AB|IJ>
     */
    pqrs_integrals<T>& rsAI = AIrs;
    rsAI.transpose();
    rsAI.collect(false);

    pqrs_integrals<T> rCAI = rsAI.transform(B, nA, cA, true);
    //SHOWIT(RCAI);
    pqrs_integrals<T> rcAI = rsAI.transform(B, na, ca, true);
    //SHOWIT(RcAI);
    pqrs_integrals<T> rJAI = rsAI.transform(B, nI, cI, true);
    //SHOWIT(RJAI);
    rsAI.free();

    pqrs_integrals<T> BCAI = rCAI.transform(A, nA, cA, false);
    //SHOWIT(BCAI);
    rCAI.free();
    BCAI.transcribe(H.getABCI()({2,0},{1,1}), true, false, NONE);
    BCAI.free();

    pqrs_integrals<T> bcAI = rcAI.transform(A, na, ca, false);
    //SHOWIT(bcAI);
    rcAI.free();
    bcAI.transcribe(H.getABCI()({1,0},{0,1}), false, false, PQ);
    bcAI.free();

    pqrs_integrals<T> BJAI = rJAI.transform(A, nA, cA, false);
    //SHOWIT(BJAI);
    rJAI.free();
    BJAI.transcribe(ABIJ__, false, false, NONE);
    BJAI.free();

    /*
     * Make <Ab|Ci>, <ab||ci>, <Ab|Ij>, and <ab|ij>
     */
    pqrs_integrals<T>& rsai = airs;
    rsai.transpose();
    rsai.collect(false);

    pqrs_integrals<T> rCai = rsai.transform(B, nA, cA, true);
    //SHOWIT(RCai);
    pqrs_integrals<T> rcai = rsai.transform(B, na, ca, true);
    //SHOWIT(Rcai);
    pqrs_integrals<T> rJai = rsai.transform(B, nI, cI, true);
    //SHOWIT(RJai);
    pqrs_integrals<T> rjai = rsai.transform(B, ni, ci, true);
    //SHOWIT(Rjai);
    rsai.free();

    pqrs_integrals<T> BCai = rCai.transform(A, nA, cA, false);
    //SHOWIT(BCai);
    rCai.free();
    BCai.transcribe(H.getABCI()({1,0},{1,0}), false, false, NONE);
    BCai.free();

    pqrs_integrals<T> bcai = rcai.transform(A, na, ca, false);
    //SHOWIT(bcai);
    rcai.free();
    bcai.transcribe(H.getABCI()({0,0},{0,0}), true, false, NONE);
    bcai.free();

    pqrs_integrals<T> BJai = rJai.transform(A, nA, cA, false);
    //SHOWIT(BJai);
    rJai.free();
    BJai.transcribe(H.getABIJ()({1,0},{0,1}), false, false, NONE);
    BJai.free();

    pqrs_integrals<T> bjai = rjai.transform(A, na, ca, false);
    //SHOWIT(bjai);
    rjai.free();
    bjai.transcribe(abij__, false, false, NONE);
    bjai.free();

    /*
     * Make <IJ||KL>, <AI||JK>, <aI|Jk>, <aI|bJ>, and <AI|BJ>
     */
    pqrs_integrals<T>& rsIJ = IJrs;
    rsIJ.transpose();
    rsIJ.collect(false);

    pqrs_integrals<T> rBIJ = rsIJ.transform(B, nA, cA, true);
    //SHOWIT(RBIJ);
    pqrs_integrals<T> rbIJ = rsIJ.transform(B, na, ca, true);
    //SHOWIT(RbIJ);
    pqrs_integrals<T> rLIJ = rsIJ.transform(B, nI, cI, true);
    //SHOWIT(RLIJ);
    pqrs_integrals<T> rlIJ = rsIJ.transform(B, ni, ci, true);
    //SHOWIT(RlIJ);
    rsIJ.free();

    pqrs_integrals<T> ABIJ = rBIJ.transform(A, nA, cA, false);
    //SHOWIT(ABIJ);
    rBIJ.free();
    ABIJ.transcribe(H.getAIBJ()({1,1},{1,1}), false, false, NONE);
    ABIJ.free();

    pqrs_integrals<T> abIJ = rbIJ.transform(A, na, ca, false);
    //SHOWIT(abIJ);
    rbIJ.free();
    abIJ.transcribe(H.getAIBJ()({0,1},{0,1}), false, false, NONE);
    abIJ.free();

    pqrs_integrals<T> akIJ = rlIJ.transform(A, na, ca, false);
    //SHOWIT(akIJ);
    rlIJ.free();
    akIJ.transcribe(H.getAIJK()({0,1},{0,1}), false, false, RS);
    akIJ.free();

    pqrs_integrals<T> AKIJ = rLIJ.transform(A, nA, cA, false);
    //SHOWIT(AKIJ);
    pqrs_integrals<T> KLIJ = rLIJ.transform(A, nI, cI, false);
    //SHOWIT(KLIJ);
    rLIJ.free();
    AKIJ.transcribe(H.getAIJK()({1,1},{0,2}), false, true, NONE);
    AKIJ.free();
    KLIJ.transcribe(H.getIJKL()({0,2},{0,2}), true, true, NONE);
    KLIJ.free();

    /*
     * Make <Ij|Kl>, <ij||kl>, <Ai|Jk>, <ai||jk>, <Ai|Bj>, and <ai|bj>
     */
    pqrs_integrals<T>& rsij = ijrs;
    rsij.transpose();
    rsij.collect(false);

    pqrs_integrals<T> rBij = rsij.transform(B, nA, cA, true);
    //SHOWIT(RBij);
    pqrs_integrals<T> rbij = rsij.transform(B, na, ca, true);
    //SHOWIT(Rbij);
    pqrs_integrals<T> rLij = rsij.transform(B, nI, cI, true);
    //SHOWIT(RLij);
    pqrs_integrals<T> rlij = rsij.transform(B, ni, ci, true);
    //SHOWIT(Rlij);
    rsij.free();

    pqrs_integrals<T> ABij = rBij.transform(A, nA, cA, false);
    //SHOWIT(ABij);
    rBij.free();
    ABij.transcribe(H.getAIBJ()({1,0},{1,0}), false, false, NONE);
    ABij.free();

    pqrs_integrals<T> abij = rbij.transform(A, na, ca, false);
    //SHOWIT(abij);
    rbij.free();
    abij.transcribe(H.getAIBJ()({0,0},{0,0}), false, false, NONE);
    abij.free();

    pqrs_integrals<T> AKij = rLij.transform(A, nA, cA, false);
    //SHOWIT(AKij);
    pqrs_integrals<T> KLij = rLij.transform(A, nI, cI, false);
    //SHOWIT(KLij);
    rLij.free();
    AKij.transcribe(H.getAIJK()({1,0},{0,1}), false, false, NONE);
    AKij.free();
    KLij.transcribe(H.getIJKL()({0,1},{0,1}), false, false, NONE);
    KLij.free();

    pqrs_integrals<T> akij = rlij.transform(A, na, ca, false);
    //SHOWIT(akij);
    pqrs_integrals<T> klij = rlij.transform(A, ni, ci, false);
    //SHOWIT(klij);
    rlij.free();
    akij.transcribe(H.getAIJK()({0,0},{0,0}), false, true, NONE);
    akij.free();
    klij.transcribe(H.getIJKL()({0,0},{0,0}), true, true, NONE);
    klij.free();

    /*
     * Make <AI||BJ> and <ai||bj>
     */
    H.getAIBJ()({1,1},{1,1})["AIBJ"] -= ABIJ__["ABJI"];
    H.getAIBJ()({0,0},{0,0})["aibj"] -= abij__["abji"];

    /*
     * Make <AB||IJ> and <ab||ij>
     */
    H.getABIJ()({2,0},{0,2})["ABIJ"] = 0.5*ABIJ__["ABIJ"];
    H.getABIJ()({0,0},{0,0})["abij"] = 0.5*abij__["abij"];

    /*
     * Make <Ai|bJ> = -<Ab|Ji> and <aI|Bj> = -<Ba|Ij>
     */
    H.getAIBJ()({1,0},{0,1})["AibJ"] = -H.getABIJ()({1,0},{0,1})["AbJi"];
    H.getAIBJ()({0,1},{1,0})["aIBj"] = -H.getABIJ()({1,0},{0,1})["BaIj"];

    /*
     * Fill in pieces which are equal by Hermicity
     */
    H.getIJAK()({0,2},{1,1})["JKAI"] = H.getAIJK()({1,1},{0,2})["AIJK"];
    H.getIJAK()({0,1},{1,0})["JkAi"] = H.getAIJK()({1,0},{0,1})["AiJk"];
    H.getIJAK()({0,1},{0,1})["JkaI"] = H.getAIJK()({0,1},{0,1})["aIJk"];
    H.getIJAK()({0,0},{0,0})["jkai"] = H.getAIJK()({0,0},{0,0})["aijk"];

    H.getAIBC()({1,1},{2,0})["AIBC"] = H.getABCI()({2,0},{1,1})["BCAI"];
    H.getAIBC()({1,0},{1,0})["AiBc"] = H.getABCI()({1,0},{1,0})["BcAi"];
    H.getAIBC()({0,1},{1,0})["aIBc"] = H.getABCI()({1,0},{0,1})["BcaI"];
    H.getAIBC()({0,0},{0,0})["aibc"] = H.getABCI()({0,0},{0,0})["bcai"];

    H.getIJAB()({0,2},{2,0})["IJAB"] = H.getABIJ()({2,0},{0,2})["ABIJ"];
    H.getIJAB()({0,1},{1,0})["IjAb"] = H.getABIJ()({1,0},{0,1})["AbIj"];
    H.getIJAB()({0,0},{0,0})["ijab"] = H.getABIJ()({0,0},{0,0})["abij"];

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

INSTANTIATE_SPECIALIZATIONS(aquarius::op::SparseAOMOIntegrals);
REGISTER_TASK(aquarius::op::SparseAOMOIntegrals<double>,"sparseaomoints");
