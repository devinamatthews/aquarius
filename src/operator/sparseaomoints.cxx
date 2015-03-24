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
SparseAOMOIntegrals<T>::SparseAOMOIntegrals(const string& name, Config& config)
: MOIntegrals<T>(name, config)
{
    this->getProduct("H").addRequirement("eri", "I");
}

template <typename T>
SparseAOMOIntegrals<T>::pqrs_integrals::pqrs_integrals(const vector<int>& norb, const ERI& aoints)
: Distributed(aoints.arena), group(aoints.group)
{
    PROFILE_FUNCTION

    ns = nr = nq = np = norb;

    auto iidx = aoints.idxs.begin();
    auto iint = aoints.ints.begin();
    for (;iidx != aoints.idxs.end();++iidx, ++iint)
    {
        idx4_t idx = *iidx;

        if (idx.i > idx.j) swap(idx.i, idx.j);
        if (idx.k > idx.l) swap(idx.k, idx.l);

        ints.push_back(*iint);
        idxs.push_back(idx);

        if (idx.i != idx.k || idx.j != idx.l)
        {
            swap(idx.i, idx.k);
            swap(idx.j, idx.l);
            ints.push_back(*iint);
            idxs.push_back(idx);
        }
    }

    PROFILE_STOP
}

template <typename T>
void SparseAOMOIntegrals<T>::pqrs_integrals::free()
{
    ints.clear();
    idxs.clear();
}

template <typename T>
void SparseAOMOIntegrals<T>::pqrs_integrals::sortInts(bool rles, size_t& nrs, vector<size_t>& rscount)
{
    size_t nrtot = sum(nr);
    size_t nstot = sum(ns);

    if (rles)
    {
        assert(nr == ns);
        nrs = nrtot*(nrtot+1)/2;
    }
    else
    {
        nrs = nrtot*nstot;
    }

    rscount.assign(nrs, 0);

    for (size_t i = 0;i < ints.size();i++)
    {
        if (rles) rscount[idxs[i].k+idxs[i].l*(idxs[i].l+1)/2]++;
        else      rscount[idxs[i].k+idxs[i].l*nrtot]++;
    }

    vector<size_t> rsoff(nrs, 0);
    for (size_t i = 1;i < nrs;i++) rsoff[i] = rsoff[i-1]+rscount[i-1];
    assert(rsoff[nrs-1]+rscount[nrs-1] == ints.size());

    vector<T> newints(ints.size());
    vector<idx4_t> newidxs(ints.size());

    auto iidx = idxs.begin();
    auto iint = ints.begin();
    for (;iidx != idxs.end();++iidx, ++iint)
    {
        size_t rs = (rles ? iidx->k+iidx->l*(iidx->l+1)/2
                          : iidx->k+iidx->l*nrtot);

        newints[rsoff[rs]] = *iint;
        newidxs[rsoff[rs]] = *iidx;

        rsoff[rs]++;
    }

    copy(newints.begin(), newints.end(), ints.begin());
    copy(newidxs.begin(), newidxs.end(), idxs.begin());
}

/*
 * Redistribute integrals such that each node has all pq for each rs pair
 */
template <typename T>
void SparseAOMOIntegrals<T>::pqrs_integrals::collect(bool rles)
{
    PROFILE_FUNCTION

    static Datatype IDX4_T_TYPE;

    if (IDX4_T_TYPE == MPI_DATATYPE_NULL)
    {
        MPI_Type_contiguous(4, MPI_TYPE_<uint16_t>::value(), IDX4_T_TYPE);
        MPI_Type_commit(IDX4_T_TYPE);
    }

    size_t nrs;
    vector<size_t> rscount;
    sortInts(rles, nrs, rscount);

    int nproc = arena.size;

    vector<MPI_Int> sendcount(nproc, 0);
    vector<MPI_Int> recvcount(nproc, 0);

    for (int i = 0;i < nproc;i++)
    {
        for (size_t rs = (nrs*i)/nproc;rs < (nrs*(i+1))/nproc;rs++)
        {
            assert(rs >= 0 && rs < nrs);
            sendcount[i] += rscount[rs];
        }
    }

    PROFILE_SECTION(collect_comm)
    this->arena.Alltoall(sendcount, recvcount);
    PROFILE_STOP

    vector<T> oldints(ints.begin(), ints.end());
    vector<idx4_t> oldidxs(idxs.begin(), idxs.end());
    ints.clear(); ints.shrink_to_fit();
    idxs.clear(); idxs.shrink_to_fit();

    size_t nnewints = sum(recvcount);
    vector<T> newints(nnewints);
    vector<idx4_t> newidxs(nnewints);

    PROFILE_SECTION(collect_comm)
    this->arena.Alltoall(oldints, sendcount, newints, recvcount);
    this->arena.Alltoall(oldidxs, sendcount, newidxs, recvcount, IDX4_T_TYPE);
    PROFILE_STOP

    oldints.clear(); oldints.shrink_to_fit();
    oldidxs.clear(); oldidxs.shrink_to_fit();

    ints.assign(newints.begin(), newints.end());
    idxs.assign(newidxs.begin(), newidxs.end());

    newints.clear(); newints.shrink_to_fit();
    newidxs.clear(); newidxs.shrink_to_fit();

    sortInts(rles, nrs, rscount);

    PROFILE_STOP
}

template <typename T>
typename SparseAOMOIntegrals<T>::pqrs_integrals
SparseAOMOIntegrals<T>::pqrs_integrals::transform(Index index, const vector<int>& nc, const vector<vector<T> >& C, bool pleq)
{
    pqrs_integrals out(arena, group);

    int n = group.getNumIrreps();
    int nptot = sum(np);
    int nqtot = sum(nq);
    int nctot = sum(nc);

    vector<int> startc(n);
    for (int i = 1;i < n;i++) startc[i] = startc[i-1]+nc[i-1];

    vector<int> startp(n);
    for (int i = 1;i < n;i++) startp[i] = startp[i-1]+np[i-1];
    vector<int> startq(n);
    for (int i = 1;i < n;i++) startq[i] = startq[i-1]+nq[i-1];

    out.nr = nr;
    out.ns = ns;

    if (index == A)
    {
        out.np = nc;
        out.nq = nq;
    }
    else
    {
        out.np = np;
        out.nq = nc;
    }

    int r = -1;
    int s = -1;

    auto ibegin = idxs.begin();
    auto iend = idxs.end();
    auto iidx = ibegin;
    auto iint = ints.begin();

    matrix<double> before(nptot, nqtot), after;

    if (index == A)
    {
        after.reset(nctot, nqtot);
    }
    else
    {
        after.reset(nptot, nctot);
    }

    while (true)
    {
        if (iidx == iend || iidx->l != s || iidx->k != r)
        {
            if (iidx != ibegin)
            {
                int irrepr = 0;
                for (int r_ = r;r_ > nr[irrepr];r_ -= nr[irrepr], irrepr++);

                int irreps = 0;
                for (int s_ = s;s_ > ns[irreps];s_ -= ns[irreps], irreps++);

                for (int irrepc = 0;irrepc < n;irrepc++)
                {
                    if (index == A)
                    {
                        int irrepp = irrepc;
                        Representation rprs = group.getIrrep(irrepp)*
                                              group.getIrrep(irrepr)*
                                              group.getIrrep(irreps);

                        for (int irrepq = 0;irrepq < n;irrepq++)
                        {
                            if (!(group.getIrrep(irrepq)*rprs).isTotallySymmetric()) continue;

                            gemm('N', 'N', nq[irrepq], nc[irrepc], np[irrepp],
                                 1.0, &before[startp[irrepp]][startq[irrepq]], before.stride(0),
                                                             C[irrepc].data(),       np[irrepp],
                                 0.0,  &after[startc[irrepc]][startq[irrepq]],  after.stride(0));

                            for (int c = startc[irrepc];c < startc[irrepc]+nc[irrepc];c++)
                            {
                                for (int q = startq[irrepq];q < startq[irrepq]+nq[irrepq];q++)
                                {
                                    double val = after[c][q];
                                    if (abs(val) > 1e-12)
                                    {
                                        out.idxs.emplace_back(c, q, r, s);
                                        out.ints.push_back(val);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        int irrepq = irrepc;
                        Representation rqrs = group.getIrrep(irrepq)*
                                              group.getIrrep(irrepr)*
                                              group.getIrrep(irreps);

                        for (int irrepp = 0;irrepp < n;irrepp++)
                        {
                            if (!(group.getIrrep(irrepp)*rqrs).isTotallySymmetric()) continue;

                            gemm('T', 'N', nc[irrepc], np[irrepp], nq[irrepq],
                                 1.0,                        C[irrepc].data(),       nq[irrepq],
                                      &before[startp[irrepp]][startq[irrepq]], before.stride(0),
                                 0.0,  &after[startp[irrepp]][startc[irrepc]],  after.stride(0));

                            for (int p = startp[irrepp];p < startp[irrepp]+np[irrepp];p++)
                            {
                                for (int c = startc[irrepc];c < startc[irrepc]+nc[irrepc];c++)
                                {
                                    double val = after[p][c];
                                    if (abs(val) > 1e-12)
                                    {
                                        out.idxs.emplace_back(p, c, r, s);
                                        out.ints.push_back(val);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (iidx == iend)
            {
                break;
            }
            else
            {
                r = iidx->k;
                s = iidx->l;
                before = 0.0;
            }
        }

        before[iidx->i][iidx->j] = *iint;
        if (pleq && iidx->i != iidx->j)
            before[iidx->j][iidx->i] = *iint;

        ++iidx;
        ++iint;
    }

    return out;
}

template <typename T>
void SparseAOMOIntegrals<T>::pqrs_integrals::transcribe(SymmetryBlockedTensor<T>& symtensor, bool assympr, bool assymqs, Side side)
{
    PROFILE_FUNCTION

    if (assympr) assert(np == nr);
    if (assymqs) assert(nq == ns);

    int n = group.getNumIrreps();

    vector<int> irrepp;
    for (int i = 0;i < n;i++) irrepp += vector<int>(np[i],i);
    vector<int> irrepq;
    for (int i = 0;i < n;i++) irrepq += vector<int>(nq[i],i);

    vector<int> startr(n);
    for (int i = 1;i < n;i++) startr[i] = startr[i-1]+nr[i-1];
    vector<int> starts(n);
    for (int i = 1;i < n;i++) starts[i] = starts[i-1]+ns[i-1];

    vector<int> startp(n);
    for (int i = 1;i < n;i++) startp[i] = startp[i-1]+np[i-1];
    vector<int> startq(n);
    for (int i = 1;i < n;i++) startq[i] = startq[i-1]+nq[i-1];

    vector<size_t> offpq(n*n);
    vector<int> symirrs(4);
    Representation irrrs(group), irrqrs(group), irrpqrs(group);

    int nrtot = sum(nr);
    int nstot = sum(ns);

    vector<size_t> numrs(nrtot*nstot);
    for (auto& idx : idxs)
    {
        size_t rs = idx.k+idx.l*nrtot;
        numrs[rs]++;
    }

    vector<size_t> offrs(nrtot*nstot);
    for (int i = 1;i < nrtot*nstot;i++) offrs[i] = offrs[i-1]+numrs[i-1];

    /*
     * (pq|rs) -> <pr|qs>
     */
    for (int irrs = 0;irrs < n;irrs++)
    {
        for (int irrr = 0;irrr < n;irrr++)
        {
            irrrs = group.getIrrep(irrr)*group.getIrrep(irrs);

            matrix<vector<tkv_pair<T>>> pairs(n,n);

            for (int s = starts[irrs];s < starts[irrs]+ns[irrs];s++)
            {
                for (int r = startr[irrr];r < startr[irrr]+nr[irrr];r++)
                {
                    auto iidx = idxs.begin()+offrs[r+s*nrtot];
                    auto iend =         iidx+numrs[r+s*nrtot];
                    auto iint = ints.begin()+offrs[r+s*nrtot];
                    for (;iidx != iend;++iidx, ++iint)
                    {
                        int p = iidx->i;
                        int q = iidx->j;
                        int irrp = irrepp[p];
                        int irrq = irrepp[q];

                        if (assympr && p == r) continue;
                        if (assymqs && q == s) continue;
                        if (assympr && assymqs && q > s) continue;

                        if (side == PQ || (assympr && p > r))
                        {
                            if (side == RS || (assymqs && q > s))
                            {
                                pairs[irrp][irrq].emplace_back(((((int64_t)q)*ns[irrs]+s)*np[irrp]+p)*nr[irrr]+r, *iint);
                            }
                            else
                            {
                                pairs[irrp][irrq].emplace_back(((((int64_t)s)*nq[irrq]+q)*np[irrp]+p)*nr[irrr]+r, -*iint);
                            }
                        }
                        else
                        {
                            if (side == RS || (assymqs && q > s))
                            {
                                pairs[irrp][irrq].emplace_back(((((int64_t)q)*ns[irrs]+s)*nr[irrr]+r)*np[irrp]+p, -*iint);
                            }
                            else
                            {
                                pairs[irrp][irrq].emplace_back(((((int64_t)s)*nq[irrq]+q)*nr[irrr]+r)*np[irrp]+p, *iint);
                            }
                        }
                    }
                }
            }

            for (int irrq = 0;irrq < n;irrq++)
            {
                irrqrs = irrrs;
                irrqrs *= group.getIrrep(irrq);

                for (int irrp = 0;irrp < n;irrp++)
                {
                    irrpqrs = irrqrs;
                    irrpqrs *= group.getIrrep(irrp);

                    if (!irrpqrs.isTotallySymmetric()) continue;
                    if (assympr && assymqs && irrq > irrs) continue;

                    symirrs[0] = irrp;
                    symirrs[1] = irrr;
                    symirrs[2] = irrq;
                    symirrs[3] = irrs;
                    if (side == PQ || (assympr && irrp > irrr)) swap(symirrs[0], symirrs[1]);
                    if (side == RS || (assymqs && irrq > irrs)) swap(symirrs[2], symirrs[3]);

                    symtensor.writeRemoteData(symirrs, 1, 1, pairs[irrp][irrq]);
                }
            }
        }
    }

    PROFILE_STOP
}

template <typename T>
void SparseAOMOIntegrals<T>::pqrs_integrals::transpose()
{
    for (auto& idx : idxs)
    {
        swap(idx.i, idx.k);
        swap(idx.j, idx.l);
    }
    swap(np, nr);
    swap(nq, ns);
}

template <typename T>
T absmax(const vector<T>& c)
{
    //return c[iamax(c.size(), c.data(), 1)];
    return nrm2(c.size(), c.data(), 1);
}

template <typename T>
bool SparseAOMOIntegrals<T>::run(TaskDAG& dag, const Arena& arena)
{
    CTF_Timer_epoch ep("SparseAOMOIntegrals");
    ep.begin();
    const MOSpace<T>& occ = this->template get<MOSpace<T> >("occ");
    const MOSpace<T>& vrt = this->template get<MOSpace<T> >("vrt");

    const ERI& ints = this->template get<ERI>("I");

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

    vector<vector<typename real_type<T>::type> >& Ea =
        this->template get<vector<vector<typename real_type<T>::type> > >("Ea");
    vector<vector<typename real_type<T>::type> >& Eb =
        this->template get<vector<vector<typename real_type<T>::type> > >("Eb");

    SymmetryBlockedTensor<T>& Fa = this->template get<SymmetryBlockedTensor<T> >("Fa");
    SymmetryBlockedTensor<T>& Fb = this->template get<SymmetryBlockedTensor<T> >("Fb");

    //this->put("H", new TwoElectronOperator<T>("V", OneElectronOperator<T>("f", arena, occ, vrt)));
    this->put("H", new TwoElectronOperator<T>("V", OneElectronOperator<T>("f", occ, vrt, Fa, Fb)));
    TwoElectronOperator<T>& H = this->template get<TwoElectronOperator<T> >("H");

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

    vector<vector<T> > cA(n), ca(n), cI(n), ci(n);

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
    pqrs_integrals pqrs(N, ints);
    pqrs.collect(true);

    /*
     * First quarter-transformation
     */
    //SHOWIT(PQrs);
    pqrs_integrals pArs = pqrs.transform(B, nA, cA, true);
    //SHOWIT(PArs);
    pqrs_integrals pars = pqrs.transform(B, na, ca, true);
    //SHOWIT(Pars);
    pqrs_integrals pIrs = pqrs.transform(B, nI, cI, true);
    //SHOWIT(PIrs);
    pqrs_integrals pirs = pqrs.transform(B, ni, ci, true);
    //SHOWIT(Pirs);
    pqrs.free();

    /*
     * Second quarter-transformation
     */
    pqrs_integrals ABrs = pArs.transform(A, nA, cA, false);
    //SHOWIT(ABrs);
    pArs.free();
    pqrs_integrals abrs = pars.transform(A, na, ca, false);
    //SHOWIT(abrs);
    pars.free();
    pqrs_integrals AIrs = pIrs.transform(A, nA, cA, false);
    //SHOWIT(AIrs);
    pqrs_integrals IJrs = pIrs.transform(A, nI, cI, false);
    //SHOWIT(IJrs);
    pIrs.free();
    pqrs_integrals airs = pirs.transform(A, na, ca, false);
    //SHOWIT(airs);
    pqrs_integrals ijrs = pirs.transform(A, ni, ci, false);
    //SHOWIT(ijrs);
    pirs.free();

    /*
     * Make <AB||CD>
     */
    pqrs_integrals& rsAB = ABrs;
    rsAB.transpose();
    rsAB.collect(false);

    pqrs_integrals rDAB = rsAB.transform(B, nA, cA, true);
    rsAB.free();

    pqrs_integrals CDAB = rDAB.transform(A, nA, cA, false);
    //SHOWIT(CDAB);
    rDAB.free();
    CDAB.transcribe(H.getABCD()({2,0},{2,0}), true, true, NONE);
    CDAB.free();

    /*
     * Make <Ab|Cd> and <ab||cd>
     */
    pqrs_integrals& rsab = abrs;
    rsab.transpose();
    rsab.collect(false);

    pqrs_integrals rDab = rsab.transform(B, nA, cA, true);
    //SHOWIT(RDab);
    pqrs_integrals rdab = rsab.transform(B, na, ca, true);
    //SHOWIT(Rdab);
    rsab.free();

    pqrs_integrals CDab = rDab.transform(A, nA, cA, false);
    //SHOWIT(CDab);
    rDab.free();
    CDab.transcribe(H.getABCD()({1,0},{1,0}), false, false, NONE);
    CDab.free();

    pqrs_integrals cdab = rdab.transform(A, na, ca, false);
    //SHOWIT(cdab);
    rdab.free();
    cdab.transcribe(H.getABCD()({0,0},{0,0}), true, true, NONE);
    cdab.free();

    /*
     * Make <AB||CI>, <Ab|cI>, and <AB|IJ>
     */
    pqrs_integrals& rsAI = AIrs;
    rsAI.transpose();
    rsAI.collect(false);

    pqrs_integrals rCAI = rsAI.transform(B, nA, cA, true);
    //SHOWIT(RCAI);
    pqrs_integrals rcAI = rsAI.transform(B, na, ca, true);
    //SHOWIT(RcAI);
    pqrs_integrals rJAI = rsAI.transform(B, nI, cI, true);
    //SHOWIT(RJAI);
    rsAI.free();

    pqrs_integrals BCAI = rCAI.transform(A, nA, cA, false);
    //SHOWIT(BCAI);
    rCAI.free();
    BCAI.transcribe(H.getABCI()({2,0},{1,1}), true, false, NONE);
    BCAI.free();

    pqrs_integrals bcAI = rcAI.transform(A, na, ca, false);
    //SHOWIT(bcAI);
    rcAI.free();
    bcAI.transcribe(H.getABCI()({1,0},{0,1}), false, false, PQ);
    bcAI.free();

    pqrs_integrals BJAI = rJAI.transform(A, nA, cA, false);
    //SHOWIT(BJAI);
    rJAI.free();
    BJAI.transcribe(ABIJ__, false, false, NONE);
    BJAI.free();

    /*
     * Make <Ab|Ci>, <ab||ci>, <Ab|Ij>, and <ab|ij>
     */
    pqrs_integrals& rsai = airs;
    rsai.transpose();
    rsai.collect(false);

    pqrs_integrals rCai = rsai.transform(B, nA, cA, true);
    //SHOWIT(RCai);
    pqrs_integrals rcai = rsai.transform(B, na, ca, true);
    //SHOWIT(Rcai);
    pqrs_integrals rJai = rsai.transform(B, nI, cI, true);
    //SHOWIT(RJai);
    pqrs_integrals rjai = rsai.transform(B, ni, ci, true);
    //SHOWIT(Rjai);
    rsai.free();

    pqrs_integrals BCai = rCai.transform(A, nA, cA, false);
    //SHOWIT(BCai);
    rCai.free();
    BCai.transcribe(H.getABCI()({1,0},{1,0}), false, false, NONE);
    BCai.free();

    pqrs_integrals bcai = rcai.transform(A, na, ca, false);
    //SHOWIT(bcai);
    rcai.free();
    bcai.transcribe(H.getABCI()({0,0},{0,0}), true, false, NONE);
    bcai.free();

    pqrs_integrals BJai = rJai.transform(A, nA, cA, false);
    //SHOWIT(BJai);
    rJai.free();
    BJai.transcribe(H.getABIJ()({1,0},{0,1}), false, false, NONE);
    BJai.free();

    pqrs_integrals bjai = rjai.transform(A, na, ca, false);
    //SHOWIT(bjai);
    rjai.free();
    bjai.transcribe(abij__, false, false, NONE);
    bjai.free();

    /*
     * Make <IJ||KL>, <AI||JK>, <aI|Jk>, <aI|bJ>, and <AI|BJ>
     */
    pqrs_integrals& rsIJ = IJrs;
    rsIJ.transpose();
    rsIJ.collect(false);

    pqrs_integrals rBIJ = rsIJ.transform(B, nA, cA, true);
    //SHOWIT(RBIJ);
    pqrs_integrals rbIJ = rsIJ.transform(B, na, ca, true);
    //SHOWIT(RbIJ);
    pqrs_integrals rLIJ = rsIJ.transform(B, nI, cI, true);
    //SHOWIT(RLIJ);
    pqrs_integrals rlIJ = rsIJ.transform(B, ni, ci, true);
    //SHOWIT(RlIJ);
    rsIJ.free();

    pqrs_integrals ABIJ = rBIJ.transform(A, nA, cA, false);
    //SHOWIT(ABIJ);
    rBIJ.free();
    ABIJ.transcribe(H.getAIBJ()({1,1},{1,1}), false, false, NONE);
    ABIJ.free();

    pqrs_integrals abIJ = rbIJ.transform(A, na, ca, false);
    //SHOWIT(abIJ);
    rbIJ.free();
    abIJ.transcribe(H.getAIBJ()({0,1},{0,1}), false, false, NONE);
    abIJ.free();

    pqrs_integrals akIJ = rlIJ.transform(A, na, ca, false);
    //SHOWIT(akIJ);
    rlIJ.free();
    akIJ.transcribe(H.getAIJK()({0,1},{0,1}), false, false, RS);
    akIJ.free();

    pqrs_integrals AKIJ = rLIJ.transform(A, nA, cA, false);
    //SHOWIT(AKIJ);
    pqrs_integrals KLIJ = rLIJ.transform(A, nI, cI, false);
    //SHOWIT(KLIJ);
    rLIJ.free();
    AKIJ.transcribe(H.getAIJK()({1,1},{0,2}), false, true, NONE);
    AKIJ.free();
    KLIJ.transcribe(H.getIJKL()({0,2},{0,2}), true, true, NONE);
    KLIJ.free();

    /*
     * Make <Ij|Kl>, <ij||kl>, <Ai|Jk>, <ai||jk>, <Ai|Bj>, and <ai|bj>
     */
    pqrs_integrals& rsij = ijrs;
    rsij.transpose();
    rsij.collect(false);

    pqrs_integrals rBij = rsij.transform(B, nA, cA, true);
    //SHOWIT(RBij);
    pqrs_integrals rbij = rsij.transform(B, na, ca, true);
    //SHOWIT(Rbij);
    pqrs_integrals rLij = rsij.transform(B, nI, cI, true);
    //SHOWIT(RLij);
    pqrs_integrals rlij = rsij.transform(B, ni, ci, true);
    //SHOWIT(Rlij);
    rsij.free();

    pqrs_integrals ABij = rBij.transform(A, nA, cA, false);
    //SHOWIT(ABij);
    rBij.free();
    ABij.transcribe(H.getAIBJ()({1,0},{1,0}), false, false, NONE);
    ABij.free();

    pqrs_integrals abij = rbij.transform(A, na, ca, false);
    //SHOWIT(abij);
    rbij.free();
    abij.transcribe(H.getAIBJ()({0,0},{0,0}), false, false, NONE);
    abij.free();

    pqrs_integrals AKij = rLij.transform(A, nA, cA, false);
    //SHOWIT(AKij);
    pqrs_integrals KLij = rLij.transform(A, nI, cI, false);
    //SHOWIT(KLij);
    rLij.free();
    AKij.transcribe(H.getAIJK()({1,0},{0,1}), false, false, NONE);
    AKij.free();
    KLij.transcribe(H.getIJKL()({0,1},{0,1}), false, false, NONE);
    KLij.free();

    pqrs_integrals akij = rlij.transform(A, na, ca, false);
    //SHOWIT(akij);
    pqrs_integrals klij = rlij.transform(A, ni, ci, false);
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
