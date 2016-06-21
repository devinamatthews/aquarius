#include "aomoints.hpp"

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
pqrs_integrals<T>::pqrs_integrals(const vector<int>& norb, const ERI& aoints)
: Distributed(aoints.arena), group(aoints.group)
{
    PROFILE_FUNCTION

    ns = nr = nq = np = norb;

    const deque<T>& oldints = aoints.ints;
    const deque<idx4_t>& oldidxs = aoints.idxs;
    size_t noldints = oldints.size();
    size_t nints = noldints;

    for (auto& idx : oldidxs)
    {
        if (!((idx.i == idx.k && idx.j == idx.l) ||
              (idx.i == idx.l && idx.j == idx.k))) nints++;
    }
    ints.reserve(nints);
    idxs.reserve(nints);

    auto iidx = oldidxs.begin();
    auto iint = oldints.begin();
    int j = 0;
    for (;iidx != oldidxs.end();++iidx, ++iint)
    {
        idx4_t idx = *iidx;

        if (idx.i > idx.j) swap(idx.i, idx.j);
        if (idx.k > idx.l) swap(idx.k, idx.l);

        ints.push_back(*iint);
        idxs.push_back(idx);
        j++;

        if (idx.i != idx.k || idx.j != idx.l)
        {
            swap(idx.i, idx.k);
            swap(idx.j, idx.l);
            ints.push_back(*iint);
            idxs.push_back(idx);
            j++;
        }
    }
    assert(j == nints);

    PROFILE_STOP
}

template <typename T>
pqrs_integrals<T>::pqrs_integrals(abrs_integrals<T>& abrs)
: Distributed(abrs.arena), group(abrs.group)
{
    PROFILE_FUNCTION

    np = abrs.nr;
    nq = abrs.ns;
    nr = abrs.na;
    ns = abrs.nb;

    int n = group.getNumIrreps();

    vector<int> irrepp;
    for (int i = 0;i < n;i++) irrepp += vector<int>(np[i],i);
    vector<int> irrepq;
    for (int i = 0;i < n;i++) irrepq += vector<int>(nq[i],i);

    vector<int> startr(n,0);
    for (int i = 1;i < n;i++) startr[i] = startr[i-1]+nr[i-1];
    vector<int> starts(n,0);
    for (int i = 1;i < n;i++) starts[i] = starts[i-1]+ns[i-1];

    Representation irrpq(group), irrpqs(group), irrpqrs(group);

    ints.swap(abrs.ints);
    idxs.resize(ints.size());
    size_t ipqrs, irs;
    for (ipqrs = 0, irs = 0;irs < abrs.rs.size();irs++)
    {
        int p = abrs.rs[irs].i;
        int q = abrs.rs[irs].j;
        int irrp = irrepp[p];
        int irrq = irrepq[q];
        irrpq = group.getIrrep(irrp)*group.getIrrep(irrq);

        for (int irrs = 0;irrs < n;irrs++)
        {
            irrpqs = irrpq;
            irrpqs *= group.getIrrep(irrs);

            for (int irrr = 0;irrr < n;irrr++)
            {
                irrpqrs = irrpqs;
                irrpqrs *= group.getIrrep(irrr);
                if (!irrpqrs.isTotallySymmetric()) continue;

                for (int s = starts[irrs];s < starts[irrs]+ns[irrs];s++)
                {
                    for (int r = startr[irrr];r < startr[irrr]+nr[irrr];r++)
                    {
                        idxs[ipqrs].i = p;
                        idxs[ipqrs].j = q;
                        idxs[ipqrs].k = r;
                        idxs[ipqrs].l = s;
                        ipqrs++;
                    }
                }
            }
        }
    }
    assert(ipqrs == ints.size());

    abrs.free();

    PROFILE_STOP
}

template <typename T>
void pqrs_integrals<T>::free()
{
    ints.clear();
    idxs.clear();
}

template <typename T>
void pqrs_integrals<T>::sortInts(bool rles, size_t& nrs, vector<size_t>& rscount)
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

    for (size_t i = 0;i < ints.size();i++)
    {
        size_t rs = (rles ? idxs[i].k+idxs[i].l*(idxs[i].l+1)/2
                          : idxs[i].k+idxs[i].l*nrtot);

        newints[rsoff[rs]] = ints[i];
        newidxs[rsoff[rs]] = idxs[i];

        rsoff[rs]++;
    }

    swap(ints, newints);
    swap(idxs, newidxs);
}


template <typename T>
void pqrs_integrals<T>::transcribe(SymmetryBlockedTensor<T>& symtensor, bool assympr, bool assymqs, Side side)
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
pqrs_integrals<T> pqrs_integrals<T>::transform(Index index, const vector<int>& nc,
                                               const vector<vector<T>>& C, bool pleq)
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
                                    if (aquarius::abs(val) > 1e-12)
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
                                    if (aquarius::abs(val) > 1e-12)
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
void pqrs_integrals<T>::transpose()
{
    for (auto& idx : idxs)
    {
        swap(idx.i, idx.k);
        swap(idx.j, idx.l);
    }
    swap(np, nr);
    swap(nq, ns);
}

/*
 * Redistribute integrals such that each node has all pq for each rs pair
 */
template <typename T>
void pqrs_integrals<T>::collect(bool rles)
{
    PROFILE_FUNCTION

    Datatype IDX4_T_TYPE = MPI_TYPE_<uint16_t>::value()*4;

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
    this->arena.comm().Alltoall(const_cast<const vector<MPI_Int>&>(sendcount), recvcount);
    PROFILE_STOP

    size_t nnewints = sum(recvcount);
    vector<T> newints(nnewints);
    vector<idx4_t> newidxs(nnewints);

    PROFILE_SECTION(collect_comm)
    this->arena.comm().Alltoall(ints, sendcount, newints, recvcount);
    this->arena.comm().Alltoall(idxs, sendcount, newidxs, recvcount, IDX4_T_TYPE);
    PROFILE_STOP

    swap(ints, newints);
    swap(idxs, newidxs);
    newints.clear();
    newidxs.clear();

    sortInts(rles, nrs, rscount);

    PROFILE_STOP
}

template <typename T>
abrs_integrals<T>::abrs_integrals(pqrs_integrals<T>& pqrs, const bool pleq)
: Distributed(pqrs.arena), group(pqrs.group)
{
    PROFILE_FUNCTION

    int n = group.getNumIrreps();
    Representation irrrs(group), irrbrs(group), irrabrs(group);

    na = pqrs.np;
    nb = pqrs.nq;
    nr = pqrs.nr;
    ns = pqrs.ns;

    vector<int> irrepa;
    for (int i = 0;i < n;i++) irrepa += vector<int>(na[i],i);
    vector<int> irrepb;
    for (int i = 0;i < n;i++) irrepb += vector<int>(nb[i],i);

    vector<int> starta(n,0);
    for (int i = 1;i < n;i++) starta[i] = starta[i-1]+na[i-1];
    vector<int> startb(n,0);
    for (int i = 1;i < n;i++) startb[i] = startb[i-1]+nb[i-1];

    size_t nrs = 0;
    size_t nints = 0;
    for (size_t ipqrs = 0;ipqrs < pqrs.ints.size();ipqrs++)
    {
        assert( ipqrs == 0 ||
                pqrs.idxs[ipqrs].l  > pqrs.idxs[ipqrs-1].l ||
               (pqrs.idxs[ipqrs].l == pqrs.idxs[ipqrs-1].l &&
                pqrs.idxs[ipqrs].k >= pqrs.idxs[ipqrs-1].k));
        if (ipqrs == 0 ||
            pqrs.idxs[ipqrs].k != pqrs.idxs[ipqrs-1].k ||
            pqrs.idxs[ipqrs].l != pqrs.idxs[ipqrs-1].l)
        {
            nints += getNumAB(idx2_t(pqrs.idxs[ipqrs].k, pqrs.idxs[ipqrs].l));
            nrs++;
        }
    }

    rs.resize(nrs);
    ints.assign(nints, (T)0);

    size_t iabrs = 0;
    size_t irs = -1;
    size_t nab = 0;
    vector<size_t> offab(n*n);
    for (size_t ipqrs = 0;ipqrs < pqrs.ints.size();ipqrs++)
    {
        if (ipqrs == 0 ||
            pqrs.idxs[ipqrs].k != pqrs.idxs[ipqrs-1].k ||
            pqrs.idxs[ipqrs].l != pqrs.idxs[ipqrs-1].l)
        {
            irs++;
            iabrs += nab;
            assert(irs >= 0 && irs < nrs);
            rs[irs].i = pqrs.idxs[ipqrs].k;
            rs[irs].j = pqrs.idxs[ipqrs].l;
            fill(offab.begin(), offab.end(), SIZE_MAX);
            nab = getNumAB(rs[irs], offab);
        }

        int irra = irrepa[pqrs.idxs[ipqrs].i];
        int irrb = irrepb[pqrs.idxs[ipqrs].j];
        int a = pqrs.idxs[ipqrs].i-starta[irra];
        int b = pqrs.idxs[ipqrs].j-startb[irrb];

        assert(offab[irra+irrb*n] != SIZE_MAX);
        assert(a >= 0 && a < na[irra]);
        assert(b >= 0 && b < nb[irrb]);
                  ints[iabrs+offab[irra+irrb*n]+a+b*na[irra]] = pqrs.ints[ipqrs];
        if (pleq) ints[iabrs+offab[irrb+irra*n]+b+a*nb[irrb]] = pqrs.ints[ipqrs];
    }
    assert(irs+1 == nrs);
    assert(iabrs+nab == ints.size());

    pqrs.free();

    PROFILE_STOP
}

template <typename T>
abrs_integrals<T> abrs_integrals<T>::transform(Index index, const vector<int>& nc, const vector<vector<T>>& C)
{
    abrs_integrals out(arena, group);

    out.nr = nr;
    out.ns = ns;
    out.rs = rs;

    int n = group.getNumIrreps();

    if (index == A)
    {
        out.na = nc;
        out.nb = nb;

        size_t nints = 0;
        vector<size_t> offabrs(rs.size());
        vector<size_t> offcbrs(rs.size());
        for (size_t irs = 0;irs < rs.size();irs++)
        {
            size_t nab =     getNumAB(rs[irs]);
            size_t ncb = out.getNumAB(rs[irs]);
            if (irs < rs.size()-1)
            {
                offabrs[irs+1] = offabrs[irs]+nab;
                offcbrs[irs+1] = offcbrs[irs]+ncb;
            }
            nints += ncb;
        }
        out.ints.clear();
        out.ints.resize(nints, 0);

        PROFILE_SECTION(transform)
        long_int flops = 0;
        #pragma omp parallel
        {
            vector<size_t> offab(n*n);
            vector<size_t> offcb(n*n);

            #pragma omp for schedule(static), reduction(+:flops)
            for (size_t irs = 0;irs < rs.size();irs++)
            {
                fill(offab.begin(), offab.end(), SIZE_MAX);
                fill(offcb.begin(), offcb.end(), SIZE_MAX);
                    getNumAB(rs[irs], offab);
                out.getNumAB(rs[irs], offcb);

                for (int irrb = 0;irrb < n;irrb++)
                {
                    for (int irra = 0;irra < n;irra++)
                    {
                        if (offab[irra+irrb*n] == SIZE_MAX || offcb[irra+irrb*n] == SIZE_MAX) continue;
                        if (nc[irra] == 0 || nb[irrb] == 0 || na[irra] == 0) continue;

                        assert(offabrs[irs]+offab[irra+irrb*n] >= 0);
                        assert(offabrs[irs]+offab[irra+irrb*n] < ints.size());
                        assert(offabrs[irs]+offab[irra+irrb*n]+na[irra]*nb[irrb] > 0);
                        assert(offabrs[irs]+offab[irra+irrb*n]+na[irra]*nb[irrb] <= ints.size());
                        assert(offcbrs[irs]+offcb[irra+irrb*n] >= 0);
                        assert(offcbrs[irs]+offcb[irra+irrb*n] < out.ints.size());
                        assert(offcbrs[irs]+offcb[irra+irrb*n]+nc[irra]*nb[irrb] > 0);
                        assert(offcbrs[irs]+offcb[irra+irrb*n]+nc[irra]*nb[irrb] <= out.ints.size());
                        gemm('T', 'N', nc[irra], nb[irrb], na[irra],
                             1.0,                                  C[irra].data(), na[irra],
                                      ints.data()+offabrs[irs]+offab[irra+irrb*n], na[irra],
                             1.0, out.ints.data()+offcbrs[irs]+offcb[irra+irrb*n], nc[irra]);
                        flops += 2*nc[irra]*nb[irrb]*na[irra]+nc[irra]*nb[irrb];
                    }
                }
            }
        }
        //PROFILE_FLOPS(flops);
        PROFILE_STOP
    }
    else
    {
        out.na = na;
        out.nb = nc;

        size_t nints = 0;
        vector<size_t> offabrs(rs.size(), 0);
        vector<size_t> offacrs(rs.size(), 0);
        for (size_t irs = 0;irs < rs.size();irs++)
        {
            size_t nab =     getNumAB(rs[irs]);
            size_t nac = out.getNumAB(rs[irs]);
            if (irs < rs.size()-1)
            {
                offabrs[irs+1] = offabrs[irs]+nab;
                offacrs[irs+1] = offacrs[irs]+nac;
            }
            nints += nac;
        }
        out.ints.clear();
        out.ints.resize(nints, 0);

        PROFILE_SECTION(transform)
        long_int flops = 0;
        #pragma omp parallel
        {
            vector<size_t> offab(n*n);
            vector<size_t> offac(n*n);

            #pragma omp for schedule(static)
            for (size_t irs = 0;irs < rs.size();irs++)
            {
                fill(offab.begin(), offab.end(), SIZE_MAX);
                fill(offac.begin(), offac.end(), SIZE_MAX);
                    getNumAB(rs[irs], offab);
                out.getNumAB(rs[irs], offac);

                for (int irrb = 0;irrb < n;irrb++)
                {
                    for (int irra = 0;irra < n;irra++)
                    {
                        if (offab[irra+irrb*n] == SIZE_MAX || offac[irra+irrb*n] == SIZE_MAX) continue;
                        if (nc[irrb] == 0 || nb[irrb] == 0 || na[irra] == 0) continue;

                        assert(offabrs[irs]+offab[irra+irrb*n] >= 0);
                        assert(offabrs[irs]+offab[irra+irrb*n] < ints.size());
                        assert(offabrs[irs]+offab[irra+irrb*n]+na[irra]*nb[irrb] > 0);
                        assert(offabrs[irs]+offab[irra+irrb*n]+na[irra]*nb[irrb] <= ints.size());
                        assert(offacrs[irs]+offac[irra+irrb*n] >= 0);
                        assert(offacrs[irs]+offac[irra+irrb*n] < out.ints.size());
                        assert(offacrs[irs]+offac[irra+irrb*n]+na[irra]*nc[irrb] > 0);
                        assert(offacrs[irs]+offac[irra+irrb*n]+na[irra]*nc[irrb] <= out.ints.size());
                        gemm('N', 'N', na[irra], nc[irrb], nb[irrb],
                             1.0,     ints.data()+offabrs[irs]+offab[irra+irrb*n], na[irra],
                                                                   C[irrb].data(), nb[irrb],
                             1.0, out.ints.data()+offacrs[irs]+offac[irra+irrb*n], na[irra]);
                        flops += 2*nc[irra]*nb[irrb]*na[irra]+nc[irra]*na[irra];
                    }
                }
            }
        }
        //PROFILE_FLOPS(flops);
        PROFILE_STOP
    }

    return out;
}

template <typename T>
void abrs_integrals<T>::transcribe(SymmetryBlockedTensor<T>& symtensor, bool assymij, bool assymkl, Side side)
{
    PROFILE_FUNCTION

    if (assymij) assert(na == nr);
    if (assymkl) assert(nb == ns);

    int n = group.getNumIrreps();

    vector<int> irrepr;
    for (int i = 0;i < n;i++) irrepr += vector<int>(nr[i],i);
    vector<int> irreps;
    for (int i = 0;i < n;i++) irreps += vector<int>(ns[i],i);

    vector<int> startr(n);
    for (int i = 1;i < n;i++) startr[i] = startr[i-1]+nr[i-1];
    vector<int> starts(n);
    for (int i = 1;i < n;i++) starts[i] = starts[i-1]+ns[i-1];

    vector<int> starta(n);
    for (int i = 1;i < n;i++) starta[i] = starta[i-1]+na[i-1];
    vector<int> startb(n);
    for (int i = 1;i < n;i++) startb[i] = startb[i-1]+nb[i-1];

    vector<size_t> offab(n*n);
    vector<int> symirrs(4);
    Representation irrrs(group), irrbrs(group), irrabrs(group);

    int nrtot = sum(nr);
    int nstot = sum(ns);

    vector<size_t> offrs(nrtot*nstot, SIZE_MAX);
    {
        size_t off = 0, i = 0, j = 0;
        for (int s = 0;s < nstot;s++)
        {
            for (int r = 0;r < nrtot;r++)
            {
                if (j < rs.size() && rs[j].i == r && rs[j].j == s)
                {
                    offrs[i] = off;
                    off += getNumAB(rs[j]);
                    j++;
                }
                i++;
            }
        }
        assert(j == rs.size());
    }

    /*
     * (ab|rs) -> <ar|bs>
     */
    for (int irrs = 0;irrs < n;irrs++)
    {
        for (int irrr = 0;irrr < n;irrr++)
        {
            if (nr[irrr] == 0 || ns[irrs] == 0) continue;

            irrrs = group.getIrrep(irrr)*group.getIrrep(irrs);
            fill(offab.begin(), offab.end(), SIZE_MAX);
            getNumAB(idx2_t(startr[irrr],starts[irrs]), offab);

            for (int irrb = 0;irrb < n;irrb++)
            {
                irrbrs = irrrs;
                irrbrs *= group.getIrrep(irrb);

                for (int irra = 0;irra < n;irra++)
                {
                    irrabrs = irrbrs;
                    irrabrs *= group.getIrrep(irra);

                    if (!irrabrs.isTotallySymmetric()) continue;
                    if (assymij && assymkl && irrb > irrs) continue;

                    symirrs[0] = irra;
                    symirrs[1] = irrr;
                    symirrs[2] = irrb;
                    symirrs[3] = irrs;
                    if (side == PQ || (assymij && irra > irrr)) swap(symirrs[0], symirrs[1]);
                    if (side == RS || (assymkl && irrb > irrs)) swap(symirrs[2], symirrs[3]);

                    size_t npair = 0;
                    for (int s = 0;s < ns[irrs];s++)
                    {
                        for (int r = 0;r < nr[irrr];r++)
                        {
                            if (offrs[(r+startr[irrr])+(s+starts[irrs])*nrtot] == SIZE_MAX) continue;

                            for (int b = 0;b < nb[irrb];b++)
                            {
                                for (int a = 0;a < na[irra];a++)
                                {
                                    int i = a; int i0 = starta[irra];
                                    int j = r; int j0 = startr[irrr];
                                    int k = b; int k0 = startb[irrb];
                                    int l = s; int l0 = starts[irrs];

                                    if ((assymij && i+i0 == j+j0) ||
                                        (assymkl && k+k0 == l+l0)) continue;

                                    if (assymij && assymkl && k+k0 > l+l0) continue;

                                    npair++;
                                }
                            }
                        }
                    }

                    vector<kv_pair> pairs(npair);

                    size_t pair = 0;
                    for (int s = 0;s < ns[irrs];s++)
                    {
                        for (int r = 0;r < nr[irrr];r++)
                        {
                            if (offrs[(r+startr[irrr])+(s+starts[irrs])*nrtot] == SIZE_MAX) continue;

                            size_t idx = offrs[(r+startr[irrr])+(s+starts[irrs])*nrtot]+offab[irra+irrb*n]-1;

                            for (int b = 0;b < nb[irrb];b++)
                            {
                                for (int a = 0;a < na[irra];a++)
                                {
                                    idx++;
                                    assert(idx < ints.size());

                                    T val = ints[idx];

                                    int i = a; int i0 = starta[irra]; int ni = na[irra];
                                    int j = r; int j0 = startr[irrr]; int nj = nr[irrr];
                                    int k = b; int k0 = startb[irrb]; int nk = nb[irrb];
                                    int l = s; int l0 = starts[irrs]; int nl = ns[irrs];

                                    if ((assymij && i+i0 == j+j0) ||
                                        (assymkl && k+k0 == l+l0)) continue;

                                    if (side == PQ || (assymij && i+i0 > j+j0))
                                    {
                                        val = -val;
                                        swap(i,j);
                                        swap(ni,nj);
                                    }
                                    if (side == RS || (assymkl && k+k0 > l+l0))
                                    {
                                        if (side != RS && assymij) continue;
                                        val = -val;
                                        swap(k,l);
                                        swap(nk,nl);
                                    }

                                    pairs[pair].k = ((((int64_t)l)*nk+k)*nj+j)*ni+i;
                                    pairs[pair].d = val;

                                    pair++;
                                }
                            }
                        }
                    }
                    assert(pair == pairs.size());

                    symtensor.writeRemoteData(symirrs, 1, 1, pairs);
                }
            }
        }
    }

    PROFILE_STOP
}

template <typename T>
void abrs_integrals<T>::free()
{
    ints.clear();
    rs.clear();
}

template <typename T>
size_t abrs_integrals<T>::getNumAB(idx2_t rs)
{
    int n = group.getNumIrreps();

    Representation irrrs(group), irrbrs(group), irrabrs(group);

    vector<int> irrepr;
    for (int i = 0;i < n;i++) irrepr += vector<int>(nr[i],i);
    vector<int> irreps;
    for (int i = 0;i < n;i++) irreps += vector<int>(ns[i],i);

    irrrs = group.getIrrep(irrepr[rs.i])*
            group.getIrrep(irreps[rs.j]);

    size_t nab = 0;
    for (int irrb = 0;irrb < n;irrb++)
    {
        irrbrs = irrrs;
        irrbrs *= group.getIrrep(irrb);

        for (int irra = 0;irra < n;irra++)
        {
            irrabrs = irrbrs;
            irrabrs *= group.getIrrep(irra);

            if (irrabrs.isTotallySymmetric()) nab += na[irra]*nb[irrb];
        }
    }

    return nab;
}

template <typename T>
size_t abrs_integrals<T>::getNumAB(idx2_t rs, vector<size_t>& offab)
{
    int n = group.getNumIrreps();

    Representation irrrs(group), irrbrs(group), irrabrs(group);

    vector<int> irrepr;
    for (int i = 0;i < n;i++) irrepr += vector<int>(nr[i],i);
    vector<int> irreps;
    for (int i = 0;i < n;i++) irreps += vector<int>(ns[i],i);

    irrrs = group.getIrrep(irrepr[rs.i])*
            group.getIrrep(irreps[rs.j]);

    size_t nab = 0;
    for (int irrb = 0;irrb < n;irrb++)
    {
        irrbrs = irrrs;
        irrbrs *= group.getIrrep(irrb);

        for (int irra = 0;irra < n;irra++)
        {
            irrabrs = irrbrs;
            irrabrs *= group.getIrrep(irra);

            if (irrabrs.isTotallySymmetric())
            {
                offab[irra+irrb*n] = nab;
                nab += na[irra]*nb[irrb];
            }
        }
    }

    return nab;
}

template <typename T>
T absmax(const vector<T>& c)
{
    //return c[iamax(c.size(), c.data(), 1)];
    return nrm2(c.size(), c.data(), 1);
}

template <typename T>
AOMOIntegrals<T>::AOMOIntegrals(const string& name, Config& config)
: MOIntegrals<T>(name, config)
{
    this->getProduct("H").addRequirement("eri", "I");
}

template <typename T>
bool AOMOIntegrals<T>::run(TaskDAG& dag, const Arena& arena)
{
    CTF_Timer_epoch ep("AOMOIntegrals");
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
     * are dense blocks for each sparse rs pair
     */
    pqrs_integrals<T> pqrs(N, ints);
    pqrs.collect(true);
    abrs_integrals<T> PQrs(pqrs, true);

    /*
     * First quarter-transformation
     */
    //SHOWIT(PQrs);
    abrs_integrals<T> PArs = PQrs.transform(B, nA, cA);
    //SHOWIT(PArs);
    abrs_integrals<T> Pars = PQrs.transform(B, na, ca);
    //SHOWIT(Pars);
    abrs_integrals<T> PIrs = PQrs.transform(B, nI, cI);
    //SHOWIT(PIrs);
    abrs_integrals<T> Pirs = PQrs.transform(B, ni, ci);
    //SHOWIT(Pirs);
    PQrs.free();

    /*
     * Second quarter-transformation
     */
    abrs_integrals<T> ABrs = PArs.transform(A, nA, cA);
    //SHOWIT(ABrs);
    PArs.free();
    abrs_integrals<T> abrs = Pars.transform(A, na, ca);
    //SHOWIT(abrs);
    Pars.free();
    abrs_integrals<T> AIrs = PIrs.transform(A, nA, cA);
    //SHOWIT(AIrs);
    abrs_integrals<T> IJrs = PIrs.transform(A, nI, cI);
    //SHOWIT(IJrs);
    PIrs.free();
    abrs_integrals<T> airs = Pirs.transform(A, na, ca);
    //SHOWIT(airs);
    abrs_integrals<T> ijrs = Pirs.transform(A, ni, ci);
    //SHOWIT(ijrs);
    Pirs.free();

    /*
     * Make <AB||CD>
     */
    pqrs_integrals<T> rsAB(ABrs);
    rsAB.collect(false);

    abrs_integrals<T> RSAB(rsAB, true);
    //SHOWIT(RSAB)<T>;
    abrs_integrals<T> RDAB = RSAB.transform(B, nA, cA);
    //SHOWIT(RDAB);
    RSAB.free();

    abrs_integrals<T> CDAB = RDAB.transform(A, nA, cA);
    //SHOWIT(CDAB);
    RDAB.free();
    CDAB.transcribe(H.getABCD()({2,0},{2,0}), true, true, NONE);
    CDAB.free();

    /*
     * Make <Ab|Cd> and <ab||cd>
     */
    pqrs_integrals<T> rsab(abrs);
    rsab.collect(false);

    abrs_integrals<T> RSab(rsab, true);
    //SHOWIT(RSab);
    abrs_integrals<T> RDab = RSab.transform(B, nA, cA);
    //SHOWIT(RDab);
    abrs_integrals<T> Rdab = RSab.transform(B, na, ca);
    //SHOWIT(Rdab);
    RSab.free();

    abrs_integrals<T> CDab = RDab.transform(A, nA, cA);
    //SHOWIT(CDab);
    RDab.free();
    CDab.transcribe(H.getABCD()({1,0},{1,0}), false, false, NONE);
    CDab.free();

    abrs_integrals<T> cdab = Rdab.transform(A, na, ca);
    //SHOWIT(cdab);
    Rdab.free();
    cdab.transcribe(H.getABCD()({0,0},{0,0}), true, true, NONE);
    cdab.free();

    /*
     * Make <AB||CI>, <Ab|cI>, and <AB|IJ>
     */
    pqrs_integrals<T> rsAI(AIrs);
    rsAI.collect(false);

    abrs_integrals<T> RSAI(rsAI, true);
    //SHOWIT(RSAI);
    abrs_integrals<T> RCAI = RSAI.transform(B, nA, cA);
    //SHOWIT(RCAI);
    abrs_integrals<T> RcAI = RSAI.transform(B, na, ca);
    //SHOWIT(RcAI);
    abrs_integrals<T> RJAI = RSAI.transform(B, nI, cI);
    //SHOWIT(RJAI);
    RSAI.free();

    abrs_integrals<T> BCAI = RCAI.transform(A, nA, cA);
    //SHOWIT(BCAI);
    RCAI.free();
    BCAI.transcribe(H.getABCI()({2,0},{1,1}), true, false, NONE);
    BCAI.free();

    abrs_integrals<T> bcAI = RcAI.transform(A, na, ca);
    //SHOWIT(bcAI);
    RcAI.free();
    bcAI.transcribe(H.getABCI()({1,0},{0,1}), false, false, PQ);
    bcAI.free();

    abrs_integrals<T> BJAI = RJAI.transform(A, nA, cA);
    //SHOWIT(BJAI);
    RJAI.free();
    BJAI.transcribe(ABIJ__, false, false, NONE);
    BJAI.free();

    /*
     * Make <Ab|Ci>, <ab||ci>, <Ab|Ij>, and <ab|ij>
     */
    pqrs_integrals<T> rsai(airs);
    rsai.collect(false);

    abrs_integrals<T> RSai(rsai, true);
    //SHOWIT(RSai);
    abrs_integrals<T> RCai = RSai.transform(B, nA, cA);
    //SHOWIT(RCai);
    abrs_integrals<T> Rcai = RSai.transform(B, na, ca);
    //SHOWIT(Rcai);
    abrs_integrals<T> RJai = RSai.transform(B, nI, cI);
    //SHOWIT(RJai);
    abrs_integrals<T> Rjai = RSai.transform(B, ni, ci);
    //SHOWIT(Rjai);
    RSai.free();

    abrs_integrals<T> BCai = RCai.transform(A, nA, cA);
    //SHOWIT(BCai);
    RCai.free();
    BCai.transcribe(H.getABCI()({1,0},{1,0}), false, false, NONE);
    BCai.free();

    abrs_integrals<T> bcai = Rcai.transform(A, na, ca);
    //SHOWIT(bcai);
    Rcai.free();
    bcai.transcribe(H.getABCI()({0,0},{0,0}), true, false, NONE);
    bcai.free();

    abrs_integrals<T> BJai = RJai.transform(A, nA, cA);
    //SHOWIT(BJai);
    RJai.free();
    BJai.transcribe(H.getABIJ()({1,0},{0,1}), false, false, NONE);
    BJai.free();

    abrs_integrals<T> bjai = Rjai.transform(A, na, ca);
    //SHOWIT(bjai);
    Rjai.free();
    bjai.transcribe(abij__, false, false, NONE);
    bjai.free();

    /*
     * Make <IJ||KL>, <AI||JK>, <aI|Jk>, <aI|bJ>, and <AI|BJ>
     */
    pqrs_integrals<T> rsIJ(IJrs);
    rsIJ.collect(false);

    abrs_integrals<T> RSIJ(rsIJ, true);
    //SHOWIT(RSIJ);
    abrs_integrals<T> RBIJ = RSIJ.transform(B, nA, cA);
    //SHOWIT(RBIJ);
    abrs_integrals<T> RbIJ = RSIJ.transform(B, na, ca);
    //SHOWIT(RbIJ);
    abrs_integrals<T> RLIJ = RSIJ.transform(B, nI, cI);
    //SHOWIT(RLIJ);
    abrs_integrals<T> RlIJ = RSIJ.transform(B, ni, ci);
    //SHOWIT(RlIJ);
    RSIJ.free();

    abrs_integrals<T> ABIJ = RBIJ.transform(A, nA, cA);
    //SHOWIT(ABIJ);
    RBIJ.free();
    ABIJ.transcribe(H.getAIBJ()({1,1},{1,1}), false, false, NONE);
    ABIJ.free();

    abrs_integrals<T> abIJ = RbIJ.transform(A, na, ca);
    //SHOWIT(abIJ);
    RbIJ.free();
    abIJ.transcribe(H.getAIBJ()({0,1},{0,1}), false, false, NONE);
    abIJ.free();

    abrs_integrals<T> akIJ = RlIJ.transform(A, na, ca);
    //SHOWIT(akIJ);
    RlIJ.free();
    akIJ.transcribe(H.getAIJK()({0,1},{0,1}), false, false, RS);
    akIJ.free();

    abrs_integrals<T> AKIJ = RLIJ.transform(A, nA, cA);
    //SHOWIT(AKIJ);
    abrs_integrals<T> KLIJ = RLIJ.transform(A, nI, cI);
    //SHOWIT(KLIJ);
    RLIJ.free();
    AKIJ.transcribe(H.getAIJK()({1,1},{0,2}), false, true, NONE);
    AKIJ.free();
    KLIJ.transcribe(H.getIJKL()({0,2},{0,2}), true, true, NONE);
    KLIJ.free();

    /*
     * Make <Ij|Kl>, <ij||kl>, <Ai|Jk>, <ai||jk>, <Ai|Bj>, and <ai|bj>
     */
    pqrs_integrals<T> rsij(ijrs);
    rsij.collect(false);

    abrs_integrals<T> RSij(rsij, true);
    //SHOWIT(RSij);
    abrs_integrals<T> RBij = RSij.transform(B, nA, cA);
    //SHOWIT(RBij);
    abrs_integrals<T> Rbij = RSij.transform(B, na, ca);
    //SHOWIT(Rbij);
    abrs_integrals<T> RLij = RSij.transform(B, nI, cI);
    //SHOWIT(RLij);
    abrs_integrals<T> Rlij = RSij.transform(B, ni, ci);
    //SHOWIT(Rlij);
    RSij.free();

    abrs_integrals<T> ABij = RBij.transform(A, nA, cA);
    //SHOWIT(ABij);
    RBij.free();
    ABij.transcribe(H.getAIBJ()({1,0},{1,0}), false, false, NONE);
    ABij.free();

    abrs_integrals<T> abij = Rbij.transform(A, na, ca);
    //SHOWIT(abij);
    Rbij.free();
    abij.transcribe(H.getAIBJ()({0,0},{0,0}), false, false, NONE);
    abij.free();

    abrs_integrals<T> AKij = RLij.transform(A, nA, cA);
    //SHOWIT(AKij);
    abrs_integrals<T> KLij = RLij.transform(A, nI, cI);
    //SHOWIT(KLij);
    RLij.free();
    AKij.transcribe(H.getAIJK()({1,0},{0,1}), false, false, NONE);
    AKij.free();
    KLij.transcribe(H.getIJKL()({0,1},{0,1}), false, false, NONE);
    KLij.free();

    abrs_integrals<T> akij = Rlij.transform(A, na, ca);
    //SHOWIT(akij);
    abrs_integrals<T> klij = Rlij.transform(A, ni, ci);
    //SHOWIT(klij);
    Rlij.free();
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

    SpinorbitalTensor<double> D("D", arena, PointGroup::C1(), {vrt,occ}, {0,1}, {0,1});

    if (arena.rank == 0)
    {
        vector<kv_pair> pairsa, pairsb;
        for (int i = 0;i < 5;i++) pairsa.emplace_back(i+i*5, 1.0);
        for (int i = 0;i < 5;i++) pairsb.emplace_back(i+i*5, 1.0);
        D({0,1},{0,1})({0,0}).writeRemoteData(pairsa);
        D({0,0},{0,0})({0,0}).writeRemoteData(pairsb);
    }
    else
    {
        D({0,1},{0,1})({0,0}).writeRemoteData();
        D({0,0},{0,0})({0,0}).writeRemoteData();
    }

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
    //this->log(arena) << "AB:   " << setprecision(15) << H.getAB()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "ab:   " << setprecision(15) << H.getAB()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "AI:   " << setprecision(15) << H.getAI()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "ai:   " << setprecision(15) << H.getAI()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IA:   " << setprecision(15) << H.getIA()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "ia:   " << setprecision(15) << H.getIA()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IJ:   " << setprecision(15) << H.getIJ()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "ij:   " << setprecision(15) << H.getIJ()({0,0},{0,0}).norm(2) << endl;
    ep.end();

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::op::pqrs_integrals);
INSTANTIATE_SPECIALIZATIONS(aquarius::op::abrs_integrals);
INSTANTIATE_SPECIALIZATIONS(aquarius::op::AOMOIntegrals);
REGISTER_TASK(aquarius::op::AOMOIntegrals<double>,"aomoints");
