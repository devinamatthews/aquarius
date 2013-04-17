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

#ifndef _AQUARIUS_SCF_AOMOINTS_HPP_
#define _AQUARIUS_SCF_AOMOINTS_HPP_

#include "aoints.hpp"
#include "aoscf.hpp"
#include "moints.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class AOMOIntegrals : public MOIntegrals<T>
{
    public:
        AOMOIntegrals(const AOUHF<T>& uhf)
        : MOIntegrals<T>(uhf)
        {
            doTransformation(uhf.ints);
        }

    private:
        using Distributed<T>::nproc;
        using Distributed<T>::rank;
        enum Side {PQ, RS};
        enum Index {A, B};

        static bool sortIntsByRS(const idx4_t& i1, const idx4_t& i2)
        {
            if (i1.l < i2.l)
            {
                return true;
            }
            else if (i1.l > i2.l)
            {
                return false;
            }
            if (i1.k < i2.k)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        struct abrs_integrals;

        struct pqrs_integrals : Distributed<T>
        {
            int np, nq, nr, ns;
            size_t nints;
            T *ints;
            idx4_t *idxs;

            /*
             * Read integrals in and break (pq|rs)=(rs|pq) symmetry
             */
            pqrs_integrals(const AOIntegrals<T>& aoints)
            : Distributed<T>(aoints.ctf)
            {
                ns = nr = nq = np = aoints.molecule.getNumOrbitals();

                size_t noldints = aoints.getNumInts();
                const T *oldints = aoints.getInts();
                const idx4_t *oldidxs = aoints.getIndices();
                nints = noldints;

                for (size_t i = 0;i < noldints;i++)
                {
                    idx4_t idx = oldidxs[i];

                    if (!((idx.i == idx.k && idx.j == idx.l) ||
                          (idx.i == idx.l && idx.j == idx.k)))
                    {
                        nints++;
                    }
                }

                ints = SAFE_MALLOC(T, nints);
                idxs = SAFE_MALLOC(idx4_t, nints);

                size_t j = 0;
                for (size_t i = 0;i < noldints;i++)
                {
                    T val = oldints[i];
                    idx4_t idx = oldidxs[i];

                    if (idx.i > idx.j) std::swap(idx.i, idx.j);
                    if (idx.k > idx.l) std::swap(idx.k, idx.l);

                    assert(idx.i >= 0 && idx.i < np);
                    assert(idx.j >= 0 && idx.j < nq);
                    assert(idx.k >= 0 && idx.k < nr);
                    assert(idx.l >= 0 && idx.l < ns);

                    assert(j < nints);
                    ints[j] = val;
                    idxs[j] = idx;
                    j++;

                    if (idx.i != idx.k || idx.j != idx.l)
                    {
                        std::swap(idx.i, idx.k);
                        std::swap(idx.j, idx.l);
                        assert(j < nints);
                        ints[j] = val;
                        idxs[j] = idx;
                        j++;
                    }
                }
                assert(j == nints);
            }

            pqrs_integrals(abrs_integrals& abrs)
            : Distributed<T>(abrs.ctf)
            {
                np = abrs.nr;
                nq = abrs.ns;
                nr = abrs.na;
                ns = abrs.nb;

                nints = abrs.nints;
                ints = abrs.ints;
                idxs = SAFE_MALLOC(idx4_t, nints);
                size_t ipqrs, irs;
                for (ipqrs = 0, irs = 0;irs < abrs.nrs;irs++)
                {
                    int p = abrs.rs[irs].i;
                    int q = abrs.rs[irs].j;
                    assert(p >= 0 && p < np);
                    assert(q >= 0 && q < nq);

                    for (int s = 0;s < ns;s++)
                    {
                        for (int r = 0;r < nr;r++)
                        {
                            assert(ipqrs < nints);
                            idxs[ipqrs].i = p;
                            idxs[ipqrs].j = q;
                            idxs[ipqrs].k = r;
                            idxs[ipqrs].l = s;
                            ipqrs++;
                        }
                    }
                }
                assert(ipqrs == nints);

                FREE(abrs.rs);
            }

            void free()
            {
                FREE(ints);
            }

            /*
             * Redistribute integrals such that each node has all pq for each rs pair
             */
            void collect(const bool rles)
            {
                MPI::Datatype IDX4_T_TYPE = MPI::Datatype(MPI_INT16_T).Create_contiguous(4);
                IDX4_T_TYPE.Commit();

                int nrs;
                if (rles)
                {
                    nrs = nr*(nr+1)/2;
                }
                else
                {
                    nrs = nr*ns;
                }

                std::cosort(idxs, idxs+nints, ints, ints+nints, sortIntsByRS);

                int *rscount = SAFE_MALLOC(int, nrs);
                std::fill(rscount, rscount+nrs, 0);

                for (int i = 0;i < nints;i++)
                {
                    if (rles)
                    {
                        assert(idxs[i].k+idxs[i].l*(idxs[i].l+1)/2 < nrs);
                        rscount[idxs[i].k+idxs[i].l*(idxs[i].l+1)/2]++;
                    }
                    else
                    {
                        assert(idxs[i].k+idxs[i].l*nr < nrs);
                        rscount[idxs[i].k+idxs[i].l*nr]++;
                    }
                }

                int *rscountall = SAFE_MALLOC(int, nrs*nproc);
                this->comm.Allgather(rscount, nrs, MPI::INT, rscountall, nrs, MPI::INT);
                FREE(rscount);

                size_t nnewints = 0;
                int *sendcount = SAFE_MALLOC(int, nproc);
                int *sendoff = SAFE_MALLOC(int, nproc);
                int *recvcount = SAFE_MALLOC(int, nproc);
                int *recvoff = SAFE_MALLOC(int, nproc);
                std::fill(sendcount, sendcount+nproc, 0);
                sendoff[0] = 0;
                std::fill(recvcount, recvcount+nproc, 0);
                recvoff[0] = 0;
                for (int i = 0;i < nproc;i++)
                {
                    for (int rs = (nrs*rank)/nproc;rs < (nrs*(rank+1))/nproc;rs++)
                    {
                        assert(rs+nrs*i < nproc*nrs);
                        recvcount[i] += rscountall[rs+nrs*i];
                    }
                    if (i > 0) recvoff[i] = recvoff[i-1]+recvcount[i-1];
                    nnewints += recvcount[i];

                    for (int rs = (nrs*i)/nproc;rs < (nrs*(i+1))/nproc;rs++)
                    {
                        assert(rs+nrs*rank < nproc*nrs);
                        sendcount[i] += rscountall[rs+nrs*rank];
                    }
                    if (i > 0) sendoff[i] = sendoff[i-1]+sendcount[i-1];
                }
                assert(recvoff[nproc-1]+recvcount[nproc-1] == nnewints);
                assert(sendoff[nproc-1]+sendcount[nproc-1] == nints);

                assert(allsum((long)nints) == allsum((long)nnewints));

                FREE(rscountall);

                T* newints = SAFE_MALLOC(T, nnewints);
                idx4_t* newidxs = SAFE_MALLOC(idx4_t, nnewints);

                this->comm.Alltoallv(   ints, sendcount, sendoff, MPI::DOUBLE,
                                     newints, recvcount, recvoff, MPI::DOUBLE);

                this->comm.Alltoallv(   idxs, sendcount, sendoff, IDX4_T_TYPE,
                                     newidxs, recvcount, recvoff, IDX4_T_TYPE);

                FREE(sendcount);
                FREE(sendoff);
                FREE(recvcount);
                FREE(recvoff);

                nints = nnewints;
                std::swap(ints, newints);
                std::swap(idxs, newidxs);
                FREE(newints);
                FREE(newidxs);

                std::cosort(idxs, idxs+nints, ints, ints+nints, sortIntsByRS);
            }
        };

        struct abrs_integrals : Distributed<T>
        {
            int na, nb, nr, ns;
            size_t nints;
            T *ints;
            size_t nrs;
            idx2_t *rs;

            /*
             * Expand (p_i q_j|r_k s_l) into (pq|r_k s_l), where pleq = true indicates
             * that (pq|rs) = (qp|rs) and only p_i <= q_j is stored in the input
             */
            abrs_integrals(pqrs_integrals& pqrs, const bool pleq)
            : Distributed<T>(pqrs.ctf)
            {
                na = pqrs.np;
                nb = pqrs.nq;
                nr = pqrs.nr;
                ns = pqrs.ns;

                nrs = (pqrs.nints > 0 ? 1 : 0);
                for (size_t ipqrs = 1;ipqrs < pqrs.nints;ipqrs++)
                {
                    assert( pqrs.idxs[ipqrs].l  > pqrs.idxs[ipqrs-1].l ||
                           (pqrs.idxs[ipqrs].l == pqrs.idxs[ipqrs-1].l &&
                            pqrs.idxs[ipqrs].k >= pqrs.idxs[ipqrs-1].k));
                    if (pqrs.idxs[ipqrs].k != pqrs.idxs[ipqrs-1].k ||
                        pqrs.idxs[ipqrs].l != pqrs.idxs[ipqrs-1].l) nrs++;
                }

                assert(allsum((long)nrs) <= nr*ns);

                rs = SAFE_MALLOC(idx2_t, nrs);
                nints = nrs*na*nb;
                ints = SAFE_MALLOC(T, nints);
                std::fill(ints, ints+nints, (T)0);

                size_t ipqrs = 0, irs = 0, iabrs = 0;
                if (nrs > 0)
                {
                    rs[0].i = pqrs.idxs[0].k;
                    rs[0].j = pqrs.idxs[0].l;
                    assert(rs[0].i >= 0 && rs[0].i < nr);
                    assert(rs[0].j >= 0 && rs[0].j < ns);
                }
                for (;ipqrs < pqrs.nints;ipqrs++)
                {
                    assert(ipqrs >= 0 && ipqrs < pqrs.nints);
                    if (ipqrs > 0 &&
                        (pqrs.idxs[ipqrs].k != pqrs.idxs[ipqrs-1].k ||
                         pqrs.idxs[ipqrs].l != pqrs.idxs[ipqrs-1].l))
                    {
                        irs++;
                        iabrs += na*nb;
                        assert(irs >= 0 && irs < nrs);
                        rs[irs].i = pqrs.idxs[ipqrs].k;
                        rs[irs].j = pqrs.idxs[ipqrs].l;
                        assert(rs[irs].i >= 0 && rs[irs].i < nr);
                        assert(rs[irs].j >= 0 && rs[irs].j < ns);
                    }

                    int p = pqrs.idxs[ipqrs].i;
                    int q = pqrs.idxs[ipqrs].j;

                    assert(iabrs >= 0 && iabrs+na*nb <= nints);
                    assert(p >= 0 && p < na);
                    assert(q >= 0 && q < nb);
                    ints[iabrs+p+q*na] = pqrs.ints[ipqrs];
                    if (p != q && pleq)
                        ints[iabrs+q+p*na] = pqrs.ints[ipqrs];
                }
                assert(irs == nrs-1 || nrs == 0);
                assert(iabrs == nints-na*nb || nrs == 0);
                assert(ipqrs == pqrs.nints);

                pqrs.free();
            }

            /*
             * Transform (ab|rs) -> (cb|rs) (index = A) or (ab|rs) -> (ac|rs) (index = B)
             *
             * C is ldc*nc if trans = 'N' and ldc*[na|nb] if trans = 'T'
             */
            abrs_integrals transform(Index index, const char trans, const int nc, const double* C, const int ldc)
            {
                abrs_integrals out = *this;

                out.rs = SAFE_MALLOC(idx2_t, nrs);
                std::copy(rs, rs+nrs, out.rs);

                if (index == A)
                {
                    out.na = nc;
                    out.nints = nrs*nc*nb;
                    out.ints = SAFE_MALLOC(T, out.nints);

                    if (nc == 0) return out;

                    if (trans == 'N')
                    {
                        gemm('T', 'N', nc, nb*nrs, na, 1.0,        C, ldc,
                                                                ints,  na,
                                                       0.0, out.ints,  nc);
                    }
                    else
                    {
                        gemm('N', 'N', nc, nb*nrs, na, 1.0,        C, ldc,
                                                                ints,  na,
                                                       0.0, out.ints,  nc);
                    }
                }
                else
                {
                    out.nb = nc;
                    out.nints = nrs*na*nc;
                    out.ints = SAFE_MALLOC(T, out.nints);

                    if (nc == 0) return out;

                    size_t iin = 0;
                    size_t iout = 0;
                    for (size_t irs = 0;irs < nrs;irs++)
                    {
                        if (trans == 'N')
                        {
                            gemm('N', 'N', na, nc, nb, 1.0,      ints+iin,  na,
                                                                        C, ldc,
                                                       0.0, out.ints+iout,  na);
                        }
                        else
                        {
                            gemm('N', 'T', na, nc, nb, 1.0,      ints+iin,  na,
                                                                        C, ldc,
                                                       0.0, out.ints+iout,  na);
                        }

                        iin += na*nb;
                        iout += na*nc;
                    }
                }

                return out;
            }

            void transcribe(tensor::DistTensor<T>& tensor, bool assymij, bool assymkl, bool reverse)
            {
                assert(nr*ns == allsum((long)nrs));

                if (reverse)
                {
                    if (assymij) assert(ns == nb);
                    if (assymkl) assert(nr == na);

                    std::vector<tkv_pair<T> > pairs;

                    /*
                     * (ab|rs) -> <sb|ra>
                     */
                    for (size_t idx = 0,irs = 0;irs < nrs;irs++)
                    {
                        int r = rs[irs].i;
                        int s = rs[irs].j;

                        for (int b = 0;b < nb;b++)
                        {
                            for (int a = 0;a < na;a++)
                            {
                                if ((!assymij || s < b) && (!assymkl || r < a))
                                    pairs.push_back(tkv_pair<T>(((a*nr+r)*nb+b)*ns+s, ints[idx]));
                                idx++;
                            }
                        }
                    }

                    if (assymij) assert(tensor.getSymmetry()[0] == AS);
                    if (assymkl) assert(tensor.getSymmetry()[2] == AS);
                    assert(tensor.getLengths()[0] == ns);
                    assert(tensor.getLengths()[1] == nb);
                    assert(tensor.getLengths()[2] == nr);
                    assert(tensor.getLengths()[3] == na);
                    tensor.writeRemoteData(1, 0, pairs.size(), pairs.data());
                    pairs.clear();

                    if (assymij)
                    {
                        /*
                         * -(ab|rs) -> <bs|ra>
                         */
                        for (size_t idx = 0,irs = 0;irs < nrs;irs++)
                        {
                            int r = rs[irs].i;
                            int s = rs[irs].j;

                            for (int b = 0;b < nb;b++)
                            {
                                for (int a = 0;a < na;a++)
                                {
                                    if (b < s && (!assymkl || r < a))
                                        pairs.push_back(tkv_pair<T>(((a*nr+r)*ns+s)*nb+b, ints[idx]));
                                    idx++;
                                }
                            }
                        }

                        tensor.writeRemoteData(-1, 1, pairs.size(), pairs.data());
                    }
                    else if (assymkl)
                    {
                        /*
                         * -(ab|rs) -> <sb|ar>
                         */
                        for (size_t idx = 0,irs = 0;irs < nrs;irs++)
                        {
                            int r = rs[irs].i;
                            int s = rs[irs].j;

                            for (int b = 0;b < nb;b++)
                            {
                                for (int a = 0;a < na;a++)
                                {
                                    if (a < r)
                                        pairs.push_back(tkv_pair<T>(((r*na+a)*nb+b)*ns+s, ints[idx]));
                                    idx++;
                                }
                            }
                        }

                        tensor.writeRemoteData(-1, 1, pairs.size(), pairs.data());
                    }
                }
                else
                {
                    if (assymij) assert(na == nr);
                    if (assymkl) assert(nb == ns);

                    std::vector<tkv_pair<T> > pairs;

                    /*
                     * (ab|rs) -> <ar|bs>
                     */
                    for (size_t idx = 0,irs = 0;irs < nrs;irs++)
                    {
                        int r = rs[irs].i;
                        int s = rs[irs].j;

                        for (int b = 0;b < nb;b++)
                        {
                            for (int a = 0;a < na;a++)
                            {
                                if ((!assymij || a < r) && (!assymkl || b < s))
                                {
                                    pairs.push_back(tkv_pair<T>(((s*nb+b)*nr+r)*na+a, ints[idx]));
                                }
                                idx++;
                            }
                        }
                    }

                    if (assymij) assert(tensor.getSymmetry()[0] == AS);
                    if (assymkl) assert(tensor.getSymmetry()[2] == AS);
                    assert(tensor.getLengths()[0] == na);
                    assert(tensor.getLengths()[1] == nr);
                    assert(tensor.getLengths()[2] == nb);
                    assert(tensor.getLengths()[3] == ns);
                    tensor.writeRemoteData(1, 0, pairs.size(), pairs.data());
                    pairs.clear();

                    if (assymij)
                    {
                        /*
                         * -(ab|rs) -> <ra|bs>
                         */
                        for (size_t idx = 0,irs = 0;irs < nrs;irs++)
                        {
                            int r = rs[irs].i;
                            int s = rs[irs].j;

                            for (int b = 0;b < nb;b++)
                            {
                                for (int a = 0;a < na;a++)
                                {
                                    if (r < a && (!assymkl || b < s))
                                        pairs.push_back(tkv_pair<T>(((s*nb+b)*na+a)*nr+r, ints[idx]));
                                    idx++;
                                }
                            }
                        }

                        tensor.writeRemoteData(-1, 1, pairs.size(), pairs.data());
                    }
                    else if (assymkl)
                    {
                        /*
                         * -(ab|rs) -> <ar|sb>
                         */
                        for (size_t idx = 0,irs = 0;irs < nrs;irs++)
                        {
                            int r = rs[irs].i;
                            int s = rs[irs].j;

                            for (int b = 0;b < nb;b++)
                            {
                                for (int a = 0;a < na;a++)
                                {
                                    if (s < b)
                                        pairs.push_back(tkv_pair<T>(((b*ns+s)*nr+r)*na+a, ints[idx]));
                                    idx++;
                                }
                            }
                        }

                        tensor.writeRemoteData(-1, 1, pairs.size(), pairs.data());
                    }
                }
            }

            void free()
            {
                FREE(ints);
                FREE(rs);
            }
        };

    protected:
        void doTransformation(const AOIntegrals<T>& ints)
        {
            int N = this->uhf.getMolecule().getNumOrbitals();
            int nI = this->uhf.getMolecule().getNumAlphaElectrons();
            int ni = this->uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;

            int sizeAAII[] = {nA, nA, nI, nI};
            int sizeaaii[] = {na, na, ni, ni};
            int shapeNNNN[] = {NS, NS, NS, NS};

            tensor::DistTensor<T> ABIJ__(this->ctf, 4, sizeAAII, shapeNNNN, false);
            tensor::DistTensor<T> abij__(this->ctf, 4, sizeaaii, shapeNNNN, false);

            int64_t npair;
            T *cA, *ca, *cI, *ci;

            /*
             * Read transformation coefficients
             */
            this->uhf.getCA().getAllData(npair, cA);
            assert(npair == N*nA);
            this->uhf.getCa().getAllData(npair, ca);
            assert(npair == N*na);
            this->uhf.getCI().getAllData(npair, cI);
            assert(npair == N*nI);
            this->uhf.getCi().getAllData(npair, ci);
            assert(npair == N*ni);

            /*
             * Resort integrals so that each node has (pq|r_k s_l) where pq
             * are dense blocks for each sparse rs pair
             */
            pqrs_integrals pqrs(ints);
            pqrs.collect(true);
            abrs_integrals PQrs(pqrs, true);

            /*
             * First quarter-transformation
             */
            abrs_integrals PArs = PQrs.transform(B, 'N', nA, cA, N);
            abrs_integrals Pars = PQrs.transform(B, 'N', na, ca, N);
            abrs_integrals PIrs = PQrs.transform(B, 'N', nI, cI, N);
            abrs_integrals Pirs = PQrs.transform(B, 'N', ni, ci, N);
            PQrs.free();

            /*
             * Second quarter-transformation
             */
            abrs_integrals ABrs = PArs.transform(A, 'N', nA, cA, N);
            PArs.free();
            abrs_integrals abrs = Pars.transform(A, 'N', na, ca, N);
            Pars.free();
            abrs_integrals AIrs = PIrs.transform(A, 'N', nA, cA, N);
            abrs_integrals IJrs = PIrs.transform(A, 'N', nI, cI, N);
            PIrs.free();
            abrs_integrals airs = Pirs.transform(A, 'N', na, ca, N);
            abrs_integrals ijrs = Pirs.transform(A, 'N', ni, ci, N);
            Pirs.free();

            /*
             * Make <AB||CD>
             */
            pqrs_integrals rsAB(ABrs);
            rsAB.collect(false);

            abrs_integrals RSAB(rsAB, true);
            abrs_integrals RDAB = RSAB.transform(B, 'N', nA, cA, N);
            RSAB.free();

            abrs_integrals CDAB = RDAB.transform(A, 'N', nA, cA, N);
            RDAB.free();
            CDAB.transcribe(*this->ABCD_, true, true, false);
            CDAB.free();

            /*
             * Make <Ab|Cd> and <ab||cd>
             */
            pqrs_integrals rsab(abrs);
            rsab.collect(false);

            abrs_integrals RSab(rsab, true);
            abrs_integrals RDab = RSab.transform(B, 'N', nA, cA, N);
            abrs_integrals Rdab = RSab.transform(B, 'N', na, ca, N);
            RSab.free();

            abrs_integrals CDab = RDab.transform(A, 'N', nA, cA, N);
            RDab.free();
            CDab.transcribe(*this->AbCd_, false, false, false);
            CDab.free();

            abrs_integrals cdab = Rdab.transform(A, 'N', na, ca, N);
            Rdab.free();
            cdab.transcribe(*this->abcd_, true, true, false);
            cdab.free();

            /*
             * Make <AB||CI>, <aB|cI>, and <AB|IJ>
             */
            pqrs_integrals rsAI(AIrs);
            rsAI.collect(false);

            abrs_integrals RSAI(rsAI, true);
            abrs_integrals RCAI = RSAI.transform(B, 'N', nA, cA, N);
            abrs_integrals RcAI = RSAI.transform(B, 'N', na, ca, N);
            abrs_integrals RJAI = RSAI.transform(B, 'N', nI, cI, N);
            RSAI.free();

            abrs_integrals BCAI = RCAI.transform(A, 'N', nA, cA, N);
            RCAI.free();
            BCAI.transcribe(*this->ABCI_, true, false, false);
            BCAI.free();

            abrs_integrals bcAI = RcAI.transform(A, 'N', na, ca, N);
            RcAI.free();
            bcAI.transcribe(*this->aBcI_, false, false, false);
            bcAI.free();

            abrs_integrals BJAI = RJAI.transform(A, 'N', nA, cA, N);
            RJAI.free();
            BJAI.transcribe(ABIJ__, false, false, false);
            BJAI.free();

            /*
             * Make <Ab|Ci>, <ab||ci>, <Ab|Ij>, and <ab|ij>
             */
            pqrs_integrals rsai(airs);
            rsai.collect(false);

            abrs_integrals RSai(rsai, true);
            abrs_integrals RCai = RSai.transform(B, 'N', nA, cA, N);
            abrs_integrals Rcai = RSai.transform(B, 'N', na, ca, N);
            abrs_integrals RJai = RSai.transform(B, 'N', nI, cI, N);
            abrs_integrals Rjai = RSai.transform(B, 'N', ni, ci, N);
            RSai.free();

            abrs_integrals BCai = RCai.transform(A, 'N', nA, cA, N);
            RCai.free();
            BCai.transcribe(*this->AbCi_, false, false, false);
            BCai.free();

            abrs_integrals bcai = Rcai.transform(A, 'N', na, ca, N);
            Rcai.free();
            bcai.transcribe(*this->abci_, true, false, false);
            bcai.free();

            abrs_integrals BJai = RJai.transform(A, 'N', nA, cA, N);
            RJai.free();
            BJai.transcribe(*this->AbIj_, false, false, false);
            BJai.free();

            abrs_integrals bjai = Rjai.transform(A, 'N', na, ca, N);
            Rjai.free();
            bjai.transcribe(abij__, false, false, false);
            bjai.free();

            /*
             * Make <IJ||KL>, <IJ||KA>, <Ij|Ka>, <aI|bJ>, and <AI|BJ>
             */
            pqrs_integrals rsIJ(IJrs);
            rsIJ.collect(false);

            abrs_integrals RSIJ(rsIJ, true);
            abrs_integrals RBIJ = RSIJ.transform(B, 'N', nA, cA, N);
            abrs_integrals RbIJ = RSIJ.transform(B, 'N', na, ca, N);
            abrs_integrals RLIJ = RSIJ.transform(B, 'N', nI, cI, N);
            abrs_integrals RlIJ = RSIJ.transform(B, 'N', ni, ci, N);
            RSIJ.free();

            abrs_integrals ABIJ = RBIJ.transform(A, 'N', nA, cA, N);
            RBIJ.free();
            ABIJ.transcribe(*this->AIBJ_, false, false, false);
            ABIJ.free();

            abrs_integrals abIJ = RbIJ.transform(A, 'N', na, ca, N);
            RbIJ.free();
            abIJ.transcribe(*this->aIbJ_, false, false, false);
            abIJ.free();

            abrs_integrals akIJ = RlIJ.transform(A, 'N', na, ca, N);
            RlIJ.free();
            akIJ.transcribe(*this->IjKa_, false, false, true);
            akIJ.free();

            abrs_integrals AKIJ = RLIJ.transform(A, 'N', nA, cA, N);
            abrs_integrals KLIJ = RLIJ.transform(A, 'N', nI, cI, N);
            RLIJ.free();
            AKIJ.transcribe(*this->IJKA_, true, false, true);
            AKIJ.free();
            KLIJ.transcribe(*this->IJKL_, true, true, false);
            KLIJ.free();

            /*
             * Make <Ij|Kl>, <ij||kl>, <iJ|kA>, <ij||ka>, <Ai|Bj>, and <ai|bj>
             */
            pqrs_integrals rsij(ijrs);
            rsij.collect(false);

            abrs_integrals RSij(rsij, true);
            abrs_integrals RBij = RSij.transform(B, 'N', nA, cA, N);
            abrs_integrals Rbij = RSij.transform(B, 'N', na, ca, N);
            abrs_integrals RLij = RSij.transform(B, 'N', nI, cI, N);
            abrs_integrals Rlij = RSij.transform(B, 'N', ni, ci, N);
            RSij.free();

            abrs_integrals ABij = RBij.transform(A, 'N', nA, cA, N);
            RBij.free();
            ABij.transcribe(*this->AiBj_, false, false, false);
            ABij.free();

            abrs_integrals abij = Rbij.transform(A, 'N', na, ca, N);
            Rbij.free();
            abij.transcribe(*this->aibj_, false, false, false);
            abij.free();

            abrs_integrals AKij = RLij.transform(A, 'N', nA, cA, N);
            abrs_integrals KLij = RLij.transform(A, 'N', nI, cI, N);
            RLij.free();
            AKij.transcribe(*this->iJkA_, false, false, true);
            AKij.free();
            KLij.transcribe(*this->IjKl_, false, false, false);
            KLij.free();

            abrs_integrals akij = Rlij.transform(A, 'N', na, ca, N);
            abrs_integrals klij = Rlij.transform(A, 'N', ni, ci, N);
            Rlij.free();
            akij.transcribe(*this->ijka_, true, false, true);
            akij.free();
            klij.transcribe(*this->ijkl_, true, true, false);
            klij.free();

            /*
             * Make <AI||BJ> and <ai||bj>
             */
            (*this->AIBJ_)["AIBJ"] -= ABIJ__["ABJI"];
            (*this->aibj_)["aibj"] -= abij__["abji"];

            /*
             * Make <AB||IJ> and <ab||ij>
             */
            //(*this->ABIJ_)["ABIJ"] = 0.5*ABIJ__["ABIJ"];
            //(*this->abij_)["abij"] = 0.5*abij__["abij"];
            (*this->ABIJ_)["ABIJ"]  = ABIJ__["ABIJ"];
            (*this->ABIJ_)["ABIJ"] -= ABIJ__["ABJI"];
            (*this->abij_)["abij"]  = abij__["abij"];
            (*this->abij_)["abij"] -= abij__["abji"];

            /*
             * Make <Ai|bJ> = -<Ab|Ji> and <aI|Bj> = -<Ba|Ij>
             */
            (*this->AibJ_)["AbJi"] -= (*this->AbIj_)["AbJi"];
            (*this->aIBj_)["BaIj"] -= (*this->AbIj_)["BaIj"];

            if (nA > 0) free(cA);
            if (na > 0) free(ca);
            if (nI > 0) free(cI);
            if (ni > 0) free(ci);
        }
};

}
}

#endif
