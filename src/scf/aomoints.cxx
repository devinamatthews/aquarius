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
 * ARE DISCLAIMED. IN NO EVENT SHALL EDGAR SOLOMONIK BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "aomoints.hpp"

#include <vector>
#include <algorithm>
#include <cassert>

#include "memory/memory.h"

using namespace std;
using namespace MPI;

namespace aquarius
{
namespace scf
{

enum Side {PQ, RS};
enum Index {A, B};

static bool sortIntsByRS(const integral_t& i1, const integral_t& i2)
{
    if (i1.idx.k < i2.idx.k)
    {
        return true;
    }
    else if (i1.idx.k > i2.idx.k)
    {
        return false;
    }
    if (i1.idx.l < i2.idx.l)
    {
        return true;
    }
    else
    {
        return false;
    }
}

struct abrs_integrals;

struct pqrs_integrals
{
    int np, nq, nr, ns;
    size_t nints;
    integral_t *ints;

    /*
     * Read integrals in and break (pq|rs)=(rs|pq) symmetry
     */
    pqrs_integrals(const AOIntegrals& aoints)
    {
        ns = nr = nq = np = aoints.getMolecule().getNumOrbitals();

        size_t noldints = aoints.getNumInts();
        integral_t *oldints = aoints.getInts();
        nints = noldints;

        for (size_t i = 0;i < noldints;i++)
        {
            idx4_t idx = oldints[i].idx;

            if (idx.i != idx.k || idx.j != idx.l)
            {
                nints++;
            }
        }

        ints = SAFE_MALLOC(integral_t, nints);
        copy(oldints, oldints+noldints, ints);

        size_t j = noldints;
        for (size_t i = 0;i < noldints;i++)
        {
            idx4_t idx = oldints[i].idx;

            if (idx.i != idx.k || idx.j != idx.l)
            {
                ints[j] = oldints[i];
                swap(ints[j].idx.i, ints[j].idx.k);
                swap(ints[j].idx.j, ints[j].idx.l);
            }
        }
    }

    pqrs_integrals(const abrs_integrals& abrs)
    {
        np = abrs.na;
        nq = abrs.nb;
        nr = abrs.nr;
        ns = abrs.ns;

        nints = abrs.nints;
        ints = SAFE_MALLOC(integral_t, nints);
        for (size_t ipqrs = 0, irs = 0;irs < abrs.nrs;irs++)
        {
            int r = abrs.rs[irs].i;
            int s = abrs.rs[irs].j;

            for (int q = 0;q < nq;q++)
            {
                for (int p = 0;p < np;p++)
                {
                    ints[ipqrs].idx.i = p;
                    ints[ipqrs].idx.j = q;
                    ints[ipqrs].idx.k = r;
                    ints[ipqrs].idx.l = s;
                    ints[ipqrs].value = abrs.ints[ipqrs];
                    ipqrs++;
                }
            }
        }
    }

    void free()
    {
        FREE(ints);
    }

    void transpose()
    {
        for (size_t ipqrs = 0;ipqrs < nints;ipqrs++)
        {
            int p = ints[ipqrs].idx.i;
            int q = ints[ipqrs].idx.j;
            int r = ints[ipqrs].idx.k;
            int s = ints[ipqrs].idx.l;

            ints[ipqrs].idx.i = r;
            ints[ipqrs].idx.j = s;
            ints[ipqrs].idx.k = p;
            ints[ipqrs].idx.l = q;
        }
    }

    /*
     * Redistribute integrals such that each node has all pq for each rs pair
     */
    void collect(Intracomm& comm, const bool rles)
    {
        int nproc = comm.Get_size();
        int rank = comm.Get_rank();

        int nrs;
        if (rles)
        {
            nrs = nr*(nr+1)/2;
        }
        else
        {
            nrs = nr*ns;
        }

        sort(ints, ints+nints, sortIntsByRS);

        int *rscount = SAFE_MALLOC(int, nrs);
        fill(rscount, rscount+nrs, 0);

        for (int i = 0;i < nints;i++)
        {
            if (rles)
            {
                rscount[ints[i].idx.k+ints[i].idx.l*(ints[i].idx.l-1)/2]++;
            }
            else
            {
                rscount[ints[i].idx.k+ints[i].idx.l*nr]++;
            }
        }

        int *rscountall = SAFE_MALLOC(int, nrs*nproc);
        comm.Allgather(rscount, nrs, INT, rscountall, nrs*nproc, INT);

        FREE(rscount);

        size_t nnewints = 0;
        int *sendcount = SAFE_MALLOC(int, nproc);
        int *sendoff = SAFE_MALLOC(int, nproc);
        int *recvcount = SAFE_MALLOC(int, nproc);
        int *recvoff = SAFE_MALLOC(int, nproc);
        fill(sendcount, sendcount+nproc, 0);
        sendoff[0] = 0;
        fill(recvcount, recvcount+nproc, 0);
        recvoff[0] = 0;
        fill(rscount, rscount+nrs, 0);
        for (int i = 0;i < nproc;i++)
        {
            for (int rs = (nrs*rank)/nproc;rs < (nrs*(rank+1))/nproc;rs++)
            {
                recvcount[i] += rscountall[rs+nrs*i];
            }
            if (i > 0) recvoff[i] = recvoff[i-1]+recvcount[i];
            nnewints += recvcount[i];

            for (int rs = (nrs*i)/nproc;rs < (nrs*(i+1))/nproc;rs++)
            {
                sendcount[i] += rscountall[rs+nrs*rank];
            }
            if (i > 0) sendoff[i] = sendoff[i-1]+sendcount[i];
        }

        FREE(rscountall);

        integral_t* newints = SAFE_MALLOC(integral_t, nnewints);

        comm.Alltoallv(   ints, sendcount, sendoff, integral_t::mpi_type,
                       newints, recvcount, recvoff, integral_t::mpi_type);

        FREE(sendcount);
        FREE(sendoff);
        FREE(recvcount);
        FREE(recvoff);

        nints = nnewints;
        swap(ints, newints);
        FREE(newints);
    }
};

struct abrs_integrals
{
    int na, nb, nr, ns;
    size_t nints;
    double *ints;
    size_t nrs;
    idx2_t *rs;

    /*
     * Expand (p_i q_j|r_k s_l) into (pq|r_k s_l), where pleq = true indicates
     * that (pq|rs) = (qp|rs) and only p_i <= q_j is stored in the input
     */
    abrs_integrals(const pqrs_integrals& pqrs, const bool pleq = true)
    {
        na = pqrs.np;
        nb = pqrs.nq;
        nr = pqrs.nr;
        ns = pqrs.ns;

        nrs = 1;
        for (size_t ipqrs = 1;ipqrs < pqrs.nints;ipqrs++)
        {
            if (pqrs.ints[ipqrs].idx.k != pqrs.ints[ipqrs-1].idx.k ||
                pqrs.ints[ipqrs].idx.l != pqrs.ints[ipqrs-1].idx.l) nrs++;
        }

        rs = SAFE_MALLOC(idx2_t, nrs);
        nints = nrs*na*nb;
        ints = SAFE_MALLOC(double, nints);
        fill(ints, ints+nints, 0.0);

        rs[0].i = pqrs.ints[0].idx.k;
        rs[0].j = pqrs.ints[0].idx.l;
        for (size_t ipqrs = 1, irs = 1, iabrs = 0;ipqrs < pqrs.nints;ipqrs++)
        {
            if (pqrs.ints[ipqrs].idx.k != pqrs.ints[ipqrs-1].idx.k ||
                pqrs.ints[ipqrs].idx.l != pqrs.ints[ipqrs-1].idx.l)
            {
                rs[irs].i = pqrs.ints[ipqrs].idx.k;
                rs[irs].j = pqrs.ints[ipqrs].idx.l;
                irs++;
                iabrs += na*nb;
            }

            int p = pqrs.ints[ipqrs].idx.i;
            int q = pqrs.ints[ipqrs].idx.j;

            ints[iabrs+p+q*na] = pqrs.ints[ipqrs].value;
            if (p != q && pleq)
                ints[iabrs+q+p*na] = pqrs.ints[ipqrs].value;
        }
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
        copy(rs, rs+nrs, out.rs);

        if (index == A)
        {
            out.na = nc;
            out.nints = nrs*nc*nb;
            out.ints = SAFE_MALLOC(double, out.nints);

            size_t iin = 0;
            size_t iout = 0;
            for (size_t irs = 0;irs < nrs;irs++)
            {
                if (trans == 'N')
                {
                    dgemm('T', 'N', nc, nb, na, 1.0,             C, ldc,
                                                          ints+iin,  na,
                                                0.0, out.ints+iout,  nc);
                }
                else
                {
                    dgemm('N', 'N', nc, nb, na, 1.0,             C, ldc,
                                                          ints+iin,  na,
                                                0.0, out.ints+iout,  nc);
                }

                iin += na*nb;
                iout += nc*nb;
            }
        }
        else
        {
            out.nb = nc;
            out.nints = nrs*na*nc;
            out.ints = SAFE_MALLOC(double, out.nints);

            size_t iin = 0;
            size_t iout = 0;
            for (size_t irs = 0;irs < nrs;irs++)
            {
                if (trans == 'N')
                {
                    dgemm('N', 'N', na, nc, nb, 1.0,      ints+iin,  na,
                                                                 C, ldc,
                                                0.0, out.ints+iout,  na);
                }
                else
                {
                    dgemm('N', 'T', na, nc, nb, 1.0,      ints+iin,  na,
                                                                 C, ldc,
                                                0.0, out.ints+iout,  na);
                }

                iin += na*nb;
                iout += na*nc;
            }
        }

        return out;
    }

    void free()
    {
        FREE(ints);
        FREE(rs);
    }
};

AOMOIntegrals::AOMOIntegrals(DistWorld *dw, AOIntegrals& ints, const UHF& uhf)
: MOIntegrals(dw, uhf)
{
    doTransformation(ints);
}

void AOMOIntegrals::doTransformation(AOIntegrals& ints)
{
    int N = uhf.getMolecule().getNumOrbitals();
    int nI = uhf.getMolecule().getNumAlphaElectrons();
    int ni = uhf.getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;

    int sizeAAII[] = {nA, nA, nI, nI};
    int sizeaaii[] = {na, na, ni, ni};
    int shapeNNNN[] = {NS, NS, NS, NS};

    DistTensor ABIJ__(4, sizeAAII, shapeNNNN, dw, false);
    DistTensor abij__(4, sizeaaii, shapeNNNN, dw, false);

    int npair;
    double *cA, *ca, *cI, *ci;

    /*
     * Read transformation coefficients
     */
    uhf.getCA().getAllData(&npair, &cA);
    assert(npair == N*nA);
    uhf.getCa().getAllData(&npair, &ca);
    assert(npair == N*na);
    uhf.getCI().getAllData(&npair, &cI);
    assert(npair == N*nI);
    uhf.getCi().getAllData(&npair, &ci);
    assert(npair == N*ni);

    /*
     * Canonicalize integrals (p<=q, r<=s, pq<=rs), so that
     * indexing and sorting is predictable
     */
    ints.canonicalize();

    /*
     * Resort integrals so that each node has (pq|r_k s_l) where pq
     * are dense blocks for each sparse rs pair
     */
    pqrs_integrals pqrs(ints);
    pqrs.collect(comm, true);
    abrs_integrals PQrs(pqrs, true);
    pqrs.free();

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
    abrs_integrals abrs = PArs.transform(A, 'N', na, ca, N);
    abrs_integrals AIrs = PIrs.transform(A, 'N', nA, cA, N);
    abrs_integrals airs = Pirs.transform(A, 'N', na, ca, N);
    abrs_integrals IJrs = PIrs.transform(A, 'N', nI, cI, N);
    abrs_integrals ijrs = Pirs.transform(A, 'N', ni, ci, N);
    PArs.free();
    Pars.free();
    PIrs.free();
    Pirs.free();

    /*
     * Make <AB||CD>
     */
    pqrs_integrals rsAB(ABrs);
    ABrs.free();
    rsAB.transpose();
    rsAB.collect(comm, false);

    abrs_integrals RSAB(rsAB);
    rsAB.free();
    abrs_integrals RDAB = RSAB.transform(B, 'N', nA, cA, N);
    RSAB.free();
    abrs_integrals CDAB = RDAB.transform(A, 'N', nA, cA, N);
    RDAB.free();

    /*
     * TODO: get (AB|CD) to <AC|BD> and antisymmetrize
     */

    /*
     * Make <Ab|Cd> and <ab||cd>
     */
    pqrs_integrals rsab(abrs);
    abrs.free();
    rsab.transpose();
    rsab.collect(comm, false);

    abrs_integrals RSab(rsab);
    rsab.free();
    abrs_integrals RDab = RSab.transform(B, 'N', nA, cA, N);
    abrs_integrals Rdab = RSab.transform(B, 'N', na, ca, N);
    RSab.free();
    abrs_integrals CDab = RDab.transform(A, 'N', nA, cA, N);
    abrs_integrals cdab = Rdab.transform(A, 'N', na, ca, N);
    RDab.free();
    Rdab.free();

    /*
     * Make <AB||CI>, <aB|cI>, and <AB|IJ>
     */
    pqrs_integrals rsAI(AIrs);
    AIrs.free();
    rsAI.transpose();
    rsAI.collect(comm, false);

    abrs_integrals RSAI(rsAI);
    rsAI.free();
    abrs_integrals RCAI = RSAI.transform(B, 'N', nA, cA, N);
    abrs_integrals RcAI = RSAI.transform(B, 'N', na, ca, N);
    abrs_integrals RJAI = RSAI.transform(B, 'N', nI, cI, N);
    RSAI.free();
    abrs_integrals BCAI = RCAI.transform(A, 'N', nA, cA, N);
    abrs_integrals bcAI = RcAI.transform(A, 'N', na, ca, N);
    abrs_integrals BJAI = RJAI.transform(A, 'N', nA, cA, N);
    RCAI.free();
    RcAI.free();
    RJAI.free();

    /*
     * Make <Ab|Ci>, <ab||ci>, <Ab|Ij>, and <ab|ij>
     */
    pqrs_integrals rsai(airs);
    airs.free();
    rsai.transpose();
    rsai.collect(comm, false);

    abrs_integrals RSai(rsai);
    rsai.free();
    abrs_integrals RCai = RSai.transform(B, 'N', nA, cA, N);
    abrs_integrals Rcai = RSai.transform(B, 'N', na, ca, N);
    abrs_integrals RJai = RSai.transform(B, 'N', nI, cI, N);
    abrs_integrals Rjai = RSai.transform(B, 'N', ni, ci, N);
    RSai.free();
    abrs_integrals BCai = RCai.transform(A, 'N', nA, cA, N);
    abrs_integrals bcai = Rcai.transform(A, 'N', na, ca, N);
    abrs_integrals BJai = RJai.transform(A, 'N', nA, cA, N);
    abrs_integrals bjai = Rjai.transform(A, 'N', na, ca, N);
    RCai.free();
    Rcai.free();
    RJai.free();
    Rjai.free();

    /*
     * Make <IJ||KL>, <IJ||KA>, <Ij|Ka>, <aI|bJ>, and <AI|BJ>
     */
    pqrs_integrals rsIJ(IJrs);
    IJrs.free();
    rsIJ.transpose();
    rsIJ.collect(comm, false);

    abrs_integrals RSIJ(rsIJ);
    rsIJ.free();
    abrs_integrals RBIJ = RSIJ.transform(B, 'N', nA, cA, N);
    abrs_integrals RbIJ = RSIJ.transform(B, 'N', na, ca, N);
    abrs_integrals RLIJ = RSIJ.transform(B, 'N', nI, cI, N);
    abrs_integrals RlIJ = RSIJ.transform(B, 'N', ni, ci, N);
    RSIJ.free();
    abrs_integrals ABIJ = RBIJ.transform(A, 'N', nA, cA, N);
    abrs_integrals abIJ = RbIJ.transform(A, 'N', na, ca, N);
    abrs_integrals AKIJ = RLIJ.transform(A, 'N', nA, cA, N);
    abrs_integrals akIJ = RlIJ.transform(A, 'N', na, ca, N);
    abrs_integrals KLIJ = RLIJ.transform(A, 'N', nI, cI, N);
    RBIJ.free();
    RbIJ.free();
    RLIJ.free();
    RlIJ.free();

    /*
     * Make <Ij|Kl>, <ij||kl>, <iJ|kA>, <ij||ka>, <Ai|Bj>, and <ai|bj>
     */
    pqrs_integrals rsij(ijrs);
    ijrs.free();
    rsij.transpose();
    rsij.collect(comm, false);

    abrs_integrals RSij(rsij);
    rsij.free();
    abrs_integrals RBij = RSij.transform(B, 'N', nA, cA, N);
    abrs_integrals Rbij = RSij.transform(B, 'N', na, ca, N);
    abrs_integrals RLij = RSij.transform(B, 'N', nI, cI, N);
    abrs_integrals Rlij = RSij.transform(B, 'N', ni, ci, N);
    RSij.free();
    abrs_integrals ABij = RBij.transform(A, 'N', nA, cA, N);
    abrs_integrals abij = Rbij.transform(A, 'N', na, ca, N);
    abrs_integrals AKij = RLij.transform(A, 'N', nA, cA, N);
    abrs_integrals akij = Rlij.transform(A, 'N', na, ca, N);
    abrs_integrals KLij = RLij.transform(A, 'N', nI, cI, N);
    abrs_integrals klij = Rlij.transform(A, 'N', ni, ci, N);
    RBij.free();
    Rbij.free();
    RLij.free();
    Rlij.free();

    /*
     * Make <AI||BJ> and <ai||bj>
     */
    AIBJ_["AIBJ"] -= ABIJ__["ABJI"];
    aibj_["aibj"] -= abij__["abji"];

    /*
     * Make <AB||IJ> and <ab||ij>
     */
    ABIJ_["ABIJ"] = 0.5*ABIJ__["ABIJ"];
    abij_["abij"] = 0.5*abij__["abij"];

    free(cA);
    free(ca);
    free(cI);
    free(ci);
}

}
}
