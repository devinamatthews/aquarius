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

#include "1eints.hpp"
#include "internal.h"

/**
 * Compute the index of a function in cartesian angular momentum.
 */
#define FUNC_CART(x,y,z) ((((x)*(3+(x)+2*((y)+(z))))/2) + (y))
/**
 * Compute the index of a function in spherical harmonic angular momentum.
 *
 * Regular spherical harmonics are referenced by n=l, l>=m>=-l. Contaminants may also be referenced by
 * n>l>=0, n-l even.
 */
#define FUNC_SPHER(n,l,m) ((((n)-(l))*((n)+(l)-1))/2 + 2*(n) + ((m) > 0 ? -2*(m) : 2*(m)+1))

#define IDX_EQ(i,r,e,j,s,f) ((i) == (j) && (r) == (s) && (e) == (f))
#define IDX_GE(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >= (f)))))
#define IDX_GT(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >  (f)))))

using namespace std;
using namespace aquarius;
using namespace aquarius::integrals;
using namespace aquarius::input;
using namespace aquarius::symmetry;
using namespace aquarius::task;
using namespace aquarius::tensor;

void OVIEvaluator::operator()(int la, const double* ca, int na, const double *za,
                              int lb, const double* cb, int nb, const double *zb,
                              double *ints) const
{
    size_t nfunc = (la+1)*(la+2)*(lb+1)*(lb+2)/4;

    int nt = omp_get_max_threads();
    vector<double*> work(nt);
    for (int i = 0;i < nt;i++) work[i] = SAFE_MALLOC(double, 3*(la+1)*(lb+1));

    #pragma omp parallel for
    for (int m = 0;m < na*nb;m++)
    {
        int f = m/na;
        int e = m%na;
        oviprim(la, lb, ca, cb, za[e], zb[f], ints+nfunc*m, work[omp_get_thread_num()]);
    }

    for (int i = 0;i < nt;i++) FREE(work[i]);
}

void KEIEvaluator::operator()(int la, const double* ca, int na, const double *za,
                              int lb, const double* cb, int nb, const double *zb,
                              double *ints) const
{
    size_t nfunc = (la+1)*(la+2)*(lb+1)*(lb+2)/4;

    int nt = omp_get_max_threads();
    vector<double*> work1(nt), work2(nt);
    for (int i = 0;i < nt;i++)
    {
        work1[i] = SAFE_MALLOC(double, 3*(la+2)*(lb+2));
        work2[i] = SAFE_MALLOC(double, 3*(la+1)*(lb+1));
    }

    #pragma omp parallel for
    for (int m = 0;m < na*nb;m++)
    {
        int tid = omp_get_thread_num();
        int f = m/na;
        int e = m%na;
        keiprim(la, lb, ca, cb, za[e], zb[f], ints+nfunc*m, work1[tid], work2[tid]);
    }

    for (int i = 0;i < nt;i++)
    {
        FREE(work1[i]);
        FREE(work2[i]);
    }
}

void NAIEvaluator::operator()(int la, const double* ca, int na, const double *za,
                              int lb, const double* cb, int nb, const double *zb,
                              double *ints) const
{
    size_t nfunc = (la+1)*(la+2)*(lb+1)*(lb+2)/4;

    int nt = omp_get_max_threads();
    vector<double*> work(nt);
    for (int i = 0;i < nt;i++) work[i] = SAFE_MALLOC(double, (la+1)*(lb+1)*(la+lb+1));

    #pragma omp parallel for
    for (int m = 0;m < na*nb;m++)
    {
        int f = m/na;
        int e = m%na;

        fill(ints+m*nfunc, ints+(m+1)*nfunc, 0.0);

        for (int c = 0;c < centers.size();c++)
        {
            for (int i = 0;i < centers[c].getCenters().size();i++)
            {
                naiprim(la, lb, ca, cb, centers[c].getCenter(i), za[e], zb[f],
                        centers[c].getElement().getCharge(), ints+m*nfunc,
                        work[omp_get_thread_num()]);
            }
        }
    }

    for (int i = 0;i < nt;i++) FREE(work[i]);
}

void OneElectronHamiltonianEvaluator::operator()(int la, const double* ca, int na, const double *za,
                                                 int lb, const double* cb, int nb, const double *zb,
                                                 double *ints) const
{
    size_t nfunc = (la+1)*(la+2)*(lb+1)*(lb+2)/4;

    int nt = omp_get_max_threads();
    vector<double*> work1(nt), work2(nt), work3(nt);
    for (int i = 0;i < nt;i++)
    {
        work1[i] = SAFE_MALLOC(double, (la+1)*(lb+1)*(la+lb+1));
        work2[i] = SAFE_MALLOC(double, 3*(la+2)*(lb+2));
        work3[i] = SAFE_MALLOC(double, 3*(la+1)*(lb+1));
    }

    #pragma omp parallel for
    for (int m = 0;m < na*nb;m++)
    {
        int tid = omp_get_thread_num();
        int f = m/na;
        int e = m%na;

        keiprim(la, lb, ca, cb, za[e], zb[f], ints+nfunc*m, work2[tid], work3[tid]);

        for (int c = 0;c < centers.size();c++)
        {
            for (int i = 0;i < centers[c].getCenters().size();i++)
            {
                naiprim(la, lb, ca, cb, centers[c].getCenter(i), za[e], zb[f],
                        centers[c].getElement().getCharge(), ints+m*nfunc, work1[tid]);
            }
        }
    }

    for (int i = 0;i < nt;i++)
    {
        FREE(work1[i]);
        FREE(work2[i]);
        FREE(work3[i]);
    }
}

OneElectronIntegrals::OneElectronIntegrals(const Shell& a, const Shell& b, const OneElectronIntegralEvaluator& eval)
: a(a), b(b), eval(eval), num_processed(0)
{
    const Center& ca = a.getCenter();
    const Center& cb = b.getCenter();
    const PointGroup& group = ca.getPointGroup();

    size_t nfunccart = (a.getL()+1)*(a.getL()+2)*(b.getL()+1)*(b.getL()+2)/4;
    size_t nfuncspher = a.getNFunc()*b.getNFunc();
    size_t nprim = a.getNPrim()*b.getNPrim();
    size_t ncontr = a.getNContr()*b.getNContr();

    nints = 0;
    for (int j = 0;j < b.getNFunc();j++)
    {
        for (int i = 0;i < a.getNFunc();i++)
        {
            for (int s = 0;s < b.getDegeneracy();s++)
            {
                for (int r = 0;r < a.getDegeneracy();r++)
                {
                    const Representation& w = group.getIrrep(a.getIrrepOfFunc(i, r));
                    const Representation& x = group.getIrrep(b.getIrrepOfFunc(j, s));

                    if ((w*x).isTotallySymmetric()) nints += ncontr;
                }
            }
        }
    }

    ints = SAFE_MALLOC(double, nints);
    fill(ints, ints+nints, 0.0);

    double *aobuf1 = SAFE_MALLOC(double, nfunccart*nprim);
    double *aobuf2 = SAFE_MALLOC(double, nfunccart*nprim);

    int lambdar;
    vector<int> dcrr = group.DCR(ca.getStabilizer(), cb.getStabilizer(), lambdar);
    double coef = (double)group.getOrder()/(double)lambdar;

    for (int i = 0;i < dcrr.size();i++)
    {
        eval(a.getL(), ca.getCenter(0),                            a.getNPrim(), a.getExponents().data(),
             b.getL(), cb.getCenter(cb.getCenterAfterOp(dcrr[i])), b.getNPrim(), b.getExponents().data(),
             aobuf2);

        prim2contr2r(nfunccart, aobuf2, aobuf1);
        cart2spher2r(ncontr, aobuf1, aobuf2);

        transpose(nfuncspher, ncontr, coef, aobuf2, nfuncspher,
                                       0.0, aobuf1, ncontr);

        ao2so2(ncontr, dcrr[i], aobuf1, ints);
    }

    FREE(aobuf1);
    FREE(aobuf2);
}

OneElectronIntegrals::~OneElectronIntegrals()
{
    FREE(ints);
}

size_t OneElectronIntegrals::process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                                     size_t nprocess, double* integrals, idx2_t* indices, double cutoff)
{
    const PointGroup& group = a.getCenter().getPointGroup();

    size_t m = 0;
    size_t n = 0;
    for (int j = 0;j < b.getNFunc();j++)
    {
        for (int i = 0;i < a.getNFunc();i++)
        {
            for (int s = 0;s < b.getDegeneracy();s++)
            {
                for (int r = 0;r < a.getDegeneracy();r++)
                {
                    const Representation w = group.getIrrep(a.getIrrepOfFunc(i,r));
                    const Representation x = group.getIrrep(b.getIrrepOfFunc(j,s));

                    if (!(w*x).isTotallySymmetric()) continue;

                    for (int f = 0;f < b.getNContr();f++)
                    {
                        for (int e = 0;e < a.getNContr();e++)
                        {
                            if (num_processed > m)
                            {
                                m++;
                                continue;
                            }

                            if (abs(ints[m]) > cutoff && (&a != &b || IDX_GE(i,r,e,j,s,f)))
                            {
                                indices[n].i = a.getIndex(ctx, idxa, i, e, r);
                                indices[n].j = b.getIndex(ctx, idxb, j, f, s);
                                integrals[n++] = ints[m];
                            }

                            num_processed++;
                            m++;

                            if (n >= nprocess) return n;
                        }
                    }
                }
            }
        }
    }

    return n;
}

void OneElectronIntegrals::ao2so2(size_t nother, int r, double* aointegrals, double* sointegrals)
{
    const PointGroup& group = a.getCenter().getPointGroup();

    for (int j = 0;j < b.getNFunc();j++)
    {
        for (int i = 0;i < a.getNFunc();i++)
        {
            for (int f = 0;f < b.getDegeneracy();f++)
            {
                for (int e = 0;e < a.getDegeneracy();e++)
                {
                    int w = a.getIrrepOfFunc(i,e);
                    int x = b.getIrrepOfFunc(j,f);

                    if (!(group.getIrrep(w)*group.getIrrep(x)).isTotallySymmetric()) continue;

                    double fac = b.getParity(j,r)*group.character(x,r);
                    axpy(nother, fac, aointegrals, 1, sointegrals, 1);

                    sointegrals += nother;
                }
            }
            aointegrals += nother;
        }
    }
}

void OneElectronIntegrals::cart2spher2r(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    m = b.getNFunc();
    n = (a.getL()+1)*(a.getL()+2)/2*nother;
    k = (b.getL()+1)*(b.getL()+2)/2;

    if (b.isSpherical())
    {
        // [b,j]' x [xa,b]' = [j,xa]
        gemm('T', 'T', m, n, k, 1.0, b.getCart2Spher().data(), k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [xa,j]' = [j,xa]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = a.getNFunc();
    n = b.getNFunc()*nother;
    k = (a.getL()+1)*(a.getL()+2)/2;

    if (a.isSpherical())
    {
        // [a,i]' x [jx,a]' = [i,jx]
        gemm('T', 'T', m, n, k, 1.0, a.getCart2Spher().data(), k, buf2, n, 0.0, buf1, m);
    }
    else
    {
        // [jx,i]' = [i,jx]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::cart2spher2l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    n = a.getNFunc();
    m = (b.getL()+1)*(b.getL()+2)/2*nother;
    k = (a.getL()+1)*(a.getL()+2)/2;

    if (a.isSpherical())
    {
        // [a,bx]' x [a,i]  = [bx,i]
        gemm('T', 'N', m, n, k, 1.0, buf1, k, b.getCart2Spher().data(), k, 0.0, buf2, m);
    }
    else
    {
        // [i,bx]' = [bx,i]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    n = b.getNFunc();
    m = a.getNFunc()*nother;
    k = (b.getL()+1)*(b.getL()+2)/2;

    if (b.isSpherical())
    {
        // [b,xi]' x [b,j] = [xi,j]
        gemm('T', 'N', m, n, k, 1.0, buf2, k, a.getCart2Spher().data(), k, 0.0, buf1, m);
    }
    else
    {
        // [j,xi]' = [xi,j]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::prim2contr2r(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    // [b,j]' x [xa,b]' = [j,xa]
    m = b.getNContr();
    n = a.getNPrim()*nother;
    k = b.getNPrim();
    gemm('T', 'T', m, n, k, 1.0, b.getCoefficients().data(), k, buf1, n, 0.0, buf2, m);

    // [a,i]' x [jx,a]' = [i,jx]
    m = a.getNContr();
    n = b.getNContr()*nother;
    k = a.getNPrim();
    gemm('T', 'T', m, n, k, 1.0, a.getCoefficients().data(), k, buf2, n, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::prim2contr2l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    // [a,bx]' x [a,i] = [bx,i]
    m = a.getNContr();
    n = b.getNPrim()*nother;
    k = a.getNPrim();
    gemm('T', 'N', m, n, k, 1.0, buf1, k, b.getCoefficients().data(), k, 0.0, buf2, m);

    // [b,xi]' x [b,j] = [xi,j]
    m = b.getNContr();
    n = a.getNContr()*nother;
    k = b.getNPrim();
    gemm('T', 'N', m, n, k, 1.0, buf2, k, a.getCoefficients().data(), k, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

OneElectronIntegralsTask::OneElectronIntegralsTask(const string& name, const Config& config)
: Task("1eints", name)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("molecule", "molecule"));
    addProduct(Product("ovi", "S", reqs));
    addProduct(Product("kei", "T", reqs));
    addProduct(Product("nai", "G", reqs));
    addProduct(Product("1ehamiltonian", "H", reqs));
}

void OneElectronIntegralsTask::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = get<Molecule>("molecule");

    Context ctx(Context::ISCF);

    const vector<int>& N = molecule.getNumOrbitals();
    int n = molecule.getGroup().getNumIrreps();

    vector<int> irrep;
    for (int i = 0;i < n;i++) irrep += vector<int>(N[i],i);

    vector<uint16_t> start(n,0);
    for (int i = 1;i < n;i++) start[i] = start[i-1]+N[i-1];

    vector<vector<int> > idx = Shell::setupIndices(ctx, molecule);
    vector<Shell> shells(molecule.getShellsBegin(), molecule.getShellsEnd());
    vector<vector<tkv_pair<double> > > ovi_pairs(n), nai_pairs(n), kei_pairs(n);
    vector<Center> centers;

    for (vector<Atom>::const_iterator i = molecule.getAtomsBegin();i != molecule.getAtomsEnd();++i)
    {
        centers.push_back(i->getCenter());
    }

    int block = 0;
    for (int a = 0;a < shells.size();++a)
    {
        for (int b = 0;b <= a;++b)
        {
            if (block%arena.nproc == arena.rank)
            {
                OneElectronIntegrals s(shells[a], shells[b], OVIEvaluator());
                OneElectronIntegrals t(shells[a], shells[b], KEIEvaluator());
                OneElectronIntegrals g(shells[a], shells[b], NAIEvaluator(centers));

                size_t nint = s.getNumInts();
                vector<double> ints(nint);
                vector<idx2_t> idxs(nint);
                size_t nproc;

                nproc = s.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                for (int k = 0;k < nproc;k++)
                {
                    int irr = irrep[idxs[k].i];
                    assert(irr == irrep[idxs[k].j]);

                    uint16_t i = idxs[k].i-start[irr];
                    uint16_t j = idxs[k].j-start[irr];

                                ovi_pairs[irr].push_back(tkv_pair<double>(i*N[irr]+j, ints[k]));
                    if (i != j) ovi_pairs[irr].push_back(tkv_pair<double>(j*N[irr]+i, ints[k]));
                }

                nproc = t.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                for (int k = 0;k < nproc;k++)
                {
                    int irr = irrep[idxs[k].i];
                    assert(irr == irrep[idxs[k].j]);

                    uint16_t i = idxs[k].i-start[irr];
                    uint16_t j = idxs[k].j-start[irr];

                                kei_pairs[irr].push_back(tkv_pair<double>(i*N[irr]+j, ints[k]));
                    if (i != j) kei_pairs[irr].push_back(tkv_pair<double>(j*N[irr]+i, ints[k]));
                }

                nproc = g.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                for (int k = 0;k < nproc;k++)
                {
                    int irr = irrep[idxs[k].i];
                    assert(irr == irrep[idxs[k].j]);

                    uint16_t i = idxs[k].i-start[irr];
                    uint16_t j = idxs[k].j-start[irr];

                                nai_pairs[irr].push_back(tkv_pair<double>(i*N[irr]+j, ints[k]));
                    if (i != j) nai_pairs[irr].push_back(tkv_pair<double>(j*N[irr]+i, ints[k]));
                }
            }

            block++;
        }
    }

    OVI *ovi = new OVI(arena, molecule.getGroup(), N);
    KEI *kei = new KEI(arena, molecule.getGroup(), N);
    NAI *nai = new NAI(arena, molecule.getGroup(), N);
    OneElectronHamiltonian *oeh = new OneElectronHamiltonian(arena, molecule.getGroup(), N);

    for (int i = 0;i < n;i++)
    {
        vector<int> irreps(2,i);
        (*ovi).writeRemoteData(irreps, ovi_pairs[i]);
        (*kei).writeRemoteData(irreps, kei_pairs[i]);
        (*nai).writeRemoteData(irreps, nai_pairs[i]);
        (*oeh).writeRemoteData(irreps, kei_pairs[i]);
        (*oeh).writeRemoteData(irreps, 1.0, 1.0, nai_pairs[i]);
    }

    put("S", ovi);
    put("T", kei);
    put("G", nai);
    put("H", oeh);
}

REGISTER_TASK(OneElectronIntegralsTask,"1eints");
