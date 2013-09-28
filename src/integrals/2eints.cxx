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

#include "2eints.hpp"
#include "internal.h"

#define TMP_BUFSIZE 65536
#define INTEGRAL_CUTOFF 1e-13

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

void ERIEvaluator::operator()(int la, const double* ca, int na, const double *za,
                              int lb, const double* cb, int nb, const double *zb,
                              int lc, const double* cc, int nc, const double *zc,
                              int ld, const double* cd, int nd, const double *zd,
                              double *ints) const
{
    size_t nfunc = (la+1)*(la+2)*(lb+1)*(lb+2)*(lc+1)*(lc+2)*(ld+1)*(ld+2)/16;

    int nt = omp_get_max_threads();
    vector<double*> work(nt);
    for (int i = 0;i < nt;i++)
        work[i] = SAFE_MALLOC(double, (la+1)*(lb+1)*(lc+1)*(ld+1)*(la+lb+lc+ld+1));

    #pragma omp parallel for
    for (int m = 0;m < na*nb*nc*nd;m++)
    {
        int tid = omp_get_thread_num();
        int h = m/(na*nb*nc);
        int r = m%(na*nb*nc);
        int g = r/(na*nb);
        int s = r%(na*nb);
        int f = s/na;
        int e = s%na;
        osprim(la, lb, lc, ld, ca, cb, cc, cd, za[e], zb[f], zc[g], zd[h],
               ints+nfunc*m, work[tid]);
    }

    for (int i = 0;i < nt;i++) FREE(work[i]);
}

TwoElectronIntegrals::TwoElectronIntegrals(const Shell& a, const Shell& b, const Shell& c, const Shell& d,
                                           const TwoElectronIntegralEvaluator& eval)
: a(a), b(b), c(c), d(d), eval(eval), num_processed(0)
{
    const Center& ca = a.getCenter();
    const Center& cb = b.getCenter();
    const Center& cc = c.getCenter();
    const Center& cd = d.getCenter();
    const PointGroup& group = ca.getPointGroup();

    size_t nfunccart = (a.getL()+1)*(a.getL()+2)*(b.getL()+1)*(b.getL()+2)*
                       (c.getL()+1)*(c.getL()+2)*(d.getL()+1)*(d.getL()+2)/16;
    size_t nfuncspher = a.getNFunc()*b.getNFunc()*c.getNFunc()*d.getNFunc();
    size_t nprim = a.getNPrim()*b.getNPrim()*c.getNPrim()*d.getNPrim();
    size_t ncontr = a.getNContr()*b.getNContr()*c.getNContr()*d.getNContr();

    nints = 0;
    for (int l = 0;l < d.getNFunc();l++)
    {
        for (int k = 0;k < c.getNFunc();k++)
        {
            for (int j = 0;j < b.getNFunc();j++)
            {
                for (int i = 0;i < a.getNFunc();i++)
                {
                    for (int u = 0;u < d.getDegeneracy();u++)
                    {
                        for (int t = 0;t < c.getDegeneracy();t++)
                        {
                            for (int s = 0;s < b.getDegeneracy();s++)
                            {
                                for (int r = 0;r < a.getDegeneracy();r++)
                                {
                                    const Representation& w = group.getIrrep(a.getIrrepOfFunc(i, r));
                                    const Representation& x = group.getIrrep(b.getIrrepOfFunc(j, s));
                                    const Representation& y = group.getIrrep(c.getIrrepOfFunc(k, t));
                                    const Representation& z = group.getIrrep(d.getIrrepOfFunc(l, u));

                                    if ((w*x*y*z).isTotallySymmetric()) nints += ncontr;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    ints = SAFE_MALLOC(double, nints);
    fill(ints, ints+nints, 0.0);

    double *aobuf1 = SAFE_MALLOC(double, nfunccart*nprim);
    double *aobuf2 = SAFE_MALLOC(double, nfunccart*nprim);

    int lambdar, lambdas, lambdat;
    vector<int> dcrr = group.DCR(ca.getStabilizer(), cb.getStabilizer(), lambdar);
    vector<int> dcrs = group.DCR(cc.getStabilizer(), cd.getStabilizer(), lambdas);
    vector<int> dcrt = group.DCR(intersection(ca.getStabilizer(), cb.getStabilizer()),
                                 intersection(cc.getStabilizer(), cd.getStabilizer()), lambdat);
    double coef = (double)group.getOrder()/(double)lambdat;

    for (int i = 0;i < dcrr.size();i++)
    {
        for (int j = 0;j < dcrs.size();j++)
        {
            for (int k = 0;k < dcrt.size();k++)
            {
                int r = dcrr[i];
                int s = dcrs[j];
                int t = dcrt[k];
                int st = group.getOpProduct(s,t);

                eval(a.getL(), ca.getCenter(0),                       a.getNPrim(), a.getExponents(),
                     b.getL(), cb.getCenter(cb.getCenterAfterOp(r)),  b.getNPrim(), b.getExponents(),
                     c.getL(), cc.getCenter(cc.getCenterAfterOp(t)),  c.getNPrim(), c.getExponents(),
                     d.getL(), cd.getCenter(cd.getCenterAfterOp(st)), d.getNPrim(), d.getExponents(),
                     aobuf2);

                prim2contr4r(nfunccart, aobuf2, aobuf1);
                cart2spher4r(ncontr, aobuf1, aobuf2);

                transpose(nfuncspher, ncontr, coef, aobuf2, nfuncspher,
                                               0.0, aobuf1, ncontr);

                ao2so4(ncontr, r, t, st, aobuf1, ints);
            }
        }
    }

    FREE(aobuf1);
    FREE(aobuf2);
}

TwoElectronIntegrals::~TwoElectronIntegrals()
{
    FREE(ints);
}

size_t TwoElectronIntegrals::process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                                     const vector<int>& idxc, const vector<int>& idxd,
                                     size_t nprocess, double* integrals, idx4_t* indices, double cutoff)
{
    const PointGroup& group = a.getCenter().getPointGroup();
    Representation z = group.getIrrep(0);
    Representation yz = group.getIrrep(0);
    Representation xyz = group.getIrrep(0);
    Representation wxyz = group.getIrrep(0);

    size_t m = 0;
    size_t n = 0;
    for (int l = 0;l < d.getNFunc();l++)
    {
        for (int k = 0;k < c.getNFunc();k++)
        {
            for (int j = 0;j < b.getNFunc();j++)
            {
                for (int i = 0;i < a.getNFunc();i++)
                {
                    for (int u = 0;u < d.getDegeneracy();u++)
                    {
                        z = group.getIrrep(d.getIrrepOfFunc(l,u));
                        for (int t = 0;t < c.getDegeneracy();t++)
                        {
                            yz = z;
                            yz *= group.getIrrep(c.getIrrepOfFunc(k,t));
                            for (int s = 0;s < b.getDegeneracy();s++)
                            {
                                xyz = yz;
                                xyz *= group.getIrrep(b.getIrrepOfFunc(j,s));
                                for (int r = 0;r < a.getDegeneracy();r++)
                                {
                                    wxyz = xyz;
                                    wxyz *= group.getIrrep(a.getIrrepOfFunc(i,r));

                                    if (!wxyz.isTotallySymmetric()) continue;

                                    for (int h = 0;h < d.getNContr();h++)
                                    {
                                        for (int g = 0;g < c.getNContr();g++)
                                        {
                                            for (int f = 0;f < b.getNContr();f++)
                                            {
                                                for (int e = 0;e < a.getNContr();e++)
                                                {
                                                    if (num_processed > m)
                                                    {
                                                        m++;
                                                        continue;
                                                    }

                                                    bool bad = false;

                                                    if (&a == &b && !IDX_GE(i,r,e,j,s,f)) bad = true;
                                                    if (&c == &d && !IDX_GE(k,t,g,l,u,h)) bad = true;
                                                    if (&a == &c && &b == &d && !(IDX_GT(i,r,e,k,t,g) ||
                                                        (IDX_EQ(i,r,e,k,t,g) && IDX_GE(j,s,f,l,u,h)))) bad = true;

                                                    if (!bad && abs(ints[m]) > cutoff)
                                                    {
                                                        indices[n].i = a.getIndex(ctx, idxa, i, e, r);
                                                        indices[n].j = b.getIndex(ctx, idxb, j, f, s);
                                                        indices[n].k = c.getIndex(ctx, idxc, k, g, t);
                                                        indices[n].l = d.getIndex(ctx, idxd, l, h, u);
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
                        }
                    }
                }
            }
        }
    }

    return n;
}

void TwoElectronIntegrals::ao2so4(size_t nother, int r, int t, int st, double* aointegrals, double* sointegrals)
{
    const PointGroup& group = a.getCenter().getPointGroup();
    Representation yz = group.getIrrep(0);
    Representation xyz = group.getIrrep(0);
    Representation wxyz = group.getIrrep(0);

    for (int l = 0;l < d.getNFunc();l++)
    {
        for (int k = 0;k < c.getNFunc();k++)
        {
            for (int j = 0;j < b.getNFunc();j++)
            {
                for (int i = 0;i < a.getNFunc();i++)
                {
                    for (int h = 0;h < d.getDegeneracy();h++)
                    {
                        int z = d.getIrrepOfFunc(l,h);
                        for (int g = 0;g < c.getDegeneracy();g++)
                        {
                            int y = c.getIrrepOfFunc(k,g);
                            yz = group.getIrrep(z);
                            yz *= group.getIrrep(y);
                            for (int f = 0;f < b.getDegeneracy();f++)
                            {
                                int x = b.getIrrepOfFunc(j,f);
                                xyz = yz;
                                xyz *= group.getIrrep(x);
                                for (int e = 0;e < a.getDegeneracy();e++)
                                {
                                    int w = a.getIrrepOfFunc(i,e);
                                    wxyz = xyz;
                                    wxyz *= group.getIrrep(w);

                                    if (!wxyz.isTotallySymmetric()) continue;

                                    double fac = b.getParity(j,r) *group.character(x,r)*
                                                 c.getParity(k,t) *group.character(y,t)*
                                                 d.getParity(l,st)*group.character(z,st);
                                    axpy(nother, fac, aointegrals, 1, sointegrals, 1);

                                    sointegrals += nother;
                                }
                            }
                        }
                    }
                    aointegrals += nother;
                }
            }
        }
    }
}

void TwoElectronIntegrals::cart2spher4r(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    m = d.getNFunc();
    n = (a.getL()+1)*(a.getL()+2)/2*(b.getL()+1)*(b.getL()+2)/2*(c.getL()+1)*(c.getL()+2)/2*nother;
    k = (d.getL()+1)*(d.getL()+2)/2;

    if (d.isSpherical())
    {
        // [d,l]' x [xabc,d]' = [l,xabc]
        gemm('T', 'T', m, n, k, 1.0, d.getCart2Spher(), k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [xabc,l]' = [l,xabc]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = c.getNFunc();
    n = d.getNFunc()*(a.getL()+1)*(a.getL()+2)/2*(b.getL()+1)*(b.getL()+2)/2*nother;
    k = (c.getL()+1)*(c.getL()+2)/2;

    if (c.isSpherical())
    {
        // [c,k]' x [lxab,c]' = [k,lxab]
        gemm('T', 'T', m, n, k, 1.0, c.getCart2Spher(), k, buf2, n, 0.0, buf1, m);
    }
    else
    {
        // [lxab,k]' = [k,lxab]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    m = b.getNFunc();
    n = c.getNFunc()*d.getNFunc()*(a.getL()+1)*(a.getL()+2)/2*nother;
    k = (b.getL()+1)*(b.getL()+2)/2;

    if (b.isSpherical())
    {
        // [b,j]' x [klxa,b]' = [j,klxa]
        gemm('T', 'T', m, n, k, 1.0, b.getCart2Spher(), k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [klxa,j]' = [j,klxa]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = a.getNFunc();
    n = b.getNFunc()*c.getNFunc()*d.getNFunc()*nother;
    k = (a.getL()+1)*(a.getL()+2)/2;

    if (a.isSpherical())
    {
        // [a,i]' x [jklx,a]' = [i,jklx]
        gemm('T', 'T', m, n, k, 1.0, a.getCart2Spher(), k, buf2, n, 0.0, buf1, m);
    }
    else
    {
        // [jklx,i]' = [i,jklx]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    copy(m*n, buf1, 1, buf2, 1);
}

void TwoElectronIntegrals::cart2spher4l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    m = a.getNFunc();
    n = (b.getL()+1)*(b.getL()+2)/2*(c.getL()+1)*(c.getL()+2)/2*(d.getL()+1)*(d.getL()+2)/2*nother;
    k = (a.getL()+1)*(a.getL()+2)/2;

    if (a.isSpherical())
    {
        // [a,bcdx]' x [a,i]' = [bcdx,i]
        gemm('T', 'N', m, n, k, 1.0, buf1, k, d.getCart2Spher(), k, 0.0, buf2, m);
    }
    else
    {
        // [i,bcdx]' = [bcdx,i]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = b.getNFunc();
    n = a.getNFunc()*(c.getL()+1)*(c.getL()+2)/2*(d.getL()+1)*(d.getL()+2)/2*nother;
    k = (b.getL()+1)*(b.getL()+2)/2;

    if (b.isSpherical())
    {
        // [b,cdxi]' x [b,j]' = [cdxi,j]
        gemm('T', 'N', m, n, k, 1.0, buf2, k, c.getCart2Spher(), k, 0.0, buf1, m);
    }
    else
    {
        // [j,cdxi]' = [cdxi,j]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    m = c.getNFunc();
    n = a.getNFunc()*b.getNFunc()*(d.getL()+1)*(d.getL()+2)/2*nother;
    k = (c.getL()+1)*(c.getL()+2)/2;

    if (c.isSpherical())
    {
        // [c,dxij]' x [c,k]' = [dxij,k]
        gemm('T', 'N', m, n, k, 1.0, buf1, k, b.getCart2Spher(), k, 0.0, buf2, m);
    }
    else
    {
        // [k,dxij]' = [dxij,k]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = d.getNFunc();
    n = a.getNFunc()*b.getNFunc()*c.getNFunc()*nother;
    k = (d.getL()+1)*(d.getL()+2)/2;

    if (d.isSpherical())
    {
        // [d,xijk]' x [d,l]' = [xijk,l]
        gemm('T', 'N', m, n, k, 1.0, buf2, k, a.getCart2Spher(), k, 0.0, buf1, m);
    }
    else
    {
        // [l,xijk]' = [xijk,l]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    copy(m*n, buf1, 1, buf2, 1);
}

void TwoElectronIntegrals::prim2contr4r(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    // [d,l]' x [xabc,d]' = [l,xabc]
    m = d.getNContr();
    n = a.getNPrim()*b.getNPrim()*c.getNPrim()*nother;
    k = d.getNPrim();
    gemm('T', 'T', m, n, k, 1.0, d.getCoefficients(), k, buf1, n, 0.0, buf2, m);

    // [c,k]' x [lxab,c]' = [k,lxab]
    m = c.getNContr();
    n = d.getNContr()*a.getNPrim()*b.getNPrim()*nother;
    k = c.getNPrim();
    gemm('T', 'T', m, n, k, 1.0, c.getCoefficients(), k, buf2, n, 0.0, buf1, m);

    // [b,j]' x [klxa,b]' = [j,klxa]
    m = b.getNContr();
    n = c.getNContr()*d.getNContr()*a.getNPrim()*nother;
    k = b.getNPrim();
    gemm('T', 'T', m, n, k, 1.0, b.getCoefficients(), k, buf1, n, 0.0, buf2, m);

    // [a,i]' x [jkl,xa]' = [i,jklx]
    m = a.getNContr();
    n = b.getNContr()*c.getNContr()*d.getNContr()*nother;
    k = a.getNPrim();
    gemm('T', 'T', m, n, k, 1.0, a.getCoefficients(), k, buf2, n, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

void TwoElectronIntegrals::prim2contr4l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    // [a,bcdx]' x [a,i] = [bcdx,i]
    m = a.getNContr();
    n = b.getNPrim()*c.getNPrim()*d.getNPrim()*nother;
    k = a.getNPrim();
    gemm('T', 'N', m, n, k, 1.0, buf1, k, a.getCoefficients(), k, 0.0, buf2, m);

    // [b,cdxi]' x [b,j] = [cdxi,j]
    m = b.getNContr();
    n = a.getNContr()*c.getNPrim()*d.getNPrim()*nother;
    k = b.getNPrim();
    gemm('T', 'N', m, n, k, 1.0, buf2, k, b.getCoefficients(), k, 0.0, buf1, m);

    // [c,dxij]' x [c,k] = [dxij,k]
    m = c.getNContr();
    n = a.getNContr()*b.getNContr()*d.getNPrim()*nother;
    k = c.getNPrim();
    gemm('T', 'N', m, n, k, 1.0, buf1, k, c.getCoefficients(), k, 0.0, buf2, m);

    // [d,xijk]' x [d,l] = [xijk,l]
    m = d.getNContr();
    n = a.getNContr()*b.getNContr()*c.getNContr()*nother;
    k = d.getNPrim();
    gemm('T', 'N', m, n, k, 1.0, buf2, k, d.getCoefficients(), k, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

void ERI::print(Printer& p) const
{
    //TODO
}

TwoElectronIntegralsTask::TwoElectronIntegralsTask(const string& name, const Config& config)
: Task("2eints", name)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("molecule", "molecule"));
    addProduct(Product("eri", "I", reqs));
}

void TwoElectronIntegralsTask::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = get<Molecule>("molecule");

    ERI* eri = new ERI(arena);

    Context ctx(Context::ISCF);

    vector<double> tmpval(TMP_BUFSIZE);
    vector<idx4_t> tmpidx(TMP_BUFSIZE);

    const vector<int>& N = molecule.getNumOrbitals();
    int nirrep = molecule.getGroup().getNumIrreps();

    vector<vector<int> > idx = Shell::setupIndices(Context(), molecule);
    vector<Shell> shells(molecule.getShellsBegin(), molecule.getShellsEnd());

    int abcd = 0;
    for (int a = 0;a < shells.size();++a)
    {
        for (int b = 0;b <= a;++b)
        {
            for (int c = 0;c <= a;++c)
            {
                int dmax = c;
                if (a == c) dmax = b;
                for (int d = 0;d <= dmax;++d)
                {
                    if (abcd%arena.nproc == arena.rank)
                    {
                        TwoElectronIntegrals block(shells[a], shells[b], shells[c], shells[d], ERIEvaluator());

                        size_t n;
                        while ((n = block.process(ctx, idx[a], idx[b], idx[c], idx[d],
                                                  TMP_BUFSIZE, tmpval.data(), tmpidx.data(), INTEGRAL_CUTOFF)) != 0)
                        {
                            eri->ints.insert(eri->ints.end(), tmpval.data(), tmpval.data()+n);
                            eri->idxs.insert(eri->idxs.end(), tmpidx.data(), tmpidx.data()+n);
                        }
                    }
                    abcd++;
                }
            }
        }
    }

    //TODO: load balance

    for (int i = 0;i < eri->ints.size();++i)
    {
        if (eri->idxs[i].i > eri->idxs[i].j) std::swap(eri->idxs[i].i, eri->idxs[i].j);
        if (eri->idxs[i].k > eri->idxs[i].l) std::swap(eri->idxs[i].k, eri->idxs[i].l);
        if (eri->idxs[i].i > eri->idxs[i].k ||
            (eri->idxs[i].i == eri->idxs[i].k && eri->idxs[i].j > eri->idxs[i].l))
        {
            std::swap(eri->idxs[i].i, eri->idxs[i].k);
            std::swap(eri->idxs[i].j, eri->idxs[i].l);
        }
    }

    put("I", eri);
}

REGISTER_TASK(TwoElectronIntegralsTask,"2eints");
