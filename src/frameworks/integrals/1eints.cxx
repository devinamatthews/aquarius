#include "1eints.hpp"

#define IDX_EQ(i,r,e,j,s,f) ((i) == (j) && (r) == (s) && (e) == (f))
#define IDX_GE(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >= (f)))))
#define IDX_GT(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >  (f)))))

using namespace aquarius::molecule;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace integrals
{

size_t OneElectronIntegrals::ShellBlock::process(const Context& ctx,
                                                 const vector<int>& idxa,
                                                 const vector<int>& idxb,
                                                 size_t nprocess,
                                                 double* integrals,
                                                 idx2_t* indices,
                                                 double cutoff)
{
    auto& group = a.getCenter().getPointGroup();

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

                            if (aquarius::abs(ints[m]) > cutoff && (&a != &b || IDX_GE(i,r,e,j,s,f)))
                            {
                                indices[n].i = ctx.getIndex(a, idxa, i, e, r);
                                indices[n].j = ctx.getIndex(b, idxb, j, f, s);
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

OneElectronIntegrals::ShellBlock OneElectronIntegrals::calculate(const Shell& a, const Shell& b)
{
    auto& group = a.getCenter().getPointGroup();

    int nints = 0;
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

                    if ((w*x).isTotallySymmetric())
                        nints += a.getNContr()*b.getNContr();
                }
            }
        }
    }

    vector<double> ints(nints);
    so(a.getCenter(), a.getL(), a.getCoefficients()ints.data());

    return ShellBlock(a, b, move(ints));
}

void OneElectronIntegrals::prim(const vec3& posa, int la, double za,
                                const vec3& posb, int lb, double zb,
                                double* integrals)
{
    assert(0);
}

void OneElectronIntegrals::prims(const vec3& posa, int la, const vector<double>& za,
                                 const vec3& posb, int lb, const vector<double>& zb,
                                 double* integrals)
{
    int na = za.size();
    int nb = zb.size();
    int fca = (la+1)*(la+2)/2;
    int fcb = (lb+1)*(lb+2)/2;

    #pragma omp parallel for
    for (int m = 0;m < na*nb;m++)
    {
        int f = m/na;
        int e = m%na;
        prim(posa, la, za[e],
             posb, la, zb[f],
             integrals+fca*fcb*m);
    }
}

void OneElectronIntegrals::contr(const vec3& posa, int la, const vector<double>& za, const matrix<double>& ca,
                                 const vec3& posb, int lb, const vector<double>& zb, const matrix<double>& cb,
                                 double* integrals)
{
    int na = za.size();
    int nb = zb.size();
    int fca = (la+1)*(la+2)/2;
    int fcb = (lb+1)*(lb+2)/2;

    vector<double> pintegrals(fca*fcb*na*nb);
    prims(posa, la, za,
          posb, lb, zb,
          pintegrals.data());
    prim2contr(fca*fcb, ca, cb, pintegrals.data(), integrals);
}

void OneElectronIntegrals::spher(const vec3& posa, int la, const vector<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                                 const vec3& posb, int lb, const vector<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
                                 double* integrals)
{
    int na = ca.length(0);
    int nb = cb.length(0);
    int ma = ca.length(1);
    int mb = cb.length(1);
    int fca = sa.length(0);
    int fcb = sb.length(0);
    int fsa = sa.length(1);
    int fsb = sb.length(1);

    vector<double> cintegrals(fca*fcb*na*nb);
    contr(posa, la, za, ca,
          posb, lb, zb, cb,
          integrals);
    cart2spher(ma*mb, sa, sb, integrals, cintegrals.data());
    transpose(fsa*fsb, ma*mb, 1.0, cintegrals.data(), fsa*fsb,
                              0.0,  integrals       ,   ma*mb);
}

void OneElectronIntegrals::so(const molecule::Center& posa, int la, const vector<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                              const molecule::Center& posb, int lb, const vector<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
                              double* integrals)
{
    const PointGroup& group = posa.getPointGroup();

    int na = za.size();
    int nb = zb.size();
    int fca = (la+1)*(la+2)/2;
    int fcb = (lb+1)*(lb+2)/2;

    vector<double> aointegrals(fca*fcb*na*nb);

    int lambdar;
    vector<int> dcrr = group.DCR(posa.getStabilizer(), b.getStabilizer(), lambdar);
    double coef = (double)group.getOrder()/(double)lambdar;

    for (int r : dcrr)
    {
        spher(ca.getCenter(0), cb.getCenter(cb.getCenterAfterOp(r)), aointegrals.data());
        scal(aointegrals.size(), coef, aointegrals.data(), 1);
        ao2so2(ma*mb, r, aointegrals.data(), integrals);
    }
}

void OneElectronIntegrals::ao2so2(size_t nother, int r, double* aointegrals, double* sointegrals)
{
    for (int j = 0;j < fsb;j++)
    {
        for (int i = 0;i < fsa;i++)
        {
            for (int f = 0;f < db;f++)
            {
                for (int e = 0;e < da;e++)
                {
                    int w = sa.getIrrepOfFunc(i,e);
                    int x = sb.getIrrepOfFunc(j,f);

                    if (!(group.getIrrep(w)*group.getIrrep(x)).isTotallySymmetric()) continue;

                    double fac = sb.getParity(j,r)*group.character(x,r);
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

    m = fsb;
    n = fca*nother;
    k = fcb;

    if (sb.isSpherical())
    {
        // [b,j]' x [xa,b]' = [j,xa]
        gemm('T', 'T', m, n, k, 1.0, sb.getCart2Spher().data(), k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [xa,j]' = [j,xa]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = fsa;
    n = fsb*nother;
    k = fca;

    if (sa.isSpherical())
    {
        // [a,i]' x [jx,a]' = [i,jx]
        gemm('T', 'T', m, n, k, 1.0, sa.getCart2Spher().data(), k, buf2, n, 0.0, buf1, m);
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

    n = fsa;
    m = fcb*nother;
    k = fca;

    if (sa.isSpherical())
    {
        // [a,bx]' x [a,i]  = [bx,i]
        gemm('T', 'N', m, n, k, 1.0, buf1, k, sb.getCart2Spher().data(), k, 0.0, buf2, m);
    }
    else
    {
        // [i,bx]' = [bx,i]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    n = fsb;
    m = fsa*nother;
    k = fcb;

    if (sb.isSpherical())
    {
        // [b,xi]' x [b,j] = [xi,j]
        gemm('T', 'N', m, n, k, 1.0, buf2, k, sa.getCart2Spher().data(), k, 0.0, buf1, m);
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
    m = mb;
    n = na*nother;
    k = nb;
    gemm('T', 'T', m, n, k, 1.0, sb.getCoefficients().data(), k, buf1, n, 0.0, buf2, m);

    // [a,i]' x [jx,a]' = [i,jx]
    m = ma;
    n = mb*nother;
    k = na;
    gemm('T', 'T', m, n, k, 1.0, sa.getCoefficients().data(), k, buf2, n, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::prim2contr2l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    // [a,bx]' x [a,i] = [bx,i]
    m = ma;
    n = nb*nother;
    k = na;
    gemm('T', 'N', m, n, k, 1.0, buf1, k, sb.getCoefficients().data(), k, 0.0, buf2, m);

    // [b,xi]' x [b,j] = [xi,j]
    m = mb;
    n = ma*nother;
    k = nb;
    gemm('T', 'N', m, n, k, 1.0, buf2, k, sa.getCoefficients().data(), k, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

}
}
