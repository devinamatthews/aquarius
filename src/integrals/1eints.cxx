#include "1eints.hpp"

#include "kei.hpp"
#include "ovi.hpp"
#include "nai.hpp"

#define IDX_EQ(i,r,e,j,s,f) ((i) == (j) && (r) == (s) && (e) == (f))
#define IDX_GE(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >= (f)))))
#define IDX_GT(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >  (f)))))

using namespace aquarius::input;
using namespace aquarius::symmetry;
using namespace aquarius::task;
using namespace aquarius::tensor;

namespace aquarius
{
namespace integrals
{

OneElectronIntegrals::OneElectronIntegrals(const Shell& a, const Shell& b)
: sa(a), sb(b), group(a.getCenter().getPointGroup()),
  ca(a.getCenter()), cb(b.getCenter()), la(a.getL()), lb(b.getL()),
  na(a.getNPrim()), nb(b.getNPrim()), ma(a.getNContr()), mb(b.getNContr()),
  da(a.getDegeneracy()), db(b.getDegeneracy()), fsa(a.getNFunc()), fsb(b.getNFunc()),
  za(a.getExponents()), zb(b.getExponents()), num_processed(0)
{
    fca = (la+1)*(la+2)/2;
    fcb = (lb+1)*(lb+2)/2;

    int nints = 0;
    for (int i = 0;i < fsa;i++)
    {
        for (int j = 0;j < fsb;j++)
        {
            for (int r = 0;r < da;r++)
            {
                for (int s = 0;s < db;s++)
                {
                    const Representation& w = group.getIrrep(a.getIrrepOfFunc(i, r));
                    const Representation& x = group.getIrrep(b.getIrrepOfFunc(j, s));

                    if ((w*x).isTotallySymmetric()) nints += ma*mb;
                }
            }
        }
    }

    ints.resize(nints);
}

void OneElectronIntegrals::run()
{
    so(ints.data());
}

size_t OneElectronIntegrals::process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                                     size_t nprocess, double* integrals, idx2_t* indices, double cutoff)
{
    const PointGroup& group = ca.getPointGroup();

    size_t m = 0;
    size_t n = 0;
    for (int i = 0;i < fsa;i++)
    {
        for (int j = 0;j < fsb;j++)
        {
            for (int r = 0;r < da;r++)
            {
                for (int s = 0;s < db;s++)
                {
                    const Representation w = group.getIrrep(sa.getIrrepOfFunc(i,r));
                    const Representation x = group.getIrrep(sb.getIrrepOfFunc(j,s));

                    if (!(w*x).isTotallySymmetric()) continue;

                    for (int e = 0;e < ma;e++)
                    {
                        for (int f = 0;f < mb;f++)
                        {
                            if (num_processed > m)
                            {
                                m++;
                                continue;
                            }

                            if (aquarius::abs(ints[m]) > cutoff && (&sa != &sb || IDX_GE(i,r,e,j,s,f)))
                            {
                                indices[n].i = sa.getIndex(ctx, idxa, i, e, r);
                                indices[n].j = sb.getIndex(ctx, idxb, j, f, s);
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

void OneElectronIntegrals::prim(const vec3& posa, int e,
                                const vec3& posb, int f, double* integrals)
{
    assert(0);
}

void OneElectronIntegrals::prims(const vec3& posa, const vec3& posb,
                                 double* integrals)
{
    #pragma omp parallel for collapse(2)
    for (int e = 0;e < na;e++)
    {
        for (int f = 0;f < nb;f++)
        {
            prim(posa, e, posb, f, integrals+fca*fcb*(e*nb+f));
        }
    }
}

void OneElectronIntegrals::contr(const vec3& posa, const vec3& posb,
                                 double* integrals)
{
    vector<double> pintegrals(fca*fcb*na*nb);
    prims(posa, posb, pintegrals.data());
    prim2contr2r(fca*fcb, pintegrals.data(), integrals);
}

void OneElectronIntegrals::spher(const vec3& posa, const vec3& posb,
                                 double* integrals)
{
    vector<double> cintegrals(fca*fcb*na*nb);
    contr(posa, posb, integrals);
    cart2spher2r(ma*mb, integrals, cintegrals.data());
    transpose(fsa*fsb, ma*mb, 1.0, cintegrals.data(), fsa*fsb,
                              0.0,  integrals       ,   ma*mb);
}

void OneElectronIntegrals::so(double* integrals)
{
    const PointGroup& group = ca.getPointGroup();

    vector<double> aointegrals(fca*fcb*na*nb);

    int lambdar;
    vector<int> dcrr = group.DCR(ca.getStabilizer(), cb.getStabilizer(), lambdar);
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
    for (int i = 0;i < fsa;i++)
    {
        for (int j = 0;j < fsb;j++)
        {
            for (int e = 0;e < da;e++)
            {
                for (int f = 0;f < db;f++)
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

    m = fsa;
    n = fcb*nother;
    k = fca;

    if (sa.isSpherical())
    {
        gemm('T', 'T', m, n, k, 1.0, sa.getCart2Spher().data(), k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = fsb;
    n = fsa*nother;
    k = fcb;

    if (sb.isSpherical())
    {
        gemm('T', 'T', m, n, k, 1.0, sb.getCart2Spher().data(), k, buf2, n, 0.0, buf1, m);
    }
    else
    {
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::cart2spher2l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    n = fsb;
    m = fca*nother;
    k = fcb;

    if (sb.isSpherical())
    {
        gemm('T', 'N', m, n, k, 1.0, buf1, k, sb.getCart2Spher().data(), k, 0.0, buf2, m);
    }
    else
    {
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    n = fsa;
    m = fsb*nother;
    k = fcb;

    if (sa.isSpherical())
    {
        gemm('T', 'N', m, n, k, 1.0, buf2, k, sa.getCart2Spher().data(), k, 0.0, buf1, m);
    }
    else
    {
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::prim2contr2r(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    m = ma;
    n = nb*nother;
    k = na;
    gemm('T', 'T', m, n, k, 1.0, sa.getCoefficients().data(), k, buf1, n, 0.0, buf2, m);

    m = mb;
    n = ma*nother;
    k = nb;
    gemm('T', 'T', m, n, k, 1.0, sb.getCoefficients().data(), k, buf2, n, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::prim2contr2l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    m = mb;
    n = na*nother;
    k = nb;
    gemm('T', 'N', m, n, k, 1.0, buf1, k, sb.getCoefficients().data(), k, 0.0, buf2, m);

    m = ma;
    n = mb*nother;
    k = na;
    gemm('T', 'N', m, n, k, 1.0, buf2, k, sa.getCoefficients().data(), k, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

}
}

REGISTER_TASK(aquarius::integrals::Ishida1eIntegralsTask,"1eints");
