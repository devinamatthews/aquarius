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

void transform(const matrix<double>& ai, const matrix<double>& bj, size_t nother,
               double* buf1, double* buf2)
{
    using transpose::T;

    int na = ai.length(0);
    int nb = bj.length(0);
    int ma = ai.length(1);
    int mb = bj.length(1);

    marray_view<double,3> xab(nother, na, nb, buf1);
    marray_view<double,3> jxa(mb, nother, na, buf2);
    marray_view<double,3> ijx(ma, mb, nother, buf1);

    gemm(1.0, bj^T, xab.lower(2)^T, 0.0, jxa.lower(1));
    gemm(1.0, ai^T, jxa.lower(2)^T, 0.0, ijx.lower(1));

    copy(ma*mb*nother, buf1, 1, buf2, 1);
}

void transform(size_t nother, const matrix<double>& ai, const matrix<double>& bj,
               double* buf1, double* buf2)
{
    using transpose::T;

    int na = ai.length(0);
    int nb = bj.length(0);
    int ma = ai.length(1);
    int mb = bj.length(1);

    marray_view<double,3> abx(na, nb, nother, buf1);
    marray_view<double,3> bxi(nb, nother, ma, buf2);
    marray_view<double,3> xij(nother, ma, mb, buf1);

    gemm(1.0, abx.lower(1)^T, ai, 0.0, bxi.lower(2));
    gemm(1.0, bxi.lower(1)^T, bj, 0.0, xij.lower(2));

    copy(ma*mb*nother, buf1, 1, buf2, 1);
}

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

                            if (std::abs(ints[m]) > cutoff && (&a != &b || IDX_GE(i,r,e,j,s,f)))
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
    so(a, b, ints.data());

    return {a, b, Context::SFIC, move(ints)};
}

void OneElectronIntegrals::prim(const vec3& posa, int la, double za,
                                const vec3& posb, int lb, double zb,
                                double* integrals)
{
    assert(0);
}

void OneElectronIntegrals::prims(const vec3& posa, int la, const row<double>& za,
                                 const vec3& posb, int lb, const row<double>& zb,
                                 double* integrals)
{
    int na = za.length();
    int nb = zb.length();
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

void OneElectronIntegrals::contr(const vec3& posa, int la, const row<double>& za, const matrix<double>& ca,
                                 const vec3& posb, int lb, const row<double>& zb, const matrix<double>& cb,
                                 double* integrals)
{
    int na = za.length();
    int nb = zb.length();
    int fca = (la+1)*(la+2)/2;
    int fcb = (lb+1)*(lb+2)/2;

    vector<double> pintegrals(fca*fcb*na*nb);
    prims(posa, la, za,
          posb, lb, zb,
          pintegrals.data());
    transform(fca*fcb, ca, cb, pintegrals.data(), integrals);
}

void OneElectronIntegrals::spher(const vec3& posa, int la, const row<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                                 const vec3& posb, int lb, const row<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
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
    transform(ma*mb, sa, sb, integrals, cintegrals.data());
    transpose(fsa*fsb, ma*mb, 1.0, cintegrals.data(), fsa*fsb,
                              0.0,  integrals       ,   ma*mb);
}

void OneElectronIntegrals::so(const Shell& a, const Shell& b, double* integrals)
{
    auto& posa = a.getCenter();
    auto& posb = b.getCenter();
    auto& group = posa.getPointGroup();
    int la = a.getL();
    int lb = b.getL();
    auto& za = a.getExponents();
    auto& zb = b.getExponents();
    auto& ca = a.getCoefficients();
    auto& cb = b.getCoefficients();
    auto& sa = a.getCart2Spher();
    auto& sb = b.getCart2Spher();

    int na = ca.length(0);
    int nb = cb.length(0);
    int ma = ca.length(1);
    int mb = cb.length(1);
    int fca = sa.length(0);
    int fcb = sb.length(0);
    int fsa = sa.length(1);
    int fsb = sb.length(1);
    int da = a.getDegeneracy();
    int db = b.getDegeneracy();

    vector<double> aointegrals(fca*fcb*na*nb);

    int lambdar;
    vector<int> dcrr = group.DCR(posa.getStabilizer(), posb.getStabilizer(), lambdar);
    double coef = (double)group.getOrder()/(double)lambdar;

    for (int r : dcrr)
    {
        spher(posa.getCenterAfterOp(0), la, za, ca, sa,
              posb.getCenterAfterOp(r), lb, zb, cb, sb,
              aointegrals.data());
        scal(aointegrals.size(), coef, aointegrals.data(), 1);

        double* aoints = aointegrals.data();
        double* soints = integrals;

        for (int j = 0;j < fsb;j++)
        {
            for (int i = 0;i < fsa;i++)
            {
                for (int f = 0;f < db;f++)
                {
                    for (int e = 0;e < da;e++)
                    {
                        int w = a.getIrrepOfFunc(i,e);
                        int x = b.getIrrepOfFunc(j,f);

                        if (!(group.getIrrep(w)*group.getIrrep(x)).isTotallySymmetric()) continue;

                        double fac = b.getParity(j,r)*group.character(x,r);
                        axpy(ma*mb, fac, aoints, 1, soints, 1);

                        soints += ma*mb;
                    }
                }
                aoints += ma*mb;
            }
        }
    }
}

}
}
