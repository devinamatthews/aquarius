#include "2eints.hpp"

using namespace aquarius::molecule;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace integrals
{

void transform(const matrix<double>& ai, const matrix<double>& bj,
               const matrix<double>& ck, const matrix<double>& dl,
               size_t nother, double* buf1, double* buf2)
{
    using transpose::T;

    int na = ai.length(0);
    int nb = bj.length(0);
    int nc = ck.length(0);
    int nd = dl.length(0);
    int ma = ai.length(1);
    int mb = bj.length(1);
    int mc = ck.length(1);
    int md = dl.length(1);

    marray_view<double,5> xabcd(nother, na, nb, nc, nd, buf1);
    marray_view<double,5> lxabc(md, nother, na, nc, nc, buf2);
    marray_view<double,5> klxab(mc, md, nother, na, nb, buf1);
    marray_view<double,5> jklxa(mb, mc, md, nother, na, buf2);
    marray_view<double,5> ijklx(ma, mb, mc, md, nother, buf1);

    gemm(1.0, dl^T, xabcd.lower(4)^T, 0.0, lxabc.lower(1));
    gemm(1.0, ck^T, lxabc.lower(4)^T, 0.0, klxab.lower(1));
    gemm(1.0, bj^T, klxab.lower(4)^T, 0.0, jklxa.lower(1));
    gemm(1.0, ai^T, jklxa.lower(4)^T, 0.0, ijklx.lower(1));

    copy(ma*mb*mc*md*nother, buf1, 1, buf2, 1);
}

void transform(size_t nother,
               const matrix<double>& ai, const matrix<double>& bj,
               const matrix<double>& ck, const matrix<double>& dl,
               double* buf1, double* buf2)
{
    using transpose::T;

    int na = ai.length(0);
    int nb = bj.length(0);
    int nc = ck.length(0);
    int nd = dl.length(0);
    int ma = ai.length(1);
    int mb = bj.length(1);
    int mc = ck.length(1);
    int md = dl.length(1);

    marray_view<double,5> abcdx(na, nb, nc, nd, nother, buf1);
    marray_view<double,5> bcdxi(nb, nc, nd, nother, ma, buf1);
    marray_view<double,5> cdxij(nc, nd, nother, ma, mb, buf1);
    marray_view<double,5> dxijk(nd, nother, ma, mb, mc, buf1);
    marray_view<double,5> xijkl(nother, ma, mb, mc, md, buf1);

    gemm(1.0, abcdx.lower(1)^T, ai, 0.0, bcdxi.lower(4));
    gemm(1.0, bcdxi.lower(1)^T, bj, 0.0, cdxij.lower(4));
    gemm(1.0, cdxij.lower(1)^T, ck, 0.0, dxijk.lower(4));
    gemm(1.0, dxijk.lower(1)^T, dl, 0.0, xijkl.lower(4));

    copy(ma*mb*mc*md*nother, buf1, 1, buf2, 1);
}

size_t TwoElectronIntegrals::ShellBlock::process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                                                 const vector<int>& idxc, const vector<int>& idxd,
                                                 size_t nprocess, double* integrals, idx4_t* indices, double cutoff)
{
    auto& group = a.getCenter().getPointGroup();

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

                                                    if (!bad && std::abs(ints[m]) > cutoff)
                                                    {
                                                        indices[n].i = ctx.getIndex(a, idxa, i, e, r);
                                                        indices[n].j = ctx.getIndex(b, idxb, j, f, s);
                                                        indices[n].k = ctx.getIndex(c, idxc, k, g, t);
                                                        indices[n].l = ctx.getIndex(d, idxd, l, h, u);
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

TwoElectronIntegrals::ShellBlock TwoElectronIntegrals::calculate(const Shell& a, const Shell& b,
                                                                 const Shell& c, const Shell& d)
{
    auto& group = a.getCenter().getPointGroup();

    int nints = 0;
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
                        auto& z = group.getIrrep(d.getIrrepOfFunc(l, u));

                        for (int t = 0;t < c.getDegeneracy();t++)
                        {
                            auto& y = group.getIrrep(c.getIrrepOfFunc(k, t));

                            for (int s = 0;s < b.getDegeneracy();s++)
                            {
                                auto& x = group.getIrrep(b.getIrrepOfFunc(j, s));

                                for (int r = 0;r < a.getDegeneracy();r++)
                                {
                                    auto& w = group.getIrrep(a.getIrrepOfFunc(i, r));

                                    if ((w*x*y*z).isTotallySymmetric()) nints +=
                                        a.getNContr()*b.getNContr()*
                                        c.getNContr()*d.getNContr();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    vector<double> ints(nints);
    so(a, b, c, d, ints.data());

    return {a, b, c, d, Context::SFIC, move(ints)};
}

void TwoElectronIntegrals::prim(const vec3& posa, int la, double za,
                                const vec3& posb, int lb, double zb,
                                const vec3& posc, int lc, double zc,
                                const vec3& posd, int ld, double zd,
                                double* integrals)
{
    assert(0);
}

void TwoElectronIntegrals::prims(const vec3& posa, int la, const row<double>& za,
                                 const vec3& posb, int lb, const row<double>& zb,
                                 const vec3& posc, int lc, const row<double>& zc,
                                 const vec3& posd, int ld, const row<double>& zd,
                                 double* integrals)
{
    int na = za.length();
    int nb = zb.length();
    int nc = zc.length();
    int nd = zd.length();
    int fca = (la+1)*(la+2)/2;
    int fcb = (lb+1)*(lb+2)/2;
    int fcc = (lc+1)*(lc+2)/2;
    int fcd = (ld+1)*(ld+2)/2;

    #pragma omp parallel for
    for (int m = 0;m < na*nb*nc*nd;m++)
    {
        int e = m%na;
        int x = m/na;
        int f = x%nb;
        int y = x/nb;
        int g = y%nc;
        int z = y/nc;
        int h = z%nd;

        prim(posa, la, za[e],
             posb, la, zb[f],
             posc, lc, zc[g],
             posd, ld, zd[h],
             integrals+fca*fcb*fcc*fcd*m);
    }
}

void TwoElectronIntegrals::contr(const vec3& posa, int la, const row<double>& za, const matrix<double>& ca,
                                 const vec3& posb, int lb, const row<double>& zb, const matrix<double>& cb,
                                 const vec3& posc, int lc, const row<double>& zc, const matrix<double>& cc,
                                 const vec3& posd, int ld, const row<double>& zd, const matrix<double>& cd,
                                 double* integrals)
{
    int na = za.length();
    int nb = zb.length();
    int nc = zc.length();
    int nd = zd.length();
    int fca = (la+1)*(la+2)/2;
    int fcb = (lb+1)*(lb+2)/2;
    int fcc = (lc+1)*(lc+2)/2;
    int fcd = (ld+1)*(ld+2)/2;

    vector<double> pintegrals(fca*fcb*fcc*fcd*na*nb*nc*nd);
    prims(posa, la, za,
          posb, lb, zb,
          posc, lc, zc,
          posd, ld, zd,
          pintegrals.data());
    transform(fca*fcb*fcc*fcd, ca, cb, cc, cd, pintegrals.data(), integrals);
}

void TwoElectronIntegrals::spher(const vec3& posa, int la, const row<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                                 const vec3& posb, int lb, const row<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
                                 const vec3& posc, int lc, const row<double>& zc, const matrix<double>& cc, const matrix<double>& sc,
                                 const vec3& posd, int ld, const row<double>& zd, const matrix<double>& cd, const matrix<double>& sd,
                                 double* integrals)
{
    int na = ca.length(0);
    int nb = cb.length(0);
    int nc = cc.length(0);
    int nd = cd.length(0);
    int ma = ca.length(1);
    int mb = cb.length(1);
    int mc = cc.length(1);
    int md = cd.length(1);
    int fca = sa.length(0);
    int fcb = sb.length(0);
    int fcc = sc.length(0);
    int fcd = sd.length(0);
    int fsa = sa.length(1);
    int fsb = sb.length(1);
    int fsc = sc.length(1);
    int fsd = sd.length(1);

    vector<double> cintegrals(fca*fcb*fca*fcd*na*nb*nc*nd);
    contr(posa, la, za, ca,
          posb, lb, zb, cb,
          posc, lc, zc, cc,
          posd, ld, zd, cd,
          integrals);
    transform(ma*mb*mc*md, sa, sb, sc, sd, integrals, cintegrals.data());
    transpose(fsa*fsb*fsc*fsd, ma*mb*mc*md, 1.0, cintegrals.data(), fsa*fsb*fsc*fsd,
                                            0.0,  integrals       ,     ma*mb*mc*md);
}

void TwoElectronIntegrals::so(const Shell& a, const Shell& b, const Shell& c, const Shell& d,
                              double* integrals)
{
    auto& posa = a.getCenter();
    auto& posb = b.getCenter();
    auto& posc = c.getCenter();
    auto& posd = d.getCenter();
    auto& group = posa.getPointGroup();
    int la = a.getL();
    int lb = b.getL();
    int lc = c.getL();
    int ld = d.getL();
    auto& za = a.getExponents();
    auto& zb = b.getExponents();
    auto& zc = c.getExponents();
    auto& zd = d.getExponents();
    auto& ca = a.getCoefficients();
    auto& cb = b.getCoefficients();
    auto& cc = c.getCoefficients();
    auto& cd = d.getCoefficients();
    auto& sa = a.getCart2Spher();
    auto& sb = b.getCart2Spher();
    auto& sc = c.getCart2Spher();
    auto& sd = d.getCart2Spher();

    int na = ca.length(0);
    int nb = cb.length(0);
    int nc = cc.length(0);
    int nd = cd.length(0);
    int ma = ca.length(1);
    int mb = cb.length(1);
    int mc = cc.length(1);
    int md = cd.length(1);
    int fca = sa.length(0);
    int fcb = sb.length(0);
    int fcc = sc.length(0);
    int fcd = sd.length(0);
    int fsa = sa.length(1);
    int fsb = sb.length(1);
    int fsc = sc.length(1);
    int fsd = sd.length(1);
    int da = a.getDegeneracy();
    int db = b.getDegeneracy();
    int dc = c.getDegeneracy();
    int dd = d.getDegeneracy();

    vector<double> aointegrals(fca*fcb*fcc*fcd*na*nb*nc*nd);

    int lambdar, lambdas, lambdat;
    vector<int> dcrr = group.DCR(posa.getStabilizer(), posb.getStabilizer(), lambdar);
    vector<int> dcrs = group.DCR(posc.getStabilizer(), posd.getStabilizer(), lambdas);
    vector<int> dcrt = group.DCR(intersection(posa.getStabilizer(), posb.getStabilizer()),
                                 intersection(posc.getStabilizer(), posd.getStabilizer()), lambdat);
    double coef = (double)group.getOrder()/(double)lambdat;

    Representation yz = group.getIrrep(0);
    Representation xyz = group.getIrrep(0);
    Representation wxyz = group.getIrrep(0);

    for (int r : dcrr)
    {
        for (int s : dcrs)
        {
            for (int t : dcrt)
            {
                int st = group.getOpProduct(s,t);

                spher(posa.getCenterAfterOp( 0), la, za, ca, sa,
                      posb.getCenterAfterOp( r), lb, zb, cb, sb,
                      posc.getCenterAfterOp( t), lc, zc, cc, sc,
                      posd.getCenterAfterOp(st), ld, zd, cd, sd,
                      aointegrals.data());
                scal(aointegrals.size(), coef, aointegrals.data(), 1);

                double* aoints = aointegrals.data();
                double* soints = integrals;

                for (int l = 0;l < fsd;l++)
                {
                    for (int k = 0;k < fsc;k++)
                    {
                        for (int j = 0;j < fsb;j++)
                        {
                            for (int i = 0;i < fsa;i++)
                            {
                                for (int h = 0;h < dd;h++)
                                {
                                    int z = d.getIrrepOfFunc(l,h);
                                    for (int g = 0;g < dc;g++)
                                    {
                                        int y = c.getIrrepOfFunc(k,g);
                                        yz = group.getIrrep(z);
                                        yz *= group.getIrrep(y);
                                        for (int f = 0;f < db;f++)
                                        {
                                            int x = b.getIrrepOfFunc(j,f);
                                            xyz = yz;
                                            xyz *= group.getIrrep(x);
                                            for (int e = 0;e < da;e++)
                                            {
                                                int w = a.getIrrepOfFunc(i,e);
                                                wxyz = xyz;
                                                wxyz *= group.getIrrep(w);

                                                if (!wxyz.isTotallySymmetric()) continue;

                                                double fac = b.getParity(j, r)*group.character(x, r)*
                                                             c.getParity(k, t)*group.character(y, t)*
                                                             d.getParity(l,st)*group.character(z,st);
                                                axpy(ma*mb*mc*md, fac, aoints, 1, soints, 1);

                                                soints += ma*mb*mc*md;
                                            }
                                        }
                                    }
                                }
                                aoints += ma*mb*mc*md;
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
