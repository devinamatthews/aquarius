#include "2eints.hpp"
#include "os.hpp"

using namespace aquarius::input;
using namespace aquarius::symmetry;
using namespace aquarius::task;

namespace aquarius
{
namespace integrals
{

TwoElectronIntegrals::TwoElectronIntegrals(const Shell& a, const Shell& b, const Shell& c, const Shell& d)
: sa(a), sb(b), sc(c), sd(d), group(a.getCenter().getPointGroup()),
  ca(a.getCenter()), cb(b.getCenter()), cc(c.getCenter()), cd(d.getCenter()),
  la(a.getL()), lb(b.getL()), lc(c.getL()), ld(d.getL()),
  na(a.getNPrim()), nb(b.getNPrim()), nc(c.getNPrim()), nd(d.getNPrim()),
  ma(a.getNContr()), mb(b.getNContr()), mc(c.getNContr()), md(d.getNContr()),
  da(a.getDegeneracy()), db(b.getDegeneracy()), dc(c.getDegeneracy()), dd(d.getDegeneracy()),
  fsa(a.getNFunc()), fsb(b.getNFunc()), fsc(c.getNFunc()), fsd(d.getNFunc()),
  za(a.getExponents()), zb(b.getExponents()), zc(c.getExponents()), zd(d.getExponents()),
  num_processed(0), accuracy_(0)
{
    fca = (la+1)*(la+2)/2;
    fcb = (lb+1)*(lb+2)/2;
    fcc = (lc+1)*(lc+2)/2;
    fcd = (ld+1)*(ld+2)/2;

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

                                    if ((w*x*y*z).isTotallySymmetric()) nints += ma*mb*mc*md;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    ints.resize(nints);
}

void TwoElectronIntegrals::run()
{
    so(ints.data());
}

size_t TwoElectronIntegrals::process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                                     const vector<int>& idxc, const vector<int>& idxd,
                                     size_t nprocess, double* integrals, idx4_t* indices, double cutoff)
{
    Representation z = group.getIrrep(0);
    Representation yz = group.getIrrep(0);
    Representation xyz = group.getIrrep(0);
    Representation wxyz = group.getIrrep(0);

    size_t m = 0;
    size_t n = 0;
    for (int l = 0;l < fsd;l++)
    {
        for (int k = 0;k < fsc;k++)
        {
            for (int j = 0;j < fsb;j++)
            {
                for (int i = 0;i < fsa;i++)
                {
                    for (int u = 0;u < dd;u++)
                    {
                        z = group.getIrrep(sd.getIrrepOfFunc(l,u));
                        for (int t = 0;t < dc;t++)
                        {
                            yz = z;
                            yz *= group.getIrrep(sc.getIrrepOfFunc(k,t));
                            for (int s = 0;s < db;s++)
                            {
                                xyz = yz;
                                xyz *= group.getIrrep(sb.getIrrepOfFunc(j,s));
                                for (int r = 0;r < da;r++)
                                {
                                    wxyz = xyz;
                                    wxyz *= group.getIrrep(sa.getIrrepOfFunc(i,r));

                                    if (!wxyz.isTotallySymmetric()) continue;

                                    for (int h = 0;h < md;h++)
                                    {
                                        for (int g = 0;g < mc;g++)
                                        {
                                            for (int f = 0;f < mb;f++)
                                            {
                                                for (int e = 0;e < ma;e++)
                                                {
                                                    if (num_processed > m)
                                                    {
                                                        m++;
                                                        continue;
                                                    }

                                                    bool bad = false;

                                                    if (&sa == &sb && !IDX_GE(i,r,e,j,s,f)) bad = true;
                                                    if (&sc == &sd && !IDX_GE(k,t,g,l,u,h)) bad = true;
                                                    if (&sa == &sc && &sb == &sd && !(IDX_GT(i,r,e,k,t,g) ||
                                                        (IDX_EQ(i,r,e,k,t,g) && IDX_GE(j,s,f,l,u,h)))) bad = true;

                                                    if (!bad && aquarius::abs(ints[m]) > cutoff)
                                                    {
                                                        indices[n].i = sa.getIndex(ctx, idxa, i, e, r);
                                                        indices[n].j = sb.getIndex(ctx, idxb, j, f, s);
                                                        indices[n].k = sc.getIndex(ctx, idxc, k, g, t);
                                                        indices[n].l = sd.getIndex(ctx, idxd, l, h, u);
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

void TwoElectronIntegrals::prim(const vec3& posa, int e, const vec3& posb, int f,
                                const vec3& posc, int g, const vec3& posd, int h,
                                double* integrals)
{
    assert(0);
}

void TwoElectronIntegrals::prims(const vec3& posa, const vec3& posb, const vec3& posc, const vec3& posd,
                                 double* integrals)
{
    constexpr double TWO_PI_52 = 34.98683665524972497; // 2*pi^(5/2)

    matrix<double> Kab(na, nb), Kcd(nc, nd);

    #pragma omp parallel for
    for (int64_t j = 0;j < na*nb;j++)
    {
        int f = j/na;
        int e = j%na;

        double zp = za[e]+zb[f];
        Kab[e][f] = exp(-za[e]*zb[f]*norm2(posa-posb)/zp)/zp;
    }

    #pragma omp parallel for
    for (int64_t j = 0;j < nc*nd;j++)
    {
        int h = j/nc;
        int g = j%nc;

        double zq = zc[g]+zd[h];
        Kcd[g][h] = exp(-zc[g]*zd[h]*norm2(posc-posd)/zq)/zq;
    }

    #pragma omp parallel
    {
        int nt = omp_get_num_threads();
        int i = omp_get_thread_num();
        int len = fca*fcb*fcc*fcd;
        for (int64_t j = i;j < na*nb*nc*nd;j += nt)
        {
            int h = j/(na*nb*nc);
            int r = j%(na*nb*nc);
            int g = r/(na*nb);
            int s = r%(na*nb);
            int f = s/na;
            int e = s%na;

            double A0 = TWO_PI_52*Kab[e][f]*Kcd[g][h]/sqrt(za[e]+zb[f]+zc[g]+zd[h]);

            if (Kab[e][f] < accuracy_ ||
                Kcd[g][h] < accuracy_ ||
                A0 < accuracy_)
            {
                fill_n(integrals+j*len, len, 0.0);
                continue;
            }

            prim(posa, e, posb, f, posc, g, posd, h, integrals+j*len);
        }
    }
}

void TwoElectronIntegrals::contr(const vec3& posa, const vec3& posb, const vec3& posc, const vec3& posd,
                                 double* integrals)
{
    vector<double> pintegrals(fca*fcb*fcc*fcd*na*nb*nc*nd);
    prims(posa, posb, posc, posd, pintegrals.data());
    prim2contr4r(fca*fcb*fcc*fcd, pintegrals.data(), integrals);
}

void TwoElectronIntegrals::spher(const vec3& posa, const vec3& posb, const vec3& posc, const vec3& posd,
                                 double* integrals)
{
    vector<double> cintegrals(fca*fcb*fcc*fcd*na*nb*nc*nd);
    contr(posa, posb, posc, posd, integrals);
    cart2spher4r(ma*mb*mc*md, integrals, cintegrals.data());
    transpose(fsa*fsb*fsc*fsd, ma*mb*mc*md, 1.0, cintegrals.data(), fsa*fsb*fsc*fsd,
                                            0.0,  integrals       ,     ma*mb*mc*md);
}

void TwoElectronIntegrals::so(double* integrals)
{
    vector<double> aointegrals(fca*fcb*fcc*fcd*na*nb*nc*nd);

    int lambdar, lambdas, lambdat;
    vector<int> dcrr = group.DCR(ca.getStabilizer(), cb.getStabilizer(), lambdar);
    vector<int> dcrs = group.DCR(cc.getStabilizer(), cd.getStabilizer(), lambdas);
    vector<int> dcrt = group.DCR(intersection(ca.getStabilizer(), cb.getStabilizer()),
                                 intersection(cc.getStabilizer(), cd.getStabilizer()), lambdat);
    double coef = (double)group.getOrder()/(double)lambdat;

    for (int r : dcrr)
    {
        for (int s : dcrs)
        {
            for (int t : dcrt)
            {
                int st = group.getOpProduct(s,t);

                spher(ca.getCenter(0),
                      cb.getCenter(cb.getCenterAfterOp(r)),
                      cc.getCenter(cc.getCenterAfterOp(t)),
                      cd.getCenter(cd.getCenterAfterOp(st)),
                      aointegrals.data());
                scal(aointegrals.size(), coef, aointegrals.data(), 1);
                ao2so4(ma*mb*mc*md, r, t, st, aointegrals.data(), ints.data());
            }
        }
    }
}

void TwoElectronIntegrals::ao2so4(size_t nother, int r, int t, int st, double* aointegrals, double* sointegrals)
{
    Representation yz = group.getIrrep(0);
    Representation xyz = group.getIrrep(0);
    Representation wxyz = group.getIrrep(0);

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
                        int z = sd.getIrrepOfFunc(l,h);
                        for (int g = 0;g < dc;g++)
                        {
                            int y = sc.getIrrepOfFunc(k,g);
                            yz = group.getIrrep(z);
                            yz *= group.getIrrep(y);
                            for (int f = 0;f < db;f++)
                            {
                                int x = sb.getIrrepOfFunc(j,f);
                                xyz = yz;
                                xyz *= group.getIrrep(x);
                                for (int e = 0;e < da;e++)
                                {
                                    int w = sa.getIrrepOfFunc(i,e);
                                    wxyz = xyz;
                                    wxyz *= group.getIrrep(w);

                                    if (!wxyz.isTotallySymmetric()) continue;

                                    double fac = sb.getParity(j,r) *group.character(x,r)*
                                                 sc.getParity(k,t) *group.character(y,t)*
                                                 sd.getParity(l,st)*group.character(z,st);
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

    m = fsd;
    n = fca*fcb*fcc*nother;
    k = fcd;

    if (sd.isSpherical())
    {
        // [d,l]' x [xabc,d]' = [l,xabc]
        gemm('T', 'T', m, n, k, 1.0, sd.getCart2Spher().data(), k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [xabc,l]' = [l,xabc]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = fsc;
    n = fca*fcb*fsd*nother;
    k = fcc;

    if (sc.isSpherical())
    {
        // [c,k]' x [lxab,c]' = [k,lxab]
        gemm('T', 'T', m, n, k, 1.0, sc.getCart2Spher().data(), k, buf2, n, 0.0, buf1, m);
    }
    else
    {
        // [lxab,k]' = [k,lxab]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    m = fsb;
    n = fca*fsc*fsd*nother;
    k = fcb;

    if (sb.isSpherical())
    {
        // [b,j]' x [klxa,b]' = [j,klxa]
        gemm('T', 'T', m, n, k, 1.0, sb.getCart2Spher().data(), k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [klxa,j]' = [j,klxa]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = fsa;
    n = fsb*fsc*fsd*nother;
    k = fca;

    if (sa.isSpherical())
    {
        // [a,i]' x [jklx,a]' = [i,jklx]
        gemm('T', 'T', m, n, k, 1.0, sa.getCart2Spher().data(), k, buf2, n, 0.0, buf1, m);
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

    m = fsa;
    n = fcb*fcc*fcd*nother;
    k = fca;

    if (sa.isSpherical())
    {
        // [a,bcdx]' x [a,i]' = [bcdx,i]
        gemm('T', 'N', m, n, k, 1.0, buf1, k, sd.getCart2Spher().data(), k, 0.0, buf2, m);
    }
    else
    {
        // [i,bcdx]' = [bcdx,i]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = fsb;
    n = fsa*fcc*fcd*nother;
    k = fcb;

    if (sb.isSpherical())
    {
        // [b,cdxi]' x [b,j]' = [cdxi,j]
        gemm('T', 'N', m, n, k, 1.0, buf2, k, sc.getCart2Spher().data(), k, 0.0, buf1, m);
    }
    else
    {
        // [j,cdxi]' = [cdxi,j]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    m = fsc;
    n = fsa*fsb*fcd*nother;
    k = fcc;

    if (sc.isSpherical())
    {
        // [c,dxij]' x [c,k]' = [dxij,k]
        gemm('T', 'N', m, n, k, 1.0, buf1, k, sb.getCart2Spher().data(), k, 0.0, buf2, m);
    }
    else
    {
        // [k,dxij]' = [dxij,k]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = fsd;
    n = fsa*fsb*fsc*nother;
    k = fcd;

    if (sd.isSpherical())
    {
        // [d,xijk]' x [d,l]' = [xijk,l]
        gemm('T', 'N', m, n, k, 1.0, buf2, k, sa.getCart2Spher().data(), k, 0.0, buf1, m);
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
    m = md;
    n = na*nb*nc*nother;
    k = nd;
    gemm('T', 'T', m, n, k, 1.0, sd.getCoefficients().data(), k, buf1, n, 0.0, buf2, m);

    // [c,k]' x [lxab,c]' = [k,lxab]
    m = mc;
    n = na*nb*md*nother;
    k = nc;
    gemm('T', 'T', m, n, k, 1.0, sc.getCoefficients().data(), k, buf2, n, 0.0, buf1, m);

    // [b,j]' x [klxa,b]' = [j,klxa]
    m = mb;
    n = na*mc*md*nother;
    k = nb;
    gemm('T', 'T', m, n, k, 1.0, sb.getCoefficients().data(), k, buf1, n, 0.0, buf2, m);

    // [a,i]' x [jkl,xa]' = [i,jklx]
    m = ma;
    n = mb*mc*md*nother;
    k = na;
    gemm('T', 'T', m, n, k, 1.0, sa.getCoefficients().data(), k, buf2, n, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

void TwoElectronIntegrals::prim2contr4l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    // [a,bcdx]' x [a,i] = [bcdx,i]
    m = ma;
    n = nb*nc*nd*nother;
    k = na;
    gemm('T', 'N', m, n, k, 1.0, buf1, k, sa.getCoefficients().data(), k, 0.0, buf2, m);

    // [b,cdxi]' x [b,j] = [cdxi,j]
    m = mb;
    n = ma*nc*nd*nother;
    k = nb;
    gemm('T', 'N', m, n, k, 1.0, buf2, k, sb.getCoefficients().data(), k, 0.0, buf1, m);

    // [c,dxij]' x [c,k] = [dxij,k]
    m = mc;
    n = ma*mb*nd*nother;
    k = nc;
    gemm('T', 'N', m, n, k, 1.0, buf1, k, sc.getCoefficients().data(), k, 0.0, buf2, m);

    // [d,xijk]' x [d,l] = [xijk,l]
    m = md;
    n = ma*mb*mc*nother;
    k = nd;
    gemm('T', 'N', m, n, k, 1.0, buf2, k, sd.getCoefficients().data(), k, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

void ERI::print(Printer& p) const
{
    //TODO
}

}
}

static const char* spec = R"(

storage_cutoff?
    double 1e-14,
calc_cutoff?
    double 1e-15

)";

REGISTER_TASK(aquarius::integrals::OS2eIntegralsTask,"2eints",spec);
