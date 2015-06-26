#include "libint2eints.hpp"

#include "fmgamma.hpp"

using namespace aquarius::input;
using namespace aquarius::symmetry;
using namespace aquarius::task;

namespace aquarius
{
namespace integrals
{

Libint2eIntegrals::Libint2eIntegrals(const Shell& a, const Shell& b, const Shell& c, const Shell& d)
: TwoElectronIntegrals(a, b, c, d)
{
    int nt = omp_get_max_threads();
    inteval.resize(nt);
    for (int i = 0;i < nt;i++)
    {
        libint2_init_eri(&inteval[i], max(max(la,lb),max(lc,ld)), NULL);
        inteval[i].contrdepth = 1;
        assert(LIBINT2_MAX_VECLEN == 1);
    }
}

void Libint2eIntegrals::prims(const vec3& posa, const vec3& posb, const vec3& posc, const vec3& posd,
                              double* integrals)
{
    constexpr double TWO_PI_52 = 34.98683665524972497; // 2*pi^(5/2)
    Fm fm;

    matrix<double> Kab(na, nb), Kcd(nc, nd);

    accuracy_ = 1e-13;

    #pragma omp parallel
    {
        #pragma omp for
        for (int64_t j = 0;j < na*nb;j++)
        {
            int f = j/na;
            int e = j%na;

            double zp = za[e]+zb[f];
            Kab[e][f] = exp(-za[e]*zb[f]*norm2(posa-posb)/zp)/zp;
        }

        #pragma omp for
        for (int64_t j = 0;j < nc*nd;j++)
        {
            int h = j/nc;
            int g = j%nc;

            double zq = zc[g]+zd[h];
            Kcd[g][h] = exp(-zc[g]*zd[h]*norm2(posc-posd)/zq)/zq;
        }

        int nt = omp_get_num_threads();
        int i = omp_get_thread_num();
        int len = fca*fcb*fcc*fcd;
        for (int64_t j = i;j < len;j += nt)
        {
            int h = j/(na*nb*nc);
            int r = j%(na*nb*nc);
            int g = r/(na*nb);
            int s = r%(na*nb);
            int f = s/na;
            int e = s%na;

            double zp = za[e]+zb[f];
            double zq = zc[g]+zd[h];

            vec3 posp = (posa*za[e] + posb*zb[f])/zp;
            vec3 posq = (posc*zc[g] + posd*zd[h])/zq;
            vec3 posw = (posp*zp    + posq*zq   )/(zp+zq);

            double A0 = TWO_PI_52*Kab[e][f]*Kcd[g][h]/sqrt(zp+zq);
            double Z = norm2(posp-posq)*zp*zq/(zp+zq);

            if (Kab[e][f] < accuracy_ ||
                Kcd[g][h] < accuracy_ ||
                A0 < accuracy_)
            {
                fill_n(integrals+j*len, len, 0.0);
                continue;
            }

            bool swapab = la < lb;
            bool swapcd = lc < ld;
            bool swappq = la+lb > lc+ld;

            vec3 afac = posp - (swapab ? posb : posa);
            vec3 cfac = posq - (swapcd ? posd : posc);
            vec3 pfac = posw - posp;
            vec3 qfac = posw - posq;
            vec3 abfac = posa - posb;
            vec3 cdfac = posc - posd;

            if (swapab) abfac *= -1;
            if (swapcd) cdfac *= -1;
            if (swappq)
            {
                swap(afac, cfac);
                swap(pfac, qfac);
                swap(abfac, cdfac);
            }

            double s1fac = 0.5/zp;
            double s2fac = 0.5/zq;
            double gfac = 0.5/(zp+zq);
            double rho = zp*zq/(zp+zq);

            if (swappq) swap(s1fac, s2fac);

            copy_n( afac.data(), 3, inteval[i].PA_x);
            copy_n( cfac.data(), 3, inteval[i].QC_x);
            copy_n( pfac.data(), 3, inteval[i].WP_x);
            copy_n( qfac.data(), 3, inteval[i].WQ_x);
            copy_n(abfac.data(), 3, inteval[i].AB_x);
            copy_n(cdfac.data(), 3, inteval[i].CD_x);
            inteval[i].oo2z[0] = s1fac;
            inteval[i].oo2e[0] = s2fac;
            inteval[i].oo2ze[0] = gfac;
            inteval[i].roe[0] = gfac/s1fac;
            inteval[i].roz[0] = gfac/s2fac;

            fm(Z, la+lb+lc+ld, inteval[i]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0);
            scal(la+lb+lc+ld+1, A0, inteval[i]._aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0, 1);

            if (swappq)
                if (swapab)
                    if (swapcd) libint2_build_eri[ld][lc][lb][la](&inteval[i]);
                    else        libint2_build_eri[lc][ld][lb][la](&inteval[i]);
                else
                    if (swapcd) libint2_build_eri[ld][lc][la][lb](&inteval[i]);
                    else        libint2_build_eri[lc][ld][la][lb](&inteval[i]);
            else
                if (swapab)
                    if (swapcd) libint2_build_eri[lb][la][ld][lc](&inteval[i]);
                    else        libint2_build_eri[lb][la][lc][ld](&inteval[i]);
                else
                    if (swapcd) libint2_build_eri[la][lb][ld][lc](&inteval[i]);
                    else        libint2_build_eri[la][lb][lc][ld](&inteval[i]);

            array<unsigned, 4> perm = {3,2,1,0};
            if (swappq)
            {
                swap(perm[0], perm[2]);
                swap(perm[1], perm[3]);
            }
            if (swapab) swap(perm[2], perm[3]);
            if (swapcd) swap(perm[0], perm[1]);

            array<marray<double, 4>::idx_type, 4> to_len = {fcd, fcc, fcb, fca};
            array<marray<double, 4>::idx_type, 4> from_len;
            from_len[perm[0]] = to_len[0];
            from_len[perm[1]] = to_len[1];
            from_len[perm[2]] = to_len[2];
            from_len[perm[3]] = to_len[3];

            marray<double, 4> from(from_len, inteval[i].targets[0]);
            marray<double, 4> to(to_len, integrals+j*len);
            copy(from.permute(perm), to);

            PROFILE_FLOPS(inteval[i].nflops[0]);
        }
    }
}

}
}

static const char* spec = R"(

storage_cutoff?
    double 1e-14,
calc_cutoff?
    double 1e-15

)";

REGISTER_TASK(aquarius::integrals::Libint2eIntegralsTask,"libint2eints",spec);
