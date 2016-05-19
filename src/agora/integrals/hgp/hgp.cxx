#include "util/global.hpp"

namespace aquarius
{
namespace integrals
{

using slice::all;

class HGP
{
    /*
     * Redistribute angular momentum: (e=a...a+b,0|X) -> (a,b|X)
     */
    void hrr(int la, int lb,
             const vec3& posa, const vec3& posb,
             double* integrals1, double* integrals2)
    {
        assert(la >= lb);

        double* ints1 = integrals1;
        double* ints2 = integrals2;

        for (int la1 = la+lb;la1 > la;la1--)
        {
            hrr1(la, la1, la+lb-la1, posa, posb, ints1, ints2);
            swap(ints1, ints2);
        }

        if (ints1 != integrals2)
        {
            size_t sz = ((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2);
            copy(sz, ints1, 1, integrals2, 1);
        }
    }

    /*
     * Transfer angular momentum from A to B:
     *
     * (a,b+1_l|X) = (a+1_l,b|X) + (A-B)_l (a,b|X)
     *
     * where l = x,y,z
     */
    void hrr1(int la0, int la1, int lb,
              const vec3& posa, const vec3& posb,
              double* restrict integrals1, double* restrict integrals2)
    {
        vec3 fac = posa-posb;

        size_t inca1 = 1;
        size_t incb1 = inca1*((la1+1)*(la1+2)*(la1+3)/6-la0*(la0+1)*(la0+2)/6);

        size_t inca2 = 1;
        size_t incb2 = inca2*(la1*(la1+1)*(la1+2)/6-la0*(la0+1)*(la0+2)/6);

        double *ap1b, *abp1;
        double *ab = integrals1+incb1*(lb+1)*(lb+2)/2-1-inca1*(la1+1)*(la1+2)/2;
        for (int la = la1-1;la >= la0;la--)
        {
            for (int ax = 0;ax <= la;ax++)
            {
                for (int ay = 0;ay <= la-ax;ay++)
                {
                    for (int bx = 0;bx <= lb;bx++)
                    {
                        for (int by = 0;by <= lb-bx;by++)
                        {
                            ap1b = ab+inca1*(la+1)*(la+2)/2+inca1*ax;
                            abp1 = ab+incb2*(lb+1)*(lb+2)/2+incb2*bx;
                            *abp1 = *ap1b + fac[2]*(*ab);
                            ab -= incb1;
                        }
                        ap1b = ab+inca1*(la+1)*(la+2)/2+inca1*(ax+1);
                        abp1 = ab+incb2*(lb+1)*(lb+2)/2+incb2*(bx+1);
                        *abp1 = *ap1b + fac[1]*(*ab);
                    }
                    ap1b = ab+inca1*(la+1)*(la+2)/2+inca1*(la+2);
                    abp1 = ab+incb2*(lb+1)*(lb+2)/2+incb2*(lb+2);
                    *abp1 = *ap1b + fac[0]*(*ab);
                    ab -= inca1-incb1*(lb+1)*(lb+2)/2;
                }
            }
        }
    }

    void vrr(const Shell& a, const Shell& b, const Shell& c, const Shell& d,
             const vec3& posa, const vec3& posb, const vec3& posc, const vec3& posd,
             double* restrict integrals)
    {
        constexpr double PI_52 = 17.493418327624862846262821679872; // pi^(5/2)

        int na = a.getNPrim();
        int nb = b.getNPrim();
        int nc = c.getNPrim();
        int nd = d.getNPrim();

        int la = a.getL();
        int lb = b.getL();
        int lc = c.getL();
        int ld = d.getL();

        const vector<double>& za = a.getExponents();
        const vector<double>& zb = b.getExponents();
        const vector<double>& zc = c.getExponents();
        const vector<double>& zd = d.getExponents();

        int le0 = la;
        int le1 = la+lb;
        int lf0 = lc;
        int lf1 = lc+ld;

        int vmax = le1+lf1;

        int nabcd = na*nb*nc*nd;

        marray<double,4> xtable({le1+1,lf1+1,vmax+1,nabcd});

        matrix<double> efac(3, nabcd);
        matrix<double> ffac(3, nabcd);
        matrix<double> pfac(3, nabcd);
        matrix<double> qfac(3, nabcd);
        row<double> s1fac(nabcd);
        row<double> s2fac(nabcd);
        row<double> t1fac(nabcd);
        row<double> t2fac(nabcd);
        row<double>  gfac(nabcd);

        Fm fm;
        vector<double> ssssm(vmax+1);

        for (int e = 0, i = 0;e < na;e++)
        {
            for (int f = 0;f < nb;f++)
            {
                for (int g = 0;g < nc;g++)
                {
                    for (int h = 0;h < nd;h++, i++)
                    {
                        double zp = za[e]+zb[f];
                        double zq = zc[g]+zd[h];

                        vec3 posp = (posa*za[e] + posb*zb[f])/zp;
                        vec3 posq = (posc*zc[g] + posd*zd[h])/zq;
                        vec3 posw = (posp*zp    + posq*zq   )/(zp+zq);

                        efac[all][i] = posp - posa;
                        ffac[all][i] = posq - posc;
                        pfac[all][i] = posw - posp;
                        qfac[all][i] = posw - posq;

                        s1fac[i] = 1.0/(2*zp);
                        s2fac[i] = 1.0/(2*zq);
                         gfac[i] = 1.0/(2*(zp+zq));
                        t1fac[i] = -gfac[i]*zq/zp;
                        t2fac[i] = -gfac[i]*zp/zq;

                        double A0 = 2*PI_52*exp(-za[e]*zb[f]*norm2(posa-posb)/zp
                                                -zc[g]*zd[h]*norm2(posc-posd)/zq)/(sqrt(zp+zq)*zp*zq);
                        double Z = norm2(posp-posq)*zp*zq/(zp+zq);

                        /*
                         * Get (00|00)[m] and scale by A0
                         */
                        fm(Z, ssssm);
                        for (int v = 0;v <= vmax;v++) xtable[0][0][v][i] = A0*ssssm[v];
                    }
                }
            }
        }

        // fill table with x
        filltable(efac[0], ffac[0], pfac[0], qfac[0],
                  s1fac, t1fac, s2fac, t2fac, gfac,
                  xtable);

        // loop over all possible distributions of x momenta
        for (int fx = lf1;fx >= 0;fx--)
        {
            for (int ex = le1;ex >= 0;ex--)
            {
                // and fill remainder with y from that point
                filltable(efac[1], ffac[1], pfac[1], qfac[1],
                          s1fac, t1fac, s2fac, t2fac, gfac,
                          xtable[range(ex,le1+1)][range(fx,lf1+1)]);

                // loop over all possible distirubtions of y momenta given x
                for (int fy = lf1-fx;fy >= 0;fy--)
                {
                    for (int ey = le1-ex;ey >= 0;ey--)
                    {
                        // and fill remainder with z from that point
                        filltable(efac[2], ffac[2], pfac[2], qfac[2],
                                  s1fac, t1fac, s2fac, t2fac, gfac,
                                  xtable[range(ex+ey,le1+1)][range(fx+fy,lf1+1)]);

                        int fxyz = ???;
                        for (int fz = lf1-fx-fy;fz >= max(0,lf0-fx-fy);fz--)
                        {
                            int exyz = ???;
                            for (int ez = le1-ex-ey;ez >= max(0,le0-ex-ey);ez--)
                            {
                                integral[exyz][0][fxyz][0] = xtable[ex+ey+ez][fx+fy+fz][0];
                            }
                        }
                    }
                }
            }
        }
    }

    void filltable(marray<double,1>&& efac, marray<double,1>&& ffac,
                   marray<double,1>&& pfac, marray<double,1>&& qfac,
                   const row<double>& s1fac, const row<double>& t1fac,
                   const row<double>& s2fac, const row<double>& t2fac,
                   const row<double>& gfac,
                   marray<double,4>&& table);
    {
        filltable(move(efac), move(ffac), move(pfac), move(qfac),
                  s1fac, t1fac, s2fac, t2fac, gfac, table);
    }

    void filltable(marray<double,1>&& efac, marray<double,1>&& ffac,
                   marray<double,1>&& pfac, marray<double,1>&& qfac,
                   const row<double>& s1fac, const row<double>& t1fac,
                   const row<double>& s2fac, const row<double>& t2fac,
                   const row<double>& gfac,
                   marray<double,4>& table);
    {
        int le = table.length(0)-1;
        int lf = table.length(1)-1;
        int np = table.length(2);
        int vmax = le+lf;

        int e, f;

        f = 0;
        if (f < lf)
        {
            for (int v = 0;v < vmax-f;v++)
            {
                for (int p = 0;p < np;p++)
                {
                    table[e][f+1][v][p] = ffac[p]*table[e][f][  v][p] +
                                          qfac[p]*table[e][f][v+1][p];
                }
            }
        }

        e = 0;
        if (e < le)
        {
            for (int v = 0;v < vmax-e-f;v++)
            {
                for (int p = 0;p < np;p++)
                {
                    table[e+1][f][v][p] = efac[p]*table[e][f][  v][p] +
                                          pfac[p]*table[e][f][v+1][p];
                }
            }
        }

        for (e = 1;e <= le;e++)
        {
            if (e < le)
            {
                for (int v = 0;v < vmax-e-f;v++)
                {
                    for (int p = 0;p < np;p++)
                    {
                        table[e+1][f][v][p] =     efac[p]*table[  e][f][  v][p] +
                                                  pfac[p]*table[  e][f][v+1][p] +
                                              e*(s1fac[p]*table[e-1][f][  v][p] +
                                                 t1fac[p]*table[e-1][f][v+1][p]);
                    }
                }
            }
        }

        for (f = 1;f <= lf;f++)
        {
            if (f < lf)
            {
                for (int v = 0;v < vmax-f;v++)
                {
                    for (int p = 0;p < np;p++)
                    {
                        table[e][f+1][v][p] =     ffac[p]*table[e][  f][  v][p] +
                                                  qfac[p]*table[e][  f][v+1][p] +
                                              f*(s2fac[p]*table[e][f-1][  v][p] +
                                                 t2fac[p]*table[e][f-1][v+1][p]);
                    }
                }
            }

            e = 0;
            if (e < le)
            {
                for (int v = 0;v < vmax-e-f;v++)
                {
                    for (int p = 0;p < np;p++)
                    {
                        table[e+1][f][v][p] =   efac[p]*table[e][  f][  v][p] +
                                                pfac[p]*table[e][  f][v+1][p] +
                                              f*gfac[p]*table[e][f-1][  v][p];
                    }
                }
            }

            for (e = 1;e <= le;e++)
            {
                if (e < le)
                {
                    for (int v = 0;v < vmax-e-f;v++)
                    {
                        for (int p = 0;p < np;p++)
                        {
                            table[e+1][f][v][p] =     efac[p]*table[  e][  f][  v][p] +
                                                      pfac[p]*table[  e][  f][v+1][p] +
                                                  f*  gfac[p]*table[  e][f-1][  v][p] +
                                                  e*(s1fac[p]*table[e-1][  f][  v][p] +
                                                     t1fac[p]*table[e-1][  f][v+1][p]);
                        }
                    }
                }
            }
        }
    }
}

}
}
