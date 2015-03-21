#include "ishida.hpp"

namespace aquarius
{
namespace integrals
{

void IshidaERI::prim(const vec3& posa, int e, const vec3& posb, int f,
                     const vec3& posc, int g, const vec3& posd, int h, double* restrict integrals)
{
    constexpr double PI_52 = 17.493418327624862846262821679872;
    int nrys = (la+lb+lc+ld)/2 + 1;

    marray<double,6> xtable(3, ld+1, lc+1, lb+1, la+1, nrys);

    double zp = za[e] + zb[f];
    double zq = zc[g] + zd[h];

    vec3 posp = (posa*za[e] + posb*zb[f])/zp;
    vec3 posq = (posc*zc[g] + posd*zd[h])/zq;

    double A0 = 2*PI_52*exp(-za[e]*zb[f]*norm2(posa-posb)/zp
                            -zc[g]*zd[h]*norm2(posc-posd)/zq)/(sqrt(zp+zq)*zp*zq);
    double Z = norm2(posp-posq)*zp*zq/(zp+zq);

    Rys rys;
    row<double> rts(nrys), wts(nrys);
    rys(Z, nrys, rts.data(), wts.data());

    for (int xyz = 0;xyz < 3;xyz++)
    {
        double pfac = zq*(posq[xyz] - posp[xyz])/(zp+zq);
        double qfac = zp*(posp[xyz] - posq[xyz])/(zp+zq);

        double afac = posp[xyz] - posa[xyz];
        double bfac = posp[xyz] - posb[xyz];
        double cfac = posq[xyz] - posc[xyz];
        double dfac = posq[xyz] - posd[xyz];

        row<double> aafac(nrys);
        row<double> bbfac(nrys);
        row<double> ccfac(nrys);
        row<double> ddfac(nrys);

        /*
         * there is a typo in Ishida (JCP v98), the definition of G after Eq. 5 should read G = \xi s_i^2
         */
        row<double> gfac(nrys);

        for (int v = 0;v < nrys;v++)
        {
            aafac[v] = afac + pfac*rts[v];
            bbfac[v] = bfac + pfac*rts[v];
            ccfac[v] = cfac + qfac*rts[v];
            ddfac[v] = dfac + qfac*rts[v];

            /*
             * there is a typo in Ishida (JCP v98), the definition of G after Eq. 5 should read G = \xi s_i^2
             */
            gfac[v] = rts[v]/(2*(zp+zq));
        }

        double s1fac = (0.5 - zq*gfac)/zp;
        double s2fac = (0.5 - zp*gfac)/zq;

        filltable((xyz == 0 ? A0 : 1.0), aafac, bbfac, ccfac, ddfac, s1fac, s2fac, gfac, xtable[xyz]);
    }

    for (int dx = 0;dx <= ld;dx++)
    {
        for (int dy = 0;dy <= ld-dx;dy++)
        {
            int dz = ld-dx-dy;
            for (int cx = 0;cx <= lc;cx++)
            {
                for (int cy = 0;cy <= lc-cx;cy++)
                {
                    int cz = lc-cx-cy;
                    for (int bx = 0;bx <= lb;bx++)
                    {
                        for (int by = 0;by <= lb-bx;by++)
                        {
                            int bz = lb-bx-by;
                            for (int ax = 0;ax <= la;ax++)
                            {
                                for (int ay = 0;ay <= la-ax;ay++)
                                {
                                    int az = la-ax-ay;

                                    *integrals = 0.0;
                                    for (int v = 0;v < nrys;v++)
                                    {
                                        *integrals += xtable[0][dx][cx][bx][ax][v] *
                                                      xtable[1][dy][cy][by][ay][v] *
                                                      xtable[2][dz][cz][bz][az][v];
                                    }
                                    integrals++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void IshidaERI::filltable(double factor,
                          row<double>& aafac, row<double>& bbfac, row<double>& ccfac, row<double>& ddfac,
                          double s1fac, double s2fac, row<double>& gfac, marray<double,5>&& xtable)
{
    int nrys = (la+lb+lc+ld)/2 + 1;

    for (int v = 0;v < nrys;v++)
    {
        xtable[0][0][0][0][v] = factor;
        xtable[0][0][0][1][v] = aafac[v]*factor;
        xtable[0][0][1][0][v] = bbfac[v]*factor;
        xtable[0][1][0][0][v] = ccfac[v]*factor;
        xtable[1][0][0][0][v] = ddfac[v]*factor;

        for (int a = 1;a < la;a++)
        {
            xtable[0][0][0][a+1][v] =   aafac[v]*xtable[0][0][0][  a][v] +
                                      a*s1fac   *xtable[0][0][0][a-1][v];
        }

        for (int b = 1;b < lb;b++)
        {
            xtable[0][0][b+1][0][v] =   bbfac[v]*xtable[0][0][  b][0][v] +
                                      b*s1fac   *xtable[0][0][b-1][0][v];
        }

        for (int c = 1;c <= lc;c++)
        {
            xtable[0][c+1][0][0][v] =   ccfac[v]*xtable[0][  c][0][0][v] +
                                      c*s2fac   *xtable[0][c-1][0][0][v];
        }

        for (int d = 1;d < ld;d++)
        {
            xtable[d+1][0][0][0][v] =   ddfac[v]*xtable[  d][0][0][0][v] +
                                      d*s2fac   *xtable[d-1][0][0][0][v];
        }

        for (int c = 1;c <= lc;c++)
        {
            xtable[1][c][0][0][v] =   ddfac[v]*xtable[0][  c][0][0][v] +
                                    c*s2fac   *xtable[0][c-1][0][0][v];

            for (int d = 1;d < ld;d++)
            {
                xtable[d+1][c][0][0][v] =   ddfac[v]*xtable[  d][  c][0][0][v] +
                                          c*s2fac   *xtable[  d][c-1][0][0][v] +
                                          d*s2fac   *xtable[d-1][  c][0][0][v];
            }
        }

        for (int b = 1;b <= lb;b++)
        {
            xtable[0][1][b][0][v] =   ccfac[v]*xtable[0][0][  b][0][v] +
                                    b* gfac[v]*xtable[0][0][b-1][0][v];

            for (int c = 1;c < lc;c++)
            {
                xtable[0][c+1][b][0][v] =   ccfac[v]*xtable[0][  c][  b][0][v] +
                                          b* gfac[v]*xtable[0][  c][b-1][0][v] +
                                          c*s2fac   *xtable[0][c-1][  b][0][v];
            }
        }

        for (int b = 1;b <= lb;b++)
        {
            xtable[1][0][b][0][v] =   ddfac[v]*xtable[0][0][  b][0][v] +
                                    b* gfac[v]*xtable[0][0][b-1][0][v];

            for (int d = 1;d < ld;d++)
            {
                xtable[d+1][0][b][0][v] =   ddfac[v]*xtable[  d][0][  b][0][v] +
                                          b* gfac[v]*xtable[  d][0][b-1][0][v] +
                                          d*s2fac   *xtable[d-1][0][  b][0][v];
            }
        }

        for (int a = 1;a <= la;a++)
        {
            xtable[0][0][1][a][v] =   bbfac[v]*xtable[0][0][0][  a][v] +
                                    a*s1fac   *xtable[0][0][0][a-1][v];

            for (int b = 1;b < lb;b++)
            {
                xtable[0][0][b+1][a][v] =   bbfac[v]*xtable[0][0][  b][  a][v] +
                                          a*s1fac   *xtable[0][0][  b][a-1][v] +
                                          b*s1fac   *xtable[0][0][b-1][  a][v];
            }
        }

        for (int a = 1;a <= la;a++)
        {
            xtable[0][1][0][a][v] =   ccfac[v]*xtable[0][0][0][  a][v] +
                                    a* gfac[v]*xtable[0][0][0][a-1][v];

            for (int c = 1;c < lc;c++)
            {
                xtable[0][c+1][0][a][v] =   ccfac[v]*xtable[0][  c][0][  a][v] +
                                          a* gfac[v]*xtable[0][  c][0][a-1][v] +
                                          c*s2fac   *xtable[0][c-1][0][  a][v];
            }
        }

        for (int a = 1;a <= la;a++)
        {
            xtable[1][0][0][a][v] =   ddfac[v]*xtable[0][0][0][  a][v] +
                                    a* gfac[v]*xtable[0][0][0][a-1][v];

            for (int d = 1;d < ld;d++)
            {
                xtable[d+1][0][0][a][v] =   ddfac[v]*xtable[  d][0][0][  a][v] +
                                          a* gfac[v]*xtable[  d][0][0][a-1][v] +
                                          d*s2fac   *xtable[d-1][0][0][  a][v];
            }
        }

        for (int b = 1;b <= lb;b++)
        {
            for (int c = 1;c <= lc;c++)
            {
                xtable[1][c][b][0][v] =   ddfac[v]*xtable[0][  c][  b][0][v] +
                                        b* gfac[v]*xtable[0][  c][b-1][0][v] +
                                        c*s2fac   *xtable[0][c-1][  b][0][v];

                for (int d = 1;d < ld;d++)
                {
                    xtable[d+1][c][b][0][v] =   ddfac[v]*xtable[  d][  c][  b][0][v] +
                                              b* gfac[v]*xtable[  d][  c][b-1][0][v] +
                                              c*s2fac   *xtable[  d][c-1][  b][0][v] +
                                              d*s2fac   *xtable[d-1][  c][  b][0][v];
                }
            }
        }

        for (int a = 1;a <= la;a++)
        {
            for (int c = 1;c <= lc;c++)
            {
                xtable[1][c][0][a][v] =   ddfac[v]*xtable[0][  c][0][  a][v] +
                                        a* gfac[v]*xtable[0][  c][0][a-1][v] +
                                        c*s2fac   *xtable[0][c-1][0][  a][v];

                for (int d = 1;d < ld;d++)
                {
                    xtable[d+1][c][0][a][v] =   ddfac[v]*xtable[  d][  c][0][  a][v] +
                                              a* gfac[v]*xtable[  d][  c][0][a-1][v] +
                                              c*s2fac   *xtable[  d][c-1][0][  a][v] +
                                              d*s2fac   *xtable[d-1][  c][0][  a][v];
                }
            }
        }

        for (int a = 1;a <= la;a++)
        {
            for (int b = 1;b <= lb;b++)
            {
                xtable[0][1][b][a][v] =   ccfac[v]*xtable[0][0][  b][  a][v] +
                                        a* gfac[v]*xtable[0][0][  b][a-1][v] +
                                        b* gfac[v]*xtable[0][0][b-1][  a][v];

                for (int c = 1;c < lc;c++)
                {
                    xtable[0][c+1][b][a][v] =   ccfac[v]*xtable[0][  c][  b][  a][v] +
                                              a* gfac[v]*xtable[0][  c][  b][a-1][v] +
                                              b* gfac[v]*xtable[0][  c][b-1][  a][v] +
                                              c*s2fac   *xtable[0][c-1][  b][  a][v];
                }
            }
        }

        for (int a = 1;a <= la;a++)
        {
            for (int b = 1;b <= lb;b++)
            {
                xtable[1][0][b][a][v] =   ddfac[v]*xtable[0][0][  b][  a][v] +
                                        a* gfac[v]*xtable[0][0][  b][a-1][v] +
                                        b* gfac[v]*xtable[0][0][b-1][  a][v];

                for (int d = 1;d < ld;d++)
                {
                    xtable[d+1][0][b][a][v] =   ddfac[v]*xtable[  d][0][  b][  a][v] +
                                              a* gfac[v]*xtable[  d][0][  b][a-1][v] +
                                              b* gfac[v]*xtable[  d][0][b-1][  a][v] +
                                              d*s2fac   *xtable[d-1][0][  b][  a][v];
                }
            }
        }

        for (int a = 1;a <= la;a++)
        {
            for (int b = 1;b <= lb;b++)
            {
                for (int c = 1;c <= lc;c++)
                {
                    xtable[1][c][b][a][v] =   ddfac[v]*xtable[0][  c][  b][  a][v] +
                                            a* gfac[v]*xtable[0][  c][  b][a-1][v] +
                                            b* gfac[v]*xtable[0][  c][b-1][  a][v] +
                                            c*s2fac   *xtable[0][c-1][  b][  a][v];

                    for (int d = 1;d < ld;d++)
                    {
                        xtable[d+1][c][b][a][v] =   ddfac[v]*xtable[  d][  c][  b][  a][v] +
                                                  a* gfac[v]*xtable[  d][  c][  b][a-1][v] +
                                                  b* gfac[v]*xtable[  d][  c][b-1][  a][v] +
                                                  c*s2fac   *xtable[  d][c-1][  b][  a][v] +
                                                  d*s2fac   *xtable[d-1][  c][  b][  a][v];
                    }
                }
            }
        }
    }
}

}
}
