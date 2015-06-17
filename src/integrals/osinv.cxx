#include "osinv.hpp"

#include "fmgamma.hpp"

/*
 * marray flags (bit field in class w/ accessors) -- alignment, contiguity (single stride), vectorizability (unit stride)
 * flags for expressions (template parameter?)
 * make flags trinary (yes aligned, no aligned, check)
 * functions to force/prohibit alignment etc.
 * add stride to expression
 * different algorithms for assignment (vector, no vector etc.)
 */

using MArray::slice::all;

namespace aquarius
{
namespace integrals
{

void OSInvERI::prims(const vec3& posa, const vec3& posb,
                     const vec3& posc, const vec3& posd, double* restrict integrals)
{
    constexpr double TWO_PI_52 = 34.98683665524972497; // 2*pi^(5/2)
    int vmax = la+lb+lc+ld;

    matrix<double> Kab(na, nb), Kcd(nc, nd);
    matrix<double> zp(na, nb), zq(nc, nd);
    matrix<vec3> posp(na, nb, 3), posq(na, nb, 3);

    for (int e = 0;e < na;e++)
    {
        for (int f = 0;f < nb;f++)
        {
            zp[e][f] = za[e]+zb[f];
            Kab[e][f] = exp(-za[e]*zb[f]*norm2(posa-posb)/zp[e][f])/zp[e][f];
            posp[e][f] = (posa*za[e] + posb*zb[f])/zp[e][f];
        }
    }

    for (int g = 0;g < nc;g++)
    {
        for (int h = 0;h < nd;h++)
        {
            zq[g][h] = zc[g]+zd[h];
            Kcd[g][h] = exp(-zc[g]*zd[h]*norm2(posc-posd)/zq[g][h])/zq[g][h];
            posq[g][h] = (posc*zc[g] + posd*zd[h])/zq[g][h];
        }
    }

    marray<double,5> afac(3, na, nb, nc, nd);
    marray<double,5> bfac(3, na, nb, nc, nd);
    marray<double,5> cfac(3, na, nb, nc, nd);
    marray<double,5> dfac(3, na, nb, nc, nd);
    marray<double,5> pfac(3, na, nb, nc, nd);
    marray<double,5> qfac(3, na, nb, nc, nd);
    marray<double,4> s1fac(na, nb, nc, nd);
    marray<double,4> s2fac(na, nb, nc, nd);
    marray<double,4> t1fac(na, nb, nc, nd);
    marray<double,4> t2fac(na, nb, nc, nd);
    marray<double,4>  gfac(na, nb, nc, nd);

    marray<double,9> xtable(ld+1, lc+1, lb+1, la+1, vmax+1, na, nb, nc, nd);

    Fm fm;
    for (int e = 0;e < na;e++)
    {
        for (int f = 0;f < nb;f++)
        {
            for (int g = 0;g < nc;g++)
            {
                for (int h = 0;h < nd;h++)
                {
                    double A0 = TWO_PI_52*Kab[e][f]*Kcd[g][h]/sqrt(zp[e][f]+zq[g][h]);
                    double Z = norm2(posp[e][f]-posq[g][h])*zp[e][f]*zq[g][h]/(zp[e][f]+zq[g][h]);

                    vec3 posw = (posp[e][f]*zp[e][f] + posq[g][h]*zq[g][h])/(zp[e][f]+zq[g][h]);
                    afac[all][e][f][g][h] = posp[e][f]-posa;
                    bfac[all][e][f][g][h] = posp[e][f]-posb;
                    cfac[all][e][f][g][h] = posq[g][h]-posc;
                    dfac[all][e][f][g][h] = posq[g][h]-posd;
                    pfac[all][e][f][g][h] = posw-posp[e][f];
                    qfac[all][e][f][g][h] = posw-posq[g][h];
                    s1fac[e][f][g][h] = 0.5/zp[e][f];
                    s2fac[e][f][g][h] = 0.5/zq[g][h];
                     gfac[e][f][g][h] = 0.5/(zp[e][f]+zq[g][h]);
                    t1fac[e][f][g][h] = -gfac[e][f]*zq[g][h]/zp[e][f];
                    t2fac[e][f][g][h] = -gfac[e][f]*zp[e][f]/zq[g][h];

                    fm(A0, Z, xtable[0][0][0][0][all][e][f][g][h]);
                }
            }
        }
    }

    marray<double,8> integral(fcd, fcc, fcb, fca, na, nb, nc, nd, integrals);

    // fill table with x
    filltable(afac[0], bfac[0], cfac[0], dfac[0], pfac[0], qfac[0],
              s1fac, t1fac, s2fac, t2fac, gfac,
              xtable);

    // loop over all possible distributions of x momenta
    for (int dx = ld;dx >= 0;dx--)
    {
        for (int cx = lc;cx >= 0;cx--)
        {
            for (int bx = lb;bx >= 0;bx--)
            {
                for (int ax = la;ax >= 0;ax--)
                {
                    // and fill remainder with y from that point
                    filltable(afac[1], bfac[1], cfac[1], dfac[1], pfac[1], qfac[1],
                              s1fac, t1fac, s2fac, t2fac, gfac,
                              xtable[range(dx,ld+1)][range(cx,lc+1)][range(bx,lb+1)][range(ax,la+1)]);

                    // loop over all possible distirubtions of y momenta given x
                    for (int dy = ld-dx;dy >= 0;dy--)
                    {
                        for (int cy = lc-cx;cy >= 0;cy--)
                        {
                            for (int by = lb-bx;by >= 0;by--)
                            {
                                for (int ay = la-ax;ay >= 0;ay--)
                                {
                                    int az = la-ax-ay;
                                    int bz = lb-bx-by;
                                    int cz = lc-cx-cy;
                                    int dz = ld-dx-dy;

                                    // and fill remainder with z from that point
                                    filltable(afac[2], bfac[2], cfac[2], dfac[2], pfac[2], qfac[2],
                                              s1fac, t1fac, s2fac, t2fac, gfac,
                                              xtable[range(dx+dy,ld+1)][range(cx+cy,lc+1)][range(bx+by,lb+1)][range(ax+ay,la+1)]);

                                    integral[XYZ(dx,dy,dz)][XYZ(cx,cy,cz)][XYZ(bx,by,bz)][XYZ(ax,ay,az)] =
                                        xtable[ld][lc][lb][la][0];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void OSInvERI::filltable(const marray<double,4>&  afac, const marray<double,4>&  bfac,
                         const marray<double,4>&  cfac, const marray<double,4>&  dfac,
                         const marray<double,4>&  pfac, const marray<double,4>&  qfac,
                         const marray<double,4>& s1fac, const marray<double,4>& s2fac,
                         const marray<double,4>& t1fac, const marray<double,4>& t2fac,
                         const marray<double,4>&  gfac, marray<double,9>& table)
{
    int ld = table.length(0)-1;
    int lc = table.length(1)-1;
    int lb = table.length(2)-1;
    int la = table.length(3)-1;
    int vmax = la+lb+lc+ld;

    for (int d = 0;d <= ld;d++)
    {
        if (d < ld)
        {
            for (int v = 0;v < vmax-d;v++)
            {
                table[d+1][0][0][0][v] = dfac*table[d][0][0][0][  v] +
                                         qfac*table[d][0][0][0][v+1];

                if (d > 0)
                {
                    table[d+1][0][0][0][v] += d*s2fac*table[d-1][0][0][0][  v] +
                                              d*t2fac*table[d-1][0][0][0][v+1];
                }
            }
        }

        for (int c = 0;c <= lc;c++)
        {
            if (c < lc)
            {
                for (int v = 0;v < vmax-c-d;v++)
                {
                    table[d][c+1][0][0][v] = cfac*table[d][c][0][0][  v] +
                                             qfac*table[d][c][0][0][v+1];

                    if (c > 0)
                    {
                        table[d][c+1][0][0][v] += c*s2fac*table[d][c-1][0][0][  v] +
                                                  c*t2fac*table[d][c-1][0][0][v+1];
                    }

                    if (d > 0)
                    {
                        table[d][c+1][0][0][v] += d*s2fac*table[d-1][c][0][0][  v] +
                                                  d*t2fac*table[d-1][c][0][0][v+1];
                    }
                }
            }

            for (int b = 0;b <= lb;b++)
            {
                if (b < lb)
                {
                    for (int v = 0;v < vmax-b-c-d;v++)
                    {
                        table[d][c][b+1][0][v] = bfac*table[d][c][b][0][  v] +
                                                 pfac*table[d][c][b][0][v+1];

                        if (b > 0)
                        {
                            table[d][c][b+1][0][v] += b*s1fac*table[d][c][b-1][0][  v] +
                                                      b*t1fac*table[d][c][b-1][0][v+1];
                        }

                        if (c > 0)
                        {
                            table[d][c][b+1][0][v] += c*gfac*table[d][c-1][b][0][v+1];
                        }

                        if (d > 0)
                        {
                            table[d][c][b+1][0][v] += d*gfac*table[d-1][c][b][0][v+1];
                        }
                    }
                }

                for (int a = 0;a < la;a++)
                {
                    for (int v = 0;v < vmax-a-b-c-d;v++)
                    {
                        table[d][c][b][a+1][v] = afac*table[d][c][b][a][  v] +
                                                 pfac*table[d][c][b][a][v+1];

                        if (a > 0)
                        {
                            table[d][c][b][a+1][v] += a*s1fac*table[d][c][b][a-1][  v] +
                                                      a*t1fac*table[d][c][b][a-1][v+1];
                        }

                        if (b > 0)
                        {
                            table[d][c][b][a+1][v] += b*s1fac*table[d][c][b-1][a][  v] +
                                                      b*t1fac*table[d][c][b-1][a][v+1];
                        }

                        if (c > 0)
                        {
                            table[d][c][b][a+1][v] += c*gfac*table[d][c-1][b][a][v+1];
                        }

                        if (d > 0)
                        {
                            table[d][c][b][a+1][v] += d*gfac*table[d-1][c][b][a][v+1];
                        }
                    }
                }
            }
        }
    }
}

}
}
