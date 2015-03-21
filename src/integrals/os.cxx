#include "os.hpp"

#include "fmgamma.hpp"

namespace aquarius
{
namespace integrals
{

void OSERI::prim(const vec3& posa, int e, const vec3& posb, int f,
                 const vec3& posc, int g, const vec3& posd, int h, double* restrict integrals)
{
    constexpr double TWO_PI_52 = 34.98683665524972497; // 2*pi^(5/2)
    int vmax = la+lb+lc+ld;

    double zp = za[e]+zb[f];
    double zq = zc[g]+zd[h];

    double A0 = TWO_PI_52*exp(-za[e]*zb[f]*norm2(posa-posb)/zp
    		                  -zc[g]*zd[h]*norm2(posc-posd)/zq)/(zp*zq*sqrt(zp+zq));

    /*
    if (fabs(A0) < 1e-14)
    {
    	dzero(((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)*
    		  ((lc+1)*(lc+2)/2)*((ld+1)*(ld+2)/2), integrals, 1);
    	return;
    }
    */

    vec3 posp = (posa*za[e] + posb*zb[f])/zp;
    vec3 posq = (posc*zc[g] + posd*zd[h])/zq;
    vec3 posw = (posp*zp    + posq*zq   )/(zp+zq);

    vec3 afac = posp - posa;
    vec3 bfac = posp - posb;
    vec3 cfac = posq - posc;
    vec3 dfac = posq - posd;
    vec3 pfac = posw - posp;
    vec3 qfac = posw - posq;

    double s1fac = 1.0/(2*zp);
    double s2fac = 1.0/(2*zq);
    double gfac = 1.0/(2*(zp+zq));
    double t1fac = -gfac*zq/zp;
    double t2fac = -gfac*zp/zq;

    double Z = norm2(posp-posq)*zp*zq/(zp+zq);

    marray<double,5> xtable(ld+1, lc+1, lb+1, la+1, vmax+1);

    Fm fm;
    fm(Z, vmax, xtable[0][0][0][0].data());
    for (int v = 0;v <= vmax;v++)
    {
        xtable[0][0][0][0][v] *= A0;
    }

    marray<double,4> integral((ld+1)*(ld+2)/2, (lc+1)*(lc+2)/2,
                              (lb+1)*(lb+2)/2, (la+1)*(la+2)/2, integrals);

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

void OSERI::filltable(double afac, double bfac, double cfac, double dfac, double pfac, double qfac,
                      double s1fac, double t1fac, double s2fac, double t2fac, double gfac,
                      marray<double,5>& table)
{
    int ld = table.length(0)-1;
    int lc = table.length(1)-1;
    int lb = table.length(2)-1;
    int la = table.length(3)-1;
    int vmax = la+lb+lc+ld;

    /*
    if (ld > 0)
    {
        for (int v = 0;v < vmax;v++)
        {
            table[1][0][0][0][v] = dfac*table[0][0][0][0][  v] +
                                   qfac*table[0][0][0][0][v+1];
        }
    }

    for (int d = 1;d < ld;d++)
    {
        for (int v = 0;v < vmax-d;v++)
        {
            table[d+1][0][0][0][v] =    dfac*table[  d][0][0][0][  v] +
                                        qfac*table[  d][0][0][0][v+1] +
                                     d*s2fac*table[d-1][0][0][0][  v] +
                                     d*t2fac*table[d-1][0][0][0][v+1];
        }
    }

    if (lc > 0)
    {
        for (int v = 0;v < vmax;v++)
        {
            table[0][1][0][0][v] = cfac*table[0][0][0][0][  v] +
                                   qfac*table[0][0][0][0][v+1];
        }
    }

    for (int c = 1;c < lc;c++)
    {
        for (int v = 0;v < vmax-c;v++)
        {
            table[0][c+1][0][0][v] =    cfac*table[0][  c][0][0][  v] +
                                        qfac*table[0][  c][0][0][v+1] +
                                     c*s2fac*table[0][c-1][0][0][  v] +
                                     c*t2fac*table[0][c-1][0][0][v+1];
        }
    }

    if (lb > 0)
    {
        for (int v = 0;v < vmax;v++)
        {
            table[0][0][1][0][v] = bfac*table[0][0][0][0][  v] +
                                   pfac*table[0][0][0][0][v+1];
        }
    }

    for (int b = 1;b < lb;b++)
    {
        for (int v = 0;v < vmax-b;v++)
        {
            table[0][0][b+1][0][v] =    bfac*table[0][0][  b][0][  v] +
                                        pfac*table[0][0][  b][0][v+1] +
                                     b*s1fac*table[0][0][b-1][0][  v] +
                                     b*t1fac*table[0][0][b-1][0][v+1];
        }
    }

    if (la > 0)
    {
        for (int v = 0;v < vmax;v++)
        {
            table[0][0][0][1][v] = afac*table[0][0][0][0][  v] +
                                   pfac*table[0][0][0][0][v+1];
        }
    }

    for (int a = 1;a < la;a++)
    {
        for (int v = 0;v < vmax-a;v++)
        {
            table[0][0][0][a+1][v] =    afac*table[0][0][0][  a][  v] +
                                        pfac*table[0][0][0][  a][v+1] +
                                     a*s1fac*table[0][0][0][a-1][  v] +
                                     a*t1fac*table[0][0][0][a-1][v+1];
        }
    }

    for (int b = 1;b <= lb;b++)
    {
        if (la > 0)
        {
            for (int v = 0;v < vmax-b;v++)
            {
                table[0][0][b][1][v] =    afac*table[0][0][  b][0][  v] +
                                          pfac*table[0][0][  b][0][v+1] +
                                       b*s1fac*table[0][0][b-1][0][  v] +
                                       b*t1fac*table[0][0][b-1][0][v+1];
            }
        }

        for (int a = 1;a < la;a++)
        {
            for (int v = 0;v < vmax-a-b;v++)
            {
                table[0][0][b][a+1][v] =    afac*table[0][0][  b][  a][  v] +
                                            pfac*table[0][0][  b][  a][v+1] +
                                         a*s1fac*table[0][0][  b][a-1][  v] +
                                         a*t1fac*table[0][0][  b][a-1][v+1] +
                                         b*s1fac*table[0][0][b-1][  a][  v] +
                                         b*t1fac*table[0][0][b-1][  a][v+1];
            }
        }
    }

    for (int c = 1;c <= lc;c++)
    {
        if (lb > 0)
        {
            for (int v = 0;v < vmax-c;v++)
            {
                table[0][c][1][0][v] =   bfac*table[0][  c][0][0][  v] +
                                         pfac*table[0][  c][0][0][v+1] +
                                       c*gfac*table[0][c-1][0][0][v+1];
            }
        }

        for (int b = 1;b < lb;b++)
        {
            for (int v = 0;v < vmax-b-c;v++)
            {
                table[0][c][b+1][0][v] =    bfac*table[0][  c][  b][0][  v] +
                                            pfac*table[0][  c][  b][0][v+1] +
                                         b*s1fac*table[0][  c][b-1][0][  v] +
                                         b*t1fac*table[0][  c][b-1][0][v+1] +
                                         c* gfac*table[0][c-1][  b][0][v+1];
            }
        }
    }

    for (int c = 1;c <= lc;c++)
    {
        if (la > 0)
        {
            for (int v = 0;v < vmax-c;v++)
            {
                table[0][c][0][1][v] =   afac*table[0][  c][0][0][  v] +
                                         pfac*table[0][  c][0][0][v+1] +
                                       c*gfac*table[0][c-1][0][0][v+1];
            }
        }

        for (int a = 1;a < la;a++)
        {
            for (int v = 0;v < vmax-a-c;v++)
            {
                table[0][c][0][a+1][v] =    afac*table[0][  c][0][  a][  v] +
                                            pfac*table[0][  c][0][  a][v+1] +
                                         a*s1fac*table[0][  c][0][a-1][  v] +
                                         a*t1fac*table[0][  c][0][a-1][v+1] +
                                         c* gfac*table[0][c-1][0][  a][v+1];
            }
        }
    }

    for (int d = 1;d <= ld;d++)
    {
        if (lc > 0)
        {
            for (int v = 0;v < vmax-d;v++)
            {
                table[d][1][0][0][v] =    cfac*table[  d][0][0][0][  v] +
                                          qfac*table[  d][0][0][0][v+1] +
                                       d*s2fac*table[d-1][0][0][0][  v] +
                                       d*t2fac*table[d-1][0][0][0][v+1];
            }
        }

        for (int c = 1;c < lc;c++)
        {
            for (int v = 0;v < vmax-c-d;v++)
            {
                table[d][c+1][0][0][v] =    cfac*table[  d][  c][0][0][  v] +
                                            qfac*table[  d][  c][0][0][v+1] +
                                         c*s2fac*table[  d][c-1][0][0][  v] +
                                         c*t2fac*table[  d][c-1][0][0][v+1] +
                                         d*s2fac*table[d-1][  c][0][0][  v] +
                                         d*t2fac*table[d-1][  c][0][0][v+1];
            }
        }
    }

    for (int d = 1;d <= ld;d++)
    {
        if (lb > 0)
        {
            for (int v = 0;v < vmax-d;v++)
            {
                table[d][0][1][0][v] =   bfac*table[  d][0][0][0][  v] +
                                         pfac*table[  d][0][0][0][v+1] +
                                       d*gfac*table[d-1][0][0][0][v+1];
            }
        }

        for (int b = 1;b < lb;b++)
        {
            for (int v = 0;v < vmax-b-d;v++)
            {
                table[d][0][b+1][0][v] =    bfac*table[  d][0][  b][0][  v] +
                                            pfac*table[  d][0][  b][0][v+1] +
                                         b*s1fac*table[  d][0][b-1][0][  v] +
                                         b*t1fac*table[  d][0][b-1][0][v+1] +
                                         d* gfac*table[d-1][0][  b][0][v+1];
            }
        }
    }

    for (int d = 1;d <= ld;d++)
    {
        if (la > 0)
        {
            for (int v = 0;v < vmax-d;v++)
            {
                table[d][0][0][1][v] =   afac*table[  d][0][0][0][  v] +
                                         pfac*table[  d][0][0][0][v+1] +
                                       d*gfac*table[d-1][0][0][0][v+1];
            }
        }

        for (int a = 1;a < la;a++)
        {
            for (int v = 0;v < vmax-a-d;v++)
            {
                table[d][0][0][a+1][v] =    afac*table[  d][0][0][  a][  v] +
                                            pfac*table[  d][0][0][  a][v+1] +
                                         a*s1fac*table[  d][0][0][a-1][  v] +
                                         a*t1fac*table[  d][0][0][a-1][v+1] +
                                         d* gfac*table[d-1][0][0][  a][v+1];
            }
        }
    }

    for (int c = 1;c <= lc;c++)
    {
        for (int b = 1;b <= lb;b++)
        {
            if (la > 0)
            {
                for (int v = 0;v < vmax-b-c;v++)
                {
                    table[0][c][b][1][v] =    afac*table[0][  c][  b][0][  v] +
                                              pfac*table[0][  c][  b][0][v+1] +
                                           b*s1fac*table[0][  c][b-1][0][  v] +
                                           b*t1fac*table[0][  c][b-1][0][v+1] +
                                           c* gfac*table[0][c-1][  b][0][v+1];
                }
            }

            for (int a = 1;a < la;a++)
            {
                for (int v = 0;v < vmax-a-b-c;v++)
                {
                    table[0][c][b][a+1][v] =    afac*table[0][  c][  b][  a][  v] +
                                                pfac*table[0][  c][  b][  a][v+1] +
                                             a*s1fac*table[0][  c][  b][a-1][  v] +
                                             a*t1fac*table[0][  c][  b][a-1][v+1] +
                                             b*s1fac*table[0][  c][b-1][  a][  v] +
                                             b*t1fac*table[0][  c][b-1][  a][v+1] +
                                             c* gfac*table[0][c-1][  b][  a][v+1];
                }
            }
        }
    }

    for (int d = 1;d <= ld;d++)
    {
        for (int b = 1;b <= lb;b++)
        {
            if (la > 0)
            {
                for (int v = 0;v < vmax-b-d;v++)
                {
                    table[d][0][b][1][v] =    afac*table[  d][0][  b][0][  v] +
                                              pfac*table[  d][0][  b][0][v+1] +
                                           b*s1fac*table[  d][0][b-1][0][  v] +
                                           b*t1fac*table[  d][0][b-1][0][v+1] +
                                           d* gfac*table[d-1][0][  b][0][v+1];
                }
            }

            for (int a = 1;a < la;a++)
            {
                for (int v = 0;v < vmax-a-b-d;v++)
                {
                    table[d][0][b][a+1][v] =    afac*table[  d][0][  b][  a][  v] +
                                                pfac*table[  d][0][  b][  a][v+1] +
                                             a*s1fac*table[  d][0][  b][a-1][  v] +
                                             a*t1fac*table[  d][0][  b][a-1][v+1] +
                                             b*s1fac*table[  d][0][b-1][  a][  v] +
                                             b*t1fac*table[  d][0][b-1][  a][v+1] +
                                             d* gfac*table[d-1][0][  b][  a][v+1];
                }
            }
        }
    }

    for (int d = 1;d <= ld;d++)
    {
        for (int c = 1;c <= lc;c++)
        {
            if (lb > 0)
            {
                for (int v = 0;v < vmax-c-d;v++)
                {
                    table[d][c][1][0][v] =   bfac*table[  d][  c][0][0][  v] +
                                             pfac*table[  d][  c][0][0][v+1] +
                                           c*gfac*table[  d][c-1][0][0][v+1] +
                                           d*gfac*table[d-1][  c][0][0][v+1];
                }
            }

            for (int b = 1;b < lb;b++)
            {
                for (int v = 0;v < vmax-b-c-d;v++)
                {
                    table[d][c][b+1][0][v] =    bfac*table[  d][  c][  b][0][  v] +
                                                pfac*table[  d][  c][  b][0][v+1] +
                                             b*s1fac*table[  d][  c][b-1][0][  v] +
                                             b*t1fac*table[  d][  c][b-1][0][v+1] +
                                             c* gfac*table[  d][c-1][  b][0][v+1] +
                                             d* gfac*table[d-1][  c][  b][0][v+1];
                }
            }
        }
    }

    for (int d = 1;d <= ld;d++)
    {
        for (int c = 1;c <= lc;c++)
        {
            if (la > 0)
            {
                for (int v = 0;v < vmax-c-d;v++)
                {
                    table[d][c][0][1][v] =   afac*table[  d][  c][0][0][  v] +
                                             pfac*table[  d][  c][0][0][v+1] +
                                           c*gfac*table[  d][c-1][0][0][v+1] +
                                           d*gfac*table[d-1][  c][0][0][v+1];
                }
            }

            for (int a = 1;a < la;a++)
            {
                for (int v = 0;v < vmax-a-c-d;v++)
                {
                    table[d][c][0][a+1][v] =    afac*table[  d][  c][0][  a][  v] +
                                                pfac*table[  d][  c][0][  a][v+1] +
                                             a*s1fac*table[  d][  c][0][a-1][  v] +
                                             a*t1fac*table[  d][  c][0][a-1][v+1] +
                                             c* gfac*table[  d][c-1][0][  a][v+1] +
                                             d* gfac*table[d-1][  c][0][  a][v+1];
                }
            }
        }
    }

    for (int d = 1;d <= ld;d++)
    {
        for (int c = 1;c <= lc;c++)
        {
            for (int b = 1;b <= lb;b++)
            {
                if (la > 0)
                {
                    for (int v = 0;v < vmax-b-c-d;v++)
                    {
                        table[d][c][b][1][v] =    afac*table[  d][  c][  b][0][  v] +
                                                  pfac*table[  d][  c][  b][0][v+1] +
                                               b*s1fac*table[  d][  c][b-1][0][  v] +
                                               b*t1fac*table[  d][  c][b-1][0][v+1] +
                                               c* gfac*table[  d][c-1][  b][0][v+1] +
                                               d* gfac*table[d-1][  c][  b][0][v+1];
                    }
                }

                for (int a = 1;a < la;a++)
                {
                    for (int v = 0;v < vmax-a-b-c-d;v++)
                    {
                        table[d][c][b][a+1][v] =    afac*table[  d][  c][  b][  a][  v] +
                                                    pfac*table[  d][  c][  b][  a][v+1] +
                                                 a*s1fac*table[  d][  c][  b][a-1][  v] +
                                                 a*t1fac*table[  d][  c][  b][a-1][v+1] +
                                                 b*s1fac*table[  d][  c][b-1][  a][  v] +
                                                 b*t1fac*table[  d][  c][b-1][  a][v+1] +
                                                 c* gfac*table[  d][c-1][  b][  a][v+1] +
                                                 d* gfac*table[d-1][  c][  b][  a][v+1];
                    }
                }
            }
        }
    }
    */

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
