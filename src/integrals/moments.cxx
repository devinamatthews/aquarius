#include "moments.hpp"

namespace aquarius
{
namespace integrals
{

/*
 * Calculate moment integrals with the algorithm of Obara and Saika
 *  S. Obara; A. Saika, J. Chem. Phys. 84, 3963 (1986)
 */
void OSMoments::prim(const vec3& posa, int e,
                     const vec3& posb, int f, double* integrals)
{
    constexpr double PI_32 = 5.5683279968317078452848179821188;

    double zp = za[e] + zb[f];
    double A0 = PI_32*exp(-za[e]*zb[f]*norm2(posa-posb)/zp)/pow(zp, 1.5);

    vec3 posp = (posa*za[e] + posb*zb[f])/zp;

    vec3 afac = posp - posa;
    vec3 bfac = posp - posb;
    vec3 cfac = posp - posc;
    double sfac = 1.0/(2*zp);

    marray<double,3> table(la+1, lb+1, lc+1);
    table[0][0][0] = A0;

    marray<double,3> integral((la+1)*(la+2)/2, (lb+1)*(lb+2)/2, (lc+1)*(lc+2)/2, integrals);

    // fill table with x
    filltable(afac[0], bfac[0], cfac[0], sfac, table);

    // loop over all possible distributions of x momenta
    for (int cx = lc;cx >= 0;cx--)
    {
        for (int bx = lb;bx >= 0;bx--)
        {
            for (int ax = la;ax >= 0;ax--)
            {
                // and fill remainder with y from that point
                filltable(afac[1], bfac[1], cfac[1], sfac,
                          table[range(ax,la+1)][range(bx,lb+1)][range(cx,lc+1)]);

                // loop over all possible distirubtions of y momenta given x
                for (int cy = lc-cx;cy >= 0;cy--)
                {
                    for (int by = lb-bx;by >= 0;by--)
                    {
                        for (int ay = la-ax;ay >= 0;ay--)
                        {
                            int az = la-ax-ay;
                            int bz = lb-bx-by;
                            int cz = lc-cx-cy;

                            // and fill remainder with z from that point
                            filltable(afac[2], bfac[2], cfac[2], sfac,
                                      table[range(ax+ay,la+1)][range(bx+by,lb+1)][range(cx+cy,lc+1)]);

                            integral[XYZ(ax,ay,az)][XYZ(bx,by,bz)][XYZ(cx,cy,cz)] = table[la][lb][lc];
                        }
                    }
                }
            }
        }
    }
}

void OSMoments::filltable(double afac, double bfac, double cfac, double sfac,
                          marray<double,3>& table)
{
    int la = table.length(0)-1;
    int lb = table.length(1)-1;
    int lc = table.length(2)-1;

    table[0][0][1] = cfac*table[0][0][0];
    table[0][1][0] = bfac*table[0][0][0];
    table[1][0][0] = afac*table[0][0][0];

    for (int c = 1;c < lc;c++)
    {
        table[0][0][c+1] =   cfac*table[0][0][  c] +
                           c*sfac*table[0][0][c-1];
    }

    for (int b = 1;b < lb;b++)
    {
        table[0][b+1][0] =   bfac*table[0][  b][0] +
                           b*sfac*table[0][b-1][0];
    }

    for (int a = 1;a < la;a++)
    {
        table[a+1][0][0] =   afac*table[  a][0][0] +
                           a*sfac*table[a-1][0][0];
    }

    for (int b = 1;b <= lb;b++)
    {
        table[1][b][0] =   afac*table[0][  b][0] +
                         b*sfac*table[0][b-1][0];

        for (int a = 1;a < la;a++)
        {
            table[a+1][b][0] =   afac*table[  a][  b][0] +
                               a*sfac*table[a-1][  b][0] +
                               b*sfac*table[  a][b-1][0];
        }
    }

    for (int c = 1;c <= lc;c++)
    {
        table[0][1][c] =   bfac*table[0][0][  c] +
                         c*sfac*table[0][0][c-1];

        for (int b = 1;b < lb;b++)
        {
            table[0][b+1][c] =   bfac*table[0][  b][  c] +
                               b*sfac*table[0][b-1][  c] +
                               c*sfac*table[0][  b][c-1];
        }
    }

    for (int c = 1;c <= lc;c++)
    {
        table[1][0][c] =   afac*table[0][0][  c] +
                         c*sfac*table[0][0][c-1];

        for (int a = 1;a < la;a++)
        {
            table[a+1][0][c] =   afac*table[  a][0][  c] +
                               a*sfac*table[a-1][0][  c] +
                               c*sfac*table[  a][0][c-1];
        }
    }

    for (int c = 1;c <= lc;c++)
    {
        for (int b = 1;b <= lb;b++)
        {
            table[1][b][c] =   afac*table[0][  b][  c] +
                             b*sfac*table[0][b-1][  c] +
                             c*sfac*table[0][  b][c-1];

            for (int a = 1;a < la;a++)
            {
                table[a+1][b][c] =   afac*table[  a][  b][  c] +
                                   a*sfac*table[a-1][  b][  c] +
                                   b*sfac*table[  a][b-1][  c] +
                                   c*sfac*table[  a][  b][c-1];
            }
        }
    }
}

}
}
