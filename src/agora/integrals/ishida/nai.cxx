#include "nai.hpp"

#include "fmgamma.hpp"

namespace aquarius
{
namespace integrals
{

/*
 * Calculate NAIs with the Rys Polynomial algorithm of Ishida
 *  K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 *  Ishida, K., J. Chem. Phys., 98, 2176 (1993)
 */
void IshidaNAI::prim(const vec3& posa, int e,
                     const vec3& posb, int f, double* restrict integrals)
{
    Fm fm;

    int vmax = la+lb;

    double zp = za[e] + zb[f];
    vec3 posp = (posa*za[e]+posb*zb[f])/zp;

    marray<double,3> gtable(lb+1, la+1, vmax+1);
    matrix<double> integral((lb+1)*(lb+2)/2, (la+1)*(la+2)/2, integrals);

    for (auto& center : centers)
    {
        double charge = center.getElement().getCharge();

        for (auto& posc : center.getCenters())
        {
            vec3 afac = posp-posa;
            vec3 bfac = posp-posb;
            vec3 cfac = posp-posc;
            double sfac = 0.5/zp;

            double A0 = -charge*2*M_PI*exp(-za[e]*zb[f]*norm2(posa-posb)/zp)/zp;
            double Z = norm2(posp-posc)*zp;

            fm(Z, vmax, gtable[0][0].data());
            for (int i = 0;i <= vmax;i++) gtable[0][0][i] *= A0;

            // fill table with x
            filltable(afac[0], bfac[0], cfac[0], sfac, gtable);

            // loop over all possible distributions of x momenta
            for (int bx = lb;bx >= 0;bx--)
            {
                for (int ax = la;ax >= 0;ax--)
                {
                    // and fill remainder with y from that point
                    filltable(afac[1], bfac[1], cfac[1], sfac, gtable[range(bx,lb+1)][range(ax,la+1)]);

                    // loop over all possible distributions of y momenta given x
                    for (int by = lb-bx;by >= 0;by--)
                    {
                        for (int ay = la-ax;ay >= 0;ay--)
                        {
                            int az = la-ax-ay;
                            int bz = lb-bx-by;

                            // and fill remainder with z from that point
                            filltable(afac[2], bfac[2], cfac[2], sfac, gtable[range(bx+by,lb+1)][range(ax+ay,la+1)]);

                            integral[XYZ(bx,by,bz)][XYZ(ax,ay,az)] += gtable[lb][la][0];
                        }
                    }
                }
            }
        }
    }
}

void IshidaNAI::filltable(double afac, double bfac, double cfac, double sfac, marray<double,3>& gtable)
{
    int lb = gtable.length(0)-1;
    int la = gtable.length(1)-1;
    int vmax = la + lb;

    if (lb > 0)
    {
        for (int v = 0;v < vmax;v++)
        {
            gtable[1][0][v] = bfac*gtable[0][0][  v] -
                              cfac*gtable[0][0][v+1];
        }
    }

    if (la > 0)
    {
        for (int v = 0;v < vmax;v++)
        {
            gtable[0][1][v] = afac*gtable[0][0][  v] -
                              cfac*gtable[0][0][v+1];
        }
    }

    for (int b = 1;b < lb;b++)
    {
        for (int v = 0;v < vmax-b;v++)
        {
            gtable[b+1][0][v] =   bfac*gtable[  b][0][  v] -
                                  cfac*gtable[  b][0][v+1] +
                                b*sfac*gtable[b-1][0][  v] -
                                b*sfac*gtable[b-1][0][v+1];
        }
    }

    for (int a = 1;a < la;a++)
    {
        for (int v = 0;v < vmax-a;v++)
        {
            gtable[0][a+1][v] =   afac*gtable[0][  a][  v] -
                                  cfac*gtable[0][  a][v+1] +
                                a*sfac*gtable[0][a-1][  v] -
                                a*sfac*gtable[0][a-1][v+1];
        }
    }

    for (int b = 1;b <= lb;b++)
    {
        if (la > 0)
        {
            for (int v = 0;v < vmax-b;v++)
            {
                gtable[b][1][v] =   afac*gtable[  b][0][  v] -
                                    cfac*gtable[  b][0][v+1] +
                                  b*sfac*gtable[b-1][0][  v] -
                                  b*sfac*gtable[b-1][0][v+1];
            }
        }

        for (int a = 1;a < la;a++)
        {
            for (int v = 0;v < vmax-a-b;v++)
            {
                gtable[b][a+1][v] =   afac*gtable[  b][  a][  v] -
                                      cfac*gtable[  b][  a][v+1] +
                                    a*sfac*gtable[  b][a-1][  v] -
                                    a*sfac*gtable[  b][a-1][v+1] +
                                    b*sfac*gtable[b-1][  a][  v] -
                                    b*sfac*gtable[b-1][  a][v+1];
            }
        }
    }
}

}
}
