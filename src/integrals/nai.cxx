#include "nai.hpp"

#include "fmgamma.hpp"

using MArray::slice::all;

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
                     const vec3& posb, int f, double* integrals)
{
    Fm fm;

    int vmax = la+lb;

    double zp = za[e] + zb[f];
    vec3 posp = (posa*za[e]+posb*zb[f])/zp;

    marray<double,3> gtable(vmax+1, la+1, lb+1);
    matrix<double> integral(fca, fcb, integrals);

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

            fm(A0, Z, gtable[all][0][0]);

            // fill table with x
            filltable(afac[0], bfac[0], cfac[0], sfac, gtable);

            // loop over all possible distributions of x momenta
            for (int ax = la;ax >= 0;ax--)
            {
                for (int bx = lb;bx >= 0;bx--)
                {
                    // and fill remainder with y from that point
                    filltable(afac[1], bfac[1], cfac[1], sfac, gtable[all][range(ax,la+1)][range(bx,lb+1)]);

                    // loop over all possible distributions of y momenta given x
                    for (int ay = la-ax;ay >= 0;ay--)
                    {
                        for (int by = lb-bx;by >= 0;by--)
                        {
                            int az = la-ax-ay;
                            int bz = lb-bx-by;

                            // and fill remainder with z from that point
                            filltable(afac[2], bfac[2], cfac[2], sfac, gtable[all][range(ax+ay,la+1)][range(bx+by,lb+1)]);

                            integral[XYZ(ax,ay,az)][XYZ(bx,by,bz)] += gtable[0][la][lb];
                        }
                    }
                }
            }
        }
    }
}

void IshidaNAI::filltable(double afac, double bfac, double cfac, double sfac, marray<double,3>& gtable)
{
    int la = gtable.length(1)-1;
    int lb = gtable.length(2)-1;
    int vmax = la + lb;

    for (int a = 0;a <= la;a++)
    {
        if (a < la)
        {
            for (int v = 0;v < vmax-a;v++)
            {
                gtable[v][a+1][0] = afac*gtable[  v][a][0] -
                                    cfac*gtable[v+1][a][0];

                if (a > 0)
                {
                    gtable[v][a+1][0] += a*sfac*gtable[  v][a-1][0] -
                                         a*sfac*gtable[v+1][a-1][0];
                }
            }
        }

        for (int b = 0;b < lb;b++)
        {
            for (int v = 0;v < vmax-a-b;v++)
            {
                gtable[v][a][b+1] = bfac*gtable[  v][a][b] -
                                    cfac*gtable[v+1][a][b];

                if (a > 0)
                {
                    gtable[v][a][b+1] += a*sfac*gtable[  v][a-1][b] -
                                         a*sfac*gtable[v+1][a-1][b];
                }

                if (b > 0)
                {
                    gtable[v][a][b+1] += b*sfac*gtable[  v][a][b-1] -
                                         b*sfac*gtable[v+1][a][b-1];
                }
            }
        }
    }
}

}
}
