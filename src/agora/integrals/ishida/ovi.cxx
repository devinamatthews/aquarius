#include "ovi.hpp"

namespace aquarius
{
namespace integrals
{

void IshidaOVI::prim(const vec3& posa, int e,
                     const vec3& posb, int f, double* restrict integrals)
{
    constexpr double PI_32 = 5.5683279968317078452848179821188; // pi^(3/2)

    int vmax = la+lb;
    int nab = na*nb;

    marray<double,3> stable(3, lb+1, la+1);

    double zp = za[e] + zb[f];
    double A0 = PI_32*exp(-za[e]*zb[f]*norm2(posa-posb)/zp)/pow(zp,1.5);

    vec3 posp = (posa*za[e]+posb*zb[f])/zp;
    vec3 afac = posp-posa;
    vec3 bfac = posp-posb;
    double gfac = 0.5/zp;

    for (int xyz = 0;xyz < 3;xyz++)
    {
        stable[xyz][0][0] = 1.0;
        if (la > 0) stable[xyz][0][1] = afac[xyz];
        if (lb > 0) stable[xyz][1][0] = bfac[xyz];

        for (int b = 1;b < lb;b++)
        {
            stable[xyz][b+1][0] =   bfac[xyz]*stable[xyz][  b][0] +
                                  b*gfac     *stable[xyz][b-1][0];
        }

        for (int a = 1;a < la;a++)
        {
            stable[xyz][0][a+1] =   afac[xyz]*stable[xyz][0][  a] +
                                  a*gfac     *stable[xyz][0][a-1];
        }

        for (int a = 1;a <= la;a++)
        {
            if (lb > 0)
            {
                stable[xyz][1][a] =   bfac[xyz]*stable[xyz][0][  a] +
                                    a*gfac     *stable[xyz][0][a-1];
            }

            for (int b = 1;b < lb;b++)
            {
                stable[xyz][b+1][a] =   bfac[xyz]*stable[xyz][  b][  a] +
                                      a*gfac     *stable[xyz][  b][a-1] +
                                      b*gfac     *stable[xyz][b-1][  a];
            }
        }
    }

    for (int b = 0;b <= lb;b++)
    {
        for (int a = 0;a <= la;a++)
        {
            stable[0][b][a] *= A0;
        }
    }

    matrix<double> integral((lb+1)*(lb+2)/2, (la+1)*(la+2)/2, integrals);

    for (int bx = lb;bx >= 0;bx--)
    {
        for (int by = lb-bx;by >= 0;by--)
        {
            int bz = lb-bx-by;
            for (int ax = la;ax >= 0;ax--)
            {
                for (int ay = la-ax;ay >= 0;ay--)
                {
                    int az = la-ax-ay;
                    integral[XYZ(bx,by,bz)][XYZ(ax,ay,az)] = stable[0][bx][ax]*stable[1][by][ay]*stable[2][bz][az];
                }
            }
        }
    }
}

}
}
