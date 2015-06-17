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

    marray<double,3> stable(3, la+1, lb+1);

    double zp = za[e] + zb[f];
    double A0 = PI_32*exp(-za[e]*zb[f]*norm2(posa-posb)/zp)/pow(zp,1.5);

    vec3 posp = (posa*za[e]+posb*zb[f])/zp;
    vec3 afac = posp-posa;
    vec3 bfac = posp-posb;
    double gfac = 0.5/zp;

    for (int xyz = 0;xyz < 3;xyz++)
    {
        stable[xyz][0][0] = (xyz == 0 ? A0 : 1);

        for (int a = 0;a <= la;a++)
        {
            if (a < la)
            {
                stable[xyz][a+1][0] = afac[xyz]*stable[xyz][a][0];

                if (a > 0)
                {
                    stable[xyz][a+1][0] += a*gfac*stable[xyz][a-1][0];
                }
            }

            for (int b = 0;b < lb;b++)
            {
                stable[xyz][a][b+1] = bfac[xyz]*stable[xyz][a][b];

                if (a > 0)
                {
                    stable[xyz][a][b+1] += a*gfac*stable[xyz][a-1][b];
                }

                if (b > 0)
                {
                    stable[xyz][a][b+1] += b*gfac*stable[xyz][a][b-1];
                }
            }
        }
    }

    matrix<double> integral(fca, fcb, integrals);

    for (int ax = la;ax >= 0;ax--)
    {
        for (int ay = la-ax;ay >= 0;ay--)
        {
            int az = la-ax-ay;
            for (int bx = lb;bx >= 0;bx--)
            {
                for (int by = lb-bx;by >= 0;by--)
                {
                    int bz = lb-bx-by;
                    integral[XYZ(ax,ay,az)][XYZ(bx,by,bz)] = stable[0][ax][bx]*stable[1][ay][by]*stable[2][az][bz];
                }
            }
        }
    }
}

}
}
