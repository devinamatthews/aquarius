#include "kei.hpp"

namespace aquarius
{
namespace integrals
{

void IshidaKEI::prim(const vec3& posa, int e,
                     const vec3& posb, int f, double* restrict integrals)
{
    constexpr double PI_32 = 5.5683279968317078452848179821188; // pi^(3/2)

    int vmax = la+lb;

    marray<double,3> stable(3, lb+2, la+2);
    marray<double,3> ttable(3, lb+1, la+1);

    double zp = za[e] + zb[f];
    double A0 = PI_32*exp(-za[e]*zb[f]*norm2(posa-posb)/zp)/pow(zp,1.5);

    vec3 posp = (posa*za[e]+posb*zb[f])/zp;
    vec3 afac = posp-posa;
    vec3 bfac = posp-posb;
    double gfac = 0.5/zp;

    for (int xyz = 0;xyz < 3;xyz++)
    {
        stable[xyz][0][0] = 1.0;
        stable[xyz][0][1] = afac[xyz];
        stable[xyz][1][0] = bfac[xyz];

        for (int b = 1;b <= lb;b++)
        {
            stable[xyz][b+1][0] =   bfac[xyz]*stable[xyz][  b][0] +
                                  b*gfac     *stable[xyz][b-1][0];
        }

        for (int a = 1;a <= la;a++)
        {
            stable[xyz][0][a+1] =   afac[xyz]*stable[xyz][0][  a] +
                                  a*gfac     *stable[xyz][0][a-1];
        }

        for (int a = 1;a <= la+1;a++)
        {
            stable[xyz][1][a] =   bfac[xyz]*stable[xyz][0][  a] +
                                a*gfac     *stable[xyz][0][a-1];

            for (int b = 1;b <= lb;b++)
            {
                stable[xyz][b+1][a] =   bfac[xyz]*stable[xyz][  b][  a] +
                                      a*gfac     *stable[xyz][  b][a-1] +
                                      b*gfac     *stable[xyz][b-1][  a];
            }
        }

        ttable[xyz][0][0] = 2*za[e]*zb[f]*stable[xyz][1][1];

        for (int b = 1;b <= lb;b++)
        {
            ttable[xyz][b][0] = 2*za[e]*zb[f]*stable[xyz][b+1][1] -
                                  za[e]*    b*stable[xyz][b-1][1];
        }

        for (int a = 1;a <= la;a++)
        {
            ttable[xyz][0][a] = 2*za[e]*zb[f]*stable[xyz][1][a+1] -
                                      a*zb[f]*stable[xyz][1][a-1];
        }

        for (int a = 1;a <= la;a++)
        {
            for (int b = 1;b <= lb;b++)
            {
                ttable[xyz][b][a] = 2*za[e]*zb[f]*stable[xyz][b+1][a+1] -
                                          a*zb[f]*stable[xyz][b+1][a-1] -
                                      za[e]*    b*stable[xyz][b-1][a+1] +
                                          a*    b*stable[xyz][b-1][a-1]/2;
            }
        }
    }

    for (int b = 0;b <= lb;b++)
    {
        for (int a = 0;a <= la;a++)
        {
            stable[0][b][a] *= A0;
            ttable[0][b][a] *= A0;
        }
    }

    matrix<double> integral((lb+1)*(lb+2)/2, (la+1)*(la+2)/2, integrals);

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
                    integral[XYZ(bx,by,bz)][XYZ(ax,ay,az)] =
                        ttable[0][bx][ax]*stable[1][by][ay]*stable[2][bz][az] +
                        stable[0][bx][ax]*ttable[1][by][ay]*stable[2][bz][az] +
                        stable[0][bx][ax]*stable[1][by][ay]*ttable[2][bz][az];
                }
            }
        }
    }
}

}
}
