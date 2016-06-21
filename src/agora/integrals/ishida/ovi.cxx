#include "frameworks/util.hpp"
#include "frameworks/integrals.hpp"
#include "frameworks/molecule.hpp"

#include "agora/integrals/integrals.hpp"

namespace aquarius
{
namespace integrals
{

/*
 * Calculate OVIs with the algorithm of Ishida
 *  K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 */
class IshidaOVI : public OVI
{
    protected:
        void prim(const vec3& posa, int la, double za,
                  const vec3& posb, int lb, double zb,
                  double* integrals) override
        {
            constexpr double PI_32 = 5.5683279968317078452848179821188; // pi^(3/2)

            int vmax = la+lb;

            marray<double,3> stable(3, lb+1, la+1);

            double zp = za + zb;
            double A0 = PI_32*exp(-za*zb*norm2(posa-posb)/zp)/pow(zp,1.5);

            vec3 posp = (posa*za+posb*zb)/zp;
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

            matrix_view<double> integral((lb+1)*(lb+2)/2, (la+1)*(la+2)/2, integrals);

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
};

}
}

REGISTER_VENDOR(aquarius::integrals::Integrals,
                aquarius::integrals::OVI,
                aquarius::integrals::IshidaOVI);
