#ifndef _AQUARIUS_AGORA_INTEGRALS_ISHIDA_HPP_
#define _AQUARIUS_AGORA_INTEGRALS_ISHIDA_HPP_

#include "frameworks/util.hpp"
#include "frameworks/integrals.hpp"
#include "frameworks/molecule.hpp"
#include "agora/integrals/integrals.hpp"

namespace aquarius
{
namespace integrals
{

class IshidaERI : public ERI
{
    protected:
        void filltable(double factor,
                       row<double>& aafac, row<double>& bbfac, row<double>& ccfac, row<double>& ddfac,
                       double s1fac, double s2fac, row<double>& gfac, marray<double,5>&& xtable);

        /*
         * Calculate ERIs with the Rys Polynomial algorithm of Ishida
         *  Ishida, K. J. Chem. Phys. 95, 5198-205 (1991)
         *  Ishida, K. J. Chem. Phys. 98, 2176 (1993)
         */
        void prim(const vec3& posa, int la, double za,
                  const vec3& posb, int lb, double zb,
                  const vec3& posc, int lc, double zc,
                  const vec3& posd, int ld, double zd,
                  double* integrals) override;
};

REGISTER_VENDOR(Integrals, ERI, IshidaERI);

}
}
