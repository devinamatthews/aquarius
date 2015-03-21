#ifndef _AQUARIUS_INTEGRALS_ISHIDA_HPP_
#define _AQUARIUS_INTEGRALS_ISHIDA_HPP_

#include "util/global.hpp"

namespace aquarius
{
namespace integrals
{

class IshidaERI : public TwoElectronIntegrals
{
    protected:
        void filltable(double factor,
                       row<double>& aafac, row<double>& bbfac, row<double>& ccfac, row<double>& ddfac,
                       double s1fac, double s2fac, row<double>& gfac, marray<double,5>&& xtable);

    public:
        /*
         * Calculate ERIs with the Rys Polynomial algorithm of Ishida
         *  Ishida, K. J. Chem. Phys. 95, 5198-205 (1991)
         *  Ishida, K. J. Chem. Phys. 98, 2176 (1993)
         */
        void prim(const vec3& posa, int e, const vec3& posb, int f,
                  const vec3& posc, int g, const vec3& posd, int h, double* restrict integrals);
};

}
}
