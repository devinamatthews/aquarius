#ifndef _AQUARIUS_INTEGRALS_OS_HPP_
#define _AQUARIUS_INTEGRALS_OS_HPP_

#include "util/global.hpp"

#include "2eints.hpp"

namespace aquarius
{
namespace integrals
{

class OSERI : public TwoElectronIntegrals
{
    protected:
        void filltable(double afac, double bfac, double cfac, double dfac, double pfac, double qfac,
                       double s1fac, double s2fac, double t1fac, double t2fac, double gfac,
                       marray<double,5>&& table)
        {
            filltable(afac, bfac, cfac, dfac, pfac, qfac,
                      s1fac, s2fac, t1fac, t2fac, gfac, table);
        }

        void filltable(double afac, double bfac, double cfac, double dfac, double pfac, double qfac,
                       double s1fac, double s2fac, double t1fac, double t2fac, double gfac,
                       marray<double,5>& table);

    public:
        OSERI(const Shell& a, const Shell& b, const Shell& c, const Shell& d)
        : TwoElectronIntegrals(a, b, c, d) {}

        /**
         * Calculate ERIs with the recursive algorithm of Obara and Saika
         *  S. Obara; A. Saika, J. Chem. Phys. 84, 3963 (1986)
         */
        void prim(const vec3& posa, int e, const vec3& posb, int f,
                  const vec3& posc, int g, const vec3& posd, int h, double* restrict integrals);
};

}
}

#endif
