#ifndef _AQUARIUS_AGORA_INTEGRALS_OS_HPP_
#define _AQUARIUS_AGORA_INTEGRALS_OS_HPP_

#include "frameworks/util.hpp"
#include "frameworks/molecule.hpp"
#include "frameworks/integrals.hpp"
#include "agora/integrals/integrals.hpp"

namespace aquarius
{
namespace integrals
{

class OSERI : public ERI
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

        /**
         * Calculate ERIs with the recursive algorithm of Obara and Saika
         *  S. Obara; A. Saika, J. Chem. Phys. 84, 3963 (1986)
         */
        void prim(const vec3& posa, int la, double za,
                  const vec3& posb, int lb, double zb,
                  const vec3& posc, int lc, double zc,
                  const vec3& posd, int ld, double zd,
                  double* integrals) override;
};

REGISTER_VENDOR(Integrals, OSERI);

}
}

#endif
