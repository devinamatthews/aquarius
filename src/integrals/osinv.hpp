#ifndef _AQUARIUS_INTEGRALS_OSINV_HPP_
#define _AQUARIUS_INTEGRALS_OSINV_HPP_

#include "util/global.hpp"

#include "2eints.hpp"

namespace aquarius
{
namespace integrals
{

class OSInvERI : public TwoElectronIntegrals
{
    protected:
        void filltable(const marray<double,4>&  afac, const marray<double,4>&  bfac,
                       const marray<double,4>&  cfac, const marray<double,4>&  dfac,
                       const marray<double,4>&  pfac, const marray<double,4>&  qfac,
                       const marray<double,4>& s1fac, const marray<double,4>& s2fac,
                       const marray<double,4>& t1fac, const marray<double,4>& t2fac,
                       const marray<double,4>&  gfac, marray<double,9>&& table)
        {
            filltable(afac, bfac, cfac, dfac, pfac, qfac,
                      s1fac, s2fac, t1fac, t2fac, gfac, table);
        }

        void filltable(const marray<double,4>&  afac, const marray<double,4>&  bfac,
                       const marray<double,4>&  cfac, const marray<double,4>&  dfac,
                       const marray<double,4>&  pfac, const marray<double,4>&  qfac,
                       const marray<double,4>& s1fac, const marray<double,4>& s2fac,
                       const marray<double,4>& t1fac, const marray<double,4>& t2fac,
                       const marray<double,4>&  gfac, marray<double,9>& table);

    public:
        OSInvERI(const Shell& a, const Shell& b, const Shell& c, const Shell& d)
        : TwoElectronIntegrals(a, b, c, d) {}

        /**
         * Calculate ERIs with the recursive algorithm of Obara and Saika
         *  S. Obara; A. Saika, J. Chem. Phys. 84, 3963 (1986)
         */
        void prims(const vec3& posa, const vec3& posb,
                   const vec3& posc, const vec3& posdh, double* integrals);
};

}
}

#endif
