#ifndef _AQUARIUS_INTEGRALS_NAI_HPP_
#define _AQUARIUS_INTEGRALS_NAI_HPP_

#include "util/global.hpp"

#include "1eints.hpp"

namespace aquarius
{
namespace integrals
{

/*
 * Calculate NAIs with the Rys Polynomial algorithm of Ishida
 *  K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 *  Ishida, K., J. Chem. Phys., 98, 2176 (1993)
 */
class IshidaNAI : public OneElectronIntegrals
{
    protected:
        vector<Center> centers;

        void filltable(double afac, double bfac, double cfac, double sfac, marray<double,3>&& gtable)
        {
            filltable(afac, bfac, cfac, sfac, gtable);
        }

        void filltable(double afac, double bfac, double cfac, double sfac, marray<double,3>& gtable);

    public:
        IshidaNAI(const Shell& a, const Shell& b, const vector<Center>& centers)
        : OneElectronIntegrals(a, b), centers(centers) {}

        void prim(const vec3& posa, int e,
                  const vec3& posb, int f, double* integrals);
};

}
}

#endif
