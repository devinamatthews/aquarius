#ifndef _AQUARIUS_AGORA_INTEGRALS_NAI_HPP_
#define _AQUARIUS_AGORA_INTEGRALS_NAI_HPP_

#include "frameworks/util.hpp"
#include "frameworks/integrals.hpp"
#include "frameworks/molecule.hpp"
#include "agora/integrals/integrals.hpp"

namespace aquarius
{
namespace integrals
{

/*
 * Calculate NAIs with the Rys Polynomial algorithm of Ishida
 *  K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 *  Ishida, K., J. Chem. Phys., 98, 2176 (1993)
 */
class IshidaNAI : public NAI
{
    protected:
        void filltable(double afac, double bfac, double cfac, double sfac, marray<double,3>&& gtable)
        {
            filltable(afac, bfac, cfac, sfac, gtable);
        }

        void filltable(double afac, double bfac, double cfac, double sfac, marray<double,3>& gtable);

        void prim(const vec3& posa, int la, double za,
                  const vec3& posb, int lb, double zb,
                  const vec3& posc, double charge, double* integrals) override;
};

REGISTER_VENDOR(Integrals, NAI, IshidaNAI);

}
}

#endif
