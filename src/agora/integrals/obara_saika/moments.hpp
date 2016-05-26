#ifndef _AQUARIUS_AGORA_INTEGRALS_MOMENTS_HPP_
#define _AQUARIUS_AGORA_INTEGRALS_MOMENTS_HPP_

#include "frameworks/util.hpp"
#include "frameworks/integrals.hpp"
#include "frameworks/molecule.hpp"
#include "agora/integrals/integrals.hpp"

namespace aquarius
{
namespace integrals
{

class OSMoments : public MOM
{
    protected:
        void filltable(double afac, double bfac, double cfac, double sfac,
                       marray<double,3>&& table)
        {
            filltable(afac, bfac, cfac, sfac, table);
        }

        void filltable(double afac, double bfac, double cfac, double sfac,
                       marray<double,3>& table);

        /*
         * Calculate moment integrals with the algorithm of Obara and Saika
         *  S. Obara; A. Saika, J. Chem. Phys. 84, 3963 (1986)
         */
        void prim(const vec3& posa, int la, double za,
                  const vec3& posb, int lb, double zb,
                  const vec3& posc, int lc, double* integrals) override;
};

REGISTER_VENDOR(Integrals, MOM, OSMoments);

}
}
