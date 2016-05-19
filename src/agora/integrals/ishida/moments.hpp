#ifndef _AQUARIUS_INTEGRALS_MOMENTS_HPP_
#define _AQUARIUS_INTEGRALS_MOMENTS_HPP_

#include "util/global.hpp"

#include "shell.hpp"

namespace aquarius
{
namespace integrals
{

class OSMoments : public OneElectronIntegrals
{
    protected:
        int lc;
        vec3 posc;

        void filltable(double afac, double bfac, double cfac, double sfac,
                       marray<double,3>&& table)
        {
            filltable(afac, bfac, cfac, sfac, table);
        }

        void filltable(double afac, double bfac, double cfac, double sfac,
                       marray<double,3>& table);

    public:
        /*
         * Calculate moment integrals with the algorithm of Obara and Saika
         *  S. Obara; A. Saika, J. Chem. Phys. 84, 3963 (1986)
         */
        void prim(const vec3& posa, int e,
                  const vec3& posb, int f, double* integrals);
};

}
}
