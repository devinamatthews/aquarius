#ifndef _AQUARIUS_INTEGRALS_OVI_HPP_
#define _AQUARIUS_INTEGRALS_OVI_HPP_

#include "util/global.hpp"

#include "1eints.hpp"

namespace aquarius
{
namespace integrals
{

/*
 * Calculate OVIs with the algorithm of Ishida
 *  K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 */
class IshidaOVI : public OneElectronIntegrals
{
    public:
        IshidaOVI(const Shell& a, const Shell& b)
        : OneElectronIntegrals(a, b) {}

        void prim(const vec3& posa, int e,
                  const vec3& posb, int f, double* integrals);
};

}
}

#endif
