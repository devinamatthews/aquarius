#ifndef _AQUARIUS_AGORA_INTEGRALS_OVI_HPP_
#define _AQUARIUS_AGORA_INTEGRALS_OVI_HPP_

#include "frameworks/util.hpp"
#include "frameworks/integrals.hpp"
#include "frameworks/molecule.hpp"
#include "agora/integrals/integrals.hpp"

namespace aquarius
{
namespace integrals
{

/*
 * Calculate OVIs with the algorithm of Ishida
 *  K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 */
class IshidaOVI : public OVI
{
    protected:
        void prim(const vec3& posa, int la, double za,
                  const vec3& posb, int lb, double zb,
                  double* integrals) override;
};

REGISTER_VENDOR(Integrals, OVI, IshidaOVI);

}
}

#endif
