#ifndef _AQUARIUS_AGORA_INTEGRALS_INTEGRALS_HPP_
#define _AQUARIUS_AGORA_INTEGRALS_INTEGRALS_HPP_

#include "frameworks/agora.hpp"
#include "frameworks/integrals.hpp"

namespace aquarius
{
namespace integrals
{

class Integrals
{
    public:
        enum IntegralType {OVI, KEI, NAI, ERI};

        template <IntegralType type> class Vendor;
};

template <>
class Integrals::Vendor<Integrals::OVI> : OneElectronIntegrals {};

}
}

REGISTER_MARKETPLACE(aquarius::integrals::Integrals, "integrals");

#endif
