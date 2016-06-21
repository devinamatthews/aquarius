#ifndef _AQUARIUS_AGORA_INTEGRALS_INTEGRALS_HPP_
#define _AQUARIUS_AGORA_INTEGRALS_INTEGRALS_HPP_

#include "frameworks/agora.hpp"
#include "frameworks/integrals.hpp"

namespace aquarius
{
namespace integrals
{

class OVI : OneElectronIntegrals { void check() const {} };
class KEI : OneElectronIntegrals { void check() const {} };
class NAI :     NuclearIntegrals { void check() const {} };
class ERI : TwoElectronIntegrals { void check() const {} };
class MOM :      MomentIntegrals { void check() const {} };

class Integrals : public agora::Marketplace<Integrals>
{
    protected:
        OVI& getBestVendor(ptr_vector<OVI>& matches)
        {
            return matches[0];
        }

        KEI& getBestVendor(ptr_vector<KEI>& matches)
        {
            return matches[0];
        }

        NAI& getBestVendor(ptr_vector<NAI>& matches)
        {
            return matches[0];
        }

        ERI& getBestVendor(ptr_vector<ERI>& matches)
        {
            return matches[0];
        }

        MOM& getBestVendor(ptr_vector<MOM>& matches)
        {
            return matches[0];
        }
};

}
}

#endif
