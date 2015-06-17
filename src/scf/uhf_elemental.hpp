#ifndef _AQUARIUS_SCF_UHF_ELEMENTAL_HPP_
#define _AQUARIUS_SCF_UHF_ELEMENTAL_HPP_

#include "uhf.hpp"
#include "util/global.hpp"

namespace aquarius
{
namespace scf
{

class ElementalUHF : public UHF
{
    public:
        ElementalUHF(const string& name, input::Config& config);

    protected:
        void calcSMinusHalf();

        void diagonalizeFock();
};

}
}

#endif
