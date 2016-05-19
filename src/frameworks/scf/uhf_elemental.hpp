#ifndef _AQUARIUS_SCF_UHF_ELEMENTAL_HPP_
#define _AQUARIUS_SCF_UHF_ELEMENTAL_HPP_

#include "uhf.hpp"
#include "util/global.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class ElementalUHF : public UHF<T>
{
    public:
        ElementalUHF(const string& name, input::Config& config);

    protected:
        using UHF<T>::E_alpha;
        using UHF<T>::E_beta;

        void calcSMinusHalf();

        void diagonalizeFock();
};

}
}

#endif
