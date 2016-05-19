#ifndef _AQUARIUS_SCF_AOUHF_HPP_
#define _AQUARIUS_SCF_AOUHF_HPP_

#include "util/global.hpp"

#include "integrals/2eints.hpp"

#include "uhf_local.hpp"
#include "uhf_elemental.hpp"

namespace aquarius
{
namespace scf
{

template <typename T, template <typename T_> class WhichUHF>
class AOUHF : public WhichUHF<T>
{
    protected:
        void buildFock();

    public:
        AOUHF(const string& name, input::Config& config);
};

}
}

#endif
