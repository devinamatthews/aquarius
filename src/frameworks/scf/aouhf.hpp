#ifndef _AQUARIUS_SCF_AOUHF_HPP_
#define _AQUARIUS_SCF_AOUHF_HPP_

#include "../../frameworks/integrals/2eints.hpp"
#include "../../frameworks/scf/uhf_elemental.hpp"
#include "../../frameworks/scf/uhf_local.hpp"
#include "../../frameworks/util/global.hpp"

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
