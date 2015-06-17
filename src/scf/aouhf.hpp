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

template <class WhichUHF>
class AOUHF : public WhichUHF
{
    protected:
        void buildFock();

    public:
        AOUHF(const string& name, input::Config& config);
};

}
}

#endif
