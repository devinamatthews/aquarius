#ifndef _AQUARIUS_CC_TDA_ELEMENTAL_HPP_
#define _AQUARIUS_CC_TDA_ELEMENTAL_HPP_

#include "../../frameworks/cc/ccsd.hpp"
#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/operator/denominator.hpp"
#include "../../frameworks/operator/excitationoperator.hpp"
#include "../../frameworks/operator/st2eoperator.hpp"
#include "../../frameworks/util/global.hpp"
#include "../molecule/molecule.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class ElementalTDA : public task::Task
{
    protected:

    public:
        ElementalTDA(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
