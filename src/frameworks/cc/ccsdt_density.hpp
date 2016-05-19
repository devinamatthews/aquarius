#ifndef _AQUARIUS_CC_CCSDT_DENSITY_HPP_
#define _AQUARIUS_CC_CCSDT_DENSITY_HPP_

#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/operator/deexcitationoperator.hpp"
#include "../../frameworks/operator/excitationoperator.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/time/time.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSDTDensity : public task::Task
{
    public:
        CCSDTDensity(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
