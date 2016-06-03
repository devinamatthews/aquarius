#ifndef _AQUARIUS_CC_CCSDTQ_1B_DENSITY_HPP_
#define _AQUARIUS_CC_CCSDTQ_1B_DENSITY_HPP_

#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/operator/deexcitationoperator.hpp"
#include "../../frameworks/operator/excitationoperator.hpp"
#include "../../frameworks/operator/st2eoperator.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/time/time.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSDTQ_1b_Density : public task::Task
{
    public:
        CCSDTQ_1b_Density(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
