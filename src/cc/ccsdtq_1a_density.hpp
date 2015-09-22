#ifndef _AQUARIUS_CC_CCSDTQ_1A_DENSITY_HPP_
#define _AQUARIUS_CC_CCSDTQ_1A_DENSITY_HPP_

#include "util/global.hpp"

#include "time/time.hpp"
#include "task/task.hpp"
#include "operator/2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/deexcitationoperator.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSDTQ_1a_Density : public task::Task
{
    public:
        CCSDTQ_1a_Density(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
