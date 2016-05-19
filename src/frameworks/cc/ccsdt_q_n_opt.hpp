#ifndef _AQUARIUS_CC_CCSDT_Q_N_OPT_HPP_
#define _AQUARIUS_CC_CCSDT_Q_N_OPT_HPP_

#include "../../frameworks/operator/2eoperator.hpp"
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
class CCSDT_Q_N_opt : public task::Task
{
    public:
        CCSDT_Q_N_opt(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
