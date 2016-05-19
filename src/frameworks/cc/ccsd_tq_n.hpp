#ifndef _AQUARIUS_CC_CCSD_TQ_N_HPP_
#define _AQUARIUS_CC_CCSD_TQ_N_HPP_

#include "util/global.hpp"

#include "task/task.hpp"
#include "time/time.hpp"
#include "operator/2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/st2eoperator.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSD_TQ_N : public task::Task
{
    public:
        CCSD_TQ_N(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
