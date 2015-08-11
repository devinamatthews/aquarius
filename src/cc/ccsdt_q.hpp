#ifndef _AQUARIUS_CC_CCSDT_Q_HPP_
#define _AQUARIUS_CC_CCSDT_Q_HPP_

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
class CCSDT_Q : public task::Task
{
    public:
        CCSDT_Q(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
