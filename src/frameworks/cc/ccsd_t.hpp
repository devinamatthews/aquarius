#ifndef _AQUARIUS_CC_CCSD_T_HPP_
#define _AQUARIUS_CC_CCSD_T_HPP_

#include "util/global.hpp"

#include "task/task.hpp"
#include "time/time.hpp"
#include "operator/2eoperator.hpp"
#include "operator/excitationoperator.hpp"

#include "ccsd.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSD_T : public task::Task
{
    public:
        CCSD_T(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
