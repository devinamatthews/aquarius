#ifndef _AQUARIUS_CC_E_CCSD_T_N_HPP_
#define _AQUARIUS_CC_E_CCSD_T_N_HPP_

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
class E_CCSD_T_N : public task::Task
{
    public:
        E_CCSD_T_N(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
