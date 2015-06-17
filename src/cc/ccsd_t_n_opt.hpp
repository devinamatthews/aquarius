#ifndef _AQUARIUS_CC_CCSD_T_N_OPT_HPP_
#define _AQUARIUS_CC_CCSD_T_N_OPT_HPP_

#include "util/global.hpp"

#include "task/task.hpp"
#include "time/time.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "convergence/diis.hpp"

#include "ccsd.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSD_T_N_opt : public task::Task
{
    public:
        CCSD_T_N_opt(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
