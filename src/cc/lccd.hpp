#ifndef _AQUARIUS_CC_LCCD_HPP_
#define _AQUARIUS_CC_LCCD_HPP_

#include "util/global.hpp"

#include "time/time.hpp"
#include "task/task.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/st2eoperator.hpp"
#include "operator/denominator.hpp"
#include "convergence/diis.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class LCCD : public Iterative<U>
{
    protected:
        convergence::DIIS< op::ExcitationOperator<U,2> > diis;

    public:
        LCCD(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
