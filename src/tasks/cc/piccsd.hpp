#ifndef _AQUARIUS_CC_PICCSD_HPP_
#define _AQUARIUS_CC_PICCSD_HPP_

#include "../../frameworks/convergence/diis.hpp"
#include "../../frameworks/operator/deexcitationoperator.hpp"
#include "../../frameworks/operator/excitationoperator.hpp"
#include "../../frameworks/operator/st2eoperator.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/util/global.hpp"
#include "../../frameworks/util/iterative.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class PiCCSD : public Iterative<U>
{
    protected:
        convergence::DIIS<op::DeexcitationOperator<U,2>> diis;

    public:
        PiCCSD(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
