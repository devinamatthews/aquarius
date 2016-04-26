#ifndef _AQUARIUS_CC_PICCSD_HPP_
#define _AQUARIUS_CC_PICCSD_HPP_

#include "util/global.hpp"

#include "operator/st2eoperator.hpp"
#include "operator/deexcitationoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "convergence/diis.hpp"
#include "util/iterative.hpp"
#include "task/task.hpp"

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
