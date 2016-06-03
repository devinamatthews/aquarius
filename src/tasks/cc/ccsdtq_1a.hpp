#ifndef _AQUARIUS_CC_CCSDTQ_1A_HPP_
#define _AQUARIUS_CC_CCSDTQ_1A_HPP_

#include "../../frameworks/cc/ccsd.hpp"
#include "../../frameworks/convergence/diis.hpp"
#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/operator/excitationoperator.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/time/time.hpp"
#include "../../frameworks/util/global.hpp"
#include "../../frameworks/util/iterative.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSDTQ_1a : public Iterative<U>
{
    protected:
        convergence::DIIS<op::ExcitationOperator<U,3>> diis;

    public:
        CCSDTQ_1a(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
