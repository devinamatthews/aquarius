#ifndef _AQUARIUS_CC_CCSDTQ_HPP_
#define _AQUARIUS_CC_CCSDTQ_HPP_

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
class CCSDTQ : public Iterative<U>
{
    protected:
        convergence::DIIS<op::ExcitationOperator<U,4>> diis;
        string guess;

    public:
        CCSDTQ(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);

        /*
        double getProjectedS2() const;

        double getProjectedMultiplicity() const;
        */
};

}
}

#endif
