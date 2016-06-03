#ifndef _AQUARIUS_CC_LAMBDACCSDT_Q_HPP_
#define _AQUARIUS_CC_LAMBDACCSDT_Q_HPP_

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

/*
 * Solve the left-hand coupled cluster eigenvalue equation:
 *
 *               _
 * <0|L|Phi><Phi|H    |Phi> = 0
 *                open
 *
 *       _    -T   T       T
 * where X = e  X e  = (X e )
 *                           c
 */
template <typename U>
class LambdaCCSDT_Q : public Iterative<U>
{
    protected:
        convergence::DIIS<op::DeexcitationOperator<U,3>> diis;

    public:
        LambdaCCSDT_Q(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
