#ifndef _AQUARIUS_CC_PERTURBEDLAMBDACCSD_HPP_
#define _AQUARIUS_CC_PERTURBEDLAMBDACCSD_HPP_

#include "util/global.hpp"

#include "operator/2eoperator.hpp"
#include "operator/st2eoperator.hpp"
#include "operator/perturbedst2eoperator.hpp"
#include "operator/deexcitationoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "util/iterative.hpp"
#include "convergence/diis.hpp"

namespace aquarius
{
namespace cc
{

/*
 * Solve the frequency-dependent left-hand eigenfunction response equations:
 *
 *     A             _                    _A
 * <0|L (w)|Phi><Phi|H  + w|Phi> = - <0|L H    |Phi> = - <0|N|Phi>
 *                    open                 open
 *
 *       _    -T   T       T       _A   _    _  A
 * where X = e  X e  = (X e )  and H  = A + (H T )
 *                           c                    c
 *
 * As in LambdaCCSD, the full left-hand eigenfunction is used instead of lambda.
 */
template <typename U>
class PerturbedLambdaCCSD : public Iterative<U>
{
    protected:
        convergence::DIIS<op::DeexcitationOperator<U,2>> diis;

    public:
        PerturbedLambdaCCSD(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
