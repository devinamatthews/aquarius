#ifndef _AQUARIUS_CC_PERTURBEDCCSD_HPP_
#define _AQUARIUS_CC_PERTURBEDCCSD_HPP_

#include "../../frameworks/cc/ccsd.hpp"
#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/operator/excitationoperator.hpp"
#include "../../frameworks/operator/st2eoperator.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{
namespace cc
{

/*
 * Solve the frequency-dependent amplitude response equations:
 *
 *      _                A                _
 * <Phi|H  - w|Phi><Phi|T (w)|0> = - <Phi|A|0> = - <Phi|X|0>
 *       open
 *
 *       _    -T   T       T
 * where X = e  X e  = (X e )
 *                           c
 */
template <typename U>
class PerturbedCCSD : public Iterative<U>
{
    protected:
        U omega;
        convergence::DIIS<op::ExcitationOperator<U,2>> diis;

    public:
        PerturbedCCSD(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
