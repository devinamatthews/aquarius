#ifndef _AQUARIUS_CC_CCSD_HPP_
#define _AQUARIUS_CC_CCSD_HPP_

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
class CCSD : public Iterative<U>
{
    protected:
        convergence::DIIS< op::ExcitationOperator<U,2> > diis;

    public:
        CCSD(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);

        /*
        static double getProjectedS2(const op::MOSpace<U>& occ, const op::MOSpace<U>& vrt,
                                     const tensor::SpinorbitalTensor<U>& T1,
                                     const tensor::SpinorbitalTensor<U>& T2);
         */
};

}
}

#endif
