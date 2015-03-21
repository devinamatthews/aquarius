#ifndef _AQUARIUS_CC_TDA_HPP_
#define _AQUARIUS_CC_TDA_HPP_

#include "util/global.hpp"

#include "convergence/davidson.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/st2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/denominator.hpp"
#include "input/molecule.hpp"

#include "ccsd.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class TDA : public task::Task
{
    protected:

    public:
        TDA(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
