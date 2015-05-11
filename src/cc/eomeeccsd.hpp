#ifndef _AQUARIUS_CC_EOMEECCSD_HPP_
#define _AQUARIUS_CC_EOMEECCSD_HPP_

#include "util/global.hpp"

#include "convergence/davidson.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/st2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/denominator.hpp"

#include "ccsd.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class EOMEECCSD : public Iterative<U>
{
    protected:
        input::Config davidson_config;
        int nroot;
        int ntriplet;
        bool multiroot;

    public:
        EOMEECCSD(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
