#ifndef _AQUARIUS_CC_EOMEECCSD_HPP_
#define _AQUARIUS_CC_EOMEECCSD_HPP_

#include "../../frameworks/convergence/davidson.hpp"
#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/operator/denominator.hpp"
#include "../../frameworks/operator/excitationoperator.hpp"
#include "../../frameworks/operator/st2eoperator.hpp"
#include "../../frameworks/util/global.hpp"
#include "../../frameworks/util/iterative.hpp"
#include "../molecule/molecule.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class RHFEOMEECCSD : public Iterative<U>
{
    protected:
        input::Config davidson_config;
        int nroot;
        int nsinglet;
        int ntriplet;
        bool triplet;
        bool quintet;

    public:
        RHFEOMEECCSD(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
