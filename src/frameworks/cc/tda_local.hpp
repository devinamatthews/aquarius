#ifndef _AQUARIUS_CC_TDA_LOCAL_HPP_
#define _AQUARIUS_CC_TDA_LOCAL_HPP_

#include "util/global.hpp"

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
class LocalTDA : public task::Task
{
    protected:

    public:
        LocalTDA(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
