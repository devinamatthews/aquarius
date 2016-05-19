#ifndef _AQUARIUS_CFOURSCF_HPP_
#define _AQUARIUS_CFOURSCF_HPP_

#include "../../frameworks/convergence/diis.hpp"
#include "../../frameworks/input/config.hpp"
#include "../../frameworks/molecule/molecule.hpp"
#include "../../frameworks/operator/space.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/util/global.hpp"
#include "../../frameworks/util/iterative.hpp"
#include "tensor/symblocked_tensor.hpp"
#include "integrals/1eints.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class CFOURSCF : public task::Task
{
    protected:
        bool frozen_core;
        bool semicanonical;

    public:
        CFOURSCF(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
