#ifndef _AQUARIUS_CFOURSCF_HPP_
#define _AQUARIUS_CFOURSCF_HPP_

#include "util/global.hpp"

#include "input/config.hpp"
#include "task/task.hpp"
#include "tensor/symblocked_tensor.hpp"
#include "integrals/1eints.hpp"
#include "input/molecule.hpp"
#include "util/iterative.hpp"
#include "convergence/diis.hpp"
#include "operator/space.hpp"

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
