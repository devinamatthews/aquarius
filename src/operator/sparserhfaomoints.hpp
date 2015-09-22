#ifndef _AQUARIUS_OPERATOR_SPARSERHFAOMOINTS_HPP_
#define _AQUARIUS_OPERATOR_SPARSERHFAOMOINTS_HPP_

#include "util/global.hpp"

#include "scf/aouhf.hpp"
#include "integrals/2eints.hpp"
#include "task/task.hpp"

#include "aomoints.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class SparseRHFAOMOIntegrals : public task::Task
{
    public:
        SparseRHFAOMOIntegrals(const string& name, input::Config& config);

    protected:
        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
