#ifndef _AQUARIUS_OPERATOR_SPARSERHFAOMOINTS_HPP_
#define _AQUARIUS_OPERATOR_SPARSERHFAOMOINTS_HPP_

#include "../../frameworks/integrals/2eints.hpp"
#include "../../frameworks/operator/aomoints.hpp"
#include "../../frameworks/scf/aouhf.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/util/global.hpp"

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
