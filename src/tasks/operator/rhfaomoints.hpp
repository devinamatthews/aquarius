#ifndef _AQUARIUS_OPERATOR_RHFAOMOINTS_HPP_
#define _AQUARIUS_OPERATOR_RHFAOMOINTS_HPP_

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
class RHFAOMOIntegrals : public task::Task
{
    public:
        RHFAOMOIntegrals(const string& name, input::Config& config);

    protected:
        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
