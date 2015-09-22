#ifndef _AQUARIUS_OPERATOR_RHFAOMOINTS_HPP_
#define _AQUARIUS_OPERATOR_RHFAOMOINTS_HPP_

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
