#ifndef _AQUARIUS_OPERATOR_MOINTS_HPP_
#define _AQUARIUS_OPERATOR_MOINTS_HPP_

#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class MOIntegrals : public task::Task
{
    protected:
        MOIntegrals(const string& name, input::Config& config);
};

}
}

#endif
