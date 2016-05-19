#ifndef _AQUARIUS_OPERATOR_MOINTS_HPP_
#define _AQUARIUS_OPERATOR_MOINTS_HPP_

#include "util/global.hpp"

#include "task/task.hpp"

#include "2eoperator.hpp"

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
