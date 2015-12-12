#ifndef _AQUARIUS_OPERATOR_FCIDUMP_HPP_
#define _AQUARIUS_OPERATOR_FCIDUMP_HPP_

#include "util/global.hpp"

#include "symmetry/symmetry.hpp"
#include "task/task.hpp"

#include "2eoperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class FCIDUMP : public task::Task
{
    protected:
        string path;

    public:
        FCIDUMP(const string& name, input::Config& config);

    protected:
        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
