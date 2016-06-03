#ifndef _AQUARIUS_OPERATOR_FCIDUMP_HPP_
#define _AQUARIUS_OPERATOR_FCIDUMP_HPP_

#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/symmetry/symmetry.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class FCIDUMP : public task::Task
{
    protected:
        string path;
        bool semi;
        bool full_fock;

    public:
        FCIDUMP(const string& name, input::Config& config);

    protected:
        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
