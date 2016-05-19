#ifndef _AQUARIUS_SCF_HPP_
#define _AQUARIUS_SCF_HPP_

#include "../../frameworks/task/task.hpp"
#include "../../frameworks/util/global.hpp"
#include "../task/config.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class SCF : public task::Task
{
    public:
        SCF(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
