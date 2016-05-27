#ifndef _AQUARIUS_CC_CFOURGRAD_HPP_
#define _AQUARIUS_CC_CFOURGRAD_HPP_

#include "util/global.hpp"

#include "task/task.hpp"
#include "time/time.hpp"
#include "operator/2eoperator.hpp"

namespace aquarius
{
namespace cc
{

class CFOURGradient : public task::Task
{
    public:
        CFOURGradient(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
