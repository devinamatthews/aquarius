#ifndef _AQUARIUS_CC_MP4DQ_HPP_
#define _AQUARIUS_CC_MP4DQ_HPP_

#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/operator/denominator.hpp"
#include "../../frameworks/operator/excitationoperator.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/time/time.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class MP4DQ : public task::Task
{
    protected:

    public:
        MP4DQ(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
