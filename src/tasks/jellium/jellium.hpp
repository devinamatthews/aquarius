#ifndef _AQUARIUS_TASKS_JELLIUM_HPP_
#define _AQUARIUS_TASKS_JELLIUM_HPP_

#include "frameworks/util.hpp"
#include "frameworks/task.hpp"
#include "frameworks/tensor.hpp"

namespace aquarius
{
namespace jellium
{

class Jellium : public task::Task
{
    protected:
        int nelec;
        int norb;
        double radius;
        vector<vec3> gvecs;
        int nocc;
        double V;
        double L;
        double PotVm;

        void writeIntegrals(bool pvirt, bool qvirt, bool rvirt, bool svirt,
                            tensor::Tensor<PGSYMMETRIC>& tensor);

    public:
        Jellium(const string& name, task::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
