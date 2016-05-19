#ifndef _AQUARIUS_JELLIUM_JELLIUM_HPP_
#define _AQUARIUS_JELLIUM_JELLIUM_HPP_

#include "util/global.hpp"

#include "task/task.hpp"
#include "operator/2eoperator.hpp"
#include "tensor/symblocked_tensor.hpp"

namespace aquarius
{
namespace jellium
{

template <typename U>
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
                            tensor::SymmetryBlockedTensor<U>& tensor);

    public:
        Jellium(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
