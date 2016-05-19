#ifndef _AQUARIUS_OPERATOR_CHOLESKYMOINTS_HPP_
#define _AQUARIUS_OPERATOR_CHOLESKYMOINTS_HPP_

#include "../../frameworks/integrals/cholesky.hpp"
#include "../../frameworks/operator/moints.hpp"
#include "../../frameworks/scf/choleskyuhf.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class CholeskyMOIntegrals : public MOIntegrals<T>
{
    public:
        CholeskyMOIntegrals(const string& name, input::Config& config);

    protected:
        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
