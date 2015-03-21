#ifndef _AQUARIUS_OPERATOR_CHOLESKYMOINTS_HPP_
#define _AQUARIUS_OPERATOR_CHOLESKYMOINTS_HPP_

#include "util/global.hpp"

#include "scf/choleskyuhf.hpp"
#include "integrals/cholesky.hpp"

#include "moints.hpp"

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
