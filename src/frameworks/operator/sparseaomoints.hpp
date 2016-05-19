#ifndef _AQUARIUS_OPERATOR_SPARSEAOMOINTS_HPP_
#define _AQUARIUS_OPERATOR_SPARSEAOMOINTS_HPP_

#include "../../frameworks/integrals/2eints.hpp"
#include "../../frameworks/operator/aomoints.hpp"
#include "../../frameworks/operator/moints.hpp"
#include "../../frameworks/scf/aouhf.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class SparseAOMOIntegrals : public MOIntegrals<T>
{
    public:
        SparseAOMOIntegrals(const string& name, input::Config& config);

    protected:
        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
