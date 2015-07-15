#ifndef _AQUARIUS_CC_CCSDIPGF_HPP_
#define _AQUARIUS_CC_CCSDIPGF_HPP_

#include "util/global.hpp"

#include "convergence/complex_linear_krylov.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/st2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/denominator.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSDIPGF : public Iterative<complex_type_t<U>>
{
    protected:
        typedef complex_type_t<U> CU;

        input::Config krylov_config;
        int orbital;
        vector<CU> omegas;
        CU omega;

    public:
        CCSDIPGF(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
