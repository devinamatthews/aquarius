#ifndef _AQUARIUS_CC_CCSDIPGF_HPP_
#define _AQUARIUS_CC_CCSDIPGF_HPP_

#include "../../frameworks/convergence/complex_linear_krylov.hpp"
#include "../../frameworks/operator/2eoperator.hpp"
#include "../../frameworks/operator/denominator.hpp"
#include "../../frameworks/operator/excitationoperator.hpp"
#include "../../frameworks/operator/st2eoperator.hpp"
#include "../../frameworks/util/global.hpp"
#include "../../frameworks/util/iterative.hpp"

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
