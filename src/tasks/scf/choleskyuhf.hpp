#ifndef _AQUARIUS_SCF_CHOLESKYUHF_HPP_
#define _AQUARIUS_SCF_CHOLESKYUHF_HPP_

#include "../../frameworks/integrals/cholesky.hpp"
#include "../../frameworks/scf/uhf.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{
namespace scf
{

template <typename T, template <typename T_> class WhichUHF>
class CholeskyUHF : public WhichUHF<T>
{
    public:
        CholeskyUHF(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

    protected:
        void buildFock();
};

}
}

#endif
