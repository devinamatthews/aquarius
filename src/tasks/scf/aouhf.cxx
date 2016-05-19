#include "aouhf.hpp"

#include "scf/local_uhf.hpp"

using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;
using namespace aquarius::scf;

namespace aquarius
{

using LocalAOUHF = AOUHF<LocalUHF>;

}

static const char* spec = R"(

    frozen_core?
        bool false,
    convergence?
        double 1e-12,
    max_iterations?
        int 150,
    conv_type?
        enum { MAXE, RMSE, MAE },
    diis?
    {
        damping?
            double 0.0,
        start?
            int 8,
        order?
            int 6,
        jacobi?
            bool false
    }

)";

INSTANTIATE_SPECIALIZATIONS_2(aquarius::scf::AOUHF, aquarius::scf::LocalUHF);
REGISTER_TASK(CONCAT(aquarius::scf::AOUHF<double,aquarius::scf::LocalUHF>), "localaoscf",spec);

#if HAVE_ELEMENTAL
INSTANTIATE_SPECIALIZATIONS_2(aquarius::scf::AOUHF, aquarius::scf::ElementalUHF);
REGISTER_TASK(CONCAT(aquarius::scf::AOUHF<double,aquarius::scf::ElementalUHF>), "elementalaoscf",spec);
#endif
