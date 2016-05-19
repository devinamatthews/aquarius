#include "2eints.hpp"

using namespace aquarius::integrals;

namespace aquarius
{

using OS2eIntegralsTask = TwoElectronIntegralsTask<OSERI>;

}

static const char* spec = R"(

storage_cutoff?
    double 1e-14,
calc_cutoff?
    double 1e-15

)";

REGISTER_TASK(aquarius::OS2eIntegralsTask,"2eints",spec);
