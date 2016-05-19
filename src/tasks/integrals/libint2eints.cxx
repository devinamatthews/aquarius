#include "../../frameworks/integrals/libint2eints.hpp"

#include "2eints.hpp"


using namespace aquarius::integrals;

namespace aquarius
{

using Libint2eIntegralsTask = TwoElectronIntegralsTask<Libint2eIntegrals>;

}

static const char* spec = R"(

storage_cutoff?
    double 1e-14,
calc_cutoff?
    double 1e-15

)";

REGISTER_TASK(aquarius::Libint2eIntegralsTask,"libint2eints",spec);
