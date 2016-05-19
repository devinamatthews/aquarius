#include "1eints.hpp"

using namespace aquarius::integrals;

namespace aquarius
{

using Ishida1eIntegralsTask = OneElectronIntegralsTask<IshidaOVI, IshidaKEI, IshidaNAI>;

}

REGISTER_TASK(aquarius::Ishida1eIntegralsTask,"1eints");
