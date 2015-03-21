#ifndef _AQUARIUS_INTEGRALS_CFOUR1EINTS_HPP_
#define _AQUARIUS_INTEGRALS_CFOUR1EINTS_HPP_

#include "util/global.hpp"

#include "1eints.hpp"

namespace aquarius
{
namespace integrals
{

class CFOUROneElectronIntegralsTask : public task::Task
{
    public:
        CFOUROneElectronIntegralsTask(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
