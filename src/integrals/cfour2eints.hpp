#ifndef _AQUARIUS_INTEGRALS_CFOUR2EINTS_HPP_
#define _AQUARIUS_INTEGRALS_CFOUR2EINTS_HPP_

#include "util/global.hpp"

#include "2eints.hpp"

namespace aquarius
{
namespace integrals
{

class CFOURTwoElectronIntegralsTask : public task::Task
{
    public:
        CFOURTwoElectronIntegralsTask(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
