#ifndef _AQUARIUS_SCF_UHF_LOCAL_HPP_
#define _AQUARIUS_SCF_UHF_LOCAL_HPP_

#include "util/global.hpp"

#include "uhf.hpp"

namespace aquarius
{
namespace scf
{

class LocalUHF : public UHF
{
    public:
        LocalUHF(const string& name, input::Config& config);

    protected:
        void calcSMinusHalf();

        void diagonalizeFock();
};

}
}

#endif
