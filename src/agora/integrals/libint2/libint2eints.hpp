#ifndef _AQUARIUS_INTEGRALS_LIBINT2EINTS_HPP_
#define _AQUARIUS_INTEGRALS_LIBINT2EINTS_HPP_

#include "util/global.hpp"

#include "2eints.hpp"

#include "libint2.h"

namespace aquarius
{
namespace integrals
{

class Libint2eIntegrals : public TwoElectronIntegrals
{
    protected:
        vector<Libint_t> inteval;

    public:
        Libint2eIntegrals(const Shell& a, const Shell& b, const Shell& c, const Shell& d);

    protected:
        void prims(const vec3& posa, const vec3& posb, const vec3& posc, const vec3& posd,
                   double* integrals);
};

using Libint2eIntegralsTask = TwoElectronIntegralsTask<Libint2eIntegrals>;

}
}

#endif
