#ifndef _AQUARIUS_INTEGRALS_LIBINT2EINTS_HPP_
#define _AQUARIUS_INTEGRALS_LIBINT2EINTS_HPP_

#include "frameworks/util.hpp"
#include "frameworks/integrals.hpp"
#include "frameworks/molecule.hpp"
#include "agora/integrals/integrals.hpp"

#include "libint2.h"

namespace aquarius
{
namespace integrals
{

class Libint2eIntegrals : public ERI
{
    protected:
        vector<Libint_t> inteval;

        void prims(const vec3& posa, int la, const vector<double>& za,
                   const vec3& posb, int lb, const vector<double>& zb,
                   const vec3& posc, int lc, const vector<double>& zc,
                   const vec3& posd, int ld, const vector<double>& zd,
                   double* integrals) override;
};

REGISTER_VENDOR(Integrals, ERI, Libint2eIntegrals);

}
}

#endif
