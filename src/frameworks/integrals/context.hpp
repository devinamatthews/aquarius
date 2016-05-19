#ifndef _AQUARIUS_INTEGRALS_CONTEXT_HPP_
#define _AQUARIUS_INTEGRALS_CONTEXT_HPP_

#include "util/global.hpp"

#include "symmetry/symmetry.hpp"

namespace aquarius
{
namespace integrals
{

class Context
{
    public:
        enum Ordering {ISCF,ISFC,SICF,SIFC,SCIF,SCFI,SFIC,SFCI};

    protected:
        mutable vector<vector<int>> cartesian_ordering;
        mutable vector<vector<int>> spherical_ordering;
        Ordering ordering;

    public:
        Context(Ordering ordering = ISCF)
        : ordering(ordering) {}

        void setCartesianOrdering(int L, const vector<int>& ordering);

        const vector<int>& getCartesianOrdering(int L) const;

        void setSphericalOrdering(int L, const vector<int>& ordering);

        const vector<int>& getSphericalOrdering(int L) const;

        void setOrdering(Ordering ordering) { this->ordering = ordering; }

        Ordering getOrdering() const { return ordering; }
};

}
}

#endif
