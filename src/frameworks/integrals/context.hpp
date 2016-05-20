#ifndef _AQUARIUS_FRAMEWORKS_INTEGRALS_CONTEXT_HPP_
#define _AQUARIUS_FRAMEWORKS_INTEGRALS_CONTEXT_HPP_

#include "frameworks/util.hpp"
#include "frameworks/molecule.hpp"

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

        vector<vector<int>> setupIndices(const molecule::Molecule& m);

        int getIndex(const molecule::Shell& shell, const vector<int>& idx,
                     int func, int contr, int degen) const;

        //void aoToSo(Context::Ordering primitive_ordering, double* aoso, int ld) const;
};

}
}

#endif