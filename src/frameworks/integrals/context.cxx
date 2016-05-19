#include "context.hpp"

using namespace aquarius::symmetry;

namespace aquarius
{
namespace integrals
{

void Context::setCartesianOrdering(int L, const vector<int>& ordering)
{
    assert(L >= 0);
    assert(ordering.size() == (L+1)*(L+2)/2);

    if (L >= cartesian_ordering.size())
    {
        int Lmin = cartesian_ordering.size();
        cartesian_ordering.resize(L+1);
        for (int l = Lmin;l < L;l++)
        {
            cartesian_ordering[l].resize((l+1)*(l+2)/2);
            for (int c = 0;c < (l+1)*(l+2);c++)
            {
                cartesian_ordering[l][c] = c;
            }
        }
    }

    cartesian_ordering[L] = ordering;
}

const vector<int>& Context::getCartesianOrdering(int L) const
{
    assert(L >= 0);

    if (L >= cartesian_ordering.size())
    {
        int Lmin = cartesian_ordering.size();
        cartesian_ordering.resize(L+1);
        for (int l = Lmin;l <= L;l++)
        {
            cartesian_ordering[l].resize((l+1)*(l+2)/2);
            for (int c = 0;c < (l+1)*(l+2)/2;c++)
            {
                cartesian_ordering[l][c] = c;
            }
        }
    }

    return cartesian_ordering[L];
}

void Context::setSphericalOrdering(int L, const vector<int>& ordering)
{
    assert(L >= 0);
    assert(ordering.size() == 2*L+1);

    if (L >= spherical_ordering.size())
    {
        int Lmin = spherical_ordering.size();
        spherical_ordering.resize(L+1);
        for (int l = Lmin;l < L;l++)
        {
            spherical_ordering[l].resize(2*l+1);
            for (int c = 0;c < 2*l+1;c++)
            {
                spherical_ordering[l][c] = c;
            }
        }
    }

    spherical_ordering[L] = ordering;
}

const vector<int>& Context::getSphericalOrdering(int L) const
{
    assert(L >= 0);

    if (L >= spherical_ordering.size())
    {
        int Lmin = spherical_ordering.size();
        spherical_ordering.resize(L+1);
        for (int l = Lmin;l <= L;l++)
        {
            spherical_ordering[l].resize(2*l+1);
            for (int c = 0;c < 2*l+1;c++)
            {
                spherical_ordering[l][c] = c;
            }
        }
    }

    return spherical_ordering[L];
}

}
}
