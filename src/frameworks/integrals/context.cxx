#include "context.hpp"

using namespace aquarius::molecule;

/**
 * Compute the index of a function in cartesian angular momentum in the order:
 *
 * xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
 */
#define FUNC_CART(x,y,z) ((((y)+(z))*((y)+(z)+1))/2 + (z))
#define XYZ(x,y,z) FUNC_CART(x,y,z)

/**
 * Compute the index of a function in spherical harmonic angular momentum.
 *
 * Regular spherical harmonics are referenced by n=l, l>=m>=-l. Contaminants may also be referenced by
 * n>l>=0, n-l even.
 */
#define FUNC_SPHER(n,l,m) ((((n)-(l))*((n)+(l)-1))/2 + 2*(n) + ((m) > 0 ? -2*(m) : 2*(m)+1))

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

vector<vector<int>> Context::setupIndices(const Molecule& m)
{
    int nirrep = m.getGroup().getNumIrreps();

    vector<vector<int>> idx;
    vector<int> nfunc(nirrep);

    for (auto s = m.getShellsBegin();s != m.getShellsEnd();++s)
    {
        idx.push_back(vector<int>(nirrep));
        vector<int>& index = idx.back();

        if (ordering == ISCF || ordering == ISFC)
        {
            for (int i = 0;i < nirrep;i++) index[i] = nfunc[i];
        }
        else
        {
            if (ordering == SICF || ordering == SIFC)
            {
                for (int i = 0;i < nirrep-1;i++)
                {
                    index[i+1] = index[i] + s->getNFuncInIrrep(i)*s->getNContr();
                }
            }
            else if (ordering == SFIC ||
                     ordering == SFCI ||
                     ordering == SCFI)
            {
                for (int i = 0;i < nirrep-1;i++)
                {
                    index[i+1] = index[i];
                }
            }
            else if (ordering == SCIF)
            {
                for (int i = 0;i < nirrep-1;i++)
                {
                    index[i+1] = index[i] + s->getNFuncInIrrep(i);
                }
            }

            int nfunctot = sum(nfunc);
            for (int i = 0;i < nirrep;i++) index[i] += nfunctot;
        }

        for (int i = 0;i < nirrep;i++) nfunc[i] += s->getNFuncInIrrep(i)*s->getNContr();
    }

    for (int i = 1;i < nirrep;i++) nfunc[i] += nfunc[i-1];

    if (ordering == ISCF || ordering == ISFC)
    {
        for (auto& index : idx)
        {
            for (int i = 1;i < nirrep;i++) index[i] += nfunc[i-1];
        }
    }

    return idx;
}

int Context::getIndex(const Shell& shell, const vector<int>& idx,
                      int func, int contr, int degen) const
{
    int irrep = shell.getIrrepOfFunc(func, degen);

    int ffunc;
    if (shell.isSpherical())
    {
        ffunc = getSphericalOrdering(shell.getL())[func];
    }
    else
    {
        ffunc = getCartesianOrdering(shell.getL())[func];
    }

    switch (ordering)
    {
        case ISCF:
        case SICF:
            return idx[irrep] + contr*shell.getNFuncInIrrep(irrep) + shell.getIrrepOfFunc(func, degen);
        case ISFC:
        case SIFC:
            return idx[irrep] + contr + shell.getIrrepOfFunc(func, degen)*shell.getNContr();
        case SCIF:
            return idx[irrep] + contr*shell.getNFunc()*shell.getDegeneracy() + shell.getIrrepOfFunc(func, degen);
        case SCFI:
            return idx[irrep] + (contr*shell.getNFunc() + ffunc)*shell.getDegeneracy() + shell.getPosInIrrep(func, irrep);
        case SFIC:
            return idx[irrep] + contr + (ffunc*shell.getDegeneracy() + shell.getPosInIrrep(func, irrep))*shell.getNContr();
        case SFCI:
            return idx[irrep] + (contr + ffunc*shell.getNContr())*shell.getDegeneracy() + shell.getPosInIrrep(func, irrep);
    }

    return 0;
}

}
}
