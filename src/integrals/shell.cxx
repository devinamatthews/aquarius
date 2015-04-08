#include "shell.hpp"

#include "input/molecule.hpp"

using namespace aquarius::input;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace integrals
{

Shell::Shell(const Center& pos, int L, int nprim, int ncontr, bool spherical, bool keep_contaminants,
             const vector<double>& exponents, const vector<double>& coefficients)
: center(pos), L(L), nprim(nprim), ncontr(ncontr), nfunc(spherical && !keep_contaminants ? 2*L+1 : (L+1)*(L+2)/2),
  ndegen(pos.getCenters().size()), spherical(spherical), keep_contaminants(keep_contaminants),
  exponents(exponents), coefficients(coefficients)
{
    const PointGroup& group = pos.getPointGroup();
    int nirrep = group.getNumIrreps();
    int order = group.getOrder();

    vector<int> proj(order);

    irreps.resize(nfunc, vector<int>(nirrep));
    func_irrep.resize(nfunc, vector<int>(nirrep));
    irrep_pos.resize(nfunc, vector<int>(nirrep));

    nfunc_per_irrep.resize(nirrep, 0);

    parity.resize(nfunc, vector<int>(order));

    for (int op = 0;op < order;op++)
    {
        if (spherical)
        {
            int f = 0;
            for (int l = L;l >= (keep_contaminants ? 0 : L);l -= 2)
            {
                for (int m = l;m > 0;m--)
                {
                    parity[f++][op] = (group.sphericalParity(l,  m, op) < 0 ? -1 : 1);
                    parity[f++][op] = (group.sphericalParity(l, -m, op) < 0 ? -1 : 1);
                }
                parity[f++][op] = (group.sphericalParity(l, 0, op) < 0 ? -1 : 1);
            }
        }
        else
        {
            int f = 0;
            for (int x = L;x >= 0;x--)
            {
                for (int y = L-x;y >= 0;y--)
                {
                    int z = L-x-y;
                    parity[f++][op] = (group.cartesianParity(x, y, z, op) < 0 ? -1 : 1);
                }
            }
        }
    }

    /*
     * the tricky part: determine the irrep of each final SO function
     * each AO function will give n SO functions, where n is the number of symmetry-equivalent atoms associated to this shell
     *
     * loop through the functions and for each irrep, determine the projection of this function onto the degenerate centers
     * if this is non-zero, then this is one of the irreps for this function
     *
     * this could potentially be simplified by using the DCR of the atom's stabilizer, but this works
     */
    for (int func = 0;func < nfunc;func++)
    {
        int i = 0;

        for (int irrep = 0;irrep < order;irrep++)
        {
            for (int j = 0;j < ndegen;j++)
            {
                proj[j] = 0;
            }

            /*
             * do the projection, using the characters of the irrep and the parity of either the cartesian or spherical functions
             */
            for (int op = 0;op < order;op++)
            {
                proj[pos.getCenterAfterOp(op)] += (group.character(irrep, op)*parity[func][op] < 0 ? -1 : 1);
            }

            for (int j = 0;j < ndegen;j++)
            {
                proj[j] /= (order/ndegen);
            }

            /*
             * check if the projection is non-zero
             */
            int nonzero = 0;
            for (int j = 0;j < ndegen;j++)
            {
                nonzero += abs(proj[j]);
            }
            if (nonzero > 0)
            {
                irrep_pos[func][irrep] = i;
                func_irrep[func][i] = nfunc_per_irrep[irrep];
                irreps[func][i] = irrep;
                nfunc_per_irrep[irrep]++;
                i++;
            }
        }
    }

    /*
     * Normalize the shell
     */
    const double PI2_N34 = 0.25197943553838073034791409490358;

    for (int i = 0;i < ncontr;i++)
    {
        double norm = 0.0;
        for (int j = 0;j < nprim;j++)
        {
            for (int k = 0;k < nprim;k++)
            {
                double zeta = sqrt(exponents[j]*exponents[k])/(exponents[j]+exponents[k]);
                norm += coefficients[i*nprim+j]*coefficients[i*nprim+k]*pow(2*zeta,(double)L+1.5);
            }
        }

        for (int j = 0;j < nprim;j++)
        {
            this->coefficients[i*nprim+j] *= PI2_N34*pow(4*exponents[j],((double)L+1.5)/2)/sqrt(norm);
        }
    }

    /*
     * Generate the cartesian -> spherical harmonic transformation
     */
    if (spherical)
    {
        cart2spher.resize(nfunc*(L+1)*(L+2)/2);
        int k = 0;
        for (int l = L;l >= (keep_contaminants ? 0 : L);l -= 2)
        {
            for (int m = l;m > 0;m--)
            {
                for (int x = L;x >= 0;x--)
                    for (int y = L-x;y >= 0;y--)
                        cart2spher[k++] = cartcoef(l, m, x, y, L-x-y);
                for (int x = L;x >= 0;x--)
                    for (int y = L-x;y >= 0;y--)
                        cart2spher[k++] = cartcoef(l, -m, x, y, L-x-y);
            }
            for (int x = L;x >= 0;x--)
                for (int y = L-x;y >= 0;y--)
                    cart2spher[k++] = cartcoef(l, 0, x, y, L-x-y);
        }
        assert(k == nfunc*(L+1)*(L+2)/2);
    }
}

vector<vector<int>> Shell::setupIndices(const Context& ctx, const Molecule& m)
{
    int nirrep = m.getGroup().getNumIrreps();

    vector<vector<int>> idx;
    vector<int> nfunc(nirrep, (int)0);

    for (Molecule::const_shell_iterator s = m.getShellsBegin();s != m.getShellsEnd();++s)
    {
        idx.push_back(vector<int>(nirrep));
        vector<int>& index = idx.back();

        if (ctx.getOrdering() == Context::ISCF || ctx.getOrdering() == Context::ISFC)
        {
            for (int i = 0;i < nirrep;i++) index[i] = nfunc[i];
        }
        else
        {
            if (ctx.getOrdering() == Context::SICF || ctx.getOrdering() == Context::SIFC)
            {
                for (int i = 0;i < nirrep-1;i++)
                {
                    index[i+1] = index[i] + s->getNFuncInIrrep(i)*s->getNContr();
                }
            }
            else if (ctx.getOrdering() == Context::SFIC ||
                     ctx.getOrdering() == Context::SFCI ||
                     ctx.getOrdering() == Context::SCFI)
            {
                for (int i = 0;i < nirrep-1;i++)
                {
                    index[i+1] = index[i];
                }
            }
            else if (ctx.getOrdering() == Context::SCIF)
            {
                for (int i = 0;i < nirrep-1;i++)
                {
                    index[i+1] = index[i] + s->getNFuncInIrrep(i);
                }
            }

            int nfunctot = 0;
            for (int i = 0;i < nirrep;i++) nfunctot += nfunc[i];
            for (int i = 0;i < nirrep;i++) index[i] += nfunctot;
        }

        for (int i = 0;i < nirrep;i++) nfunc[i] += s->getNFuncInIrrep(i)*s->getNContr();
    }

    for (int i = 1;i < nirrep;i++) nfunc[i] += nfunc[i-1];

    if (ctx.getOrdering() == Context::ISCF || ctx.getOrdering() == Context::ISFC)
    {
        for (vector<vector<int>>::iterator s = idx.begin();s != idx.end();++s)
        {
            vector<int>& index = *s;
            for (int i = 1;i < nirrep;i++) index[i] += nfunc[i-1];
        }
    }

    return idx;
}

int Shell::getIndex(const Context& ctx, vector<int> idx, int func, int contr, int degen) const
{
    int irrep = irreps[func][degen];

    int ffunc;
    if (spherical)
    {
        ffunc = ctx.getSphericalOrdering(L)[func];
    }
    else
    {
        ffunc = ctx.getCartesianOrdering(L)[func];
    }

    switch (ctx.getOrdering())
    {
        case Context::ISCF:
        case Context::SICF:
            return idx[irrep] + contr*nfunc_per_irrep[irrep] + func_irrep[func][degen];
        case Context::ISFC:
        case Context::SIFC:
            return idx[irrep] + contr + func_irrep[func][degen]*ncontr;
        case Context::SCIF:
            return idx[irrep] + contr*nfunc*ndegen + func_irrep[func][degen];
        case Context::SCFI:
            return idx[irrep] + (contr*nfunc + ffunc)*ndegen + irrep_pos[func][irrep];
        case Context::SFIC:
            return idx[irrep] + contr + (ffunc*ndegen + irrep_pos[func][irrep])*ncontr;
        case Context::SFCI:
            return idx[irrep] + (contr + ffunc*ncontr)*ndegen + irrep_pos[func][irrep];
    }

    return 0;
}

double Shell::cartcoef(int l, int m, int lx, int ly, int lz)
{
    int am = abs(m);
    int j = lx+ly-am;
    if ((j&1) == 1) return 0.0;
    j /= 2;

    double c = sqrt((double)(binom(2*lx,lx)*binom(2*ly,ly)*binom(2*lz,lz)*binom(l+am,am))/
                    (double)(binom(2*l,l)*binom(l,am)*binom(lx+ly+lz,lx)*binom(ly+lz,ly))/
                    (double)(dfact(2*lx-1)*dfact(2*ly-1)*dfact(2*lz-1)))/pow(2.0,l);
    if (m != 0) c *= sqrt(2.0);

    if (m >= 0)
    {
        if (((am-lx)&1) == 1) return 0.0;
        if (((am-lx)&3) == 2) c = -c;
    }
    else
    {
        if (((am-lx)&1) == 0) return 0.0;
        if (((am-lx)&3) == 3) c = -c;
    }

    double sum = 0.0;
    for (int i = 0;i <= (l-am)/2;i++)
    {
        for (int k = 0;k <= j;k++)
        {
            double tmp = binom(2*l-2*i,l+am)*binom(l,i)*binom(i,j)*binom(j,k)*binom(am,lx-2*k);
            if (((i+k)&1) == 1) tmp = -tmp;
            sum += tmp;
        }
    }

    return sum*c;
}

}
}
