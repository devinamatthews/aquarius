#ifndef _AQUARIUS_INTEGRALS_SHELL_HPP_
#define _AQUARIUS_INTEGRALS_SHELL_HPP_

#include "util/global.hpp"

#include "symmetry/symmetry.hpp"

#include "center.hpp"
#include "element.hpp"
#include "context.hpp"

#define XYZ(x,y,z) ((x)*((x)+1)/2 + (x)*(1+(y)+(z)) + (y))

namespace aquarius
{
namespace input
{
class Molecule;
}

namespace integrals
{

class Shell
{
    friend class Context;

    protected:
        Center center;
        int L;
        int nprim;
        int ncontr;
        int nfunc;
        int ndegen;
        vector<int> nfunc_per_irrep;
        bool spherical;
        bool keep_contaminants;
        vector<vector<int> > func_irrep;
        vector<vector<int> > irrep_pos;
        vector<vector<int> > irreps;
        vector<double> exponents;
        vector<double> coefficients;
        vector<vector<int> > parity;
        vector<double> cart2spher;

    public:
        Shell(const Center& pos, int L, int nprim, int ncontr, bool spherical, bool keep_contaminants,
              const vector<double>& exponents, const vector<double>& coefficients);

        static vector<vector<int> > setupIndices(const Context& ctx, const input::Molecule& m);

        int getIndex(const Context& ctx, vector<int> idx, int func, int contr, int degen) const;

        //void aoToSo(Context::Ordering primitive_ordering, double* aoso, int ld) const;

        const Center& getCenter() const { return center; }

        int getL() const { return L; }

        int getDegeneracy() const { return ndegen; }

        int getNPrim() const { return nprim; }

        int getNContr() const { return ncontr; }

        int getNFunc() const { return nfunc; }

        int getNFuncInIrrep(int irrep) const { return nfunc_per_irrep[irrep]; }

        const vector<int>& getNFuncInEachIrrep() const { return nfunc_per_irrep; }

        bool isSpherical() const { return spherical; }

        bool getContaminants() const { return keep_contaminants; }

        int getIrrepOfFunc(int func, int degen) const { return irreps[func][degen]; }

        const vector<int>& getIrrepsOfFunc(int func) const { return irreps[func]; }

        const vector<double>& getExponents() const { return exponents; }

        const vector<double>& getCoefficients() const { return coefficients; }

        int getParity(int func, int op) const { return parity[func][op]; }

        const vector<double>& getCart2Spher() const { return cart2spher; }

    protected:
        static double cartcoef(int l, int m, int lx, int ly, int lz);
};

}
}

#endif
