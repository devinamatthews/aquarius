#ifndef _AQUARIUS_FRAMEWORKS_MOLECULE_SHELL_HPP_
#define _AQUARIUS_FRAMEWORKS_MOLECULE_SHELL_HPP_

#include "center.hpp"

namespace aquarius
{
namespace molecule
{

class Shell
{
    protected:
        Center center;
        int L;
        vector<int> nfunc_per_irrep;
        bool spherical;
        bool keep_contaminants;
        matrix<int> func_irrep;
        matrix<int> irrep_pos;
        matrix<int> irreps;
        vector<double> exponents;
        matrix<double> coefficients;
        matrix<int> parity;
        matrix<double> cart2spher;

    public:
        Shell(const Center& pos, int L, bool spherical, bool keep_contaminants,
              const vector<double>& exponents, const matrix<double>& coefficients);

        const Center& getCenter() const { return center; }

        int getL() const { return L; }

        int getDegeneracy() const { return center.getCenters().size(); }

        int getNPrim() const { return coefficients.length(0); }

        int getNContr() const { return coefficients.length(1); }

        int getNFunc() const { return (spherical && !keep_contaminants ? (2*L+1) : (L+1)*(L+2)/2); }

        int getNFuncInIrrep(int irrep) const { return nfunc_per_irrep[irrep]; }

        int getPosInIrrep(int func, int irrep) const { return irrep_pos[func][irrep]; }

        const vector<int>& getNFuncInEachIrrep() const { return nfunc_per_irrep; }

        bool isSpherical() const { return spherical; }

        bool getContaminants() const { return keep_contaminants; }

        int getIrrepOfFunc(int func, int degen) const { return irreps[func][degen]; }

        const vector<int>& getIrrepsOfFunc(int func) const { return irreps[func]; }

        const vector<double>& getExponents() const { return exponents; }

        const matrix<double>& getCoefficients() const { return coefficients; }

        int getParity(int func, int op) const { return parity[func][op]; }

        const matrix<double>& getCart2Spher() const { return cart2spher; }

    protected:
        static double cartcoef(int l, int m, int lx, int ly, int lz);
};

}
}

#endif
