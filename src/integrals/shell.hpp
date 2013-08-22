/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_INTEGRALS_SHELL_HPP_
#define _AQUARIUS_INTEGRALS_SHELL_HPP_

#include <cstddef>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstring>

#include "memory/memory.h"
#include "util/math_ext.h"
#include "util/blas.h"
#include "symmetry/symmetry.hpp"

#include "center.hpp"
#include "element.hpp"
#include "context.hpp"

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

    friend void swap(Shell& a, Shell& b)
    {
        using std::swap;
        swap(a.center, b.center);
        swap(a.L, b.L);
        swap(a.nprim, b.nprim);
        swap(a.ncontr, b.ncontr);
        swap(a.nfunc, b.nfunc);
        swap(a.ndegen, b.ndegen);
        swap(a.nfunc_per_irrep, b.nfunc_per_irrep);
        swap(a.spherical, b.spherical);
        swap(a.keep_contaminants, b.keep_contaminants);
        swap(a.func_irrep, b.func_irrep);
        swap(a.irrep_pos, b.irrep_pos);
        swap(a.irreps, b.irreps);
        swap(a.exponents, b.exponents);
        swap(a.coefficients, b.coefficients);
        swap(a.parity, b.parity);
        swap(a.cart2spher, b.cart2spher);
    }

    protected:
        Center center;
        int L;
        int nprim;
        int ncontr;
        int nfunc;
        int ndegen;
        std::vector<int> nfunc_per_irrep;
        bool spherical;
        bool keep_contaminants;
        std::vector<std::vector<int> > func_irrep;
        std::vector<std::vector<int> > irrep_pos;
        std::vector<std::vector<int> > irreps;
        double* exponents;
        double* coefficients;
        std::vector<std::vector<int> > parity;
        double *cart2spher;

    public:
        Shell(const Shell& other);

        Shell(const Shell& other, bool spherical, bool keep_contaminants);

        Shell(const Center& pos, int L, int nprim, int ncontr, bool spherical, bool keep_contaminants,
              const double* exponents, const double* coefficients);

        ~Shell();

        Shell& operator=(Shell other) { swap(*this, other); return *this; }

        static std::vector<std::vector<int> > setupIndices(const Context& ctx, const input::Molecule& m);

        int getIndex(const Context& ctx, std::vector<int> idx, int func, int contr, int degen) const;

        //void aoToSo(Context::Ordering primitive_ordering, double* aoso, int ld) const;

        const Center& getCenter() const { return center; }

        int getL() const { return L; }

        int getDegeneracy() const { return ndegen; }

        int getNPrim() const { return nprim; }

        int getNContr() const { return ncontr; }

        int getNFunc() const { return nfunc; }

        int getNFuncInIrrep(int irrep) const { return nfunc_per_irrep[irrep]; }

        const std::vector<int>& getNFuncInEachIrrep() const { return nfunc_per_irrep; }

        bool isSpherical() const { return spherical; }

        bool getContaminants() const { return keep_contaminants; }

        int getIrrepOfFunc(int func, int degen) const { return irreps[func][degen]; }

        const std::vector<int>& getIrrepsOfFunc(int func) const { return irreps[func]; }

        const double* getExponents() const { return exponents; }

        const double* getCoefficients() const { return coefficients; }

        int getParity(int func, int op) const { return parity[func][op]; }

        const double* getCart2Spher() const { return cart2spher; }

    protected:
        static double cartcoef(int l, int m, int lx, int ly, int lz);
};

}
}

#endif
