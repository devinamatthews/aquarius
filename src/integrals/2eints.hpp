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

#ifndef _AQUARIUS_INTEGRALS_2EINTS_HPP_
#define _AQUARIUS_INTEGRALS_2EINTS_HPP_

#include <cstddef>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <algorithm>

#include "memory/memory.h"
#include "util/math_ext.h"
#include "util/blas.h"
#include "util/stl_ext.hpp"
#include "symmetry/symmetry.hpp"
#include "task/task.hpp"
#include "input/molecule.hpp"
#include "input/config.hpp"

#include "shell.hpp"

namespace aquarius
{

struct idx4_t
{
    uint16_t i;
    uint16_t j;
    uint16_t k;
    uint16_t l;

    idx4_t() : i(0), j(0), k(0), l(0) {}

    idx4_t(uint16_t i, uint16_t j, uint16_t k, uint16_t l) : i(i), j(j), k(k), l(l) {}
};

namespace integrals
{

class TwoElectronIntegralEvaluator
{
    public:
        virtual void operator()(int la, const double* ca, int na, const double *za,
                                int lb, const double* cb, int nb, const double *zb,
                                int lc, const double* cc, int nc, const double *zc,
                                int ld, const double* cd, int nd, const double *zd,
                                double *ints) const = 0;
};

class ERIEvaluator : public TwoElectronIntegralEvaluator
{
    public:
        void operator()(int la, const double* ca, int na, const double *za,
                        int lb, const double* cb, int nb, const double *zb,
                        int lc, const double* cc, int nc, const double *zc,
                        int ld, const double* cd, int nd, const double *zd,
                        double *ints) const;
};

class TwoElectronIntegrals
{
    protected:
        const Shell &a, &b, &c, &d;
        double *ints;
        const TwoElectronIntegralEvaluator& eval;
        size_t nints;
        size_t num_processed;

    public:
        TwoElectronIntegrals(const Shell& a, const Shell& b, const Shell& c, const Shell& d,
                             const TwoElectronIntegralEvaluator& eval);

        ~TwoElectronIntegrals();

        size_t getNumInts() const { return nints; }

        const double* getIntegrals() const { return ints; }

        size_t process(const Context& ctx, const std::vector<int>& idxa, const std::vector<int>& idxb,
                       const std::vector<int>& idxc, const std::vector<int>& idxd,
                       size_t nprocess, double* integrals, idx4_t* indices, double cutoff = -1);

    protected:
        void ao2so4(size_t nother, int r, int t, int st, double* aointegrals, double* sointegrals);

        void cart2spher4r(size_t nother, double* buf1, double* buf2);

        void cart2spher4l(size_t nother, double* buf1, double* buf2);

        void prim2contr4r(size_t nother, double* buf1, double* buf2);

        void prim2contr4l(size_t nother, double* buf1, double* buf2);
};

class ERI : public task::Resource
{
    public:
        std::vector<double> ints;
        std::vector<idx4_t> idxs;

        ERI(const Arena& arena) : Resource(arena) {}

        void print(task::Printer& p) const;
};

class TwoElectronIntegralsTask : public task::Task
{
    public:
        TwoElectronIntegralsTask(const std::string& name, const input::Config& config);

        void run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
