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

#ifndef _AQUARIUS_INTEGRALS_1EINTS_HPP_
#define _AQUARIUS_INTEGRALS_1EINTS_HPP_

#include <cstddef>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstring>

#include "memory/memory.h"
#include "util/math_ext.h"
#include "util/blas.h"
#include "util/stl_ext.hpp"
#include "symmetry/symmetry.hpp"
#include "tensor/symblocked_tensor.hpp"
#include "task/task.hpp"
#include "input/molecule.hpp"
#include "input/config.hpp"

#include "shell.hpp"

namespace aquarius
{

struct idx2_t
{
    uint16_t i;
    uint16_t j;

    idx2_t() : i(0), j(0) {}

    idx2_t(uint16_t i, uint16_t j) : i(i), j(j) {}
};

namespace integrals
{

class OneElectronIntegralEvaluator
{
    public:
				virtual ~OneElectronIntegralEvaluator() {}

        virtual void operator()(int la, const double* ca, int na, const double *za,
                                int lb, const double* cb, int nb, const double *zb,
                                double *ints) const = 0;
};

class OVIEvaluator : public OneElectronIntegralEvaluator
{
    public:
        void operator()(int la, const double* ca, int na, const double *za,
                        int lb, const double* cb, int nb, const double *zb,
                        double *ints) const;
};

class KEIEvaluator : public OneElectronIntegralEvaluator
{
    public:
        void operator()(int la, const double* ca, int na, const double *za,
                        int lb, const double* cb, int nb, const double *zb,
                        double *ints) const;
};

class NAIEvaluator : public OneElectronIntegralEvaluator
{
    protected:
        const std::vector<Center>& centers;

    public:
        NAIEvaluator(const std::vector<Center>& centers) : centers(centers) {}

        void operator()(int la, const double* ca, int na, const double *za,
                        int lb, const double* cb, int nb, const double *zb,
                        double *ints) const;
};

class OneElectronHamiltonianEvaluator : public OneElectronIntegralEvaluator
{
    protected:
        const std::vector<Center>& centers;

    public:
        OneElectronHamiltonianEvaluator(const std::vector<Center>& centers) : centers(centers) {}

        void operator()(int la, const double* ca, int na, const double *za,
                        int lb, const double* cb, int nb, const double *zb,
                        double *ints) const;
};

class OneElectronIntegrals
{
    protected:
        const Shell &a, &b;
        double *ints;
        const OneElectronIntegralEvaluator& eval;
        size_t nints;
        size_t num_processed;

    public:
        OneElectronIntegrals(const Shell& a, const Shell& b, const OneElectronIntegralEvaluator& eval);

        ~OneElectronIntegrals();

        size_t getNumInts() const { return nints; }

        const double* getIntegrals() const { return ints; }

        size_t process(const Context& ctx, const std::vector<int>& idxa, const std::vector<int>& idxb,
                       size_t nprocess, double* integrals, idx2_t* indices, double cutoff = -1);

    protected:
        void ao2so2(size_t nother, int r, double* aointegrals, double* sointegrals);

        void cart2spher2r(size_t nother, double* buf1, double* buf2);

        void cart2spher2l(size_t nother, double* buf1, double* buf2);

        void prim2contr2r(size_t nother, double* buf1, double* buf2);

        void prim2contr2l(size_t nother, double* buf1, double* buf2);
};

class OneElectronIntegralsTask : public task::Task
{
    public:
        class OneElectronIntegral : public tensor::SymmetryBlockedTensor<double>
        {
            protected:
                std::string name;

            public:
                OneElectronIntegral(const Arena& arena, const std::string& name,
                                    const symmetry::PointGroup& group, const std::vector<int>& norb)
                : tensor::SymmetryBlockedTensor<double>(name, arena, group, 2, {norb,norb}, {NS,NS}, true),
                  name(name) {}
        };

        OneElectronIntegralsTask(const std::string& name, const input::Config& config);

        void run(task::TaskDAG& dag, const Arena& arena);
};

struct OVI : public OneElectronIntegralsTask::OneElectronIntegral
{
    OVI(const Arena& arena, const symmetry::PointGroup& group, const std::vector<int>& norb)
    : OneElectronIntegral(arena, "S", group, norb) {}
};

struct NAI : public OneElectronIntegralsTask::OneElectronIntegral
{
    NAI(const Arena& arena, const symmetry::PointGroup& group, const std::vector<int>& norb)
    : OneElectronIntegral(arena, "G", group, norb) {}
};

struct KEI : public OneElectronIntegralsTask::OneElectronIntegral
{
    KEI(const Arena& arena, const symmetry::PointGroup& group, const std::vector<int>& norb)
    : OneElectronIntegral(arena, "T", group, norb) {}
};

struct OneElectronHamiltonian : public OneElectronIntegralsTask::OneElectronIntegral
{
    OneElectronHamiltonian(const Arena& arena, const symmetry::PointGroup& group, const std::vector<int>& norb)
    : OneElectronIntegral(arena, "H", group, norb) {}
};

}
}

namespace std
{
    template<> struct is_pod<aquarius::idx2_t> { static const bool value = true; };
}

#endif
