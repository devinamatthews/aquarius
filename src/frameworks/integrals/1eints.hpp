#ifndef _AQUARIUS_FRAMEWORKS_INTEGRALS_1EINTS_HPP_
#define _AQUARIUS_FRAMEWORKS_INTEGRALS_1EINTS_HPP_

#include "frameworks/util.hpp"
#include "frameworks/symmetry.hpp"
#include "frameworks/molecule.hpp"

#include "context.hpp"

namespace aquarius
{

struct idx2_t
{
    uint16_t i = 0, j = 0;
};

namespace integrals
{

void transform(size_t nother, const matrix<double>& xa, const matrix<double>& xb,
               double* buf1, double* buf2);

void transform(const matrix<double>& xa, const matrix<double>& xb, size_t nother,
               double* buf1, double* buf2);

class OneElectronIntegrals
{
    public:
        class ShellBlock
        {
            friend class OneElectronIntegrals;

            public:
                virtual ~ShellBlock() {}

                size_t process(const Context& ctx, const vector<int>& idxa,
                               const vector<int>& idxb, size_t nprocess,
                               double* integrals, idx2_t* indices,
                               double cutoff = -1);

                Context::Ordering getOrdering() const { return ordering; }

                const vector<double>& getIntegrals() const { return ints; }

            protected:
                const molecule::Shell& a;
                const molecule::Shell& b;
                Context::Ordering ordering;
                vector<double> ints;
                size_t num_processed = 0;

                ShellBlock(const molecule::Shell& a, const molecule::Shell& b,
                           Context::Ordering ordering, vector<double>&& integrals)
                : a(a), b(b), ordering(ordering), ints(move(integrals)) {}
        };

        virtual ~OneElectronIntegrals() {}

        ShellBlock calculate(const molecule::Shell& a, const molecule::Shell& b);

    protected:
        virtual void prim(const vec3& posa, int la, double za,
                          const vec3& posb, int lb, double zb,
                          double* integrals);

        virtual void prims(const vec3& posa, int la, const row<double>& za,
                           const vec3& posb, int lb, const row<double>& zb,
                           double* integrals);

        virtual void contr(const vec3& posa, int la, const row<double>& za, const matrix<double>& ca,
                           const vec3& posb, int lb, const row<double>& zb, const matrix<double>& cb,
                           double* integrals);

        virtual void spher(const vec3& posa, int la, const row<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                           const vec3& posb, int lb, const row<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
                           double* integrals);

        virtual void so(const molecule::Shell& a, const molecule::Shell& b, double* integrals);
};

}
}

#endif
