#ifndef _AQUARIUS_FRAMEWORKS_INTEGRALS_MOMINTS_HPP_
#define _AQUARIUS_FRAMEWORKS_INTEGRALS_MOMINTS_HPP_

#include "frameworks/util.hpp"
#include "frameworks/symmetry.hpp"
#include "frameworks/molecule.hpp"

namespace aquarius
{

struct idx3_t
{
    uint16_t i = 0, j = 0, k = 0;
};

namespace integrals
{

class MomentIntegrals
{
    public:
        class ShellBlock
        {
            friend class MomentIntegrals;

            public:
                virtual ~ShellBlock() {}

                size_t process(const Context& ctx, const vector<int>& idxa,
                               const vector<int>& idxb, size_t nprocess,
                               double* integrals, idx3_t* indices,
                               double cutoff = -1);

                Context::Ordering getOrdering() const { return ordering; }

                const vector<double>& getIntegrals() const { return ints; }

            protected:
                const molecule::Shell& a;
                const molecule::Shell& b;
                int L;
                Context::Ordering ordering;
                vector<double> ints;
                size_t num_processed = 0;

                ShellBlock(const molecule::Shell& a, const molecule::Shell& b,
                           int L, Context::Ordering ordering, vector<double>&& integrals)
                : a(a), b(b), L(L), ordering(ordering), ints(move(integrals)) {}
        };

        virtual ~MomentIntegrals() {}

        ShellBlock calculate(const molecule::Shell& a, const molecule::Shell& b,
                             const vec3& origin, int L);

    protected:
        virtual void prim(const vec3& posa, int la, double za,
                          const vec3& posb, int lb, double zb,
                          const vec3& posc, int lc, double* integrals);

        virtual void prims(const vec3& posa, int la, const vector<double>& za,
                           const vec3& posb, int lb, const vector<double>& zb,
                           const vec3& posc, int lc, double* integrals);

        virtual void contr(const vec3& posa, int la, const vector<double>& za, const matrix<double>& ca,
                           const vec3& posb, int lb, const vector<double>& zb, const matrix<double>& cb,
                           const vec3& posc, int lc, double* integrals);

        virtual void spher(const vec3& posa, int la, const vector<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                           const vec3& posb, int lb, const vector<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
                           const vec3& posc, int lc, double* integrals);

        virtual void so(const molecule::Shell& a, const molecule::Shell& b,
                        const vec3& posc, int lc, double* integrals);
};

}
}

#endif
