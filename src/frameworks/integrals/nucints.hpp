#ifndef _AQUARIUS_FRAMEWORKS_INTEGRALS_NUCINTS_HPP_
#define _AQUARIUS_FRAMEWORKS_INTEGRALS_NUCINTS_HPP_

#include "frameworks/util.hpp"
#include "frameworks/symmetry.hpp"
#include "frameworks/molecule.hpp"

#include "1eints.hpp"

namespace aquarius
{

namespace integrals
{

class NuclearIntegrals
{
    public:
        class ShellBlock
        {
            friend class NuclearIntegrals;

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

        virtual ~NuclearIntegrals() {}

        ShellBlock calculate(const molecule::Shell& a, const molecule::Shell& b,
                             const vector<molecule::Center>& centers);

    protected:
        virtual void prim(const vec3& posa, int la, double za,
                          const vec3& posb, int lb, double zb,
                          const vec3& posc, double charge, double* integrals);

        virtual void prims(const vec3& posa, int la, const row<double>& za,
                           const vec3& posb, int lb, const row<double>& zb,
                           const vec3& posc, double charge, double* integrals);

        virtual void contr(const vec3& posa, int la, const row<double>& za, const matrix<double>& ca,
                           const vec3& posb, int lb, const row<double>& zb, const matrix<double>& cb,
                           const vec3& posc, double charge, double* integrals);

        virtual void spher(const vec3& posa, int la, const row<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                           const vec3& posb, int lb, const row<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
                           const vec3& posc, double charge, double* integrals);

        virtual void so(const molecule::Shell& a, const molecule::Shell& b,
                        const vector<molecule::Center>& centers, double* integrals);
};

}
}

#endif
