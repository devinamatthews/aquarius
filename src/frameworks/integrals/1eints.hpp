#ifndef _AQUARIUS_INTEGRALS_1EINTS_HPP_
#define _AQUARIUS_INTEGRALS_1EINTS_HPP_

#include "frameworks/util.hpp"
#include "frameworks/symmetry.hpp"
#include "frameworks/molecule.hpp"

namespace aquarius
{

struct idx2_t
{
    uint16_t i = 0, j = 0;
};

namespace integrals
{

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
                ShellBlock(const molecule::Shell& a, const molecule::Shell& b,
                           Context::Ordering ordering, vector<double>&& integrals)
                : a(a), b(b), ordering(ordering), ints(move(integrals)) {}

                const molecule::Shell& a;
                const molecule::Shell& b;
                Context::Ordering ordering;
                vector<double> ints;
                size_t num_processed = 0;
        };

        virtual ~OneElectronIntegrals() {}

        ShellBlock calculate(const molecule::Shell& a, const molecule::Shell& b);

    protected:
        virtual void prim(const vec3& posa, int la, double za,
                          const vec3& posb, int lb, double zb,
                          double* integrals);

        virtual void prims(const vec3& posa, int la, const vector<double>& za,
                           const vec3& posb, int lb, const vector<double>& zb,
                           double* integrals);

        virtual void contr(const vec3& posa, int la, const vector<double>& za, const matrix<double>& ca,
                           const vec3& posb, int lb, const vector<double>& zb, const matrix<double>& cb,
                           double* integrals);

        virtual void spher(const vec3& posa, int la, const vector<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                           const vec3& posb, int lb, const vector<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
                           double* integrals);

        virtual void so(const molecule::Center& posa, int la, const vector<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                        const molecule::Center& posb, int lb, const vector<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
                        double* integrals);

        void ao2so2(size_t nother, int r, double* aointegrals, double* sointegrals);

        void cart2spher(size_t nother, const matrix<double>& ca, const matrix<double>& cb,
                        double* buf1, double* buf2);

        void cart2spher(const matrix<double>& ca, const matrix<double>& ca, size_t nother,
                        double* buf1, double* buf2);

        void prim2contr(size_t nother, const matrix<double>& sa, const matrix<double>& sb,
                        double* buf1, double* buf2);

        void prim2contr(const matrix<double>& sa, const matrix<double>& sb, size_t nother,
                        double* buf1, double* buf2);
};

}
}

#endif
