#ifndef _AQUARIUS_FRAMEWORKS_INTEGRALS_2EINTS_HPP_
#define _AQUARIUS_FRAMEWORKS_INTEGRALS_2EINTS_HPP_

#include "frameworks/util.hpp"
#include "frameworks/symmetry.hpp"
#include "frameworks/molecule.hpp"

#include "context.hpp"

#define IDX_EQ(i,r,e,j,s,f) ((i) == (j) && (r) == (s) && (e) == (f))
#define IDX_GE(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >= (f)))))
#define IDX_GT(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >  (f)))))

namespace aquarius
{

struct idx4_t
{
    uint16_t i = 0, j = 0, k = 0, l = 0;
};

namespace integrals
{

void transform(const matrix<double>& ai, const matrix<double>& bj,
               const matrix<double>& ck, const matrix<double>& dl,
               size_t nother, double* buf1, double* buf2);

void transform(size_t nother,
               const matrix<double>& ai, const matrix<double>& bj,
               const matrix<double>& ck, const matrix<double>& dl,
               double* buf1, double* buf2);

class TwoElectronIntegrals
{
    public:
        class ShellBlock
        {
            friend class TwoElectronIntegrals;

            public:
                const vector<double>& getIntegrals() const { return ints; }

                size_t process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                               const vector<int>& idxc, const vector<int>& idxd,
                               size_t nprocess, double* integrals, idx4_t* indices, double cutoff = -1);

                Context::Ordering getOrdering() const { return ordering; }

                double accuracy() const { return accuracy_; }

                void accuracy(double val) { accuracy_ = val; }

            protected:
                const molecule::Shell& a;
                const molecule::Shell& b;
                const molecule::Shell& c;
                const molecule::Shell& d;
                Context::Ordering ordering;
                vector<double> ints;
                size_t num_processed = 0;
                double accuracy_ = 1e-15;

                ShellBlock(const molecule::Shell& a, const molecule::Shell& b,
                           const molecule::Shell& c, const molecule::Shell& d,
                           Context::Ordering ordering, vector<double>&& ints)
                : a(a), b(b), c(c), d(d), ordering(ordering), ints(move(ints)) {}
        };

        virtual ~TwoElectronIntegrals() {}

        ShellBlock calculate(const molecule::Shell& a, const molecule::Shell& b,
                             const molecule::Shell& c, const molecule::Shell& d);

        double accuracy() const { return accuracy_; }

        void accuracy(double val) { accuracy_ = val; }

    protected:
        double accuracy_ = 1e-15;

        virtual void prim(const vec3& posa, int la, double za,
                          const vec3& posb, int lb, double zb,
                          const vec3& posc, int lc, double zc,
                          const vec3& posd, int ld, double zd,
                          double* integrals);

        virtual void prims(const vec3& posa, int la, const vector<double>& za,
                           const vec3& posb, int lb, const vector<double>& zb,
                           const vec3& posc, int lc, const vector<double>& zc,
                           const vec3& posd, int ld, const vector<double>& zd,
                           double* integrals);

        virtual void contr(const vec3& posa, int la, const vector<double>& za, const matrix<double>& ca,
                           const vec3& posb, int lb, const vector<double>& zb, const matrix<double>& cb,
                           const vec3& posc, int lc, const vector<double>& zc, const matrix<double>& cc,
                           const vec3& posd, int ld, const vector<double>& zd, const matrix<double>& cd,
                           double* integrals);

        virtual void spher(const vec3& posa, int la, const vector<double>& za, const matrix<double>& ca, const matrix<double>& sa,
                           const vec3& posb, int lb, const vector<double>& zb, const matrix<double>& cb, const matrix<double>& sb,
                           const vec3& posc, int lc, const vector<double>& zc, const matrix<double>& cc, const matrix<double>& sc,
                           const vec3& posd, int ld, const vector<double>& zd, const matrix<double>& cd, const matrix<double>& sd,
                           double* integrals);

        virtual void so(const molecule::Shell& a, const molecule::Shell& b,
                        const molecule::Shell& c, const molecule::Shell& d,
                        double* integrals);
};

class AOintegrals : public Distributed
{
    public:
        const symmetry::PointGroup& group;
        deque<double> ints;
        deque<idx4_t> idxs;

        AOintegrals(const Arena& arena, const symmetry::PointGroup& group) : Distributed(arena), group(group) {}
};

}
}

#endif
