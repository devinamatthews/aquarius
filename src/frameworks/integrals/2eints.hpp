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

class OSERI;

class TwoElectronIntegrals
{
    protected:
        const molecule::Shell& sa;
        const molecule::Shell& sb;
        const molecule::Shell& sc;
        const molecule::Shell& sd;
        const symmetry::PointGroup& group;
        const molecule::Center& ca;
        const molecule::Center& cb;
        const molecule::Center& cc;
        const molecule::Center& cd;
        int la, lb, lc, ld;
        int na, nb, nc, nd;
        int ma, mb, mc, md;
        int da, db, dc, dd;
        int fca, fcb, fcc, fcd;
        int fsa, fsb, fsc, fsd;
        const vector<double>& za;
        const vector<double>& zb;
        const vector<double>& zc;
        const vector<double>& zd;
        vector<double> ints;
        size_t num_processed;
        double accuracy_;

    public:
        TwoElectronIntegrals(const molecule::Shell& a, const molecule::Shell& b,
                             const molecule::Shell& c, const molecule::Shell& d);

        virtual ~TwoElectronIntegrals() {}

        void run();

        const vector<double>& getIntegrals() const { return ints; }

        size_t process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                       const vector<int>& idxc, const vector<int>& idxd,
                       size_t nprocess, double* integrals, idx4_t* indices, double cutoff = -1);

        double accuracy() const { return accuracy_; }

        void accuracy(double val) { accuracy_ = val; }

    protected:
        virtual void prim(const vec3& posa, int e, const vec3& posb, int f,
                          const vec3& posc, int g, const vec3& posd, int h, double* integrals);

        virtual void prims(const vec3& posa, const vec3& posb, const vec3& posc, const vec3& posd,
                           double* integrals);

        virtual void contr(const vec3& posa, const vec3& posb, const vec3& posc, const vec3& posd,
                           double* integrals);

        virtual void spher(const vec3& posa, const vec3& posb, const vec3& posc, const vec3& posd,
                           double* integrals);

        virtual void so(double* integrals);

        void ao2so4(size_t nother, int r, int t, int st, double* aointegrals, double* sointegrals);

        void cart2spher4r(size_t nother, double* buf1, double* buf2);

        void cart2spher4l(size_t nother, double* buf1, double* buf2);

        void prim2contr4r(size_t nother, double* buf1, double* buf2);

        void prim2contr4l(size_t nother, double* buf1, double* buf2);
};

class ERI : public Distributed
{
    public:
        const symmetry::PointGroup& group;
        deque<double> ints;
        deque<idx4_t> idxs;

        ERI(const Arena& arena, const symmetry::PointGroup& group) : Distributed(arena), group(group) {}
};

}
}

#endif
