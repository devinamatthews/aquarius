#ifndef _AQUARIUS_INTEGRALS_2EINTS_HPP_
#define _AQUARIUS_INTEGRALS_2EINTS_HPP_

#include "util/global.hpp"

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

class TwoElectronIntegrals
{
    protected:
        const Shell& sa;
        const Shell& sb;
        const Shell& sc;
        const Shell& sd;
        const symmetry::PointGroup& group;
        const Center& ca;
        const Center& cb;
        const Center& cc;
        const Center& cd;
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

    public:
        TwoElectronIntegrals(const Shell& a, const Shell& b, const Shell& c, const Shell& d);

        virtual ~TwoElectronIntegrals() {}

        void run();

        const vector<double>& getIntegrals() const { return ints; }

        size_t process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                       const vector<int>& idxc, const vector<int>& idxd,
                       size_t nprocess, double* integrals, idx4_t* indices, double cutoff = -1);

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

class ERI : public task::Destructible, public Distributed
{
    public:
        const symmetry::PointGroup& group;
        vector<double> ints;
        vector<idx4_t> idxs;

        ERI(const Arena& arena, const symmetry::PointGroup& group) : Distributed(arena), group(group) {}

        void print(task::Printer& p) const;
};

class TwoElectronIntegralsTask : public task::Task
{
    public:
        TwoElectronIntegralsTask(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
