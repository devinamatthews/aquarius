#ifndef _AQUARIUS_INTEGRALS_2EINTS_HPP_
#define _AQUARIUS_INTEGRALS_2EINTS_HPP_

#include "util/global.hpp"

#include "symmetry/symmetry.hpp"
#include "task/task.hpp"
#include "input/molecule.hpp"
#include "input/config.hpp"

#include "shell.hpp"

#define TMP_BUFSIZE 65536
#define INTEGRAL_CUTOFF 1e-14

#define IDX_EQ(i,r,e,j,s,f) ((i) == (j) && (r) == (s) && (e) == (f))
#define IDX_GE(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >= (f)))))
#define IDX_GT(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >  (f)))))

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
        double accuracy_;

    public:
        TwoElectronIntegrals(const Shell& a, const Shell& b, const Shell& c, const Shell& d);

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

class ERI : public task::Destructible, public Distributed
{
    public:
        const symmetry::PointGroup& group;
        deque<double> ints;
        deque<idx4_t> idxs;

        ERI(const Arena& arena, const symmetry::PointGroup& group) : Distributed(arena), group(group) {}

        void print(task::Printer& p) const;
};

template <typename ERIType>
class TwoElectronIntegralsTask : public task::Task
{
    public:
        TwoElectronIntegralsTask(const string& name, input::Config& config)
        : task::Task(name, config)
        {
            vector<task::Requirement> reqs;
            reqs.push_back(task::Requirement("molecule", "molecule"));
            addProduct(task::Product("eri", "I", reqs));
        }

        bool run(task::TaskDAG& dag, const Arena& arena)
        {
            const auto& molecule = get<input::Molecule>("molecule");

            ERI* eri = new ERI(arena, molecule.getGroup());

            Context ctx(Context::ISCF);

            vector<double> tmpval(TMP_BUFSIZE);
            vector<idx4_t> tmpidx(TMP_BUFSIZE);

            const vector<int>& N = molecule.getNumOrbitals();
            int nirrep = molecule.getGroup().getNumIrreps();

            vector<vector<int>> idx = Shell::setupIndices(Context(), molecule);
            vector<Shell> shells(molecule.getShellsBegin(), molecule.getShellsEnd());

            int abcd = 0;
            for (int a = 0;a < shells.size();++a)
            {
                for (int b = 0;b <= a;++b)
                {
                    for (int c = 0;c <= a;++c)
                    {
                        int dmax = c;
                        if (a == c) dmax = b;
                        for (int d = 0;d <= dmax;++d)
                        {
                            if (abcd%arena.size == arena.rank)
                            {
                                ERIType block(shells[a], shells[b], shells[c], shells[d]);
                                block.run();

                                size_t n;
                                while ((n = block.process(ctx, idx[a], idx[b], idx[c], idx[d],
                                                          TMP_BUFSIZE, tmpval.data(), tmpidx.data(), INTEGRAL_CUTOFF)) != 0)
                                {
                                    eri->ints.insert(eri->ints.end(), tmpval.data(), tmpval.data()+n);
                                    eri->idxs.insert(eri->idxs.end(), tmpidx.data(), tmpidx.data()+n);
                                }
                            }
                            abcd++;
                        }
                    }
                }
            }

            //TODO: load balance

            for (int i = 0;i < eri->ints.size();++i)
            {
                if (eri->idxs[i].i  > eri->idxs[i].j) swap(eri->idxs[i].i, eri->idxs[i].j);
                if (eri->idxs[i].k  > eri->idxs[i].l) swap(eri->idxs[i].k, eri->idxs[i].l);
                if (eri->idxs[i].i  > eri->idxs[i].k ||
                   (eri->idxs[i].i == eri->idxs[i].k &&
                    eri->idxs[i].j  > eri->idxs[i].l))
                {
                    swap(eri->idxs[i].i, eri->idxs[i].k);
                    swap(eri->idxs[i].j, eri->idxs[i].l);
                }
            }

            put("I", eri);

            return true;
        }
};

class OSERI;
using OS2eIntegralsTask = TwoElectronIntegralsTask<OSERI>;

}
}

#endif
