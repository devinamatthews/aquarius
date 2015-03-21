#ifndef _AQUARIUS_INTEGRALS_1EINTS_HPP_
#define _AQUARIUS_INTEGRALS_1EINTS_HPP_

#include "util/global.hpp"

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

class OneElectronIntegrals
{
    protected:
        const Shell& sa;
        const Shell& sb;
        const symmetry::PointGroup& group;
        const Center& ca;
        const Center& cb;
        int la, lb;
        int na, nb;
        int ma, mb;
        int da, db;
        int fca, fcb;
        int fsa, fsb;
        const vector<double>& za;
        const vector<double>& zb;
        vector<double> ints;
        size_t num_processed;

    public:
        OneElectronIntegrals(const Shell& a, const Shell& b);

        virtual ~OneElectronIntegrals() {}

        void run();

        const vector<double>& getIntegrals() const { return ints; }

        size_t process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                       size_t nprocess, double* integrals, idx2_t* indices, double cutoff = -1);

    protected:
        virtual void prim(const vec3& posa, int e,
                          const vec3& posb, int f, double* integrals);

        virtual void prims(const vec3& posa, const vec3& posb,
                           double* integrals);

        virtual void contr(const vec3& posa, const vec3& posb,
                           double* integrals);

        virtual void spher(const vec3& posa, const vec3& posb,
                           double* integrals);

        virtual void so(double* integrals);

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
                string name;

            public:
                OneElectronIntegral(const Arena& arena, const string& name,
                                    const symmetry::PointGroup& group, const vector<int>& norb)
                : tensor::SymmetryBlockedTensor<double>(name, arena, group, 2, {norb,norb}, {NS,NS}, true),
                  name(name) {}
        };

        OneElectronIntegralsTask(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

struct OVI : public OneElectronIntegralsTask::OneElectronIntegral
{
    OVI(const Arena& arena, const symmetry::PointGroup& group, const vector<int>& norb)
    : OneElectronIntegral(arena, "S", group, norb) {}
};

struct NAI : public OneElectronIntegralsTask::OneElectronIntegral
{
    NAI(const Arena& arena, const symmetry::PointGroup& group, const vector<int>& norb)
    : OneElectronIntegral(arena, "G", group, norb) {}
};

struct KEI : public OneElectronIntegralsTask::OneElectronIntegral
{
    KEI(const Arena& arena, const symmetry::PointGroup& group, const vector<int>& norb)
    : OneElectronIntegral(arena, "T", group, norb) {}
};

struct OneElectronHamiltonian : public OneElectronIntegralsTask::OneElectronIntegral
{
    OneElectronHamiltonian(const Arena& arena, const symmetry::PointGroup& group, const vector<int>& norb)
    : OneElectronIntegral(arena, "H", group, norb) {}
};

}
}

#endif
