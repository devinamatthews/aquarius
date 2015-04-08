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

class IshidaOVI;
class IshidaKEI;
class IshidaNAI;

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

struct OVI : public OneElectronIntegral
{
    OVI(const Arena& arena, const symmetry::PointGroup& group, const vector<int>& norb)
    : OneElectronIntegral(arena, "S", group, norb) {}
};

struct NAI : public OneElectronIntegral
{
    NAI(const Arena& arena, const symmetry::PointGroup& group, const vector<int>& norb)
    : OneElectronIntegral(arena, "G", group, norb) {}
};

struct KEI : public OneElectronIntegral
{
    KEI(const Arena& arena, const symmetry::PointGroup& group, const vector<int>& norb)
    : OneElectronIntegral(arena, "T", group, norb) {}
};

struct OneElectronHamiltonian : public OneElectronIntegral
{
    OneElectronHamiltonian(const Arena& arena, const symmetry::PointGroup& group, const vector<int>& norb)
    : OneElectronIntegral(arena, "H", group, norb) {}
};

template <typename OVIType, typename KEIType, typename NAIType>
class OneElectronIntegralsTask : public task::Task
{
    public:
        OneElectronIntegralsTask(const string& name, input::Config& config)
        : task::Task(name, config)
        {
            vector<task::Requirement> reqs;
            reqs.push_back(task::Requirement("molecule", "molecule"));
            addProduct(task::Product("ovi", "S", reqs));
            addProduct(task::Product("kei", "T", reqs));
            addProduct(task::Product("nai", "G", reqs));
            addProduct(task::Product("1ehamiltonian", "H", reqs));
        }

        bool run(task::TaskDAG& dag, const Arena& arena)
        {
            const input::Molecule& molecule = get<input::Molecule>("molecule");

            Context ctx(Context::ISCF);

            const vector<int>& N = molecule.getNumOrbitals();
            int n = molecule.getGroup().getNumIrreps();

            vector<int> irrep;
            for (int i = 0;i < n;i++) irrep += vector<int>(N[i],i);

            vector<uint16_t> start(n,0);
            for (int i = 1;i < n;i++) start[i] = start[i-1]+N[i-1];

            vector<vector<int>> idx = Shell::setupIndices(ctx, molecule);
            vector<Shell> shells(molecule.getShellsBegin(), molecule.getShellsEnd());
            vector<vector<tkv_pair<double>>> ovi_pairs(n), nai_pairs(n), kei_pairs(n);
            vector<Center> centers;

            for (auto& atom : molecule.getAtoms())
            {
                centers.push_back(atom.getCenter());
            }

            int block = 0;
            for (int a = 0;a < shells.size();++a)
            {
                for (int b = 0;b <= a;++b)
                {
                    if (block%arena.size == arena.rank)
                    {
                        OVIType s(shells[a], shells[b]);
                        KEIType t(shells[a], shells[b]);
                        NAIType g(shells[a], shells[b], centers);

                        s.run();
                        t.run();
                        g.run();

                        size_t nint = s.getIntegrals().size();
                        vector<double> ints(nint);
                        vector<idx2_t> idxs(nint);
                        size_t nproc;

                        nproc = s.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                        for (int k = 0;k < nproc;k++)
                        {
                            int irr = irrep[idxs[k].i];
                            assert(irr == irrep[idxs[k].j]);

                            uint16_t i = idxs[k].i-start[irr];
                            uint16_t j = idxs[k].j-start[irr];

                                        ovi_pairs[irr].push_back(tkv_pair<double>(i*N[irr]+j, ints[k]));
                            if (i != j) ovi_pairs[irr].push_back(tkv_pair<double>(j*N[irr]+i, ints[k]));
                        }

                        nproc = t.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                        for (int k = 0;k < nproc;k++)
                        {
                            int irr = irrep[idxs[k].i];
                            assert(irr == irrep[idxs[k].j]);

                            uint16_t i = idxs[k].i-start[irr];
                            uint16_t j = idxs[k].j-start[irr];

                                        kei_pairs[irr].push_back(tkv_pair<double>(i*N[irr]+j, ints[k]));
                            if (i != j) kei_pairs[irr].push_back(tkv_pair<double>(j*N[irr]+i, ints[k]));
                        }

                        nproc = g.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                        for (int k = 0;k < nproc;k++)
                        {
                            int irr = irrep[idxs[k].i];
                            assert(irr == irrep[idxs[k].j]);

                            uint16_t i = idxs[k].i-start[irr];
                            uint16_t j = idxs[k].j-start[irr];

                                        nai_pairs[irr].push_back(tkv_pair<double>(i*N[irr]+j, ints[k]));
                            if (i != j) nai_pairs[irr].push_back(tkv_pair<double>(j*N[irr]+i, ints[k]));
                        }
                    }

                    block++;
                }
            }

            OVI *ovi = new OVI(arena, molecule.getGroup(), N);
            KEI *kei = new KEI(arena, molecule.getGroup(), N);
            NAI *nai = new NAI(arena, molecule.getGroup(), N);
            OneElectronHamiltonian *oeh = new OneElectronHamiltonian(arena, molecule.getGroup(), N);

            for (int i = 0;i < n;i++)
            {
                vector<int> irreps(2,i);
                (*ovi).writeRemoteData(irreps, ovi_pairs[i]);
                (*kei).writeRemoteData(irreps, kei_pairs[i]);
                (*nai).writeRemoteData(irreps, nai_pairs[i]);
                (*oeh).writeRemoteData(irreps, kei_pairs[i]);
                (*oeh).writeRemoteData(irreps, 1.0, 1.0, nai_pairs[i]);
            }

            put("S", ovi);
            put("T", kei);
            put("G", nai);
            put("H", oeh);

            return true;
        }
};

using Ishida1eIntegralsTask = OneElectronIntegralsTask<IshidaOVI, IshidaKEI, IshidaNAI>;

}
}

#endif
