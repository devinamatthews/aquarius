#ifndef _AQUARIUS_INTEGRALS_1EINTS_HPP_
#define _AQUARIUS_INTEGRALS_1EINTS_HPP_

#include "util/global.hpp"

#include "symmetry/symmetry.hpp"
#include "tensor/tensor.hpp"
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
            using namespace aquarius::tensor;

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
            vector<KeyValueVector> ovi_pairs(n, KeyValueVector(Field::DOUBLE));
            vector<KeyValueVector> nai_pairs(n, KeyValueVector(Field::DOUBLE));
            vector<KeyValueVector> kei_pairs(n, KeyValueVector(Field::DOUBLE));
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

                                        ovi_pairs[irr].push_back(i*N[irr]+j, ints[k]);
                            if (i != j) ovi_pairs[irr].push_back(j*N[irr]+i, ints[k]);
                        }

                        nproc = t.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                        for (int k = 0;k < nproc;k++)
                        {
                            int irr = irrep[idxs[k].i];
                            assert(irr == irrep[idxs[k].j]);

                            uint16_t i = idxs[k].i-start[irr];
                            uint16_t j = idxs[k].j-start[irr];

                                        kei_pairs[irr].push_back(i*N[irr]+j, ints[k]);
                            if (i != j) kei_pairs[irr].push_back(j*N[irr]+i, ints[k]);
                        }

                        nproc = g.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                        for (int k = 0;k < nproc;k++)
                        {
                            int irr = irrep[idxs[k].i];
                            assert(irr == irrep[idxs[k].j]);

                            uint16_t i = idxs[k].i-start[irr];
                            uint16_t j = idxs[k].j-start[irr];

                                        nai_pairs[irr].push_back(i*N[irr]+j, ints[k]);
                            if (i != j) nai_pairs[irr].push_back(j*N[irr]+i, ints[k]);
                        }
                    }

                    block++;
                }
            }

            auto init = TensorInitializer<>("S", tensor::Field::DOUBLE) <<
                        TensorInitializer<DISTRIBUTED>(arena) <<
                        TensorInitializer<PGSYMMETRIC|BOUNDED>(molecule.getGroup(), {N,N});
            Tensor<BOUNDED|PGSYMMETRIC> S = put<Tensor<>>("S", Tensor<DISTRIBUTED|PGSYMMETRIC|BOUNDED>::construct(init));

            Tensor<BOUNDED|PGSYMMETRIC> T = put<Tensor<>>("T", S.construct("T"));
            Tensor<BOUNDED|PGSYMMETRIC> G = put<Tensor<>>("G", S.construct("G"));
            Tensor<BOUNDED|PGSYMMETRIC> H = put<Tensor<>>("H", S.construct("H"));

            for (int i = 0;i < n;i++)
            {
                vector<int> irreps = {i,i};
                S.setDataByIrrep(irreps, ovi_pairs[i]);
                T.setDataByIrrep(irreps, kei_pairs[i]);
                G.setDataByIrrep(irreps, nai_pairs[i]);
                H.setDataByIrrep(irreps, kei_pairs[i]);
                H.addDataByIrrep(irreps, 1, nai_pairs[i], 1);
            }

            return true;
        }
};

using Ishida1eIntegralsTask = OneElectronIntegralsTask<IshidaOVI, IshidaKEI, IshidaNAI>;

}
}

#endif
