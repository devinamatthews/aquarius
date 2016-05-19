#ifndef _AQUARIUS_TASK_INTEGRALS_1EINTS_HPP_
#define _AQUARIUS_TASK_INTEGRALS_1EINTS_HPP_

#include "frameworks/integrals.hpp"

#include "../../frameworks/molecule/molecule.hpp"
#include "../../frameworks/symmetry/symmetry.hpp"
#include "../../frameworks/task/config.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/tensor/tensor.hpp"
#include "../../frameworks/util/global.hpp"

namespace aquarius
{

template <typename OVIType, typename KEIType, typename NAIType>
class OneElectronIntegralsTask : public task::Task
{
    public:
        OneElectronIntegralsTask(const string& name, task::Config& config)
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

            integrals::Context ctx(integrals::Context::ISCF);

            const vector<int>& N = molecule.getNumOrbitals();
            int n = molecule.getGroup().getNumIrreps();

            vector<int> irrep;
            for (int i = 0;i < n;i++) irrep += vector<int>(N[i],i);

            vector<uint16_t> start(n,0);
            for (int i = 1;i < n;i++) start[i] = start[i-1]+N[i-1];

            vector<vector<int>> idx = integrals::Shell::setupIndices(ctx, molecule);
            vector<integrals::Shell> shells = molecule.getShells();
            vector<KeyValueVector> ovi_pairs(n), nai_pairs(n), kei_pairs(n);
            vector<integrals::Center> centers;

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

            auto init = TensorInitializer<>("S", Field::DOUBLE) <<
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
                H.addDataByIrrep(irreps, 1.0, nai_pairs[i], 1.0);
            }

            return true;
        }
};

}

#endif
