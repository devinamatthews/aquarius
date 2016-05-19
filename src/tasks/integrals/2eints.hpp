#ifndef _AQUARIUS_TASK_INTEGRALS_2EINTS_HPP_
#define _AQUARIUS_TASK_INTEGRALS_2EINTS_HPP_

#include "integrals/2eints.hpp"

#include "../../frameworks/molecule/molecule.hpp"
#include "../../frameworks/symmetry/symmetry.hpp"
#include "../../frameworks/task/config.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/util/global.hpp"

#define TMP_BUFSIZE 65536
#define INTEGRAL_CUTOFF 1e-14

namespace aquarius
{

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

            integrals::ERI* eri = new integrals::ERI(arena, molecule.getGroup());

            integrals::Context ctx(integrals::Context::ISCF);

            vector<double> tmpval(TMP_BUFSIZE);
            vector<idx4_t> tmpidx(TMP_BUFSIZE);

            const vector<int>& N = molecule.getNumOrbitals();
            int nirrep = molecule.getGroup().getNumIrreps();

            vector<vector<int>> idx = integrals::Shell::setupIndices(integrals::Context(), molecule);
            vector<integrals::Shell> shells = molecule.getShells();

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

}

#undef TMP_BUFSIZE
#undef INTEGRAL_CUTOFF

#endif
