#include "../../frameworks/input/config.hpp"
#include "../../frameworks/input/molecule.hpp"
#include "../../frameworks/molecule/molecule.hpp"
#include "../../frameworks/util/global.hpp"
using namespace aquarius::input;

#include "../../frameworks/symmetry/symmetry.hpp"
using namespace aquarius::symmetry;

#include "../../frameworks/tensor/tensor.hpp"
using namespace aquarius::tensor;

#include "../../frameworks/task/task.hpp"
using namespace aquarius::task;

#include "integrals/2eints.hpp"
using namespace aquarius::integrals;

namespace aquarius
{

class CFOURTwoElectronIntegralsTask : public task::Task
{
    public:
        CFOURTwoElectronIntegralsTask(const string& name, input::Config& config)
        : Task(name, config)
        {
            vector<Requirement> reqs;
            reqs.push_back(Requirement("molecule", "molecule"));
            addProduct(Product("eri", "I", reqs));
        }

        bool run(task::TaskDAG& dag, const Arena& arena)
        {
            const Molecule& molecule = get<Molecule>("molecule");

            ERI* eri = new ERI(arena, molecule.getGroup());

            const vector<int>& N = molecule.getNumOrbitals();
            int n = molecule.getGroup().getNumIrreps();

            ifstream ifs("IIII");

            int intsize = 8;
            int32_t junk;

            ifs.read((char*)&junk, 4);
            assert(junk == 224);

            ifs.read((char*)&junk, 4);
            if (junk != 0) intsize = 4;

            int batchsize = 600;
            vector<double> ints(batchsize);
            vector<idx4_t> idxs(batchsize);

            ifs.seekg(0);
            while (ifs)
            {
                int64_t recsize = 0;
                ifs.read((char*)&recsize, intsize);
                if (recsize != 8)
                {
                    ifs.seekg(recsize+intsize, ifstream::cur);
                    continue;
                }

                char label[9] = {};
                ifs.read(label, 8);
                ifs.seekg(intsize, ifstream::cur);
                if (strcmp(label, "TWOELSUP") != 0) continue;

                while (true)
                {
                    int64_t numints;
                    ifs.read((char*)&recsize, intsize);
                    assert(recsize == 600*16+8);
                    ifs.read((char*)ints.data(), 600*8);
                    ifs.read((char*)idxs.data(), 600*8);
                    ifs.read((char*)&numints, 8);
                    ifs.seekg(intsize, ifstream::cur);

                    if (numints < 0) break;

                    for (size_t i = 0;i < numints;i++)
                    {
                        idxs[i].i--;
                        idxs[i].j--;
                        idxs[i].k--;
                        idxs[i].l--;
                    }

                    eri->ints.insert(eri->ints.end(), ints.begin(), ints.begin()+numints);
                    eri->idxs.insert(eri->idxs.end(), idxs.begin(), idxs.begin()+numints);
                }
            }

            put("I", eri);

            return true;
        }
};

}

REGISTER_TASK(aquarius::integrals::CFOURTwoElectronIntegralsTask,"cfour2eints");
