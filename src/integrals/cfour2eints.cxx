#include "cfour2eints.hpp"

using namespace aquarius::input;
using namespace aquarius::symmetry;
using namespace aquarius::task;
using namespace aquarius::tensor;

namespace aquarius
{
namespace integrals
{

CFOURTwoElectronIntegralsTask::CFOURTwoElectronIntegralsTask(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("molecule", "molecule"));
    addProduct(Product("eri", "I", reqs));
}

bool CFOURTwoElectronIntegralsTask::run(TaskDAG& dag, const Arena& arena)
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
        SymmetryBlockedTensor<double> *tensor = NULL;
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

            eri->ints.insert(eri->ints.end(), ints.begin(), ints.begin()+numints);
            eri->idxs.insert(eri->idxs.end(), idxs.begin(), idxs.begin()+numints);
        }
    }

    for (size_t i = 0;i < eri->idxs.size();i++)
    {
        eri->idxs[i].i--;
        eri->idxs[i].j--;
        eri->idxs[i].k--;
        eri->idxs[i].l--;
    }

    put("I", eri);

    return true;
}

}
}

REGISTER_TASK(aquarius::integrals::CFOURTwoElectronIntegralsTask,"cfour2eints");
