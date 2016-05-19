#include "time/time.hpp"
#include "../../frameworks/operator/fakemoints.hpp"

#include "../../frameworks/time/time.hpp"

using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::task;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace op
{

template <typename T>
FakeMOIntegrals<T>::FakeMOIntegrals(const string& name, Config& config)
: Task(name, config),
  noa(config.get<int>("noa")),
  nob(config.get<int>("nob")),
  nva(config.get<int>("nva")),
  nvb(config.get<int>("nvb"))
{
    if (nob < 0) nob = noa;
    if (nvb < 0) nvb = nva;
    addProduct("moints", "H");
}

template <typename T>
bool FakeMOIntegrals<T>::run(TaskDAG& dag, const Arena& arena)
{
    Space occ(PointGroup::C1(), {noa}, {nob});
    Space vrt(PointGroup::C1(), {nva}, {nvb});

    this->put("H", new TwoElectronOperator<T>("V", OneElectronOperator<T>("f", arena, occ, vrt)));

    return true;
}

}
}

static const char* spec = R"!(

noa int,
nob? int -1,
nva int,
nvb? int -1

)!";

INSTANTIATE_SPECIALIZATIONS(aquarius::op::FakeMOIntegrals);
REGISTER_TASK(aquarius::op::FakeMOIntegrals<double>,"fakemoints",spec);
