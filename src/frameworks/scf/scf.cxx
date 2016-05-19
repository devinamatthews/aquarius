#include "scf.hpp"

using namespace aquarius::task;
using namespace aquarius::input;

namespace aquarius
{
namespace scf
{

template <typename T>
SCF<T>::SCF(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs += Requirement("molecule", "molecule");
    reqs += Requirement("ovi", "S");
    reqs += Requirement("1ehamiltonian", "H");
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("double", "S2", reqs));
    this->addProduct(Product("double", "multiplicity", reqs));
    this->addProduct(Product("occspace", "occ", reqs));
    this->addProduct(Product("vrtspace", "vrt", reqs));
    this->addProduct(Product("Ea", "Ea", reqs));
    this->addProduct(Product("Eb", "Eb", reqs));
    this->addProduct(Product("Fa", "Fa", reqs));
    this->addProduct(Product("Fb", "Fb", reqs));
    this->addProduct(Product("Da", "Da", reqs));
    this->addProduct(Product("Db", "Db", reqs));
}

template <typename T>
bool SCF<T>::run(TaskDAG& dag, const Arena& arena)
{
    return true;
}

INSTANTIATE_SPECIALIZATIONS(SCF);

}
}
