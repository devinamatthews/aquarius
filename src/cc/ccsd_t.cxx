#include "ccsd_t.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace cc
{

template <typename U>
CCSD_T<U>::CCSD_T(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    reqs.push_back(Requirement("ccsd.T", "T"));
    this->addProduct(Product("double", "energy", reqs));
}

template <typename U>
bool CCSD_T<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U>>("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;
    const PointGroup& group = occ.group;

    Denominator<U> D(H);
    const ExcitationOperator<U,2>& T = this->template get<ExcitationOperator<U,2>>("T");

    const SpinorbitalTensor<U>& VABIJ = H.getABIJ();
    const SpinorbitalTensor<U>& VABCI = H.getABCI();
    const SpinorbitalTensor<U>& VAIJK = H.getAIJK();

    SpinorbitalTensor<U> T3("T3", arena, group, {vrt,occ}, {3,0}, {0,3});
    SpinorbitalTensor<U> Z3("Z3", arena, group, {vrt,occ}, {3,0}, {0,3});

    Z3["abcijk"]  = VABCI["bcek"]*T(2)["aeij"];
    Z3["abcijk"] -= VAIJK["amij"]*T(2)["bcmk"];

    T3 = Z3;
    T3.weight({&D.getDA(), &D.getDI()}, {&D.getDa(), &D.getDi()});

    Z3["abcijk"] += VABIJ["abij"]*T(1)[  "ck"];

    U E_T = (1.0/36.0)*scalar(T3["efgmno"]*Z3["efgmno"]);
    this->log(arena) << strprintf("energy: %18.15f", E_T) << endl;

    this->put("energy", new U(E_T));

    return true;
}

}
}

static const char* spec = "";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSD_T);
REGISTER_TASK(aquarius::cc::CCSD_T<double>,"ccsd(t)",spec);
