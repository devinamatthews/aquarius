#include "ccsd_t_l.hpp"

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
CCSD_T_L<U>::CCSD_T_L(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    reqs.push_back(Requirement("ccsd.T", "T"));
    reqs.push_back(Requirement("ccsd.L", "L"));
    this->addProduct(Product("double", "energy", reqs));
}

template <typename U>
bool CCSD_T_L<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U>>("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;
    const PointGroup& group = occ.group;

    Denominator<U> D(H);
    const ExcitationOperator  <U,2>& T = this->template get<ExcitationOperator  <U,2>>("T");
    const DeexcitationOperator<U,2>& L = this->template get<DeexcitationOperator<U,2>>("L");

    const SpinorbitalTensor<U>& VIJAB = H.getIJAB();
    const SpinorbitalTensor<U>& VABCI = H.getABCI();
    const SpinorbitalTensor<U>& VAIBC = H.getAIBC();
    const SpinorbitalTensor<U>& VAIJK = H.getAIJK();
    const SpinorbitalTensor<U>& VIJAK = H.getIJAK();

    SpinorbitalTensor<U> T3("T3", arena, group, {vrt,occ}, {3,0}, {0,3});
    SpinorbitalTensor<U> Z3("Z3", arena, group, {vrt,occ}, {0,3}, {3,0});

    T3["abcijk"]  = VABCI["bcek"]*T(2)["aeij"];
    T3["abcijk"] -= VAIJK["amij"]*T(2)["bcmk"];
    T3.weight({&D.getDA(), &D.getDI()}, {&D.getDa(), &D.getDi()});

    Z3["ijkabc"]  = VAIBC["ekbc"]*L(2)["ijae"];
    Z3["ijkabc"] -= VIJAK["ijam"]*L(2)["mkbc"];
    Z3["ijkabc"] += VIJAB["ijab"]*L(1)[  "kc"];

    U E_T = (1.0/36.0)*scalar(T3["efgmno"]*Z3["mnoefg"]);
    this->log(arena) << printos("energy: %18.15f", E_T) << endl;

    this->put("energy", new U(E_T));

    return true;
}

}
}

static const char* spec = "";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSD_T_L);
REGISTER_TASK(aquarius::cc::CCSD_T_L<double>,"ccsd(t)_l",spec);
