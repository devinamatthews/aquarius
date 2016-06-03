#include "../../frameworks/cc/ccsdt_q_l.hpp"

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
CCSDT_Q_L<U>::CCSDT_Q_L(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    reqs.push_back(Requirement("ccsdt.T", "T"));
    reqs.push_back(Requirement("ccsdt.L", "L"));
    this->addProduct(Product("double", "energy", reqs));
}

template <typename U>
bool CCSDT_Q_L<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U>>("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;
    const PointGroup& group = occ.group;

    Denominator<U> D(H);
    const ExcitationOperator  <U,3>& T = this->template get<ExcitationOperator  <U,3>>("T");
    const DeexcitationOperator<U,3>& L = this->template get<DeexcitationOperator<U,3>>("L");

    const SpinorbitalTensor<U>& VMNEF = H.getIJAB();
    const SpinorbitalTensor<U>& VAMEF = H.getAIBC();
    const SpinorbitalTensor<U>& VABEJ = H.getABCI();
    const SpinorbitalTensor<U>& VABEF = H.getABCD();
    const SpinorbitalTensor<U>& VMNIJ = H.getIJKL();
    const SpinorbitalTensor<U>& VMNEJ = H.getIJAK();
    const SpinorbitalTensor<U>& VAMIJ = H.getAIJK();
    const SpinorbitalTensor<U>& VAMEI = H.getAIBJ();

    SpinorbitalTensor<U> WABCEJK("W(abc,ejk)", arena, group, {vrt,occ}, {3,0}, {1,2});
    SpinorbitalTensor<U> WABMIJK("W(abm,ijk)", arena, group, {vrt,occ}, {2,1}, {0,3});

    SpinorbitalTensor<U> T4("T4", arena, group, {vrt,occ}, {4,0}, {0,4});
    SpinorbitalTensor<U> L4("L4", arena, group, {vrt,occ}, {0,4}, {4,0});

    WABCEJK[  "abcejk"]  = 0.5*VABEF[  "abef"]*T(2)[  "fcjk"];

    WABMIJK[  "abmijk"]  =     VAMEI[  "amek"]*T(2)[  "ebij"];
    WABMIJK[  "abmijk"] -= 0.5*VMNIJ[  "nmjk"]*T(2)[  "abin"];

         T4["abcdijkl"]  =   WABCEJK["abcejk"]*T(2)[  "edil"];
         T4["abcdijkl"] -=   WABMIJK["abmijk"]*T(2)[  "cdml"];
         T4["abcdijkl"] +=     VABEJ[  "abej"]*T(3)["ecdikl"];
         T4["abcdijkl"] -=     VAMIJ[  "amij"]*T(3)["bcdmkl"];

         L4["ijklabcd"]  =     VMNEF[  "ijab"]*L(2)[  "klcd"];
         L4["ijklabcd"] +=     VAMEF[  "elcd"]*L(3)["ijkabe"];
         L4["ijklabcd"] -=     VMNEJ[  "ijam"]*L(3)["mklbcd"];

    T4.weight({&D.getDA(), &D.getDI()}, {&D.getDa(), &D.getDi()});

    U E = (1.0/576.0)*scalar(L4["mnopefgh"]*T4["efghmnop"]);

    this->log(arena) << printos("energy: %18.15f\n", E) << endl;

    this->put("energy", new U(E));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSDT_Q_L);
REGISTER_TASK(aquarius::cc::CCSDT_Q_L<double>,"ccsdt(q)_l");
