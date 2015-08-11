#include "ccsdt_q.hpp"

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
CCSDT_Q<U>::CCSDT_Q(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    reqs.push_back(Requirement("ccsdt.T", "T"));
    this->addProduct(Product("double", "energy", reqs));
}

template <typename U>
bool CCSDT_Q<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U>>("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;
    const PointGroup& group = occ.group;

    Denominator<U> D(H);
    const ExcitationOperator<U,3>& T = this->template get<ExcitationOperator<U,3>>("T");

    const SpinorbitalTensor<U>& VABIJ = H.getABIJ();
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
    SpinorbitalTensor<U> Z4("L4", arena, group, {vrt,occ}, {4,0}, {0,4});

    WABCEJK[  "abcejk"]  = 0.5*VABEF[  "abef"]*T(2)[  "fcjk"];

    WABMIJK[  "abmijk"]  =     VAMEI[  "amek"]*T(2)[  "ebij"];
    WABMIJK[  "abmijk"] -= 0.5*VMNIJ[  "nmjk"]*T(2)[  "abin"];

         T4["abcdijkl"]  =     VABEJ[  "abej"]*T(3)["ecdikl"];
         T4["abcdijkl"] -=     VAMIJ[  "amij"]*T(3)["bcdmkl"];

         Z4 = T4;

         T4["abcdijkl"] +=   WABCEJK["abcejk"]*T(2)[  "edil"];
         T4["abcdijkl"] -=   WABMIJK["abmijk"]*T(2)[  "cdml"];

         Z4["abcdijkl"] +=     VABIJ[  "abij"]*T(2)[  "cdkl"];

    T4.weight({&D.getDA(), &D.getDI()}, {&D.getDa(), &D.getDi()});

    U E = (1.0/576.0)*scalar(Z4*T4);

    this->log(arena) << strprintf("energy: %18.15f\n", E) << endl;

    this->put("energy", new U(E));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSDT_Q);
REGISTER_TASK(aquarius::cc::CCSDT_Q<double>,"ccsdt(q)");
