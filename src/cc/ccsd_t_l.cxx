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
    this->addProduct(Product("ccsd(t)_l.U", "U", reqs));
    this->addProduct(Product("ccsd(t)_l.S", "S", reqs));
    this->addProduct(Product("ccsd(t)_l.D", "D", reqs));
}

template <typename Type>
bool CCSD_T_L<Type>::run(task::TaskDAG& dag, const Arena& arena)
{
    const auto& H = this->template get<TwoElectronOperator<Type>>("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;
    const PointGroup& group = occ.group;

    Denominator<Type> Denom(H);
    const auto& T = this->template get<ExcitationOperator  <Type,2>>("T");
    const auto& L = this->template get<DeexcitationOperator<Type,2>>("L");

    const auto&   fIA =   H.getIA();
    const auto& VIJAB = H.getIJAB();
    const auto& VABCI = H.getABCI();
    const auto& VAIBC = H.getAIBC();
    const auto& VAIJK = H.getAIJK();
    const auto& VIJAK = H.getIJAK();

    SpinorbitalTensor<Type> T3("T3", arena, group, {vrt,occ}, {3,0}, {0,3});
    SpinorbitalTensor<Type> Z3("Z3", arena, group, {vrt,occ}, {0,3}, {3,0});

    T3["abcijk"]  = VABCI["bcek"]*T(2)["aeij"];
    T3["abcijk"] -= VAIJK["amij"]*T(2)["bcmk"];
    T3.weight({&Denom.getDA(), &Denom.getDI()}, {&Denom.getDa(), &Denom.getDi()});

    Z3["ijkabc"]  = VAIBC["ekbc"]*L(2)["ijae"];
    Z3["ijkabc"] -= VIJAK["ijam"]*L(2)["mkbc"];
    Z3["ijkabc"] += VIJAB["ijab"]*L(1)[  "kc"];
    Z3["ijkabc"] +=   fIA[  "ia"]*L(2)["jkbc"];

    Type E_T = (1.0/36.0)*scalar(T3["efgmno"]*Z3["mnoefg"]);
    this->log(arena) << printos("energy: %18.15f", E_T) << endl;

    this->put("energy", new Type(E_T));

    if (this->isUsed("U") || this->isUsed("S") || this->isUsed("D"))
    {
        auto& U = this->put("U", new ExcitationOperator<Type,2>("U", arena, occ, vrt));
        auto& S = this->put("S", new DeexcitationOperator<Type,2>("S", arena, occ, vrt));
        auto& D = this->put("D", new TwoElectronOperator<Type>("D", arena, occ, vrt));

        Z3.weight({&Denom.getDA(), &Denom.getDI()}, {&Denom.getDa(), &Denom.getDi()});

          D.getIJ()[  "ij"]  =  (1.0/12.0)*T3["efgjmn"]*   Z3["imnefg"];
          D.getAB()[  "ab"]  = -(1.0/12.0)*T3["aefmno"]*   Z3["mnobef"];
          D.getAI()[  "ai"]  =  (1.0/ 4.0)*T3["aefimn"]* L(2)[  "mnef"];

        D.getABIJ()["abij"]  =             T3["abeijm"]* L(1)[    "me"];
        D.getAIJK()["aijk"]  =  (1.0/ 2.0)*T3["aefjkm"]* L(2)[  "imef"];
        D.getABCI()["abci"]  = -(1.0/ 2.0)*T3["abemin"]* L(2)[  "mnce"];
        D.getIJAK()["ijak"]  =  (1.0/ 2.0)*Z3["ijmaef"]* T(2)[  "efkm"];
        D.getAIBC()["aibc"]  = -(1.0/ 2.0)*Z3["minbce"]* T(2)[  "aemn"];

               U(1)[  "ai"]  =  (1.0/ 4.0)*T3["aefimn"]*VIJAB[  "mnef"];
               U(2)["abij"]  =             T3["abeijm"]*  fIA[    "me"];
               U(2)["abij"] +=  (1.0/ 2.0)*T3["aefijm"]*VAIBC[  "bmef"];
               U(2)["abij"] -=  (1.0/ 2.0)*T3["abeimn"]*VIJAK[  "nmej"];

               S(2)["ijab"]  =  (1.0/ 2.0)*Z3["ijmaef"]*VABCI[  "efbm"];
               S(2)["ijab"] -=  (1.0/ 2.0)*Z3["imnabe"]*VAIJK[  "ejnm"];

        print("U1: %.15f\n", U(1)({0,0},{0,0}).norm(2));
        print("U2: %.15f\n", U(2)({1,0},{0,1}).norm(2));
        print("S2: %.15f\n", S(2)({0,1},{1,0}).norm(2));
    }

    return true;
}

}
}

static const char* spec = "";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSD_T_L);
REGISTER_TASK(aquarius::cc::CCSD_T_L<double>,"ccsd(t)_l",spec);
