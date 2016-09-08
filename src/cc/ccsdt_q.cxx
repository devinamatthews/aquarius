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
    this->addProduct(Product("ccsdt(q).G", "G", reqs));
    this->addProduct(Product("ccsdt(q).S", "S", reqs));
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

    const SpinorbitalTensor<U>&   fIA =   H.getIA();
    const SpinorbitalTensor<U>& VABIJ = H.getABIJ();
    const SpinorbitalTensor<U>& VIJAB = H.getIJAB();
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
    SpinorbitalTensor<U> Z4("Z4", arena, group, {vrt,occ}, {0,4}, {4,0});

    SpinorbitalTensor<U> T2T(SpinorbitalTensor<U>::transpose, T(2));
    SpinorbitalTensor<U> T3T(SpinorbitalTensor<U>::transpose, T(3));

    WABCEJK[  "abcejk"]  = 0.5*VABEF[  "abef"]*T(2)[    "fcjk"];

    WABMIJK[  "abmijk"]  =     VAMEI[  "amek"]*T(2)[    "ebij"];
    WABMIJK[  "abmijk"] -= 0.5*VMNIJ[  "nmjk"]*T(2)[    "abin"];

         T4["abcdijkl"]  =     VABEJ[  "abej"]*T(3)[  "ecdikl"];
         T4["abcdijkl"] -=     VAMIJ[  "amij"]*T(3)[  "bcdmkl"];

         Z4["ijklabcd"]  =                       T4["abcdijkl"];

         T4["abcdijkl"] +=   WABCEJK["abcejk"]*T(2)[    "edil"];
         T4["abcdijkl"] -=   WABMIJK["abmijk"]*T(2)[    "cdml"];

         Z4["ijklabcd"] +=       fIA[    "ia"]* T3T[  "jklbcd"];
         Z4["ijklabcd"] +=     VIJAB[  "ijab"]* T2T[    "klcd"];

    T4.weight({&D.getDA(), &D.getDI()}, {&D.getDa(), &D.getDi()});

    U E = (1.0/576.0)*scalar(Z4["ijklabcd"]*T4["abcdijkl"]);

    this->log(arena) << printos("energy: %18.15f\n", E) << endl;

    this->put("energy", new U(E));

    if (!this->isUsed("G") && !this->isUsed("S")) return true;

    printf("T4:    %.15f\n", (1.0/576.0)*scalar(T4*T4));
    printf("Z4:    %.15f\n", (1.0/576.0)*scalar(Z4*Z4));
    printf("\n");

    Z4.weight({&D.getDA(), &D.getDI()}, {&D.getDa(), &D.getDi()});

    SpinorbitalTensor<U> GAIJBCD(SpinorbitalTensor<U>::transpose, WABCEJK);
    SpinorbitalTensor<U> GIJKABL(SpinorbitalTensor<U>::transpose, WABMIJK);

    auto& G = this->put("G", new TwoElectronOperator<U>("G", arena, occ, vrt));
    auto& S = this->put("S", new DeexcitationOperator<U,3>("S", arena, occ, vrt));

    SpinorbitalTensor<U>&   DIJ =   G.getIJ();
    SpinorbitalTensor<U>&   DAB =   G.getAB();
    SpinorbitalTensor<U>&   DAI =   G.getAI();
    SpinorbitalTensor<U>& GABIJ = G.getABIJ();
    SpinorbitalTensor<U>& GABCI = G.getABCI();
    SpinorbitalTensor<U>& GABCD = G.getABCD();
    SpinorbitalTensor<U>& GIJKL = G.getIJKL();
    SpinorbitalTensor<U>& GAIJK = G.getAIJK();
    SpinorbitalTensor<U>& GAIBJ = G.getAIBJ();

    SpinorbitalTensor<U> S2T(SpinorbitalTensor<U>::transpose, S(2));
    SpinorbitalTensor<U> S3T(SpinorbitalTensor<U>::transpose, S(3));

        DIJ[      "ij"]  =  (1.0/144.0)*     Z4["imnoefgh"]*     T4["efghjmno"];
        DAB[      "ab"]  = -(1.0/144.0)*     Z4["mnopbefg"]*     T4["aefgmnop"];
        DAI[      "ai"]  =  (1.0/ 36.0)*    T3T[  "mnoefg"]*     T4["aefgimno"];
      GABIJ[    "abij"]  =  (1.0/  4.0)*    T2T[    "mnef"]*     T4["abefijmn"];

    GAIJBCD[  "aijbcd"]  = -(1.0/  2.0)*     Z4["mnijbecd"]*   T(2)[    "aemn"];
    GIJKABL[  "ijkabl"]  =  (1.0/  2.0)*     Z4["mijkeabf"]*   T(2)[    "efml"];

      GABCD[    "abcd"]  = -(1.0/  4.0)*GAIJBCD[  "amncde"]*   T(2)[    "bemn"];
      GAIBJ[    "aibj"]  = -(1.0/  2.0)*GIJKABL[  "mnibej"]*   T(2)[    "aemn"];
      GIJKL[    "ijkl"]  =  (1.0/  4.0)*GIJKABL[  "imjefl"]*   T(2)[    "efkm"];

        S2T[    "abij"]  =  (1.0/  4.0)*  VIJAB[    "mnef"]*     T4["abefijmn"];
        S3T[  "abcijk"]  =                  fIA[      "me"]*     T4["abceijkm"];

       S(2)[    "ijab"] -=  (1.0/  2.0)*  VABEF[    "efgb"]*GAIJBCD[  "gijeaf"];
       S(2)[    "ijab"] +=                VAMEI[    "eifm"]*GAIJBCD[  "fmjeab"];
       S(2)[    "ijab"] -=                VAMEI[    "eman"]*GIJKABL[  "ijnebm"];
       S(2)[    "ijab"] +=  (1.0/  2.0)*  VMNIJ[    "jmno"]*GIJKABL[  "inoabm"];

         T4["abcdijkl"] +=                                       Z4["ijklabcd"];

      GABCI[    "abci"]  = -(1.0/ 12.0)*    T3T[  "mnoefc"]*     T4["efabmnoi"];
      GAIJK[    "aijk"]  =  (1.0/ 12.0)*    T3T[  "mniefg"]*     T4["efagmnjk"];

        S3T[  "abcijk"] -=  (1.0/  2.0)*  VMNEJ[    "mnei"]*     T4["eabcmnjk"];
        S3T[  "abcijk"] +=  (1.0/  2.0)*  VAMEF[    "amef"]*     T4["efbcimjk"];

    printf("GABIJ: %.15f\n", GABIJ({1,0},{0,1}).norm(2));
    printf("GABCI: %.15f\n", GABCI({1,0},{0,1}).norm(2));
    printf("GAIJK: %.15f\n", GAIJK({1,0},{0,1}).norm(2));
    printf("GABCD: %.15f\n", GABCD({1,0},{1,0}).norm(2));
    printf("GAIBJ: %.15f\n", GAIBJ({0,1},{0,1}).norm(2));
    printf("GAIJB: %.15f\n", GAIBJ({1,0},{0,1}).norm(2));
    printf("GIJKL: %.15f\n", GIJKL({0,1},{0,1}).norm(2));
    printf("\n");
    printf("S2:    %.15f\n", S(2)({0,1},{1,0}).norm(2));
    printf("S3:    %.15f\n", (1.0/36.0)*scalar(S(3)*S(3)));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSDT_Q);
REGISTER_TASK(aquarius::cc::CCSDT_Q<double>,"ccsdt(q)");
