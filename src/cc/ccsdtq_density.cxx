#include "ccsdtq_density.hpp"

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
CCSDTQDensity<U>::CCSDTQDensity(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    reqs.push_back(Requirement("ccsdtq.T", "T"));
    reqs.push_back(Requirement("ccsdtq.L", "L"));
    this->addProduct(Product("ccsdtq.D", "D", reqs));
}

template <typename U>
bool CCSDTQDensity<U>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& H = this->template get<TwoElectronOperator <U  >>("H");
    const auto& T = this->template get<ExcitationOperator  <U,4>>("T");
    const auto& L = this->template get<DeexcitationOperator<U,4>>("L");

    const PointGroup& group = T(1).getGroup();
    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    auto& D = this->put("D", new TwoElectronOperator<U>("D", arena, occ, vrt));

    auto& DIA = D.getIA();
    auto& DAI = D.getAI();
    auto& DAB = D.getAB();
    auto& DIJ = D.getIJ();

    auto& GIJAB = D.getIJAB();
    auto& GIJAK = D.getIJAK();
    auto& GAIBC = D.getAIBC();
    auto& GABCD = D.getABCD();
    auto& GIJKL = D.getIJKL();
    auto& GAIBJ = D.getAIBJ();
    auto& GAIJK = D.getAIJK();
    auto& GABCI = D.getABCI();
    auto& GABIJ = D.getABIJ();

    SpinorbitalTensor<U> IIJAK(GIJAK);
    SpinorbitalTensor<U> IAIBC(GAIBC);

    IIJAK["ijak"]  =  L(2)["ijae"]*T(1)["ek"];
    IAIBC["aibc"]  = -L(2)["mibc"]*T(1)["am"];

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    SpinorbitalTensor<U> GAIJBCD("G(aij,bcd)", arena, group, {vrt,occ}, {1,2}, {3,0});
    SpinorbitalTensor<U> GIJKABL("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1});
    SpinorbitalTensor<U> GIJKALM("G(ijk,alm)", arena, group, {vrt,occ}, {0,3}, {1,2});
    SpinorbitalTensor<U> GAIJBCK("G(aij,bck)", arena, group, {vrt,occ}, {1,2}, {2,1});

    GAIJBCD["aijbcd"]  = -(1.0/  2.0)*   L(4)["mnijebcd"]*T(2)[    "eamn"];

    GIJKABL["ijkabl"]  =  (1.0/  2.0)*   L(4)["ijkmabef"]*T(2)[    "eflm"];

    GIJKALM["ijkalm"]  =  (1.0/  6.0)*   L(4)["ijknaefg"]*T(3)[  "efglmn"];

    GAIJBCK["aijbck"]  =  (1.0/  8.0)*   L(4)["mnijebcf"]*T(3)[  "efamkn"];

      GIJAK[  "ijak"]  =  (1.0/ 12.0)*   L(4)["ijmnaefg"]*T(3)[  "efgkmn"];

      GAIBC[  "aibc"]  = -(1.0/ 12.0)*   L(4)["mnoiefbc"]*T(3)[  "efamno"];

      GIJKL[  "ijkl"]  =  (1.0/ 48.0)*   L(4)["mnijefgh"]*T(4)["efghmnkl"];
      GIJKL[  "ijkl"] +=  (1.0/  4.0)*GIJKABL[  "mijefl"]*T(2)[    "efmk"];

      GABCD[  "abcd"]  =  (1.0/ 48.0)*   L(4)["mnopefcd"]*T(4)["efabmnop"];
      GABCD[  "abcd"] -=  (1.0/  4.0)*GAIJBCD[  "amncde"]*T(2)[    "bemn"];

      GAIBJ[  "aibj"]  = -(1.0/ 36.0)*   L(4)["imnobefg"]*T(4)["aefgjmno"];
      GAIBJ[  "aibj"] -=  (1.0/  2.0)*GIJKABL[  "mniebj"]*T(2)[    "eamn"];

      GAIJK[  "aijk"]  =  (1.0/ 12.0)*   L(3)[  "imnefg"]*T(4)["aefgjkmn"];
      GAIJK[  "aijk"] +=  (1.0/  6.0)*GAIJBCD[  "aimefg"]*T(3)[  "efgjkm"];
      GAIJK[  "aijk"] +=  (1.0/  4.0)*GIJKABL[  "mniefk"]*T(3)[  "efamnj"];

      GABCI[  "abci"]  = -(1.0/ 12.0)*   L(3)[  "mnoefc"]*T(4)["efabmnoi"];
      GABCI[  "abci"] +=  (1.0/  6.0)*GIJKABL[  "mnoeci"]*T(3)[  "eabmno"];
      GABCI[  "abci"] +=  (1.0/  4.0)*GAIJBCD[  "amncef"]*T(3)[  "efbmni"];

      GABIJ[  "abij"]  =  (1.0/  4.0)*   L(2)[    "mnef"]*T(4)["abefijmn"];
      GABIJ[  "abij"] +=  (1.0/ 12.0)*GAIJBCD[  "amnefg"]*T(4)["efgbimnj"];
      GABIJ[  "abij"] -=  (1.0/ 12.0)*GIJKABL[  "mnoefj"]*T(4)["aefbimno"];
      GABIJ[  "abij"] +=  (1.0/  6.0)*GIJKALM[  "mnoeij"]*T(3)[  "eabmno"];
      GABIJ[  "abij"] +=  (1.0/  4.0)*GAIJBCK[  "amnefj"]*T(3)[  "efbimn"];

        DIJ[    "ij"]  =  (1.0/144.0)*   L(4)["imnoefgh"]*T(4)["efghjmno"];

        DAB[    "ab"]  = -(1.0/144.0)*   L(4)["mnopbefg"]*T(4)["aefgmnop"];

        DAI[    "ai"]  =  (1.0/ 36.0)*   L(3)[  "mnoefg"]*T(4)["aefgimno"];
        DAI[    "ai"] +=  (1.0/  2.0)*  GAIBC[    "amef"]*T(2)[    "efim"];

    /*********************************************************************/

      GIJAK[  "ijak"] +=  (1.0/  2.0)*   L(3)[  "ijmaef"]*T(2)[    "efkm"];

      GAIBC[  "aibc"] -=  (1.0/  2.0)*   L(3)[  "mniebc"]*T(2)[    "eamn"];

      GIJKL[  "ijkl"] +=  (1.0/  6.0)*   L(3)[  "mijefg"]*T(3)[  "efgmkl"];
      GIJKL[  "ijkl"] +=                GIJAK[    "ijel"]*T(1)[      "ek"];

      GABCD[  "abcd"] +=  (1.0/  6.0)*   L(3)[  "mnoecd"]*T(3)[  "eabmno"];
      GABCD[  "abcd"] -=                GAIBC[    "amcd"]*T(1)[      "bm"];

      GAIBJ[  "aibj"] -=  (1.0/  4.0)*   L(3)[  "imnbef"]*T(3)[  "aefjmn"];

      GAIJK[  "aijk"] +=  (1.0/  2.0)*   L(2)[    "imef"]*T(3)[  "aefjkm"];
      GAIJK[  "aijk"] +=                GIJAK[    "miek"]*T(2)[    "aejm"];
      GAIJK[  "aijk"] +=                GAIBJ[    "aiek"]*T(1)[      "ej"];

      GABCI[  "abci"] -=  (1.0/  2.0)*   L(2)[    "mnec"]*T(3)[  "eabmni"];
      GABCI[  "abci"] +=                GAIBC[    "amce"]*T(2)[    "beim"];

      GABIJ[  "abij"] +=                 L(1)[      "me"]*T(3)[  "abeijm"];
      GABIJ[  "abij"] +=  (1.0/  2.0)*  GAIBC[    "amef"]*T(3)[  "ebfijm"];
      GABIJ[  "abij"] -=  (1.0/  2.0)*  GIJAK[    "mnej"]*T(3)[  "eabmin"];
      GABIJ[  "abij"] +=                GABCI[    "abej"]*T(1)[      "ei"];
      GABIJ[  "abij"] -=                GAIJK[    "amij"]*T(1)[      "bm"];
      GABIJ[  "abij"] +=  (1.0/  2.0)*  GABCD[    "abef"]* Tau[    "efij"];

      GAIJK[  "aijk"] -=                GIJKL[    "mijk"]*T(1)[      "am"];

      GABCI[  "abci"] -=                GAIBJ[    "amci"]*T(1)[      "bm"];
      GABCI[  "abci"] +=                GABCD[    "abce"]*T(1)[      "ei"];

        DIJ[    "ij"] +=  (1.0/ 12.0)*   L(3)[  "imnefg"]*T(3)[  "efgjmn"];

        DAB[    "ab"] -=  (1.0/ 12.0)*   L(3)[  "mnobef"]*T(3)[  "aefmno"];

        DAI[    "ai"] +=  (1.0/  4.0)*   L(2)[    "mnef"]*T(3)[  "aefimn"];
        DAI[    "ai"] -=  (1.0/  2.0)*  GIJAK[    "mnei"]*T(2)[    "eamn"];

    /*********************************************************************/

        DIA[    "ia"]  =                 L(1)[      "ia"];

        DAB[    "ab"] -=  (1.0/  2.0)*   L(2)[    "mnbe"]*T(2)[    "aemn"];

        DIJ[    "ij"] +=                 L(1)[      "ie"]*T(1)[      "ej"];
        DIJ[    "ij"] +=  (1.0/  2.0)*   L(2)[    "imef"]*T(2)[    "efjm"];

        DAI[    "ai"] +=                 L(1)[      "me"]*T(2)[    "aeim"];
        DAI[    "ai"] -=                  DIJ[      "mi"]*T(1)[      "am"];
        DAI[    "ai"] +=                  DAB[      "ae"]*T(1)[      "ei"];

        DAB[    "ab"] -=                 L(1)[      "mb"]*T(1)[      "am"];

      GIJAB[  "ijab"]  =                 L(2)[    "ijab"];

      GIJKL[  "ijkl"] +=  (1.0/  2.0)*   L(2)[    "ijef"]* Tau[    "efkl"];

      GABCD[  "abcd"] +=  (1.0/  2.0)*   L(2)[    "mncd"]* Tau[    "abmn"];

      GAIBJ[  "aibj"] -=  (1.0/  2.0)*   L(2)[    "imbe"]*T(2)[    "aejm"];
      GAIBJ[  "aibj"] -=                IIJAK[    "mibj"]*T(1)[      "am"];

      GABIJ[  "abij"] +=                                   Tau[    "abij"];
      GABIJ[  "abij"] +=                  DAI[      "ai"]*T(1)[      "bj"];
      GABIJ[  "abij"] +=  (1.0/  2.0)*  GIJKL[    "mnij"]* Tau[    "abmn"];
      GABIJ[  "abij"] -=                  DIJ[      "mi"]*T(2)[    "abmj"];
      GABIJ[  "abij"] +=                  DAB[      "ae"]*T(2)[    "ebij"];
      GABIJ[  "abij"] -=                GAIBJ[    "amei"]*T(2)[    "ebmj"];

      GAIBJ[  "aibj"] -=                GIJAK[    "mibj"]*T(1)[      "am"];
      GAIBJ[  "aibj"] +=                GAIBC[    "aibe"]*T(1)[      "ej"];
      GAIBJ[  "aibj"] -=  (1.0/  2.0)*   L(2)[    "imbe"]*T(2)[    "aejm"];
      GAIBJ[  "aibj"] -=                 L(1)[      "ib"]*T(1)[      "aj"];

      GIJAK[  "ijak"] +=                 L(2)[    "ijae"]*T(1)[      "ek"];

      GAIBC[  "aibc"] -=                 L(2)[    "mibc"]*T(1)[      "am"];

      GAIJK[  "aijk"] +=                 L(1)[      "ie"]*T(2)[    "aejk"];
      GAIJK[  "aijk"] +=  (1.0/  2.0)*  GAIBC[    "aief"]* Tau[    "efjk"];
      GAIJK[  "aijk"] +=                  DIJ[      "ik"]*T(1)[      "aj"];
      GAIJK[  "aijk"] +=                IIJAK[    "miek"]*T(2)[    "aejm"];

      GABCI[  "abci"] -=                 L(1)[      "mc"]*T(2)[    "abmi"];
      GABCI[  "abci"] +=  (1.0/  2.0)*  GIJAK[    "mnci"]* Tau[    "abmn"];
      GABCI[  "abci"] +=                  DAB[      "ac"]*T(1)[      "bi"];
      GABCI[  "abci"] +=                IAIBC[    "amce"]*T(2)[    "beim"];

        DAI[    "ai"] +=                                  T(1)[      "ai"];

        DAB *= -1;
        DIJ *= -1;
      GIJAK *= -1;
      GAIBC *= -1;
      GAIJK *= -1;
      GABCI *= -1;

    U EIA = scalar(DIA*H.getIA());
    U EAI = scalar(DAI*H.getAI());
    U EAB = scalar(DAB*H.getAB());
    U EIJ = scalar(DIJ*H.getIJ());

    U EIJAB = 0.25*scalar(GIJAB*H.getIJAB());
    U EIJAK =  0.5*scalar(GIJAK*H.getIJAK());
    U EAIBC =  0.5*scalar(GAIBC*H.getAIBC());
    U EABCD = 0.25*scalar(GABCD*H.getABCD());
    U EIJKL = 0.25*scalar(GIJKL*H.getIJKL());
    U EAIBJ =      scalar(GAIBJ*H.getAIBJ());
    U EAIJK =  0.5*scalar(GAIJK*H.getAIJK());
    U EABCI =  0.5*scalar(GABCI*H.getABCI());
    U EABIJ = 0.25*scalar(GABIJ*H.getABIJ());

    printf("|AI|: %20.15f\n", DAI({0,0},{0,0}).norm(2));
    printf("|IA|: %20.15f\n", DIA({0,0},{0,0}).norm(2));
    printf("|AB|: %20.15f\n", DAB({0,0},{0,0}).norm(2));
    printf("|IJ|: %20.15f\n", DIJ({0,0},{0,0}).norm(2));
    printf("|IJAK|: %20.15f\n", GIJAK({0,1},{0,1}).norm(2));
    printf("|AIJK|: %20.15f\n", GAIJK({0,1},{0,1}).norm(2));
    printf("|AIBC|: %20.15f\n", GAIBC({1,0},{1,0}).norm(2));
    printf("|ABCI|: %20.15f\n", GABCI({1,0},{1,0}).norm(2));

    printf("\n");
    printf("EAB: %20.15f\n", EAB);
    printf("EIJ: %20.15f\n", EIJ);
    printf("EABCD: %20.15f\n", EABCD);
    printf("EIJKL: %20.15f\n", EIJKL);
    printf("EAIBJ: %20.15f\n", EAIBJ);
    printf("EAIJK: %20.15f\n", EAIJK+EIJAK);
    printf("EABCI: %20.15f\n", EABCI+EAIBC);
    printf("EABIJ: %20.15f\n", EABIJ+EIJAB);

    U E = EIA+EAI+EAB+EIJ+
          EIJAB+EABIJ+EABCI+EAIBC+EAIJK+EIJAK+
          EABCD+EIJKL+EAIBJ;

    printf("\n");
    printf("E: %20.15f\n", E);

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSDTQDensity);
REGISTER_TASK(aquarius::cc::CCSDTQDensity<double>,"ccsdtq_density","");
