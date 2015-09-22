#include "lambdacc4.hpp"

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
LambdaCC4<U>::LambdaCC4(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    reqs.push_back(Requirement("cc4.Hbar", "Hbar"));
    reqs.push_back(Requirement("cc4.T", "T"));
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("cc4.L", "L", reqs));
}

template <typename U>
bool LambdaCC4<U>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& H    = this->template get<  TwoElectronOperator<U>>("H");
    const auto& Hbar = this->template get<STTwoElectronOperator<U>>("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;
    const PointGroup& group = occ.group;

    this->put   (      "L", new DeexcitationOperator<U,3>("L", arena, occ, vrt));
    this->puttmp(      "D", new Denominator         <U  >(H));
    this->puttmp(      "Z", new DeexcitationOperator<U,3>("Z", arena, occ, vrt));
    this->puttmp(    "DAB", new SpinorbitalTensor   <U  >(     "D(ab)", arena, group, {vrt,occ}, {1,0}, {1,0}));
    this->puttmp(    "DIJ", new SpinorbitalTensor   <U  >(     "D(ij)", arena, group, {vrt,occ}, {0,1}, {0,1}));
    this->puttmp(    "DAI", new SpinorbitalTensor   <U  >(     "D(ai)", arena, group, {vrt,occ}, {1,0}, {0,1}));
    this->puttmp(  "GABCD", new SpinorbitalTensor   <U  >(  "G(ab,cd)", arena, group, {vrt,occ}, {2,0}, {2,0}));
    this->puttmp(  "GAIBJ", new SpinorbitalTensor   <U  >(  "G(ai,bj)", arena, group, {vrt,occ}, {1,1}, {1,1}));
    this->puttmp(  "GIJKL", new SpinorbitalTensor   <U  >(  "G(ij,kl)", arena, group, {vrt,occ}, {0,2}, {0,2}));
    this->puttmp(  "GAIBC", new SpinorbitalTensor   <U  >(  "G(ai,bc)", arena, group, {vrt,occ}, {1,1}, {2,0}));
    this->puttmp(  "GIJAK", new SpinorbitalTensor   <U  >(  "G(ij,ak)", arena, group, {vrt,occ}, {0,2}, {1,1}));
    this->puttmp(  "GABCI", new SpinorbitalTensor   <U  >(  "G(ab,ci)", arena, group, {vrt,occ}, {2,0}, {1,1}));
    this->puttmp(  "GAIJK", new SpinorbitalTensor   <U  >(  "G(ai,jk)", arena, group, {vrt,occ}, {1,1}, {0,2}));
    this->puttmp("WTWABEF", new SpinorbitalTensor   <U  >( "W~(ab,ef)", arena, group, {vrt,occ}, {2,0}, {2,0}));
    this->puttmp("WTWAMEI", new SpinorbitalTensor   <U  >( "W~(am,ei)", arena, group, {vrt,occ}, {1,1}, {1,1}));
    this->puttmp("WTWMNIJ", new SpinorbitalTensor   <U  >( "W~(mn,ij)", arena, group, {vrt,occ}, {0,2}, {0,2}));
    this->puttmp("WTWABEJ", new SpinorbitalTensor   <U  >( "W~(ab,ej)", arena, group, {vrt,occ}, {2,0}, {1,1}));
    this->puttmp("WTWAMIJ", new SpinorbitalTensor   <U  >( "W~(am,ij)", arena, group, {vrt,occ}, {1,1}, {0,2}));
    this->puttmp("WABCEJK", new SpinorbitalTensor   <U  >("W(abc,ejk)", arena, group, {vrt,occ}, {3,0}, {1,2}));
    this->puttmp("WABMIJK", new SpinorbitalTensor   <U  >("W(abm,ijk)", arena, group, {vrt,occ}, {2,1}, {0,3}));
    this->puttmp("GAIJBCD", new SpinorbitalTensor   <U  >("G(aij,bcd)", arena, group, {vrt,occ}, {1,2}, {3,0}));
    this->puttmp("GIJKABL", new SpinorbitalTensor   <U  >("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1}));

    this->puttmp("T4", new SpinorbitalTensor<U>("T(abcd,ijkl)", arena, group, {vrt, occ}, {4, 0}, {0, 4}));
    this->puttmp("L4", new SpinorbitalTensor<U>("L(ijkl,abcd)", arena, group, {vrt, occ}, {0, 4}, {4, 0}));

    auto& T = this->template get   <ExcitationOperator  <U,3>>("T");
    auto& L = this->template get   <DeexcitationOperator<U,3>>("L");
    auto& Z = this->template gettmp<DeexcitationOperator<U,3>>("Z");
    auto& D = this->template gettmp<Denominator         <U  >>("D");

    Z(0) = 0;
    L(0) = 1;
    L(1)[    "ia"] = T(1)[  "ai"];
    L(2)[  "ijab"] = T(2)["abij"];
    L(3)["ijkabc"] = 0;

    const SpinorbitalTensor<U>& VMNEF = H.getIJAB();
    const SpinorbitalTensor<U>& VABEJ = H.getABCI();
    const SpinorbitalTensor<U>& VAMIJ = H.getAIJK();
    const SpinorbitalTensor<U>& VABEF = H.getABCD();
    const SpinorbitalTensor<U>& VMNIJ = H.getIJKL();
    const SpinorbitalTensor<U>& VAMEI = H.getAIBJ();
    const SpinorbitalTensor<U>& VAMEF = H.getAIBC();
    const SpinorbitalTensor<U>& VMNEJ = H.getIJAK();

    const SpinorbitalTensor<U>& WAMEF = Hbar.getAIBC();
    const SpinorbitalTensor<U>& WMNEJ = Hbar.getIJAK();

    SpinorbitalTensor<U> T12(T(2));
    T12["abij"]  = 0.5*T(1)["ai"]*T(1)["bj"];

    auto& T4 = this->template gettmp<SpinorbitalTensor<U>>("T4");

    auto& WTWABEF = this->template gettmp<SpinorbitalTensor<U>>("WTWABEF");
    auto& WTWAMEI = this->template gettmp<SpinorbitalTensor<U>>("WTWAMEI");
    auto& WTWMNIJ = this->template gettmp<SpinorbitalTensor<U>>("WTWMNIJ");
    auto& WTWABEJ = this->template gettmp<SpinorbitalTensor<U>>("WTWABEJ");
    auto& WTWAMIJ = this->template gettmp<SpinorbitalTensor<U>>("WTWAMIJ");

    auto& WABCEJK = this->template gettmp<SpinorbitalTensor<U>>("WABCEJK");
    auto& WABMIJK = this->template gettmp<SpinorbitalTensor<U>>("WABMIJK");

    WTWAMEI[  "amei"]  =       VAMEI["amei"];
    WTWAMEI[  "amei"] +=       VAMEF["amef"]*T(1)[  "fi"];
    WTWAMEI[  "amei"] -=       WMNEJ["nmei"]*T(1)[  "an"];

    WTWMNIJ[  "mnij"]  =       VMNIJ["mnij"];
    WTWMNIJ[  "mnij"] +=       VMNEJ["mnej"]*T(1)[  "ei"];
    WTWMNIJ[  "mnij"] += 0.5*  VMNEF["mnef"]* T12["efij"];

    WTWABEF[  "abef"]  =       VABEF["abef"];
    WTWABEF[  "abef"] -=       VAMEF["amef"]*T(1)[  "bm"];
    WTWABEF[  "abef"] += 0.5*  VMNEF["mnef"]* T12["abmn"];

    WTWAMIJ[  "amij"]  =       VAMIJ["amij"];
    WTWAMIJ[  "amij"] += 0.5*  VAMEF["amef"]* T12["efij"];
    WTWAMIJ[  "amij"] +=       VAMEI["amej"]*T(1)[  "ei"];
    WTWAMIJ[  "amij"] -=     WTWMNIJ["nmij"]*T(1)[  "an"];

    WTWABEJ[  "abej"]  =       VABEJ["abej"];
    WTWABEJ[  "abej"] += 0.5*  VMNEJ["mnej"]* T12["abmn"];
    WTWABEJ[  "abej"] -=       VAMEI["amej"]*T(1)[  "bm"];
    WTWABEJ[  "abej"] +=     WTWABEF["abef"]*T(1)[  "fj"];

    WABCEJK["abcejk"]  = 0.5*WTWABEF["abef"]*T(2)["fcjk"];

    WABMIJK["abmijk"]  =     WTWAMEI["amek"]*T(2)["ebij"];
    WABMIJK["abmijk"] -= 0.5*WTWMNIJ["mnkj"]*T(2)["abin"];

    T4["abcdijkl"]  = WABCEJK["abcejk"]*T(2)[  "edil"];
    T4["abcdijkl"] -= WABMIJK["abmijk"]*T(2)[  "cdml"];
    T4["abcdijkl"] += WTWABEJ[  "abej"]*T(3)["ecdikl"];
    T4["abcdijkl"] -= WTWAMIJ[  "amij"]*T(3)["bcdmkl"];

    T4.weight({&D.getDA(), &D.getDI()}, {&D.getDa(), &D.getDi()});

    WABCEJK["abcejk"]  = WTWABEF["abef"]*T(2)["fcjk"];
    WABCEJK["abcejk"] -= WTWAMEI["amej"]*T(2)["bcmk"];

    WABMIJK["abmijk"]  = WTWAMEI["amek"]*T(2)["ebij"];
    WABMIJK["abmijk"] -= WTWMNIJ["nmjk"]*T(2)["abin"];

    Iterative<U>::run(dag, arena);

    this->put("energy", new U(this->energy()));
    this->put("convergence", new U(this->conv()));

    return true;
}

template <typename U>
void LambdaCC4<U>::iterate(const Arena& arena)
{
    const auto& H    = this->template get<  TwoElectronOperator<U>>("H");
    const auto& Hbar = this->template get<STTwoElectronOperator<U>>("Hbar");

    const SpinorbitalTensor<U>&   FME =   Hbar.getIA();
    const SpinorbitalTensor<U>&   FAE =   Hbar.getAB();
    const SpinorbitalTensor<U>&   FMI =   Hbar.getIJ();
    const SpinorbitalTensor<U>& WMNEF = Hbar.getIJAB();
    const SpinorbitalTensor<U>& WAMEF = Hbar.getAIBC();
    const SpinorbitalTensor<U>& WABEJ = Hbar.getABCI();
    const SpinorbitalTensor<U>& WABEF = Hbar.getABCD();
    const SpinorbitalTensor<U>& WMNIJ = Hbar.getIJKL();
    const SpinorbitalTensor<U>& WMNEJ = Hbar.getIJAK();
    const SpinorbitalTensor<U>& WAMIJ = Hbar.getAIJK();
    const SpinorbitalTensor<U>& WAMEI = Hbar.getAIBJ();

    const SpinorbitalTensor<U>& VABEJ = H.getABCI();
    const SpinorbitalTensor<U>& VAMIJ = H.getAIJK();

    auto& T = this->template get   <ExcitationOperator  <U,3>>("T");
    auto& L = this->template get   <DeexcitationOperator<U,3>>("L");
    auto& D = this->template gettmp<Denominator         <U  >>("D");
    auto& Z = this->template gettmp<DeexcitationOperator<U,3>>("Z");

    auto& T4 = this->template gettmp<SpinorbitalTensor<U>>("T4");
    auto& L4 = this->template gettmp<SpinorbitalTensor<U>>("L4");

    auto&     DIJ = this->template gettmp<SpinorbitalTensor<U>>(    "DIJ");
    auto&     DAB = this->template gettmp<SpinorbitalTensor<U>>(    "DAB");
    auto&     DAI = this->template gettmp<SpinorbitalTensor<U>>(    "DAI");
    auto&   GABCD = this->template gettmp<SpinorbitalTensor<U>>(  "GABCD");
    auto&   GAIBJ = this->template gettmp<SpinorbitalTensor<U>>(  "GAIBJ");
    auto&   GIJKL = this->template gettmp<SpinorbitalTensor<U>>(  "GIJKL");
    auto&   GAIBC = this->template gettmp<SpinorbitalTensor<U>>(  "GAIBC");
    auto&   GIJAK = this->template gettmp<SpinorbitalTensor<U>>(  "GIJAK");
    auto&   GABCI = this->template gettmp<SpinorbitalTensor<U>>(  "GABCI");
    auto&   GAIJK = this->template gettmp<SpinorbitalTensor<U>>(  "GAIJK");
    auto& WTWABEF = this->template gettmp<SpinorbitalTensor<U>>("WTWABEF");
    auto& WTWAMEI = this->template gettmp<SpinorbitalTensor<U>>("WTWAMEI");
    auto& WTWMNIJ = this->template gettmp<SpinorbitalTensor<U>>("WTWMNIJ");
    auto& WTWABEJ = this->template gettmp<SpinorbitalTensor<U>>("WTWABEJ");
    auto& WTWAMIJ = this->template gettmp<SpinorbitalTensor<U>>("WTWAMIJ");
    auto& WABCEJK = this->template gettmp<SpinorbitalTensor<U>>("WABCEJK");
    auto& WABMIJK = this->template gettmp<SpinorbitalTensor<U>>("WABMIJK");
    auto& GAIJBCD = this->template gettmp<SpinorbitalTensor<U>>("GAIJBCD");
    auto& GIJKABL = this->template gettmp<SpinorbitalTensor<U>>("GIJKABL");

    /***************************************************************************
     *
     * Intermediates for Lambda-CCSD
     */
    DIJ["ij"]  =  0.5*T(2)["efjm"]*L(2)["imef"];
    DAB["ab"]  = -0.5*T(2)["aemn"]*L(2)["mnbe"];
    /*
     **************************************************************************/

    /***************************************************************************
     *
     * Lambda-CCSD iteration
     */
    Z(1)[  "ia"]  =       FME[  "ia"];
    Z(1)[  "ia"] +=       FAE[  "ea"]*L(1)[  "ie"];
    Z(1)[  "ia"] -=       FMI[  "im"]*L(1)[  "ma"];
    Z(1)[  "ia"] -=     WAMEI["eiam"]*L(1)[  "me"];
    Z(1)[  "ia"] += 0.5*WABEJ["efam"]*L(2)["imef"];
    Z(1)[  "ia"] -= 0.5*WAMIJ["eimn"]*L(2)["mnea"];
    Z(1)[  "ia"] -=     WMNEJ["inam"]* DIJ[  "mn"];
    Z(1)[  "ia"] -=     WAMEF["fiea"]* DAB[  "ef"];

    Z(2)["ijab"]  =     WMNEF["ijab"];
    Z(2)["ijab"] +=       FME[  "ia"]*L(1)[  "jb"];
    Z(2)["ijab"] +=     WAMEF["ejab"]*L(1)[  "ie"];
    Z(2)["ijab"] -=     WMNEJ["ijam"]*L(1)[  "mb"];
    Z(2)["ijab"] +=       FAE[  "ea"]*L(2)["ijeb"];
    Z(2)["ijab"] -=       FMI[  "im"]*L(2)["mjab"];
    Z(2)["ijab"] += 0.5*WABEF["efab"]*L(2)["ijef"];
    Z(2)["ijab"] += 0.5*WMNIJ["ijmn"]*L(2)["mnab"];
    Z(2)["ijab"] +=     WAMEI["eiam"]*L(2)["mjbe"];
    Z(2)["ijab"] -=     WMNEF["mjab"]* DIJ[  "im"];
    Z(2)["ijab"] +=     WMNEF["ijeb"]* DAB[  "ea"];
    /*
     **************************************************************************/

    /***************************************************************************
     *
     * Intermediates for Lambda-CCSDT
     */
      DIJ[  "ij"]  =  (1.0/12.0)*T(3)["efgjmn"]* L(3)["imnefg"];
      DAB[  "ab"]  = -(1.0/12.0)*T(3)["aefmno"]* L(3)["mnobef"];

    GABCD["abcd"]  =   (1.0/6.0)*T(3)["abemno"]* L(3)["mnocde"];
    GAIBJ["aibj"]  =       -0.25*T(3)["aefjmn"]* L(3)["imnbef"];
    GIJKL["ijkl"]  =   (1.0/6.0)*T(3)["efgklm"]* L(3)["ijmefg"];

    GIJAK["ijak"]  =         0.5*T(2)[  "efkm"]* L(3)["ijmaef"];
    GAIBC["aibc"]  =        -0.5*T(2)[  "aemn"]* L(3)["minbce"];

      DAI[  "ai"]  =        0.25*T(3)["aefimn"]* L(2)[  "mnef"];
      DAI[  "ai"] -=         0.5*T(2)[  "eamn"]*GIJAK[  "mnei"];
    /*
     **************************************************************************/

    /***************************************************************************
     *
     * Lambda-CCSDT iteration
     */
    Z(1)[    "ia"] +=       FME[  "ie"]*  DAB[    "ea"];
    Z(1)[    "ia"] -=       FME[  "ma"]*  DIJ[    "im"];
    Z(1)[    "ia"] -=     WMNEJ["inam"]*  DIJ[    "mn"];
    Z(1)[    "ia"] -=     WAMEF["fiea"]*  DAB[    "ef"];
    Z(1)[    "ia"] +=     WMNEF["miea"]*  DAI[    "em"];
    Z(1)[    "ia"] -= 0.5*WABEF["efga"]*GAIBC[  "gief"];
    Z(1)[    "ia"] +=     WAMEI["eifm"]*GAIBC[  "fmea"];
    Z(1)[    "ia"] -=     WAMEI["eman"]*GIJAK[  "inem"];
    Z(1)[    "ia"] += 0.5*WMNIJ["imno"]*GIJAK[  "noam"];
    Z(1)[    "ia"] -= 0.5*WAMEF["gief"]*GABCD[  "efga"];
    Z(1)[    "ia"] +=     WAMEF["fmea"]*GAIBJ[  "eifm"];
    Z(1)[    "ia"] -=     WMNEJ["inem"]*GAIBJ[  "eman"];
    Z(1)[    "ia"] += 0.5*WMNEJ["noam"]*GIJKL[  "imno"];

    Z(2)[  "ijab"] -=     WMNEF["mjab"]*  DIJ[    "im"];
    Z(2)[  "ijab"] +=     WMNEF["ijeb"]*  DAB[    "ea"];
    Z(2)[  "ijab"] += 0.5*WMNEF["ijef"]*GABCD[  "efab"];
    Z(2)[  "ijab"] +=     WMNEF["imea"]*GAIBJ[  "ejbm"];
    Z(2)[  "ijab"] += 0.5*WMNEF["mnab"]*GIJKL[  "ijmn"];
    Z(2)[  "ijab"] -=     WAMEF["fiae"]*GAIBC[  "ejbf"];
    Z(2)[  "ijab"] -=     WMNEJ["ijem"]*GAIBC[  "emab"];
    Z(2)[  "ijab"] -=     WAMEF["emab"]*GIJAK[  "ijem"];
    Z(2)[  "ijab"] -=     WMNEJ["niam"]*GIJAK[  "mjbn"];
    Z(2)[  "ijab"] += 0.5*WABEJ["efbm"]* L(3)["ijmaef"];
    Z(2)[  "ijab"] -= 0.5*WAMIJ["ejnm"]* L(3)["imnabe"];

    Z(3)["ijkabc"]  =     WMNEF["ijab"]* L(1)[    "kc"];
    Z(3)["ijkabc"] +=       FME[  "ia"]* L(2)[  "jkbc"];
    Z(3)["ijkabc"] +=     WAMEF["ekbc"]* L(2)[  "ijae"];
    Z(3)["ijkabc"] -=     WMNEJ["ijam"]* L(2)[  "mkbc"];
    Z(3)["ijkabc"] +=     WMNEF["ijae"]*GAIBC[  "ekbc"];
    Z(3)["ijkabc"] -=     WMNEF["mkbc"]*GIJAK[  "ijam"];
    Z(3)["ijkabc"] +=       FAE[  "ea"]* L(3)["ijkebc"];
    Z(3)["ijkabc"] -=       FMI[  "im"]* L(3)["mjkabc"];
    Z(3)["ijkabc"] += 0.5*WABEF["efab"]* L(3)["ijkefc"];
    Z(3)["ijkabc"] += 0.5*WMNIJ["ijmn"]* L(3)["mnkabc"];
    Z(3)["ijkabc"] +=     WAMEI["eiam"]* L(3)["mjkbec"];
    /*
     **************************************************************************/

    /***************************************************************************
     *
     * Lambda-CC4 iteration
     */
         L4["ijklabcd"]  =              WMNEF[    "ijab"]*   L(2)[    "klcd"];
         L4["ijklabcd"] +=                FME[      "ia"]*   L(3)[  "jklbcd"];
         L4["ijklabcd"] +=              WAMEF[    "ejab"]*   L(3)[  "iklecd"];
         L4["ijklabcd"] -=              WMNEJ[    "ijam"]*   L(3)[  "mklbcd"];

    L4.weight({&D.getDA(), &D.getDI()}, {&D.getDa(), &D.getDi()});

      GIJAK[    "ijak"]  =  (1.0/ 12.0)* T(3)[  "efgkmn"]*     L4["ijmnaefg"];
      GAIBC[    "aibc"]  = -(1.0/ 12.0)* T(3)[  "aefmno"]*     L4["minobcef"];

      GAIJK[    "aijk"]  =  (1.0/ 12.0)*   T4["aefgjkmn"]*   L(3)[  "imnefg"];
      GABCI[    "abci"]  = -(1.0/ 12.0)*   T4["abefmino"]*   L(3)[  "mnocef"];

        DAI[      "ai"]  =  (1.0/ 36.0)*   T4["aefgimno"]*   L(3)[  "mnoefg"];

    GIJKABL[  "ijkabl"]  =  (1.0/  2.0)* T(2)[    "eflm"]*     L4["ijkmabef"];
    GAIJBCD[  "aijbcd"]  = -(1.0/  2.0)* T(2)[    "aemn"]*     L4["mijnbcde"];

      GABCD[    "abcd"]  = -(1.0/  4.0)* T(2)[    "bemn"]*GAIJBCD[  "amncde"];
      GAIBJ[    "aibj"]  =  (1.0/  2.0)* T(2)[    "eamn"]*GIJKABL[  "imnbej"];
      GIJKL[    "ijkl"]  =  (1.0/  4.0)* T(2)[    "efmk"]*GIJKABL[  "mijefl"];

       Z(1)[      "ia"] +=              WMNEF[    "miea"]*    DAI[      "em"];
       Z(1)[      "ia"] += (1.0/ 2.0)*  WMNEF[    "imef"]*  GABCI[    "efam"];
       Z(1)[      "ia"] -= (1.0/ 2.0)*  WMNEF[    "mnea"]*  GAIJK[    "eimn"];
       Z(1)[      "ia"] -= (1.0/ 2.0)*WTWABEF[    "efga"]*  GAIBC[    "gief"];
       Z(1)[      "ia"] +=            WTWAMEI[    "eifm"]*  GAIBC[    "fmea"];
       Z(1)[      "ia"] -=            WTWAMEI[    "eman"]*  GIJAK[    "inem"];
       Z(1)[      "ia"] += (1.0/ 2.0)*WTWMNIJ[    "imno"]*  GIJAK[    "noam"];
       Z(1)[      "ia"] -= (1.0/ 2.0)*  WAMEF[    "gief"]*  GABCD[    "efga"];
       Z(1)[      "ia"] +=              WAMEF[    "fmea"]*  GAIBJ[    "eifm"];
       Z(1)[      "ia"] -=              WMNEJ[    "inem"]*  GAIBJ[    "eman"];
       Z(1)[      "ia"] += (1.0/ 2.0)*  WMNEJ[    "noam"]*  GIJKL[    "imno"];

       Z(2)[    "ijab"] += (1.0/12.0)*WABCEJK[  "efgamn"]*     L4["ijmnebfg"];
       Z(2)[    "ijab"] -= (1.0/12.0)*WABMIJK[  "efjmno"]*     L4["mnioefab"];

       Z(3)[  "ijkabc"] += (1.0/ 2.0)*WTWABEJ[    "efam"]*     L4["ijkmebcf"];
       Z(3)[  "ijkabc"] -= (1.0/ 2.0)*WTWAMIJ[    "eknm"]*     L4["ijmnabce"];
    /*
     **************************************************************************/

    Z.weight(D);
    L += Z;

    this->energy() = 0.25*real(scalar(conj(WMNEF)*L(2)));
    this->conv() = Z.norm(00);

    diis.extrapolate(L, Z);
}

}
}

static const char* spec = R"!(

convergence?
    double 1e-9,
max_iterations?
    int 50,
conv_type?
    enum { MAXE, RMSE, MAE },
diis?
{
    damping?
        double 0.0,
    start?
        int 1,
    order?
        int 5,
    jacobi?
        bool false
}

)!";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::LambdaCC4);
REGISTER_TASK(aquarius::cc::LambdaCC4<double>,"lambdacc4",spec);
