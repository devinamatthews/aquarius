#include "lambdaccsdtq.hpp"

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
LambdaCCSDTQ<U>::LambdaCCSDTQ(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsdtq.Hbar", "Hbar"));
    reqs.push_back(Requirement("ccsdtq.T", "T"));
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("ccsdtq.L", "L", reqs));
}

template <typename U>
bool LambdaCCSDTQ<U>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;
    const PointGroup& group = occ.group;

    this->put   (      "L", new DeexcitationOperator<U,4>("L", arena, occ, vrt));
    this->puttmp(      "D", new Denominator         <U  >(H));
    this->puttmp(      "Z", new DeexcitationOperator<U,4>("Z", arena, occ, vrt));
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
    this->puttmp("WABCEJK", new SpinorbitalTensor   <U  >("W(abc,ejk)", arena, group, {vrt,occ}, {3,0}, {1,2}));
    this->puttmp("WABMIJK", new SpinorbitalTensor   <U  >("W(abm,ijk)", arena, group, {vrt,occ}, {2,1}, {0,3}));
    this->puttmp("WAMNIJK", new SpinorbitalTensor   <U  >("W(amn,ijk)", arena, group, {vrt,occ}, {1,2}, {0,3}));
    this->puttmp("WABMEJI", new SpinorbitalTensor   <U  >("W(abm,ejk)", arena, group, {vrt,occ}, {2,1}, {1,2}));
    this->puttmp("GAIJBCD", new SpinorbitalTensor   <U  >("G(aij,bcd)", arena, group, {vrt,occ}, {1,2}, {3,0}));
    this->puttmp("GIJKABL", new SpinorbitalTensor   <U  >("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1}));
    this->puttmp("GIJKALM", new SpinorbitalTensor   <U  >("G(ijk,alm)", arena, group, {vrt,occ}, {0,3}, {1,2}));

    auto& T = this->template get   <ExcitationOperator  <U,4>>("T");
    auto& L = this->template get   <DeexcitationOperator<U,4>>("L");
    auto& Z = this->template gettmp<DeexcitationOperator<U,4>>("Z");
    auto& D = this->template gettmp<Denominator         <U  >>("D");

    Z(0) = 0;
    L(0) = 1;
    L(1)[      "ia"] = T(1)[      "ai"];
    L(2)[    "ijab"] = T(2)[    "abij"];
    L(3)[  "ijkabc"] = 0;//T(3)[  "abcijk"];
    L(4)["ijklabcd"] = 0;//T(4)["abcdijkl"];

    const SpinorbitalTensor<U>&   FME =   H.getIA();
    const SpinorbitalTensor<U>& WMNEF = H.getIJAB();
    const SpinorbitalTensor<U>& WAMEF = H.getAIBC();
    const SpinorbitalTensor<U>& WABEF = H.getABCD();
    const SpinorbitalTensor<U>& WMNIJ = H.getIJKL();
    const SpinorbitalTensor<U>& WMNEJ = H.getIJAK();
    const SpinorbitalTensor<U>& WAMEI = H.getAIBJ();

    auto& WABCEJK = this->template gettmp<SpinorbitalTensor<U>>("WABCEJK");
    auto& WABMIJK = this->template gettmp<SpinorbitalTensor<U>>("WABMIJK");
    auto& WAMNIJK = this->template gettmp<SpinorbitalTensor<U>>("WAMNIJK");
    auto& WABMEJI = this->template gettmp<SpinorbitalTensor<U>>("WABMEJI");

    WABCEJK["abcejk"]  =    -WAMEI["amej"]*T(2)[    "bcmk"];
    WABCEJK["abcejk"] +=     WABEF["abef"]*T(2)[    "fcjk"];
    WABCEJK["abcejk"] -=       FME[  "me"]*T(3)[  "abcmjk"];
    WABCEJK["abcejk"] += 0.5*WMNEJ["mnej"]*T(3)[  "abcmnk"];
    WABCEJK["abcejk"] +=     WAMEF["amef"]*T(3)[  "fbcmjk"];
    WABCEJK["abcejk"] -= 0.5*WMNEF["mnef"]*T(4)["abcfmjkn"];

    WABMIJK["abmijk"]  =     WAMEI["amek"]*T(2)[    "ebij"];
    WABMIJK["abmijk"] -=     WMNIJ["nmjk"]*T(2)[    "abin"];
    WABMIJK["abmijk"] +=       FME[  "me"]*T(3)[  "abeijk"];
    WABMIJK["abmijk"] += 0.5*WAMEF["bmef"]*T(3)[  "aefijk"];
    WABMIJK["abmijk"] +=     WMNEJ["nmek"]*T(3)[  "abeijn"];
    WABMIJK["abmijk"] += 0.5*WMNEF["mnef"]*T(4)["abefijkn"];

    WAMNIJK["amnijk"]  =     WMNEJ["mnek"]*T(2)[    "aeij"];
    WAMNIJK["amnijk"] += 0.5*WMNEF["mnef"]*T(3)[  "aefijk"];

    WABMEJI["abmeji"]  =     WAMEF["amef"]*T(2)[    "bfji"];
    WABMEJI["abmeji"] -=     WMNEJ["nmei"]*T(2)[    "abnj"];
    WABMEJI["abmeji"] +=     WMNEF["mnef"]*T(3)[  "abfnji"];

    Iterative<U>::run(dag, arena);

    this->put("energy", new U(this->energy()));
    this->put("convergence", new U(this->conv()));

    return true;
}

template <typename U>
void LambdaCCSDTQ<U>::iterate(const Arena& arena)
{
    const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    const SpinorbitalTensor<U>&   FME =   H.getIA();
    const SpinorbitalTensor<U>&   FAE =   H.getAB();
    const SpinorbitalTensor<U>&   FMI =   H.getIJ();
    const SpinorbitalTensor<U>& WMNEF = H.getIJAB();
    const SpinorbitalTensor<U>& WAMEF = H.getAIBC();
    const SpinorbitalTensor<U>& WABEJ = H.getABCI();
    const SpinorbitalTensor<U>& WABEF = H.getABCD();
    const SpinorbitalTensor<U>& WMNIJ = H.getIJKL();
    const SpinorbitalTensor<U>& WMNEJ = H.getIJAK();
    const SpinorbitalTensor<U>& WAMIJ = H.getAIJK();
    const SpinorbitalTensor<U>& WAMEI = H.getAIBJ();

    auto& T = this->template get   <ExcitationOperator  <U,4>>("T");
    auto& L = this->template get   <DeexcitationOperator<U,4>>("L");
    auto& D = this->template gettmp<Denominator         <U  >>("D");
    auto& Z = this->template gettmp<DeexcitationOperator<U,4>>("Z");

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
    auto& WABCEJK = this->template gettmp<SpinorbitalTensor<U>>("WABCEJK");
    auto& WABMIJK = this->template gettmp<SpinorbitalTensor<U>>("WABMIJK");
    auto& WAMNIJK = this->template gettmp<SpinorbitalTensor<U>>("WAMNIJK");
    auto& WABMEJI = this->template gettmp<SpinorbitalTensor<U>>("WABMEJI");
    auto& GAIJBCD = this->template gettmp<SpinorbitalTensor<U>>("GAIJBCD");
    auto& GIJKABL = this->template gettmp<SpinorbitalTensor<U>>("GIJKABL");
    auto& GIJKALM = this->template gettmp<SpinorbitalTensor<U>>("GIJKALM");

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
     * Intermediates for Lambda-CCSDTQ
     */
        DIJ[    "ij"]  =  (1.0/144.0)* T(4)["efghjmno"]*   L(4)["imnoefgh"];
        DAB[    "ab"]  = -(1.0/144.0)* T(4)["aefgmnop"]*   L(4)["mnopbefg"];

      GIJAK[  "ijak"]  =  (1.0/ 12.0)* T(3)[  "efgkmn"]*   L(4)["ijmnaefg"];
      GAIBC[  "aibc"]  = -(1.0/ 12.0)* T(3)[  "aefmno"]*   L(4)["minobcef"];

    GIJKABL["ijkabl"]  =  (1.0/  2.0)* T(2)[    "eflm"]*   L(4)["ijkmabef"];
    GAIJBCD["aijbcd"]  = -(1.0/  2.0)* T(2)[    "aemn"]*   L(4)["mijnbcde"];
    GIJKALM["ijkalm"]  =  (1.0/  6.0)* T(3)[  "efglmn"]*   L(4)["ijknaefg"];

      GAIJK[  "aijk"]  =  (1.0/ 12.0)* T(4)["aefgjkmn"]*   L(3)[  "imnefg"];
      GAIJK[  "aijk"] +=  (1.0/  4.0)* T(3)[  "efamnj"]*GIJKABL[  "mniefk"];
      GAIJK[  "aijk"] +=  (1.0/  6.0)* T(3)[  "efgjkm"]*GAIJBCD[  "aimefg"];

      GABCI[  "abci"]  = -(1.0/ 12.0)* T(4)["abefmino"]*   L(3)[  "mnocef"];
      GABCI[  "abci"] +=  (1.0/  4.0)* T(3)[  "efbmni"]*GAIJBCD[  "amncef"];
      GABCI[  "abci"] +=  (1.0/  6.0)* T(3)[  "eabmno"]*GIJKABL[  "mnoeci"];

      GABCD[  "abcd"]  =  (1.0/ 48.0)* T(4)["abefmnop"]*   L(4)["mnopcdef"];
      GABCD[  "abcd"] -=  (1.0/  4.0)* T(2)[    "bemn"]*GAIJBCD[  "amncde"];

      GAIBJ[  "aibj"]  = -(1.0/ 36.0)* T(4)["aefgjmno"]*   L(4)["imnobefg"];
      GAIBJ[  "aibj"] +=  (1.0/  2.0)* T(2)[    "eamn"]*GIJKABL[  "imnbej"];

      GIJKL[  "ijkl"]  =  (1.0/ 48.0)* T(4)["efghklmn"]*   L(4)["ijmnefgh"];
      GIJKL[  "ijkl"] +=  (1.0/  4.0)* T(2)[    "efmk"]*GIJKABL[  "mijefl"];

        DAI[    "ai"]  =  (1.0/ 36.0)* T(4)["aefgimno"]*   L(3)[  "mnoefg"];
        DAI[    "ai"] -=  (1.0/  2.0)* T(2)[    "eamn"]*  GIJAK[    "mnei"];
        DAI[    "ai"] +=  (1.0/  2.0)* T(2)[    "efim"]*  GAIBC[    "amef"];
    /*
     **************************************************************************/

    /***************************************************************************
     *
     * Lambda-CCSDTQ iteration
     */
    Z(1)[      "ia"] +=                FME[    "ie"]*    DAB[      "ea"];
    Z(1)[      "ia"] -=                FME[    "ma"]*    DIJ[      "im"];
    Z(1)[      "ia"] -=              WMNEJ[  "inam"]*    DIJ[      "mn"];
    Z(1)[      "ia"] -=              WAMEF[  "fiea"]*    DAB[      "ef"];
    Z(1)[      "ia"] +=              WMNEF[  "miea"]*    DAI[      "em"];
    Z(1)[      "ia"] -= (1.0/ 2.0)*  WABEF[  "efga"]*  GAIBC[    "gief"];
    Z(1)[      "ia"] +=              WAMEI[  "eifm"]*  GAIBC[    "fmea"];
    Z(1)[      "ia"] -=              WAMEI[  "eman"]*  GIJAK[    "inem"];
    Z(1)[      "ia"] += (1.0/ 2.0)*  WMNIJ[  "imno"]*  GIJAK[    "noam"];
    Z(1)[      "ia"] -= (1.0/ 2.0)*  WAMEF[  "gief"]*  GABCD[    "efga"];
    Z(1)[      "ia"] +=              WAMEF[  "fmea"]*  GAIBJ[    "eifm"];
    Z(1)[      "ia"] -=              WMNEJ[  "inem"]*  GAIBJ[    "eman"];
    Z(1)[      "ia"] += (1.0/ 2.0)*  WMNEJ[  "noam"]*  GIJKL[    "imno"];
    Z(1)[      "ia"] += (1.0/ 2.0)*  WMNEF[  "imef"]*  GABCI[    "efam"];
    Z(1)[      "ia"] -= (1.0/ 2.0)*  WMNEF[  "mnea"]*  GAIJK[    "eimn"];

    Z(2)[    "ijab"] -=              WMNEF[  "mjab"]*    DIJ[      "im"];
    Z(2)[    "ijab"] +=              WMNEF[  "ijeb"]*    DAB[      "ea"];
    Z(2)[    "ijab"] += (1.0/ 2.0)*  WMNEF[  "ijef"]*  GABCD[    "efab"];
    Z(2)[    "ijab"] +=              WMNEF[  "imea"]*  GAIBJ[    "ejbm"];
    Z(2)[    "ijab"] += (1.0/ 2.0)*  WMNEF[  "mnab"]*  GIJKL[    "ijmn"];
    Z(2)[    "ijab"] -=              WAMEF[  "fiae"]*  GAIBC[    "ejbf"];
    Z(2)[    "ijab"] -=              WMNEJ[  "ijem"]*  GAIBC[    "emab"];
    Z(2)[    "ijab"] -=              WAMEF[  "emab"]*  GIJAK[    "ijem"];
    Z(2)[    "ijab"] -=              WMNEJ[  "niam"]*  GIJAK[    "mjbn"];
    Z(2)[    "ijab"] += (1.0/12.0)*WABCEJK["efgamn"]*   L(4)["ijmnebfg"];
    Z(2)[    "ijab"] -= (1.0/12.0)*WABMIJK["efjmno"]*   L(4)["mnioefab"];

    Z(3)[  "ijkabc"] +=              WMNEF[  "ijae"]*  GAIBC[    "ekbc"];
    Z(3)[  "ijkabc"] -=              WMNEF[  "mkbc"]*  GIJAK[    "ijam"];
    Z(3)[  "ijkabc"] -=              WAMEF[  "embc"]*GIJKABL[  "ijkaem"];
    Z(3)[  "ijkabc"] += (1.0/ 2.0)*  WMNEF[  "mnbc"]*GIJKALM[  "ijkamn"];
    Z(3)[  "ijkabc"] += (1.0/ 2.0)*  WABEJ[  "efam"]*   L(4)["ijkmebcf"];
    Z(3)[  "ijkabc"] -= (1.0/ 2.0)*  WAMIJ[  "eknm"]*   L(4)["ijmnabce"];
    Z(3)[  "ijkabc"] -= (1.0/ 4.0)*WABMEJI["efkcnm"]*   L(4)["ijmnabef"];
    Z(3)[  "ijkabc"] += (1.0/ 6.0)*WAMNIJK["eijmno"]*   L(4)["mnokeabc"];

    Z(4)["ijklabcd"]  =              WMNEF[  "ijab"]*   L(2)[    "klcd"];
    Z(4)["ijklabcd"] +=                FME[    "ia"]*   L(3)[  "jklbcd"];
    Z(4)["ijklabcd"] +=              WAMEF[  "ejab"]*   L(3)[  "iklecd"];
    Z(4)["ijklabcd"] -=              WMNEJ[  "ijam"]*   L(3)[  "mklbcd"];
    Z(4)["ijklabcd"] +=              WMNEF[  "ijae"]*GAIJBCD[  "eklbcd"];
    Z(4)["ijklabcd"] -=              WMNEF[  "mlcd"]*GIJKABL[  "ijkabm"];
    Z(4)["ijklabcd"] +=                FAE[    "ea"]*   L(4)["ijklebcd"];
    Z(4)["ijklabcd"] -=                FMI[    "im"]*   L(4)["mjklabcd"];
    Z(4)["ijklabcd"] += (1.0/ 2.0)*  WABEF[  "efab"]*   L(4)["ijklefcd"];
    Z(4)["ijklabcd"] += (1.0/ 2.0)*  WMNIJ[  "ijmn"]*   L(4)["mnklabcd"];
    Z(4)["ijklabcd"] +=              WAMEI[  "eiam"]*   L(4)["mjklbecd"];
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
sub_iterations?
    int 0,
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::LambdaCCSDTQ);
REGISTER_TASK(aquarius::cc::LambdaCCSDTQ<double>,"lambdaccsdtq",spec);
