#include "ccsdt_q_n_opt.hpp"

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
CCSDT_Q_N_opt<U>::CCSDT_Q_N_opt(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    reqs.push_back(Requirement("ccsdt.Hbar", "Hbar"));
    reqs.push_back(Requirement("ccsdt.T", "T"));
    reqs.push_back(Requirement("ccsdt.L", "L"));
    this->addProduct(Product("double", "E(2)", reqs));
    this->addProduct(Product("double", "E(3)", reqs));
    this->addProduct(Product("double", "E(4)", reqs));
    this->addProduct(Product("double", "E(5)", reqs));
    this->addProduct(Product("double", "E(6)", reqs));
}

template <typename U>
bool CCSDT_Q_N_opt<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U>>("H");
    const STTwoElectronOperator<U>& Hbar = this->template get<STTwoElectronOperator<U>>("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;
    const PointGroup& group = occ.group;

    Denominator<U> D(H);
    const ExcitationOperator  <U,3>& T = this->template get<ExcitationOperator  <U,3>>("T");
    const DeexcitationOperator<U,3>& L = this->template get<DeexcitationOperator<U,3>>("L");

    SpinorbitalTensor<U> FME(Hbar.getIA());
    SpinorbitalTensor<U> FAE(Hbar.getAB());
    SpinorbitalTensor<U> FMI(Hbar.getIJ());
    FME -= H.getIA();
    FAE -= H.getAB();
    FMI -= H.getIJ();

    const SpinorbitalTensor<U>& WMNEF = Hbar.getIJAB();
    const SpinorbitalTensor<U>& WAMEF = Hbar.getAIBC();
    const SpinorbitalTensor<U>& WABEJ = Hbar.getABCI();
    const SpinorbitalTensor<U>& WABEF = Hbar.getABCD();
    const SpinorbitalTensor<U>& WMNIJ = Hbar.getIJKL();
    const SpinorbitalTensor<U>& WMNEJ = Hbar.getIJAK();
    const SpinorbitalTensor<U>& WAMIJ = Hbar.getAIJK();
    const SpinorbitalTensor<U>& WAMEI = Hbar.getAIBJ();

    SpinorbitalTensor<U> WABCEJK("W(abc,ejk)", arena, group, {vrt,occ}, {3,0}, {1,2});
    SpinorbitalTensor<U> WABMIJK("W(abm,ijk)", arena, group, {vrt,occ}, {2,1}, {0,3});
    SpinorbitalTensor<U> WAMNIJK("W(amn,ijk)", arena, group, {vrt,occ}, {1,2}, {0,3});
    SpinorbitalTensor<U> WABMEJI("W(abm,ejk)", arena, group, {vrt,occ}, {2,1}, {1,2});

    ExcitationOperator<U,4> T_1("T^(1)", arena, occ, vrt);
    ExcitationOperator<U,4> T_2("T^(2)", arena, occ, vrt);
    ExcitationOperator<U,4> T_3("T^(3)", arena, occ, vrt);
    ExcitationOperator<U,4> T_4("T^(4)", arena, occ, vrt);

    DeexcitationOperator<U,4> L_1("L^(1)", arena, occ, vrt);
    DeexcitationOperator<U,4> L_2("L^(2)", arena, occ, vrt);

    ExcitationOperator<U,4> Z("Z", arena, occ, vrt);

    SpinorbitalTensor<U> WABCEJK_1(WABCEJK);
    SpinorbitalTensor<U> WABMIJK_1(WABMIJK);

    SpinorbitalTensor<U> FMI_2(FMI);
    SpinorbitalTensor<U> FAE_2(FAE);
    SpinorbitalTensor<U> WAMEI_2(WAMEI);
    SpinorbitalTensor<U> WMNIJ_2(WMNIJ);
    SpinorbitalTensor<U> WABEF_2(WABEF);
    SpinorbitalTensor<U> WABEJ_2(WABEJ);
    SpinorbitalTensor<U> WAMIJ_2(WAMIJ);
    SpinorbitalTensor<U> WABCEJK_2(WABCEJK);
    SpinorbitalTensor<U> WABMIJK_2(WABMIJK);
    SpinorbitalTensor<U> WAMNIJK_2(WAMNIJK);
    SpinorbitalTensor<U> WABMEJI_2(WABMEJI);

    SpinorbitalTensor<U> FME_3(FME);
    SpinorbitalTensor<U> FMI_3(FMI);
    SpinorbitalTensor<U> FAE_3(FAE);
    SpinorbitalTensor<U> WAMEI_3(WAMEI);
    SpinorbitalTensor<U> WMNIJ_3(WMNIJ);
    SpinorbitalTensor<U> WABEF_3(WABEF);
    SpinorbitalTensor<U> WABEJ_3(WABEJ);
    SpinorbitalTensor<U> WAMIJ_3(WAMIJ);
    SpinorbitalTensor<U> WAMEF_3(WAMEF);
    SpinorbitalTensor<U> WMNEJ_3(WMNEJ);
    SpinorbitalTensor<U> WABMEJI_3(WABMEJI);
    SpinorbitalTensor<U> WABCEJK_3(WABCEJK);
    SpinorbitalTensor<U> WABMIJK_3(WABMIJK);
    SpinorbitalTensor<U> WAMNIJK_3(WAMNIJK);

    SpinorbitalTensor<U> WABEJ_5(WABEJ);
    SpinorbitalTensor<U> WAMIJ_5(WAMIJ);

    SpinorbitalTensor<U> DAI_1(T(1));
    SpinorbitalTensor<U> GIJAK_1(WMNEJ);
    SpinorbitalTensor<U> GAIBC_1(WAMEF);
    SpinorbitalTensor<U> GAIJK_1(WAMIJ);
    SpinorbitalTensor<U> GABCI_1(WABEJ);
    SpinorbitalTensor<U> GABCD_1(WABEF);
    SpinorbitalTensor<U> GAIBJ_1(WAMEI);
    SpinorbitalTensor<U> GIJKL_1(WMNIJ);
    SpinorbitalTensor<U> GAIJBCD_1("G(aij,bcd)", arena, group, {vrt,occ}, {1,2}, {3,0});
    SpinorbitalTensor<U> GIJKABL_1("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1});
    SpinorbitalTensor<U> GIJKALM_1("G(ijk,alm)", arena, group, {vrt,occ}, {0,3}, {1,2});

    /***************************************************************************
     *
     * L^(1)
     *
     **************************************************************************/

    L_1(4)["ijklabcd"]  = WMNEF["ijab"]*L(2)[  "klcd"];
    L_1(4)["ijklabcd"] +=   FME[  "ia"]*L(3)["jklbcd"];
    L_1(4)["ijklabcd"] += WAMEF["elcd"]*L(3)["ijkabe"];
    L_1(4)["ijklabcd"] -= WMNEJ["ijam"]*L(3)["mklbcd"];

    L_1.weight(D);

    /***************************************************************************
     *
     * T^(1)
     *
     **************************************************************************/

    WABCEJK[  "abcejk"]  =   0.5*WABEF[  "abef"]*T(2)[    "fcjk"];

    WABMIJK[  "abmijk"]  =       WAMEI[  "amek"]*T(2)[    "ebij"];
    WABMIJK[  "abmijk"] -=   0.5*WMNIJ[  "nmjk"]*T(2)[    "abin"];
    WABMIJK[  "abmijk"] +=   0.5*WAMEF[  "bmef"]*T(3)[  "aefijk"];

    WAMNIJK[  "amnijk"]  =       WMNEJ[  "mnek"]*T(2)[    "aeij"];
    WAMNIJK[  "amnijk"] +=   0.5*WMNEF[  "mnef"]*T(3)[  "aefijk"];

    WABMEJI[  "abmeji"]  =       WAMEF[  "amef"]*T(2)[    "bfji"];
    WABMEJI[  "abmeji"] -=       WMNEJ[  "nmei"]*T(2)[    "abnj"];
    WABMEJI[  "abmeji"] +=   0.5*WMNEF[  "mnef"]*T(3)[  "abfnji"];

     T_1(4)["abcdijkl"]  =     WABCEJK["abcejk"]*T(2)[    "edil"];
     T_1(4)["abcdijkl"] -=     WABMIJK["abmijk"]*T(2)[    "cdml"];
     T_1(4)["abcdijkl"] +=       WABEJ[  "abej"]*T(3)[  "ecdikl"];
     T_1(4)["abcdijkl"] -=       WAMIJ[  "amij"]*T(3)[  "bcdmkl"];
     T_1(4)["abcdijkl"] += 0.5*WAMNIJK["amnijk"]*T(3)[  "bcdmnl"];
     T_1(4)["abcdijkl"] -=     WABMEJI["abmeji"]*T(3)[  "ecdmkl"];

    U E2 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_1(4)["efghmnop"]);

    T_1.weight(D);

    /***************************************************************************
     *
     * L^(2)
     *
     **************************************************************************/

      WABCEJK[  "abcejk"]  =               WABEF[    "abef"]*     T(2)[    "fcjk"];
      WABCEJK[  "abcejk"] -=               WAMEI[    "amek"]*     T(2)[    "bcjm"];
      WABCEJK[  "abcejk"] -=                 FME[      "me"]*     T(3)[  "abcmjk"];
      WABCEJK[  "abcejk"] +=  (1.0/ 2.0)*  WMNEJ[    "mnej"]*     T(3)[  "abcmnk"];
      WABCEJK[  "abcejk"] +=               WAMEF[    "amef"]*     T(3)[  "fbcmjk"];

      WABMIJK[  "abmijk"]  =               WAMEI[    "amek"]*     T(2)[    "ebij"];
      WABMIJK[  "abmijk"] -=               WMNIJ[    "nmjk"]*     T(2)[    "abin"];
      WABMIJK[  "abmijk"] +=                 FME[      "me"]*     T(3)[  "abeijk"];
      WABMIJK[  "abmijk"] +=  (1.0/ 2.0)*  WAMEF[    "bmef"]*     T(3)[  "aefijk"];
      WABMIJK[  "abmijk"] +=               WMNEJ[    "nmek"]*     T(3)[  "abeijn"];

      WAMNIJK[  "amnijk"]  =               WMNEJ[    "mnek"]*     T(2)[    "aeij"];
      WAMNIJK[  "amnijk"] +=  (1.0/ 2.0)*  WMNEF[    "mnef"]*     T(3)[  "aefijk"];

      WABMEJI[  "abmeji"]  =               WAMEF[    "amef"]*     T(2)[    "bfji"];
      WABMEJI[  "abmeji"] -=               WMNEJ[    "nmei"]*     T(2)[    "abnj"];
      WABMEJI[  "abmeji"] +=               WMNEF[    "mnef"]*     T(3)[  "abfnji"];

    GIJKABL_1[  "ijkabl"]  =  (1.0/ 2.0)*   T(2)[    "eflm"]*   L_1(4)["ijkmabef"];
    GAIJBCD_1[  "aijbcd"]  = -(1.0/ 2.0)*   T(2)[    "aemn"]*   L_1(4)["mijnbcde"];
    GIJKALM_1[  "ijkalm"]  =  (1.0/ 6.0)*   T(3)[  "efglmn"]*   L_1(4)["ijknaefg"];

      GIJAK_1[    "ijak"]  =  (1.0/12.0)*   T(3)[  "efgkmn"]*   L_1(4)["ijmnaefg"];
      GAIBC_1[    "aibc"]  = -(1.0/12.0)*   T(3)[  "aefmno"]*   L_1(4)["minobcef"];

      GAIJK_1[    "aijk"]  =  (1.0/12.0)* T_1(4)["aefgjkmn"]*     L(3)[  "imnefg"];
      GAIJK_1[    "aijk"] +=  (1.0/ 4.0)*   T(3)[  "efamnj"]*GIJKABL_1[  "mniefk"];
      GAIJK_1[    "aijk"] +=  (1.0/ 6.0)*   T(3)[  "efgjkm"]*GAIJBCD_1[  "aimefg"];

      GABCI_1[    "abci"]  = -(1.0/12.0)* T_1(4)["abefmino"]*     L(3)[  "mnocef"];
      GABCI_1[    "abci"] +=  (1.0/ 4.0)*   T(3)[  "efbmni"]*GAIJBCD_1[  "amncef"];
      GABCI_1[    "abci"] +=  (1.0/ 6.0)*   T(3)[  "eabmno"]*GIJKABL_1[  "mnoeci"];

      GABCD_1[    "abcd"]  = -(1.0/ 4.0)*   T(2)[    "bemn"]*GAIJBCD_1[  "amncde"];
      GAIBJ_1[    "aibj"]  = -(1.0/ 2.0)*   T(2)[    "aenm"]*GIJKABL_1[  "mniebj"];
      GIJKL_1[    "ijkl"]  =  (1.0/ 4.0)*   T(2)[    "efmk"]*GIJKABL_1[  "mijefl"];

        DAI_1[      "ai"]  =  (1.0/36.0)* T_1(4)["aefgimno"]*     L(3)[  "mnoefg"];
        DAI_1[      "ai"] -=  (1.0/ 2.0)*   T(2)[    "eamn"]*  GIJAK_1[    "mnei"];
        DAI_1[      "ai"] +=  (1.0/ 2.0)*   T(2)[    "efim"]*  GAIBC_1[    "amef"];

       L_2(1)[      "ia"]  =               WMNEF[    "miea"]*    DAI_1[      "em"];
       L_2(1)[      "ia"] +=  (1.0/ 2.0)*  WMNEF[    "imef"]*  GABCI_1[    "efam"];
       L_2(1)[      "ia"] -=  (1.0/ 2.0)*  WMNEF[    "mnea"]*  GAIJK_1[    "eimn"];
       L_2(1)[      "ia"] -=  (1.0/ 2.0)*  WABEF[    "efga"]*  GAIBC_1[    "gief"];
       L_2(1)[      "ia"] +=               WAMEI[    "eifm"]*  GAIBC_1[    "fmea"];
       L_2(1)[      "ia"] -=               WAMEI[    "eman"]*  GIJAK_1[    "inem"];
       L_2(1)[      "ia"] +=  (1.0/ 2.0)*  WMNIJ[    "imno"]*  GIJAK_1[    "noam"];
       L_2(1)[      "ia"] -=  (1.0/ 2.0)*  WAMEF[    "gief"]*  GABCD_1[    "efga"];
       L_2(1)[      "ia"] +=               WAMEF[    "fmea"]*  GAIBJ_1[    "eifm"];
       L_2(1)[      "ia"] -=               WMNEJ[    "inem"]*  GAIBJ_1[    "eman"];
       L_2(1)[      "ia"] +=  (1.0/ 2.0)*  WMNEJ[    "noam"]*  GIJKL_1[    "imno"];

       L_2(2)[    "ijab"]  =              -WAMEF[    "fiae"]*  GAIBC_1[    "ejbf"];
       L_2(2)[    "ijab"] -=               WMNEJ[    "ijem"]*  GAIBC_1[    "emab"];
       L_2(2)[    "ijab"] -=               WAMEF[    "emab"]*  GIJAK_1[    "ijem"];
       L_2(2)[    "ijab"] -=               WMNEJ[    "niam"]*  GIJAK_1[    "mjbn"];
       L_2(2)[    "ijab"] +=  (1.0/ 2.0)*  WMNEF[    "ijef"]*  GABCD_1[    "efab"];
       L_2(2)[    "ijab"] +=               WMNEF[    "imea"]*  GAIBJ_1[    "ejbm"];
       L_2(2)[    "ijab"] +=  (1.0/ 2.0)*  WMNEF[    "mnab"]*  GIJKL_1[    "ijmn"];
       L_2(2)[    "ijab"] +=  (1.0/12.0)*WABCEJK[  "efgbmn"]*   L_1(4)["ijmnaefg"];
       L_2(2)[    "ijab"] -=  (1.0/12.0)*WABMIJK[  "efimno"]*   L_1(4)["mnojefab"];

       L_2(3)[  "ijkabc"]  =               WMNEF[    "ijae"]*  GAIBC_1[    "ekbc"];
       L_2(3)[  "ijkabc"] -=               WMNEF[    "mkbc"]*  GIJAK_1[    "ijam"];
       L_2(3)[  "ijkabc"] +=  (1.0/ 2.0)*  WMNEF[    "mnbc"]*GIJKALM_1[  "ijkamn"];
       L_2(3)[  "ijkabc"] -=               WAMEF[    "embc"]*GIJKABL_1[  "ijkaem"];
       L_2(3)[  "ijkabc"] +=  (1.0/ 2.0)*  WABEJ[    "efcm"]*   L_1(4)["ijkmabef"];
       L_2(3)[  "ijkabc"] -=  (1.0/ 2.0)*  WAMIJ[    "eknm"]*   L_1(4)["ijmnabce"];
       L_2(3)[  "ijkabc"] -=  (1.0/ 4.0)*WABMEJI[  "efkcnm"]*   L_1(4)["ijmnabef"];
       L_2(3)[  "ijkabc"] +=  (1.0/ 6.0)*WAMNIJK[  "eijmno"]*   L_1(4)["mnokeabc"];

       L_2(4)["ijklabcd"]  =               WMNEF[    "ijae"]*GAIJBCD_1[  "eklbcd"];
       L_2(4)["ijklabcd"] -=               WMNEF[    "mlcd"]*GIJKABL_1[  "ijkabm"];
       L_2(4)["ijklabcd"] +=                 FAE[      "ea"]*   L_1(4)["ijklebcd"];
       L_2(4)["ijklabcd"] -=                 FMI[      "im"]*   L_1(4)["mjklabcd"];
       L_2(4)["ijklabcd"] +=  (1.0/ 2.0)*  WABEF[    "efab"]*   L_1(4)["ijklefcd"];
       L_2(4)["ijklabcd"] +=  (1.0/ 2.0)*  WMNIJ[    "ijmn"]*   L_1(4)["mnklabcd"];
       L_2(4)["ijklabcd"] +=               WAMEI[    "eiam"]*   L_1(4)["mjklbecd"];

    L_2.weight(D);

    /***************************************************************************
     *
     * T^(2)
     *
     **************************************************************************/

    WABCEJK_1[  "abcejk"]  =  -0.5*WMNEF[  "mnef"]*T_1(4)["abcfmjkn"];
    WABMIJK_1[  "abmijk"]  =   0.5*WMNEF[  "mnef"]*T_1(4)["abefijkn"];

       T_2(2)[    "abij"]  =  0.25*WMNEF[  "mnef"]*T_1(4)["abefijmn"];

       T_2(3)[  "abcijk"]  =         FME[    "me"]*T_1(4)["abceijkm"];
       T_2(3)[  "abcijk"] +=   0.5*WAMEF[  "amef"]*T_1(4)["efbcimjk"];
       T_2(3)[  "abcijk"] -=   0.5*WMNEJ[  "mnek"]*T_1(4)["abecijmn"];

       T_2(4)["abcdijkl"]  =   WABCEJK_1["abcejk"]*  T(2)[    "edil"];
       T_2(4)["abcdijkl"] -=   WABMIJK_1["abmijk"]*  T(2)[    "cdml"];
       T_2(4)["abcdijkl"] +=         FAE[    "ae"]*T_1(4)["ebcdijkl"];
       T_2(4)["abcdijkl"] -=         FMI[    "mi"]*T_1(4)["abcdmjkl"];
       T_2(4)["abcdijkl"] +=   0.5*WABEF[  "abef"]*T_1(4)["efcdijkl"];
       T_2(4)["abcdijkl"] +=   0.5*WMNIJ[  "mnij"]*T_1(4)["abcdmnkl"];
       T_2(4)["abcdijkl"] +=       WAMEI[  "amei"]*T_1(4)["ebcdjmkl"];

    U E3 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_2(4)["efghmnop"]);

    T_2.weight(D);

    /***************************************************************************
     *
     * T^(3)
     *
     **************************************************************************/

        FMI_2[      "mi"]  =     0.5*WMNEF[  "mnef"]*T_2(2)[    "efin"];

        FAE_2[      "ae"]  =    -0.5*WMNEF[  "mnef"]*T_2(2)[    "afmn"];

      WAMIJ_2[    "amij"]  =         WMNEJ[  "nmej"]*T_2(2)[    "aein"];
      WAMIJ_2[    "amij"] +=     0.5*WAMEF[  "amef"]*T_2(2)[    "efij"];
      WAMIJ_2[    "amij"] +=     0.5*WMNEF[  "mnef"]*T_2(3)[  "aefijn"];

      WABEJ_2[    "abej"]  =         WAMEF[  "amef"]*T_2(2)[    "fbmj"];
      WABEJ_2[    "abej"] +=     0.5*WMNEJ[  "mnej"]*T_2(2)[    "abmn"];
      WABEJ_2[    "abej"] -=     0.5*WMNEF[  "mnef"]*T_2(3)[  "afbmnj"];

      WABEF_2[    "abef"]  =     0.5*WMNEF[  "mnef"]*T_2(2)[    "abmn"];

      WAMEI_2[    "amei"]  =        -WMNEF[  "mnef"]*T_2(2)[    "afin"];

      WMNIJ_2[    "mnij"]  =     0.5*WMNEF[  "mnef"]*T_2(2)[    "efij"];

    WABCEJK_2[  "abcejk"]  =   0.5*WABEF_2[  "abef"]*  T(2)[    "fcjk"];

    WABMIJK_2[  "abmijk"]  =       WAMEI_2[  "amek"]*  T(2)[    "ebij"];
    WABMIJK_2[  "abmijk"] -=   0.5*WMNIJ_2[  "nmjk"]*  T(2)[    "abin"];
    WABMIJK_2[  "abmijk"] +=     0.5*WAMEF[  "bmef"]*T_2(3)[  "aefijk"];

    WAMNIJK_2[  "amnijk"]  =     0.5*WMNEF[  "mnef"]*T_2(3)[  "aefijk"];

       T_3(1)[      "ai"]  =           FME[    "me"]*T_2(2)[    "aeim"];
       T_3(1)[      "ai"] +=     0.5*WAMEF[  "amef"]*T_2(2)[    "efim"];
       T_3(1)[      "ai"] -=     0.5*WMNEJ[  "mnei"]*T_2(2)[    "eamn"];
       T_3(1)[      "ai"] +=    0.25*WMNEF[  "mnef"]*T_2(3)[  "aefimn"];

       T_3(2)[    "abij"]  =         FAE_2[    "af"]*  T(2)[    "fbij"];
       T_3(2)[    "abij"] -=         FMI_2[    "ni"]*  T(2)[    "abnj"];
       T_3(2)[    "abij"] +=           FAE[    "af"]*T_2(2)[    "fbij"];
       T_3(2)[    "abij"] -=           FMI[    "ni"]*T_2(2)[    "abnj"];
       T_3(2)[    "abij"] +=     0.5*WABEF[  "abef"]*T_2(2)[    "efij"];
       T_3(2)[    "abij"] +=     0.5*WMNIJ[  "mnij"]*T_2(2)[    "abmn"];
       T_3(2)[    "abij"] +=         WAMEI[  "amei"]*T_2(2)[    "ebjm"];
       T_3(2)[    "abij"] +=           FME[    "me"]*T_2(3)[  "abeijm"];
       T_3(2)[    "abij"] +=     0.5*WAMEF[  "bmef"]*T_2(3)[  "aefijm"];
       T_3(2)[    "abij"] -=     0.5*WMNEJ[  "mnej"]*T_2(3)[  "abeinm"];

       T_3(3)[  "abcijk"]  =       WABEJ_2[  "bcek"]*  T(2)[    "aeij"];
       T_3(3)[  "abcijk"] -=       WAMIJ_2[  "bmjk"]*  T(2)[    "acim"];
       T_3(3)[  "abcijk"] +=         FAE_2[    "ce"]*  T(3)[  "abeijk"];
       T_3(3)[  "abcijk"] -=         FMI_2[    "mk"]*  T(3)[  "abcijm"];
       T_3(3)[  "abcijk"] +=   0.5*WABEF_2[  "abef"]*  T(3)[  "efcijk"];
       T_3(3)[  "abcijk"] +=   0.5*WMNIJ_2[  "mnij"]*  T(3)[  "abcmnk"];
       T_3(3)[  "abcijk"] +=       WAMEI_2[  "amei"]*  T(3)[  "ebcjmk"];
       T_3(3)[  "abcijk"] +=         WABEJ[  "bcek"]*T_2(2)[    "aeij"];
       T_3(3)[  "abcijk"] -=         WAMIJ[  "bmjk"]*T_2(2)[    "acim"];
       T_3(3)[  "abcijk"] +=           FAE[    "ce"]*T_2(3)[  "abeijk"];
       T_3(3)[  "abcijk"] -=           FMI[    "mk"]*T_2(3)[  "abcijm"];
       T_3(3)[  "abcijk"] +=     0.5*WABEF[  "abef"]*T_2(3)[  "efcijk"];
       T_3(3)[  "abcijk"] +=     0.5*WMNIJ[  "mnij"]*T_2(3)[  "abcmnk"];
       T_3(3)[  "abcijk"] +=         WAMEI[  "amei"]*T_2(3)[  "ebcjmk"];

       T_3(4)["abcdijkl"]  =     WABCEJK_2["abcejk"]*  T(2)[    "edil"];
       T_3(4)["abcdijkl"] -=     WABMIJK_2["abmijk"]*  T(2)[    "cdml"];
       T_3(4)["abcdijkl"] +=       WABEJ_2[  "abej"]*  T(3)[  "ecdikl"];
       T_3(4)["abcdijkl"] -=       WAMIJ_2[  "amij"]*  T(3)[  "bcdmkl"];
       T_3(4)["abcdijkl"] += 0.5*WAMNIJK_2["amnijk"]*  T(3)[  "bcdmnl"];
       T_3(4)["abcdijkl"] +=       WABCEJK["abcejk"]*T_2(2)[    "edil"];
       T_3(4)["abcdijkl"] -=       WABMIJK["abmijk"]*T_2(2)[    "cdml"];
       T_3(4)["abcdijkl"] +=         WABEJ[  "abej"]*T_2(3)[  "ecdikl"];
       T_3(4)["abcdijkl"] -=         WAMIJ[  "amij"]*T_2(3)[  "bcdmkl"];
       T_3(4)["abcdijkl"] +=   0.5*WAMNIJK["amnijk"]*T_2(3)[  "bcdmnl"];
       T_3(4)["abcdijkl"] -=       WABMEJI["abmeji"]*T_2(3)[  "ecdmkl"];

    U E4_23 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_3(4)["efghmnop"]);

    U E5_1_23 = (1.0/  1.0)*scalar(L_2(1)[      "me"]*T_3(1)[      "em"]);
    U E5_2_23 = (1.0/  4.0)*scalar(L_2(2)[    "mnef"]*T_3(2)[    "efmn"]);
    U E5_3_23 = (1.0/ 36.0)*scalar(L_2(3)[  "mnoefg"]*T_3(3)[  "efgmno"]);
    U E5_4_23 = (1.0/576.0)*scalar(L_2(4)["mnopefgh"]*T_3(4)["efghmnop"]);

    WABCEJK_2[  "abcejk"]  =    -0.5*WMNEF[  "mnef"]*T_2(4)["abcfmjkn"];
    WABMIJK_2[  "abmijk"]  =     0.5*WMNEF[  "mnef"]*T_2(4)["abefijkn"];

       T_3(2)[    "abij"] +=    0.25*WMNEF[  "mnef"]*T_2(4)["abefijmn"];

       T_3(3)[  "abcijk"] +=           FME[    "me"]*T_2(4)["abceijkm"];
       T_3(3)[  "abcijk"] +=     0.5*WAMEF[  "amef"]*T_2(4)["efbcimjk"];
       T_3(3)[  "abcijk"] -=     0.5*WMNEJ[  "mnek"]*T_2(4)["abecijmn"];

       T_3(4)["abcdijkl"] +=     WABCEJK_2["abcejk"]*  T(2)[    "edil"];
       T_3(4)["abcdijkl"] -=     WABMIJK_2["abmijk"]*  T(2)[    "cdml"];
       T_3(4)["abcdijkl"] +=           FAE[    "ae"]*T_2(4)["ebcdijkl"];
       T_3(4)["abcdijkl"] -=           FMI[    "mi"]*T_2(4)["abcdmjkl"];
       T_3(4)["abcdijkl"] +=     0.5*WABEF[  "abef"]*T_2(4)["efcdijkl"];
       T_3(4)["abcdijkl"] +=     0.5*WMNIJ[  "mnij"]*T_2(4)["abcdmnkl"];
       T_3(4)["abcdijkl"] +=         WAMEI[  "amei"]*T_2(4)["ebcdjmkl"];

      WAMIJ_2[    "amij"] +=           FME[    "me"]*T_2(2)[    "aeij"];
      WABEJ_2[    "abej"] -=           FME[    "me"]*T_2(2)[    "abmj"];

    U E4_4 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_3(4)["efghmnop"])-E4_23;

    U E5_1_4 = 0.0;
    U E5_2_4 = (1.0/  4.0)*scalar(L_2(2)[    "mnef"]*T_3(2)[    "efmn"])-E5_2_23;
    U E5_3_4 = (1.0/ 36.0)*scalar(L_2(3)[  "mnoefg"]*T_3(3)[  "efgmno"])-E5_3_23;
    U E5_4_4 = (1.0/576.0)*scalar(L_2(4)["mnopefgh"]*T_3(4)["efghmnop"])-E5_4_23;

    T_3.weight(D);

    /***************************************************************************
     *
     * T^(4)
     *
     **************************************************************************/

        FME_3[      "me"]  =         WMNEF[  "mnef"]*T_3(1)[      "fn"];

      WMNEJ_3[    "mnej"]  =         WMNEF[  "mnef"]*T_3(1)[      "fj"];

      WAMEF_3[    "amef"]  =        -WMNEF[  "nmef"]*T_3(1)[      "an"];

        FMI_3[      "mi"]  =     0.5*WMNEF[  "mnef"]*T_3(2)[    "efin"];
        FMI_3[      "mi"] +=         WMNEJ[  "nmfi"]*T_3(1)[      "fn"];

        FAE_3[      "ae"]  =    -0.5*WMNEF[  "mnef"]*T_3(2)[    "afmn"];
        FAE_3[      "ae"] +=         WAMEF[  "amef"]*T_3(1)[      "fm"];

      WAMIJ_3[    "amij"]  =         WMNEJ[  "nmej"]*T_3(2)[    "aein"];
      WAMIJ_3[    "amij"] +=     0.5*WAMEF[  "amef"]*T_3(2)[    "efij"];
      WAMIJ_3[    "amij"] +=     0.5*WMNEF[  "mnef"]*T_3(3)[  "aefijn"];
      WAMIJ_3[    "amij"] +=         WAMEI[  "amej"]*T_3(1)[      "ei"];
      WAMIJ_3[    "amij"] -=         WMNIJ[  "nmij"]*T_3(1)[      "an"];
      WAMIJ_3[    "amij"] +=         FME_3[    "me"]*  T(2)[    "aeij"];

      WABEJ_3[    "abej"]  =         WAMEF[  "amef"]*T_3(2)[    "fbmj"];
      WABEJ_3[    "abej"] +=     0.5*WMNEJ[  "mnej"]*T_3(2)[    "abmn"];
      WABEJ_3[    "abej"] -=     0.5*WMNEF[  "mnef"]*T_3(3)[  "afbmnj"];
      WABEJ_3[    "abej"] -=         WAMEI[  "amej"]*T_3(1)[      "bm"];
      WABEJ_3[    "abej"] +=         WABEF[  "abef"]*T_3(1)[      "fj"];

      WABEF_3[    "abef"]  =     0.5*WMNEF[  "mnef"]*T_3(2)[    "abmn"];
      WABEF_3[    "abef"] -=         WAMEF[  "amef"]*T_3(1)[      "bm"];

      WAMEI_3[    "amei"]  =        -WMNEF[  "mnef"]*T_3(2)[    "afin"];
      WAMEI_3[    "amei"] +=         WAMEF[  "amef"]*T_3(1)[      "fi"];
      WAMEI_3[    "amei"] -=         WMNEJ[  "nmei"]*T_3(1)[      "an"];

      WMNIJ_3[    "mnij"]  =     0.5*WMNEF[  "mnef"]*T_3(2)[    "efij"];
      WMNIJ_3[    "mnij"] +=         WMNEJ[  "mnej"]*T_3(1)[      "ei"];

    WABCEJK_3[  "abcejk"]  =   0.5*WABEF_3[  "abef"]*  T(2)[    "fcjk"];

    WABMIJK_3[  "abmijk"]  =       WAMEI_3[  "amek"]*  T(2)[    "ebij"];
    WABMIJK_3[  "abmijk"] -=   0.5*WMNIJ_3[  "nmjk"]*  T(2)[    "abin"];
    WABMIJK_3[  "abmijk"] +=   0.5*WAMEF_3[  "bmef"]*  T(3)[  "aefijk"];
    WABMIJK_3[  "abmijk"] +=     0.5*WAMEF[  "bmef"]*T_3(3)[  "aefijk"];

    WAMNIJK_3[  "amnijk"]  =       WMNEJ_3[  "mnek"]*  T(2)[    "aeij"];
    WAMNIJK_3[  "amnijk"] +=     0.5*WMNEF[  "mnef"]*T_3(3)[  "aefijk"];

    WABMEJI_3[  "abmeji"]  =       WAMEF_3[  "amef"]*  T(2)[    "bfji"];
    WABMEJI_3[  "abmeji"] -=       WMNEJ_3[  "nmei"]*  T(2)[    "abnj"];

       T_4(1)[      "ai"]  =           FAE[    "ae"]*T_3(1)[      "ei"];
       T_4(1)[      "ai"] -=           FMI[    "mi"]*T_3(1)[      "am"];
       T_4(1)[      "ai"] -=         WAMEI[  "amei"]*T_3(1)[      "em"];
       T_4(1)[      "ai"] +=           FME[    "me"]*T_3(2)[    "aeim"];
       T_4(1)[      "ai"] +=     0.5*WAMEF[  "amef"]*T_3(2)[    "efim"];
       T_4(1)[      "ai"] -=     0.5*WMNEJ[  "mnei"]*T_3(2)[    "eamn"];
       T_4(1)[      "ai"] +=    0.25*WMNEF[  "mnef"]*T_3(3)[  "aefimn"];

       T_4(2)[    "abij"]  =         FAE_3[    "af"]*  T(2)[    "fbij"];
       T_4(2)[    "abij"] -=         FMI_3[    "ni"]*  T(2)[    "abnj"];
       T_4(2)[    "abij"] +=         FME_3[    "me"]*  T(3)[  "abeijm"];
       T_4(2)[    "abij"] +=         WABEJ[  "abej"]*T_3(1)[      "ei"];
       T_4(2)[    "abij"] -=         WAMIJ[  "amij"]*T_3(1)[      "bm"];
       T_4(2)[    "abij"] +=           FAE[    "af"]*T_3(2)[    "fbij"];
       T_4(2)[    "abij"] -=           FMI[    "ni"]*T_3(2)[    "abnj"];
       T_4(2)[    "abij"] +=     0.5*WABEF[  "abef"]*T_3(2)[    "efij"];
       T_4(2)[    "abij"] +=     0.5*WMNIJ[  "mnij"]*T_3(2)[    "abmn"];
       T_4(2)[    "abij"] +=         WAMEI[  "amei"]*T_3(2)[    "ebjm"];
       T_4(2)[    "abij"] +=           FME[    "me"]*T_3(3)[  "abeijm"];
       T_4(2)[    "abij"] +=     0.5*WAMEF[  "bmef"]*T_3(3)[  "aefijm"];
       T_4(2)[    "abij"] -=     0.5*WMNEJ[  "mnej"]*T_3(3)[  "abeinm"];

        FMI_3[      "mi"] +=           FME[    "me"]*T_3(1)[      "ei"];
        FAE_3[      "ae"] -=           FME[    "me"]*T_3(1)[      "am"];

       T_4(3)[  "abcijk"]  =       WABEJ_3[  "bcek"]*  T(2)[    "aeij"];
       T_4(3)[  "abcijk"] -=       WAMIJ_3[  "bmjk"]*  T(2)[    "acim"];
       T_4(3)[  "abcijk"] +=         FAE_3[    "ce"]*  T(3)[  "abeijk"];
       T_4(3)[  "abcijk"] -=         FMI_3[    "mk"]*  T(3)[  "abcijm"];
       T_4(3)[  "abcijk"] +=   0.5*WABEF_3[  "abef"]*  T(3)[  "efcijk"];
       T_4(3)[  "abcijk"] +=   0.5*WMNIJ_3[  "mnij"]*  T(3)[  "abcmnk"];
       T_4(3)[  "abcijk"] +=       WAMEI_3[  "amei"]*  T(3)[  "ebcjmk"];
       T_4(3)[  "abcijk"] +=         WABEJ[  "bcek"]*T_3(2)[    "aeij"];
       T_4(3)[  "abcijk"] -=         WAMIJ[  "bmjk"]*T_3(2)[    "acim"];
       T_4(3)[  "abcijk"] +=           FAE[    "ce"]*T_3(3)[  "abeijk"];
       T_4(3)[  "abcijk"] -=           FMI[    "mk"]*T_3(3)[  "abcijm"];
       T_4(3)[  "abcijk"] +=     0.5*WABEF[  "abef"]*T_3(3)[  "efcijk"];
       T_4(3)[  "abcijk"] +=     0.5*WMNIJ[  "mnij"]*T_3(3)[  "abcmnk"];
       T_4(3)[  "abcijk"] +=         WAMEI[  "amei"]*T_3(3)[  "ebcjmk"];

      WABEJ_3[    "abej"] -=         FME_3[    "me"]*  T(2)[    "abmj"];

       T_4(4)["abcdijkl"] +=     WABCEJK_3["abcejk"]*  T(2)[    "edil"];
       T_4(4)["abcdijkl"] -=     WABMIJK_3["abmijk"]*  T(2)[    "cdml"];
       T_4(4)["abcdijkl"] +=       WABEJ_3[  "abej"]*  T(3)[  "ecdikl"];
       T_4(4)["abcdijkl"] -=       WAMIJ_3[  "amij"]*  T(3)[  "bcdmkl"];
       T_4(4)["abcdijkl"] += 0.5*WAMNIJK_3["amnijk"]*  T(3)[  "bcdmnl"];
       T_4(4)["abcdijkl"] -=     WABMEJI_3["abmeji"]*  T(3)[  "ecdmkl"];
       T_4(4)["abcdijkl"] +=       WABCEJK["abcejk"]*T_3(2)[    "edil"];
       T_4(4)["abcdijkl"] -=       WABMIJK["abmijk"]*T_3(2)[    "cdml"];
       T_4(4)["abcdijkl"] +=         WABEJ[  "abej"]*T_3(3)[  "ecdikl"];
       T_4(4)["abcdijkl"] -=         WAMIJ[  "amij"]*T_3(3)[  "bcdmkl"];
       T_4(4)["abcdijkl"] +=   0.5*WAMNIJK["amnijk"]*T_3(3)[  "bcdmnl"];
       T_4(4)["abcdijkl"] -=       WABMEJI["abmeji"]*T_3(3)[  "ecdmkl"];

    U E6_1_123 = (1.0/  1.0)*scalar(L_2(1)[      "me"]*T_4(1)[      "em"]);
    U E6_2_123 = (1.0/  4.0)*scalar(L_2(2)[    "mnef"]*T_4(2)[    "efmn"]);
    U E6_3_123 = (1.0/ 36.0)*scalar(L_2(3)[  "mnoefg"]*T_4(3)[  "efgmno"]);
    U E6_4_123 = (1.0/576.0)*scalar(L_2(4)["mnopefgh"]*T_4(4)["efghmnop"]);

    WABCEJK_3[  "abcejk"]  =    -0.5*WMNEF[  "mnef"]*T_3(4)["abcfmjkn"];
    WABMIJK_3[  "abmijk"]  =     0.5*WMNEF[  "mnef"]*T_3(4)["abefijkn"];

       T_4(2)[    "abij"]  =    0.25*WMNEF[  "mnef"]*T_3(4)["abefijmn"];

       T_4(3)[  "abcijk"]  =           FME[    "me"]*T_3(4)["abceijkm"];
       T_4(3)[  "abcijk"] +=     0.5*WAMEF[  "amef"]*T_3(4)["efbcimjk"];
       T_4(3)[  "abcijk"] -=     0.5*WMNEJ[  "mnek"]*T_3(4)["abecijmn"];

       T_4(4)["abcdijkl"]  =     WABCEJK_3["abcejk"]*  T(2)[    "edil"];
       T_4(4)["abcdijkl"] -=     WABMIJK_3["abmijk"]*  T(2)[    "cdml"];
       T_4(4)["abcdijkl"] +=           FAE[    "ae"]*T_3(4)["ebcdijkl"];
       T_4(4)["abcdijkl"] -=           FMI[    "mi"]*T_3(4)["abcdmjkl"];
       T_4(4)["abcdijkl"] +=     0.5*WABEF[  "abef"]*T_3(4)["efcdijkl"];
       T_4(4)["abcdijkl"] +=     0.5*WMNIJ[  "mnij"]*T_3(4)["abcdmnkl"];
       T_4(4)["abcdijkl"] +=         WAMEI[  "amei"]*T_3(4)["ebcdjmkl"];

      WAMIJ_3[    "amij"] +=           FME[    "me"]*T_3(2)[    "aeij"];
      WABEJ_3[    "abej"] -=           FME[    "me"]*T_3(2)[    "abmj"];

    U E6_1_4 = 0.0;
    U E6_2_4 = (1.0/  4.0)*scalar(L_2(2)[    "mnef"]*T_4(2)[    "efmn"]);
    U E6_3_4 = (1.0/ 36.0)*scalar(L_2(3)[  "mnoefg"]*T_4(3)[  "efgmno"]);
    U E6_4_4 = (1.0/576.0)*scalar(L_2(4)["mnopefgh"]*T_4(4)["efghmnop"]);

    /***************************************************************************
     *           _
     * <0|L^(1)[[H,T^(1)],T^(2)]|0>
     *           _
     * <0|L^(2)[[H,T^(1)],T^(2)]|0>
     *
     **************************************************************************/

       T_4(4)["abcdijkl"]  =     WABCEJK_1["abcejk"]*T_2(2)[    "edil"];
       T_4(4)["abcdijkl"] -=     WABMIJK_1["abmijk"]*T_2(2)[    "cdml"];
       T_4(4)["abcdijkl"] +=         FAE_2[    "ae"]*T_1(4)["ebcdijkl"];
       T_4(4)["abcdijkl"] -=         FMI_2[    "mi"]*T_1(4)["abcdmjkl"];
       T_4(4)["abcdijkl"] +=   0.5*WABEF_2[  "abef"]*T_1(4)["efcdijkl"];
       T_4(4)["abcdijkl"] +=   0.5*WMNIJ_2[  "mnij"]*T_1(4)["abcdmnkl"];
       T_4(4)["abcdijkl"] +=       WAMEI_2[  "amei"]*T_1(4)["ebcdjmkl"];

    U E1112 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_4(4)["efghmnop"]);
    U E2112 = (1.0/576.0)*scalar(L_2(4)["mnopefgh"]*T_4(4)["efghmnop"]);

    /***************************************************************************
     *           _
     * <0|(1+L)[[H,T^(2)],T^(3)]|0>
     *
     **************************************************************************/

    WAMIJ_5[  "amij"]  =     WAMEI_2["amej"]*T_3(1)[      "ei"];
    WAMIJ_5[  "amij"] -=     WMNIJ_2["nmij"]*T_3(1)[      "an"];

    WABEJ_5[  "abej"]  =    -WAMEI_2["amej"]*T_3(1)[      "bm"];
    WABEJ_5[  "abej"] +=     WABEF_2["abef"]*T_3(1)[      "fj"];

      FMI_3[    "mi"] -=         FME[  "me"]*T_3(1)[      "ei"];
      FAE_3[    "ae"] +=         FME[  "me"]*T_3(1)[      "am"];

       Z(1)[    "ai"]  =       FAE_2[  "ae"]*T_3(1)[      "ei"];
       Z(1)[    "ai"] -=       FMI_2[  "mi"]*T_3(1)[      "am"];
       Z(1)[    "ai"] +=       FME_3[  "me"]*T_2(2)[    "aeim"];

       Z(2)[  "abij"]  =     WABEJ_2["abej"]*T_3(1)[      "ei"];
       Z(2)[  "abij"] -=     WAMIJ_2["amij"]*T_3(1)[      "bm"];
       Z(2)[  "abij"] +=       FAE_3[  "af"]*T_2(2)[    "fbij"];
       Z(2)[  "abij"] -=       FMI_3[  "ni"]*T_2(2)[    "abnj"];
       Z(2)[  "abij"] +=       FAE_2[  "af"]*T_3(2)[    "fbij"];
       Z(2)[  "abij"] -=       FMI_2[  "ni"]*T_3(2)[    "abnj"];
       Z(2)[  "abij"] += 0.5*WABEF_2["abef"]*T_3(2)[    "efij"];
       Z(2)[  "abij"] += 0.5*WMNIJ_2["mnij"]*T_3(2)[    "abmn"];
       Z(2)[  "abij"] +=     WAMEI_2["amei"]*T_3(2)[    "ebjm"];
       Z(2)[  "abij"] +=       FME_3[  "me"]*T_2(3)[  "abeijm"];

      FMI_3[    "mi"] +=         FME[  "me"]*T_3(1)[      "ei"];
      FAE_3[    "ae"] -=         FME[  "me"]*T_3(1)[      "am"];

    WABEJ_2[  "abej"] +=         FME[  "me"]*T_2(2)[    "abmj"];
    WABEJ_3[  "abej"] +=         FME[  "me"]*T_3(2)[    "abmj"];

       Z(3)["abcijk"]  =     WABEJ_5["bcek"]*  T(2)[    "aeij"];
       Z(3)["abcijk"] -=     WAMIJ_5["bmjk"]*  T(2)[    "acim"];
       Z(3)["abcijk"] +=     WABEJ_3["bcek"]*T_2(2)[    "aeij"];
       Z(3)["abcijk"] -=     WAMIJ_3["bmjk"]*T_2(2)[    "acim"];
       Z(3)["abcijk"] +=     WABEJ_2["bcek"]*T_3(2)[    "aeij"];
       Z(3)["abcijk"] -=     WAMIJ_2["bmjk"]*T_3(2)[    "acim"];
       Z(3)["abcijk"] +=       FAE_3[  "ce"]*T_2(3)[  "abeijk"];
       Z(3)["abcijk"] -=       FMI_3[  "mk"]*T_2(3)[  "abcijm"];
       Z(3)["abcijk"] +=       FAE_2[  "ce"]*T_3(3)[  "abeijk"];
       Z(3)["abcijk"] -=       FMI_2[  "mk"]*T_3(3)[  "abcijm"];
       Z(3)["abcijk"] += 0.5*WABEF_3["abef"]*T_2(3)[  "efcijk"];
       Z(3)["abcijk"] += 0.5*WMNIJ_3["mnij"]*T_2(3)[  "abcmnk"];
       Z(3)["abcijk"] +=     WAMEI_3["amei"]*T_2(3)[  "ebcjmk"];
       Z(3)["abcijk"] += 0.5*WABEF_2["abef"]*T_3(3)[  "efcijk"];
       Z(3)["abcijk"] += 0.5*WMNIJ_2["mnij"]*T_3(3)[  "abcmnk"];
       Z(3)["abcijk"] +=     WAMEI_2["amei"]*T_3(3)[  "ebcjmk"];
       Z(3)["abcijk"] +=       FME_3[  "me"]*T_2(4)["abceijkm"];
       Z(3)["abcijk"] += 0.5*WAMEF_3["amef"]*T_2(4)["efbcimjk"];
       Z(3)["abcijk"] -= 0.5*WMNEJ_3["mnek"]*T_2(4)["abecijmn"];

    WABEJ_2[  "abej"] -=         FME[  "me"]*T_2(2)[    "abmj"];
    WABEJ_3[  "abej"] -=         FME[  "me"]*T_3(2)[    "abmj"];

    U E0123 = (1.0/ 1.0)*scalar(L(1)[    "me"]*Z(1)[    "em"])+
              (1.0/ 4.0)*scalar(L(2)[  "mnef"]*Z(2)[  "efmn"])+
              (1.0/36.0)*scalar(L(3)["mnoefg"]*Z(3)["efgmno"]);

    /***************************************************************************
     *           _
     * <0|L^(1)[[H,T^(1)],T^(3)]|0>
     *
     **************************************************************************/

    Z(4)["abcdijkl"]  =   WABCEJK_1["abcejk"]*T_3(2)[    "edil"];
    Z(4)["abcdijkl"] -=   WABMIJK_1["abmijk"]*T_3(2)[    "cdml"];
    Z(4)["abcdijkl"] +=       FAE_3[    "ae"]*T_1(4)["ebcdijkl"];
    Z(4)["abcdijkl"] -=       FMI_3[    "mi"]*T_1(4)["abcdmjkl"];
    Z(4)["abcdijkl"] += 0.5*WABEF_3[  "abef"]*T_1(4)["efcdijkl"];
    Z(4)["abcdijkl"] += 0.5*WMNIJ_3[  "mnij"]*T_1(4)["abcdmnkl"];
    Z(4)["abcdijkl"] +=     WAMEI_3[  "amei"]*T_1(4)["ebcdjmkl"];

    U E1113 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*Z(4)["efghmnop"]);

    /***************************************************************************
     *               _
     * 1/2 <0|(1+L)[[H,T^(2)],T^(2)]|0>
     *               _
     * 1/2 <0|L^(1)[[H,T^(2)],T^(2)]|0>
     *
     **************************************************************************/

    WABCEJK_2[  "abcejk"]  =       WABEF_2[  "abef"]*  T(2)[    "fcjk"];
    WABCEJK_2[  "abcejk"] -=       WAMEI_2[  "amek"]*  T(2)[    "bcjm"];
    WABCEJK_2[  "abcejk"] +=     0.5*WABEF[  "abef"]*T_2(2)[    "fcjk"];
    WABCEJK_2[  "abcejk"] -=     0.5*WMNEF[  "mnef"]*T_2(4)["abcfmjkn"];

    WABMIJK_2[  "abmijk"]  =       WAMEI_2[  "amek"]*  T(2)[    "ebij"];
    WABMIJK_2[  "abmijk"] -=       WMNIJ_2[  "nmjk"]*  T(2)[    "abin"];
    WABMIJK_2[  "abmijk"] +=         WAMEI[  "amek"]*T_2(2)[    "ebij"];
    WABMIJK_2[  "abmijk"] -=     0.5*WMNIJ[  "mnkj"]*T_2(2)[    "abin"];
    WABMIJK_2[  "abmijk"] +=     0.5*WAMEF[  "bmef"]*T_2(3)[  "aefijk"];
    WABMIJK_2[  "abmijk"] +=     0.5*WMNEF[  "mnef"]*T_2(4)["abefijkn"];

    WAMNIJK_2[  "amnijk"]  =     0.5*WMNEF[  "mnef"]*T_2(3)[  "aefijk"];
    WAMNIJK_2[  "amnijk"] +=         WMNEJ[  "mnek"]*T_2(2)[    "aeij"];

    WABMEJI_2[  "abmeji"]  =         WAMEF[  "amef"]*T_2(2)[    "bfji"];
    WABMEJI_2[  "abmeji"] -=         WMNEJ[  "nmei"]*T_2(2)[    "abnj"];
    WABMEJI_2[  "abmeji"] +=     0.5*WMNEF[  "mnef"]*T_2(3)[  "abfnji"];

      WABEJ_2[    "abej"] +=           FME[    "me"]*T_2(2)[    "abmj"];

         Z(2)[    "abij"]  =         FAE_2[    "af"]*T_2(2)[    "fbij"];
         Z(2)[    "abij"] -=         FMI_2[    "ni"]*T_2(2)[    "abnj"];
         Z(2)[    "abij"] +=   0.5*WMNIJ_2[  "mnij"]*T_2(2)[    "abmn"];
         Z(2)[    "abij"] +=   0.5*WAMEI_2[  "amei"]*T_2(2)[    "ebjm"];

         Z(3)[  "abcijk"]  =       WABEJ_2[  "bcek"]*T_2(2)[    "aeij"];
         Z(3)[  "abcijk"] -=       WAMIJ_2[  "bmjk"]*T_2(2)[    "acim"];
         Z(3)[  "abcijk"] +=         FAE_2[    "ce"]*T_2(3)[  "abeijk"];
         Z(3)[  "abcijk"] -=         FMI_2[    "mk"]*T_2(3)[  "abcijm"];
         Z(3)[  "abcijk"] +=   0.5*WABEF_2[  "abef"]*T_2(3)[  "efcijk"];
         Z(3)[  "abcijk"] +=   0.5*WMNIJ_2[  "mnij"]*T_2(3)[  "abcmnk"];
         Z(3)[  "abcijk"] +=       WAMEI_2[  "amei"]*T_2(3)[  "ebcjmk"];

      WABEJ_2[    "abej"] -=           FME[    "me"]*T_2(2)[    "abmj"];

         Z(4)["abcdijkl"]  =     WABCEJK_2["abcejk"]*T_2(2)[    "edil"];
         Z(4)["abcdijkl"] -=     WABMIJK_2["abmijk"]*T_2(2)[    "cdml"];
         Z(4)["abcdijkl"] +=       WABEJ_2[  "abej"]*T_2(3)[  "ecdikl"];
         Z(4)["abcdijkl"] -=       WAMIJ_2[  "amij"]*T_2(3)[  "bcdmkl"];
         Z(4)["abcdijkl"] += 0.5*WAMNIJK_2["amnijk"]*T_2(3)[  "bcdmnl"];
         Z(4)["abcdijkl"] -=     WABMEJI_2["abmeji"]*T_2(3)[  "ecdmkl"];
         Z(4)["abcdijkl"] +=         FAE_2[    "ae"]*T_2(4)["ebcdijkl"];
         Z(4)["abcdijkl"] -=         FMI_2[    "mi"]*T_2(4)["abcdmjkl"];
         Z(4)["abcdijkl"] +=   0.5*WABEF_2[  "abef"]*T_2(4)["efcdijkl"];
         Z(4)["abcdijkl"] +=   0.5*WMNIJ_2[  "mnij"]*T_2(4)["abcdmnkl"];
         Z(4)["abcdijkl"] +=       WAMEI_2[  "amei"]*T_2(4)["ebcdjmkl"];

    U E0122 = (1.0/  4.0)*scalar(  L(2)[    "mnef"]*Z(2)[    "efmn"])+
              (1.0/ 36.0)*scalar(  L(3)[  "mnoefg"]*Z(3)[  "efgmno"]);
    U E1122 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*Z(4)["efghmnop"]);

    this->log(arena) << strprintf("\n");
    this->log(arena) << strprintf("CCSDT(Q-2): %18.15f\n", E2);

    this->log(arena) << strprintf("\n");
    this->log(arena) << strprintf("CCSDT(Q-3): %18.15f\n", E3);

    U E4 = E4_23+E4_4;

    this->log(arena) << strprintf("\n");
    this->log(arena) << strprintf("<0|L4^[1][Hbar,T2^[2]+T3^[2]]|0>: %18.15f\n", E4_23);
    this->log(arena) << strprintf("<0|L4^[1][Hbar,T4^[2]]|0>:        %18.15f\n", E4_4);
    this->log(arena) << strprintf("CCSDT(Q-4) Total:                 %18.15f\n", E4);

    U E5 = E5_1_23+E5_2_23+E5_3_23+E5_4_23+
           E5_1_4+ E5_2_4 +E5_3_4 +E5_4_4 +E1112+E0122;

    this->log(arena) << strprintf("\n");
    this->log(arena) << strprintf("<0|(L1^[2]+L2^[2]+L3^[2])[Hbar,T2^[2]+T3^[2]]|0>: %18.15f\n", E5_1_23+E5_2_23+E5_3_23);
    this->log(arena) << strprintf("<0|L4^[2][Hbar,T2^[2]+T3^[2]]|0>:                 %18.15f\n", E5_4_23);
    this->log(arena) << strprintf("<0|(L1^[2]+L2^[2]+L3^[2])[Hbar,T4^[2]]|0>:        %18.15f\n", E5_1_4+E5_2_4+E5_3_4);
    this->log(arena) << strprintf("<0|L4^[2][Hbar,T4^[2]]|0>:                        %18.15f\n", E5_4_4);
    this->log(arena) << strprintf("<0|L4^[1][[Hbar,T4^[1]],T2^[2]]|0>:               %18.15f\n", E1112);
    this->log(arena) << strprintf("1/2 <0|(1+L)[[Hbar,T^[2]],T^[2]]|0>:              %18.15f\n", E0122);
    this->log(arena) << strprintf("CCSDT(Q-5) Total:                                 %18.15f\n", E5);

    U E6 = E6_1_123+E6_2_123+E6_3_123+E6_4_123+
           E6_1_4  +E6_2_4  +E6_3_4  +E6_4_4  +E2112+E1113+E1122+E0123;

    this->log(arena) << strprintf("\n");
    this->log(arena) << strprintf("<0|(L1^[2]+L2^[2]+L3^[2])[Hbar,T1^[3]+T2^[3]+T3^[3]]|0>: %18.15f\n", E6_1_123+E6_2_123+E6_3_123);
    this->log(arena) << strprintf("<0|L4^[2][Hbar,T1^[3]+T2^[3]+T3^[3]]|0>:                 %18.15f\n", E6_4_123);
    this->log(arena) << strprintf("<0|(L1^[2]+L2^[2]+L3^[2])[Hbar,T4^[3]]|0>:               %18.15f\n", E6_1_4+E6_2_4+E6_3_4);
    this->log(arena) << strprintf("<0|L4^[2][Hbar,T4^[3]]|0>:                               %18.15f\n", E6_4_4);
    this->log(arena) << strprintf("<0|L4^[2][[Hbar,T4^[1]],T2^[2]]|0>:                      %18.15f\n", E2112);
    this->log(arena) << strprintf("<0|L4^[1][[Hbar,T4^[1]],T1^[3]+T2^[3]]|0>:               %18.15f\n", E1113);
    this->log(arena) << strprintf("1/2 <0|L4^[1][[Hbar,T^[2]],T^[2]]|0>:                    %18.15f\n", E1122);
    this->log(arena) << strprintf("<0|(1+L)[[Hbar,T^[2]],T^[3]]|0>:                         %18.15f\n", E0123);
    this->log(arena) << strprintf("CCSDT(Q-6) Total:                                        %18.15f\n", E6);

    this->log(arena) << strprintf("\n");

    this->put("E(2)", new U(E2));
    this->put("E(3)", new U(E3));
    this->put("E(4)", new U(E4));
    this->put("E(5)", new U(E5));
    this->put("E(6)", new U(E6));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSDT_Q_N_opt);
REGISTER_TASK(aquarius::cc::CCSDT_Q_N_opt<double>,"ccsdt(q-n)_opt");
