#include "ccsd_tq_n.hpp"

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
CCSD_TQ_N<U>::CCSD_TQ_N(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("ccsd.T", "T"));
    reqs.push_back(Requirement("ccsd.L", "L"));
    this->addProduct(Product("double", "E(2)", reqs));
    this->addProduct(Product("double", "E(3)", reqs));
    this->addProduct(Product("double", "E(4)", reqs));
}

template <typename U>
bool CCSD_TQ_N<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U>>("H");
    const STTwoElectronOperator<U>& Hbar = this->template get<STTwoElectronOperator<U>>("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;
    const PointGroup& group = occ.group;

    Denominator<U> D(H);
    const ExcitationOperator  <U,2>& T = this->template get<ExcitationOperator  <U,2>>("T");
    const DeexcitationOperator<U,2>& L = this->template get<DeexcitationOperator<U,2>>("L");

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

    DeexcitationOperator<U,4> L_1("L^(1)", arena, occ, vrt);
    DeexcitationOperator<U,4> L_2("L^(2)", arena, occ, vrt);
    DeexcitationOperator<U,4> L_3("L^(3)", arena, occ, vrt);

    ExcitationOperator<U,4> DT_1("DT^(1)", arena, occ, vrt);
    ExcitationOperator<U,4> DT_2("DT^(2)", arena, occ, vrt);
    ExcitationOperator<U,4> DT_3("DT^(3)", arena, occ, vrt);

    DeexcitationOperator<U,4> DL_1("DL^(1)", arena, occ, vrt);
    DeexcitationOperator<U,4> DL_2("DL^(2)", arena, occ, vrt);
    DeexcitationOperator<U,4> DL_3("DL^(3)", arena, occ, vrt);

    SpinorbitalTensor<U> WTWABEJ(WABEJ);
    WTWABEJ["abej"] += FME["me"]*T(2)["abmj"];

    SpinorbitalTensor<U> WABEJ_1(WABEJ);
    SpinorbitalTensor<U> WAMIJ_1(WAMIJ);
    SpinorbitalTensor<U> WABCEJK_1(WABCEJK);
    SpinorbitalTensor<U> WABMEJI_1(WABMEJI);
    SpinorbitalTensor<U> WABMIJK_1(WABMIJK);
    SpinorbitalTensor<U> WAMNIJK_1(WAMNIJK);

    SpinorbitalTensor<U> FME_2(FME);
    SpinorbitalTensor<U> FMI_2(FMI);
    SpinorbitalTensor<U> FAE_2(FAE);
    SpinorbitalTensor<U> WAMEI_2(WAMEI);
    SpinorbitalTensor<U> WMNIJ_2(WMNIJ);
    SpinorbitalTensor<U> WABEF_2(WABEF);
    SpinorbitalTensor<U> WABEJ_2(WABEJ);
    SpinorbitalTensor<U> WAMIJ_2(WAMIJ);
    SpinorbitalTensor<U> WMNEJ_2(WMNEJ);
    SpinorbitalTensor<U> WAMEF_2(WAMEF);
    SpinorbitalTensor<U> WABCEJK_2(WABCEJK);
    SpinorbitalTensor<U> WABMIJK_2(WABMIJK);
    SpinorbitalTensor<U> WAMNIJK_2(WAMNIJK);

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

    SpinorbitalTensor<U> GIJAK_2(WMNEJ);
    SpinorbitalTensor<U> GAIBC_2(WAMEF);
    SpinorbitalTensor<U> GAIJBCD_2("G(aij,bcd)", arena, group, {vrt,occ}, {1,2}, {3,0});
    SpinorbitalTensor<U> GIJKABL_2("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1});
    SpinorbitalTensor<U> GIJKALM_2("G(ijk,alm)", arena, group, {vrt,occ}, {0,3}, {1,2});
    ExcitationOperator<U,3> Z("Z", arena, occ, vrt);

    /***************************************************************************
     *
     * T^(1)
     *
     **************************************************************************/

    WABCEJK[  "abcejk"]  = 0.5*WABEF[  "abef"]*T(2)["fcjk"];

    WABMIJK[  "abmijk"]  =     WAMEI[  "amek"]*T(2)["ebij"];
    WABMIJK[  "abmijk"] -= 0.5*WMNIJ[  "nmjk"]*T(2)["abin"];

    WAMNIJK[  "amnijk"]  =     WMNEJ[  "mnek"]*T(2)["aeij"];

    WABMEJI[  "abmeji"]  =     WAMEF[  "amef"]*T(2)["bfji"];
    WABMEJI[  "abmeji"] -=     WMNEJ[  "nmei"]*T(2)["abnj"];

    DT_1(3)[  "abcijk"]  =   WTWABEJ[  "bcek"]*T(2)["aeij"];
    DT_1(3)[  "abcijk"] -=     WAMIJ[  "bmjk"]*T(2)["acim"];

    DT_1(4)["abcdijkl"]  =   WABCEJK["abcejk"]*T(2)["edil"];
    DT_1(4)["abcdijkl"] -=   WABMIJK["abmijk"]*T(2)["cdml"];

    T_1 = DT_1;
    T_1.weight(D);

    /***************************************************************************
     *
     * T^(2)
     *
     **************************************************************************/

      WAMIJ_1[    "amij"]  =   0.5*WMNEF[  "mnef"]*T_1(3)[  "aefijn"];
      WABEJ_1[    "abej"]  =  -0.5*WMNEF[  "mnef"]*T_1(3)[  "afbmnj"];

    WABCEJK_1[  "abcejk"]  =  -0.5*WMNEF[  "mnef"]*T_1(4)["abcfmjkn"];

    WABMIJK_1[  "abmijk"]  =   0.5*WMNEF[  "mnef"]*T_1(4)["abefijkn"];
    WABMIJK_1[  "abmijk"] +=   0.5*WAMEF[  "bmef"]*T_1(3)[  "aefijk"];

      DT_2(1)[      "ai"]  =  0.25*WMNEF[  "mnef"]*T_1(3)[  "aefimn"];

      DT_2(2)[    "abij"]  =   0.5*WAMEF[  "bmef"]*T_1(3)[  "aefijm"];
      DT_2(2)[    "abij"] -=   0.5*WMNEJ[  "mnej"]*T_1(3)[  "abeinm"];
      DT_2(2)[    "abij"] +=         FME[    "me"]*T_1(3)[  "abeijm"];
      DT_2(2)[    "abij"] +=  0.25*WMNEF[  "mnef"]*T_1(4)["abefijmn"];

      DT_2(3)[  "abcijk"]  =     WABEJ_1[  "bcek"]*  T(2)[    "aeij"];
      DT_2(3)[  "abcijk"] -=     WAMIJ_1[  "bmjk"]*  T(2)[    "acim"];
      DT_2(3)[  "abcijk"] +=         FAE[    "ce"]*T_1(3)[  "abeijk"];
      DT_2(3)[  "abcijk"] -=         FMI[    "mk"]*T_1(3)[  "abcijm"];
      DT_2(3)[  "abcijk"] +=   0.5*WABEF[  "abef"]*T_1(3)[  "efcijk"];
      DT_2(3)[  "abcijk"] +=   0.5*WMNIJ[  "mnij"]*T_1(3)[  "abcmnk"];
      DT_2(3)[  "abcijk"] +=       WAMEI[  "amei"]*T_1(3)[  "ebcjmk"];
      DT_2(3)[  "abcijk"] +=         FME[    "me"]*T_1(4)["abceijkm"];
      DT_2(3)[  "abcijk"] +=   0.5*WAMEF[  "amef"]*T_1(4)["efbcimjk"];
      DT_2(3)[  "abcijk"] -=   0.5*WMNEJ[  "mnek"]*T_1(4)["abecijmn"];

      DT_2(4)["abcdijkl"]  =   WABCEJK_1["abcejk"]*  T(2)[    "edil"];
      DT_2(4)["abcdijkl"] -=   WABMIJK_1["abmijk"]*  T(2)[    "cdml"];
      DT_2(4)["abcdijkl"] +=       WABEJ[  "abej"]*T_1(3)[  "ecdikl"];
      DT_2(4)["abcdijkl"] -=       WAMIJ[  "amij"]*T_1(3)[  "bcdmkl"];
      DT_2(4)["abcdijkl"] += 0.5*WAMNIJK["amnijk"]*T_1(3)[  "bcdmnl"];
      DT_2(4)["abcdijkl"] -=     WABMEJI["abmeji"]*T_1(3)[  "ecdmkl"];
      DT_2(4)["abcdijkl"] +=         FAE[    "ae"]*T_1(4)["ebcdijkl"];
      DT_2(4)["abcdijkl"] -=         FMI[    "mi"]*T_1(4)["abcdmjkl"];
      DT_2(4)["abcdijkl"] +=   0.5*WABEF[  "abef"]*T_1(4)["efcdijkl"];
      DT_2(4)["abcdijkl"] +=   0.5*WMNIJ[  "mnij"]*T_1(4)["abcdmnkl"];
      DT_2(4)["abcdijkl"] +=       WAMEI[  "amei"]*T_1(4)["ebcdjmkl"];

    T_2 = DT_2;
    T_2.weight(D);

    /***************************************************************************
     *
     * T^(3)
     *
     **************************************************************************/

        FME_2[      "me"]  =         WMNEF[  "mnef"]*T_2(1)[      "fn"];

      WAMIJ_2[    "amij"]  =         WMNEJ[  "nmej"]*T_2(2)[    "aein"];
      WAMIJ_2[    "amij"] +=     0.5*WAMEF[  "amef"]*T_2(2)[    "efij"];
      WAMIJ_2[    "amij"] -=         WMNIJ[  "nmij"]*T_2(1)[      "an"];
      WAMIJ_2[    "amij"] +=         WAMEI[  "amej"]*T_2(1)[      "ei"];
      WAMIJ_2[    "amij"] +=     0.5*WMNEF[  "mnef"]*T_2(3)[  "aefijn"];
      WAMIJ_2[    "amij"] +=         FME_2[    "me"]*  T(2)[    "aeij"];

      WABEJ_2[    "abej"]  =         WAMEF[  "amef"]*T_2(2)[    "fbmj"];
      WABEJ_2[    "abej"] +=     0.5*WMNEJ[  "mnej"]*T_2(2)[    "abmn"];
      WABEJ_2[    "abej"] +=         WABEF[  "abef"]*T_2(1)[      "fj"];
      WABEJ_2[    "abej"] -=         WAMEI[  "amej"]*T_2(1)[      "bm"];
      WABEJ_2[    "abej"] -=     0.5*WMNEF[  "mnef"]*T_2(3)[  "afbmnj"];

      WABCEJK[  "abcejk"]  =         WABEF[  "abef"]*  T(2)[    "fcjk"];
      WABCEJK[  "abcejk"] -=         WAMEI[  "amek"]*  T(2)[    "bcjm"];

      WABMIJK[  "abmijk"]  =         WAMEI[  "amek"]*  T(2)[    "ebij"];
      WABMIJK[  "abmijk"] -=         WMNIJ[  "nmjk"]*  T(2)[    "abin"];

      WAMNIJK[  "amnijk"]  =         WMNEJ[  "mnek"]*  T(2)[    "aeij"];

      WABMEJI[  "abmeji"]  =         WAMEF[  "amef"]*  T(2)[    "bfji"];
      WABMEJI[  "abmeji"] -=         WMNEJ[  "nmei"]*  T(2)[    "abnj"];

    WAMNIJK_1[  "amnijk"]  =     0.5*WMNEF[  "mnef"]*T_1(3)[  "aefijk"];
    WABMEJI_1[  "abmeji"]  =     0.5*WMNEF[  "mnef"]*T_1(3)[  "abfnji"];

      WAMEI_2[    "amei"]  =         WMNEF[  "mnef"]*T_2(2)[    "afni"];
      WAMEI_2[    "amei"] +=         WAMEF[  "amef"]*T_2(1)[      "fi"];
      WAMEI_2[    "amei"] -=         WMNEJ[  "nmei"]*T_2(1)[      "an"];

      WMNIJ_2[    "mnij"]  =     0.5*WMNEF[  "mnef"]*T_2(2)[    "efij"];
      WMNIJ_2[    "mnij"] +=         WMNEJ[  "mnej"]*T_2(1)[      "ei"];

      WABEF_2[    "abef"]  =     0.5*WMNEF[  "mnef"]*T_2(2)[    "abmn"];
      WABEF_2[    "abef"] -=         WAMEF[  "amef"]*T_2(1)[      "bm"];

    WABCEJK_2[  "abcejk"]  =   0.5*WABEF_2[  "abef"]*  T(2)[    "fcjk"];
    WABCEJK_2[  "abcejk"] -=     0.5*WMNEF[  "mnef"]*T_2(4)["abcfmjkn"];

    WABMIJK_2[  "abmijk"]  =       WAMEI_2[  "amek"]*  T(2)[    "ebij"];
    WABMIJK_2[  "abmijk"] -=   0.5*WMNIJ_2[  "nmjk"]*  T(2)[    "abin"];
    WABMIJK_2[  "abmijk"] +=     0.5*WAMEF[  "bmef"]*T_2(3)[  "aefijk"];
    WABMIJK_2[  "abmijk"] +=     0.5*WMNEF[  "mnef"]*T_2(4)["abefijkn"];

      DT_3(3)[  "abcijk"]  =       WABEJ_2[  "bcek"]*  T(2)[    "aeij"];
      DT_3(3)[  "abcijk"] -=       WAMIJ_2[  "bmjk"]*  T(2)[    "acim"];
      DT_3(3)[  "abcijk"] +=         WABEJ[  "bcek"]*T_2(2)[    "aeij"];
      DT_3(3)[  "abcijk"] -=         WAMIJ[  "bmjk"]*T_2(2)[    "acim"];
      DT_3(3)[  "abcijk"] +=           FAE[    "ce"]*T_2(3)[  "abeijk"];
      DT_3(3)[  "abcijk"] -=           FMI[    "mk"]*T_2(3)[  "abcijm"];
      DT_3(3)[  "abcijk"] +=     0.5*WABEF[  "abef"]*T_2(3)[  "efcijk"];
      DT_3(3)[  "abcijk"] +=     0.5*WMNIJ[  "mnij"]*T_2(3)[  "abcmnk"];
      DT_3(3)[  "abcijk"] +=         WAMEI[  "amei"]*T_2(3)[  "ebcjmk"];
      DT_3(3)[  "abcijk"] +=           FME[    "me"]*T_2(4)["abceijkm"];
      DT_3(3)[  "abcijk"] +=     0.5*WAMEF[  "amef"]*T_2(4)["efbcimjk"];
      DT_3(3)[  "abcijk"] -=     0.5*WMNEJ[  "mnek"]*T_2(4)["abecijmn"];

      DT_3(4)["abcdijkl"]  =     WABCEJK_2["abcejk"]*  T(2)[    "edil"];
      DT_3(4)["abcdijkl"] -=     WABMIJK_2["abmijk"]*  T(2)[    "cdml"];
      DT_3(4)["abcdijkl"] +=       WABEJ_1[  "abej"]*T_1(3)[  "ecdikl"];
      DT_3(4)["abcdijkl"] -=       WAMIJ_1[  "amij"]*T_1(3)[  "bcdmkl"];
      DT_3(4)["abcdijkl"] += 0.5*WAMNIJK_1["amnijk"]*T_1(3)[  "bcdmnl"];
      DT_3(4)["abcdijkl"] -=     WABMEJI_1["abmeji"]*T_1(3)[  "ecdmkl"];
      DT_3(4)["abcdijkl"] +=       WABCEJK["abcejk"]*T_2(2)[    "edil"];
      DT_3(4)["abcdijkl"] -=       WABMIJK["abmijk"]*T_2(2)[    "cdml"];
      DT_3(4)["abcdijkl"] +=         WABEJ[  "abej"]*T_2(3)[  "ecdikl"];
      DT_3(4)["abcdijkl"] -=         WAMIJ[  "amij"]*T_2(3)[  "bcdmkl"];
      DT_3(4)["abcdijkl"] +=   0.5*WAMNIJK["amnijk"]*T_2(3)[  "bcdmnl"];
      DT_3(4)["abcdijkl"] -=       WABMEJI["abmeji"]*T_2(3)[  "ecdmkl"];
      DT_3(4)["abcdijkl"] +=           FAE[    "ae"]*T_2(4)["ebcdijkl"];
      DT_3(4)["abcdijkl"] -=           FMI[    "mi"]*T_2(4)["abcdmjkl"];
      DT_3(4)["abcdijkl"] +=     0.5*WABEF[  "abef"]*T_2(4)["efcdijkl"];
      DT_3(4)["abcdijkl"] +=     0.5*WMNIJ[  "mnij"]*T_2(4)["abcdmnkl"];
      DT_3(4)["abcdijkl"] +=         WAMEI[  "amei"]*T_2(4)["ebcdjmkl"];

    T_3 = DT_3;
    T_3.weight(D);

    /***************************************************************************
     *
     * L^(1)
     *
     **************************************************************************/

    DL_1(3)[  "ijkabc"]  = WMNEF["ijab"]*L(1)[  "kc"];
    DL_1(3)[  "ijkabc"] +=   FME[  "ia"]*L(2)["jkbc"];
    DL_1(3)[  "ijkabc"] += WAMEF["ekbc"]*L(2)["ijae"];
    DL_1(3)[  "ijkabc"] -= WMNEJ["ijam"]*L(2)["mkbc"];

    DL_1(4)["ijklabcd"]  = WMNEF["ijab"]*L(2)["klcd"];

    L_1 = DL_1;
    L_1.weight(D);

    /***************************************************************************
     *
     * L^(2)
     *
     **************************************************************************/

      GIJAK_1[    "ijak"]  =  (1.0/ 2.0)*   T(2)[  "efkm"]*   L_1(3)[  "ijmaef"];
      GAIBC_1[    "aibc"]  = -(1.0/ 2.0)*   T(2)[  "aemn"]*   L_1(3)[  "minbce"];

        DAI_1[      "ai"]  =  (1.0/ 4.0)* T_1(3)["efamni"]*     L(2)[    "mnef"];
        DAI_1[      "ai"] -=  (1.0/ 2.0)*   T(2)[  "eamn"]*  GIJAK_1[    "mnei"];

    GIJKABL_1[  "ijkabl"]  =  (1.0/ 2.0)*   T(2)[  "eflm"]*   L_1(4)["ijkmabef"];
    GAIJBCD_1[  "aijbcd"]  = -(1.0/ 2.0)*   T(2)[  "aemn"]*   L_1(4)["mijnbcde"];

      GABCD_1[    "abcd"]  = -(1.0/ 4.0)*   T(2)[  "bemn"]*GAIJBCD_1[  "amncde"];
      GAIBJ_1[    "aibj"]  = -(1.0/ 2.0)*   T(2)[  "aenm"]*GIJKABL_1[  "mniebj"];
      GIJKL_1[    "ijkl"]  =  (1.0/ 4.0)*   T(2)[  "efmk"]*GIJKABL_1[  "mijefl"];

      DL_2(1)[      "ia"]  =  (1.0/ 2.0)*WABEJ_1[  "efam"]*     L(2)[    "imef"];
      DL_2(1)[      "ia"] -=  (1.0/ 2.0)*WAMIJ_1[  "eimn"]*     L(2)[    "mnea"];
      DL_2(1)[      "ia"] +=               WMNEF[  "miea"]*    DAI_1[      "em"];
      DL_2(1)[      "ia"] -=  (1.0/ 2.0)*  WABEF[  "efga"]*  GAIBC_1[    "gief"];
      DL_2(1)[      "ia"] +=               WAMEI[  "eifm"]*  GAIBC_1[    "fmea"];
      DL_2(1)[      "ia"] -=               WAMEI[  "eman"]*  GIJAK_1[    "inem"];
      DL_2(1)[      "ia"] +=  (1.0/ 2.0)*  WMNIJ[  "imno"]*  GIJAK_1[    "noam"];
      DL_2(1)[      "ia"] -=  (1.0/ 2.0)*  WAMEF[  "gief"]*  GABCD_1[    "efga"];
      DL_2(1)[      "ia"] +=               WAMEF[  "fmea"]*  GAIBJ_1[    "eifm"];
      DL_2(1)[      "ia"] -=               WMNEJ[  "inem"]*  GAIBJ_1[    "eman"];
      DL_2(1)[      "ia"] +=  (1.0/ 2.0)*  WMNEJ[  "noam"]*  GIJKL_1[    "imno"];

      DL_2(2)[    "ijab"] -=               WAMEF[  "fiae"]*  GAIBC_1[    "ejbf"];
      DL_2(2)[    "ijab"] -=               WMNEJ[  "ijem"]*  GAIBC_1[    "emab"];
      DL_2(2)[    "ijab"] -=               WAMEF[  "emab"]*  GIJAK_1[    "ijem"];
      DL_2(2)[    "ijab"] -=               WMNEJ[  "niam"]*  GIJAK_1[    "mjbn"];
      DL_2(2)[    "ijab"] +=  (1.0/ 2.0)*  WMNEF[  "ijef"]*  GABCD_1[    "efab"];
      DL_2(2)[    "ijab"] +=               WMNEF[  "imea"]*  GAIBJ_1[    "ejbm"];
      DL_2(2)[    "ijab"] +=  (1.0/ 2.0)*  WMNEF[  "mnab"]*  GIJKL_1[    "ijmn"];
      DL_2(2)[    "ijab"] +=  (1.0/ 2.0)*  WABEJ[  "efbm"]*   L_1(3)[  "ijmaef"];
      DL_2(2)[    "ijab"] -=  (1.0/ 2.0)*  WAMIJ[  "ejnm"]*   L_1(3)[  "imnabe"];
      DL_2(2)[    "ijab"] +=  (1.0/12.0)*WABCEJK["efgbmn"]*   L_1(4)["ijmnaefg"];
      DL_2(2)[    "ijab"] -=  (1.0/12.0)*WABMIJK["efimno"]*   L_1(4)["mnojefab"];

      DL_2(3)[  "ijkabc"]  =               WMNEF[  "ijae"]*  GAIBC_1[    "ekbc"];
      DL_2(3)[  "ijkabc"] -=               WMNEF[  "mkbc"]*  GIJAK_1[    "ijam"];
      DL_2(3)[  "ijkabc"] -=               WAMEF[  "embc"]*GIJKABL_1[  "ijkaem"];
      DL_2(3)[  "ijkabc"] +=                 FAE[    "ea"]*   L_1(3)[  "ijkebc"];
      DL_2(3)[  "ijkabc"] -=                 FMI[    "im"]*   L_1(3)[  "mjkabc"];
      DL_2(3)[  "ijkabc"] +=  (1.0/ 2.0)*  WABEF[  "efab"]*   L_1(3)[  "ijkefc"];
      DL_2(3)[  "ijkabc"] +=  (1.0/ 2.0)*  WMNIJ[  "ijmn"]*   L_1(3)[  "mnkabc"];
      DL_2(3)[  "ijkabc"] +=               WAMEI[  "eiam"]*   L_1(3)[  "mjkbec"];
      DL_2(3)[  "ijkabc"] +=  (1.0/ 2.0)*  WABEJ[  "efcm"]*   L_1(4)["ijkmabef"];
      DL_2(3)[  "ijkabc"] -=  (1.0/ 2.0)*  WAMIJ[  "eknm"]*   L_1(4)["ijmnabce"];
      DL_2(3)[  "ijkabc"] -=  (1.0/ 4.0)*WABMEJI["efkcnm"]*   L_1(4)["ijmnabef"];
      DL_2(3)[  "ijkabc"] +=  (1.0/ 6.0)*WAMNIJK["eijmno"]*   L_1(4)["mnokeabc"];

      DL_2(4)["ijklabcd"]  =               WMNEF[  "ijae"]*GAIJBCD_1[  "eklbcd"];
      DL_2(4)["ijklabcd"] -=               WMNEF[  "mlcd"]*GIJKABL_1[  "ijkabm"];
      DL_2(4)["ijklabcd"] +=                 FME[    "ia"]*   L_1(3)[  "jklbcd"];
      DL_2(4)["ijklabcd"] +=               WAMEF[  "ejab"]*   L_1(3)[  "iklecd"];
      DL_2(4)["ijklabcd"] -=               WMNEJ[  "ijam"]*   L_1(3)[  "mklbcd"];
      DL_2(4)["ijklabcd"] +=                 FAE[    "ea"]*   L_1(4)["ijklebcd"];
      DL_2(4)["ijklabcd"] -=                 FMI[    "im"]*   L_1(4)["mjklabcd"];
      DL_2(4)["ijklabcd"] +=  (1.0/ 2.0)*  WABEF[  "efab"]*   L_1(4)["ijklefcd"];
      DL_2(4)["ijklabcd"] +=  (1.0/ 2.0)*  WMNIJ[  "ijmn"]*   L_1(4)["mnklabcd"];
      DL_2(4)["ijklabcd"] +=               WAMEI[  "eiam"]*   L_1(4)["mjklbecd"];

    L_2 = DL_2;
    L_2.weight(D);

    /***************************************************************************
     *
     * L^(3)
     *
     **************************************************************************/

      WMNEJ_2[    "mnej"]  =                 WMNEF[  "mnef"]*   T_2(1)[      "fj"];
      WAMEF_2[    "amef"]  =                -WMNEF[  "nmef"]*   T_2(1)[      "an"];

      GIJAK_2[    "ijak"]  =  (1.0/ 2.0)*     T(2)[  "efkm"]*   L_2(3)[  "ijmaef"];
      GIJAK_2[    "ijak"] +=  (1.0/12.0)*   T_1(3)["efgkmn"]*   L_1(4)["ijmnaefg"];

      GAIBC_2[    "aibc"]  = -(1.0/ 2.0)*     T(2)[  "aemn"]*   L_2(3)[  "minbce"];
      GAIBC_2[    "aibc"] -=  (1.0/12.0)*   T_1(3)["aefmno"]*   L_1(4)["minobcef"];

    GIJKABL_2[  "ijkabl"]  =  (1.0/ 2.0)*     T(2)[  "eflm"]*   L_2(4)["ijkmabef"];
    GAIJBCD_2[  "aijbcd"]  = -(1.0/ 2.0)*     T(2)[  "aemn"]*   L_2(4)["mijnbcde"];
    GIJKALM_2[  "ijkalm"]  =  (1.0/ 6.0)*   T_1(3)["efglmn"]*   L_1(4)["ijknaefg"];

      DL_3(3)[  "ijkabc"]  =                 FME_2[    "ia"]*     L(2)[    "jkbc"];
      DL_3(3)[  "ijkabc"] +=               WAMEF_2[  "ekbc"]*     L(2)[    "ijae"];
      DL_3(3)[  "ijkabc"] -=               WMNEJ_2[  "ijam"]*     L(2)[    "mkbc"];
      DL_3(3)[  "ijkabc"] +=  (1.0/ 2.0)*  WABEJ_1[  "efcm"]*   L_1(4)["ijkmabef"];
      DL_3(3)[  "ijkabc"] -=  (1.0/ 2.0)*  WAMIJ_1[  "eknm"]*   L_1(4)["ijmnabce"];
      DL_3(3)[  "ijkabc"] -=  (1.0/ 2.0)*WABMEJI_1["efkcnm"]*   L_1(4)["ijmnabef"]; // Add factor of 2 to account for 1/2 in def'n of WABMEJI_1
      DL_3(3)[  "ijkabc"] +=  (1.0/ 6.0)*WAMNIJK_1["eijmno"]*   L_1(4)["mnokeabc"];
      DL_3(3)[  "ijkabc"] +=                 WMNEF[  "ijab"]*   L_2(1)[      "kc"];
      DL_3(3)[  "ijkabc"] +=                   FME[    "ia"]*   L_2(2)[    "jkbc"];
      DL_3(3)[  "ijkabc"] +=                 WAMEF[  "ekbc"]*   L_2(2)[    "ijae"];
      DL_3(3)[  "ijkabc"] -=                 WMNEJ[  "ijam"]*   L_2(2)[    "mkbc"];
      DL_3(3)[  "ijkabc"] +=                 WMNEF[  "ijae"]*  GAIBC_2[    "ekbc"];
      DL_3(3)[  "ijkabc"] -=                 WMNEF[  "mkbc"]*  GIJAK_2[    "ijam"];
      DL_3(3)[  "ijkabc"] -=                 WAMEF[  "embc"]*GIJKABL_2[  "ijkaem"];
      DL_3(3)[  "ijkabc"] +=  (1.0/ 2.0)*    WMNEF[  "mnbc"]*GIJKALM_2[  "ijkamn"];
      DL_3(3)[  "ijkabc"] +=                   FAE[    "ea"]*   L_2(3)[  "ijkebc"];
      DL_3(3)[  "ijkabc"] -=                   FMI[    "im"]*   L_2(3)[  "mjkabc"];
      DL_3(3)[  "ijkabc"] +=  (1.0/ 2.0)*    WABEF[  "efab"]*   L_2(3)[  "ijkefc"];
      DL_3(3)[  "ijkabc"] +=  (1.0/ 2.0)*    WMNIJ[  "ijmn"]*   L_2(3)[  "mnkabc"];
      DL_3(3)[  "ijkabc"] +=                 WAMEI[  "eiam"]*   L_2(3)[  "mjkbec"];
      DL_3(3)[  "ijkabc"] +=  (1.0/ 2.0)*    WABEJ[  "efcm"]*   L_2(4)["ijkmabef"];
      DL_3(3)[  "ijkabc"] -=  (1.0/ 2.0)*    WAMIJ[  "eknm"]*   L_2(4)["ijmnabce"];
      DL_3(3)[  "ijkabc"] -=  (1.0/ 4.0)*  WABMEJI["efkcnm"]*   L_2(4)["ijmnabef"];
      DL_3(3)[  "ijkabc"] +=  (1.0/ 6.0)*  WAMNIJK["eijmno"]*   L_2(4)["mnokeabc"];

      DL_3(4)["ijklabcd"]  =                 WMNEF[  "ijab"]*   L_2(2)[    "klcd"];
      DL_3(4)["ijklabcd"] +=                   FME[    "ia"]*   L_2(3)[  "jklbcd"];
      DL_3(4)["ijklabcd"] +=                 WAMEF[  "ejab"]*   L_2(3)[  "iklecd"];
      DL_3(4)["ijklabcd"] -=                 WMNEJ[  "ijam"]*   L_2(3)[  "mklbcd"];
      DL_3(4)["ijklabcd"] +=                 WMNEF[  "ijae"]*GAIJBCD_2[  "eklbcd"];
      DL_3(4)["ijklabcd"] -=                 WMNEF[  "mlcd"]*GIJKABL_2[  "ijkabm"];
      DL_3(4)["ijklabcd"] +=                   FAE[    "ea"]*   L_2(4)["ijklebcd"];
      DL_3(4)["ijklabcd"] -=                   FMI[    "im"]*   L_2(4)["mjklabcd"];
      DL_3(4)["ijklabcd"] +=  (1.0/ 2.0)*    WABEF[  "efab"]*   L_2(4)["ijklefcd"];
      DL_3(4)["ijklabcd"] +=  (1.0/ 2.0)*    WMNIJ[  "ijmn"]*   L_2(4)["mnklabcd"];
      DL_3(4)["ijklabcd"] +=                 WAMEI[  "eiam"]*   L_2(4)["mjklbecd"];

    L_3 = DL_3;
    L_3.weight(D);

    /***************************************************************************
     *           _
     * <0|(1+L)[[H,T^(1)],T^(2)]|0>
     *
     **************************************************************************/

    Z(2)["abij"]  =       FME_2[  "me"]*T_1(3)["eabmij"];
    Z(2)["abij"] += 0.5*WAMEF_2["bmef"]*T_1(3)["aefijm"];
    Z(2)["abij"] -= 0.5*WMNEJ_2["nmej"]*T_1(3)["abeimn"];

    U E0112 = 0.25*scalar(L(2)["mnef"]*Z(2)["efmn"]);

    /***************************************************************************
     *           _
     * <0|L^(1)[[H,T^(1)],T^(1)]|0>
     *
     **************************************************************************/

    Z(3)["abcijk"]  = (1.0/ 2.0)*  WABEJ_1[  "efcm"]*L_1(4)["ijkmabef"];
    Z(3)["abcijk"] -= (1.0/ 2.0)*  WAMIJ_1[  "eknm"]*L_1(4)["ijmnabce"];
    Z(3)["abcijk"] -= (1.0/ 4.0)*WABMEJI_1["efkcnm"]*L_1(4)["ijmnabef"];
    Z(3)["abcijk"] += (1.0/ 6.0)*WAMNIJK_1["eijmno"]*L_1(4)["mnokeabc"];

    U E1111 = (1.0/36.0)*scalar(T_1(3)["efgmno"]*Z(3)["efgmno"]);

    /***************************************************************************
     *
     * CCSD(TQ-2)
     *
     **************************************************************************/

    U E101 =  (1.0/ 36.0)*scalar(L_1(3)[  "mnoefg"]*DT_1(3)[  "efgmno"])
             +(1.0/576.0)*scalar(L_1(4)["mnopefgh"]*DT_1(4)["efghmnop"]);
    U E2 = E101;

    /***************************************************************************
     *
     * CCSD(TQ-3)
     *
     **************************************************************************/

    U E201 =  (1.0/ 36.0)*scalar(L_2(3)[  "mnoefg"]*DT_1(3)[  "efgmno"])
             +(1.0/576.0)*scalar(L_2(4)["mnopefgh"]*DT_1(4)["efghmnop"]);
    U E102 =  (1.0/ 36.0)*scalar(L_1(3)[  "mnoefg"]*DT_2(3)[  "efgmno"])
             +(1.0/576.0)*scalar(L_1(4)["mnopefgh"]*DT_2(4)["efghmnop"]);

    printf("E201: %18.15f\n", E201);
    printf("E102: %18.15f\n", E102);
    printf("\n");

    U E3 = E102;

    /***************************************************************************
     *
     * CCSD(TQ-4)
     *
     **************************************************************************/

    U E301 =  (1.0/ 36.0)*scalar(L_3(3)[  "mnoefg"]*DT_1(3)[  "efgmno"])
             +(1.0/576.0)*scalar(L_3(4)["mnopefgh"]*DT_1(4)["efghmnop"]);
    U E202 =  (1.0/  1.0)*scalar(L_2(1)[      "me"]*DT_2(1)[      "em"])
             +(1.0/  4.0)*scalar(L_2(2)[    "mnef"]*DT_2(2)[    "efmn"])
             +(1.0/ 36.0)*scalar(L_2(3)[  "mnoefg"]*DT_2(3)[  "efgmno"])
             +(1.0/576.0)*scalar(L_2(4)["mnopefgh"]*DT_2(4)["efghmnop"]);
    U E103 =  (1.0/ 36.0)*scalar(L_1(3)[  "mnoefg"]*DT_3(3)[  "efgmno"])
             +(1.0/576.0)*scalar(L_1(4)["mnopefgh"]*DT_3(4)["efghmnop"]);

    printf("E301: %18.15f\n", E301);
    printf("E202: %18.15f\n", E202);
    printf("E103: %18.15f\n", E103);
    printf("\n");

    printf("E0112: %18.15f\n", E0112);
    printf("E1111: %18.15f\n", E1111);
    printf("\n");

    printf("E301:               %18.15f\n", E301);
    printf("E202+E0112+2*E1111: %18.15f\n", E202+E0112+2*E1111);
    printf("\n");

    printf("E202:             %18.15f\n", E202);
    printf("E103+E0112-E1111: %18.15f\n", E103+E0112-E1111);
    printf("\n");

    U E4 = E103+E0112;

    printf("CCSD(TQ-2): %18.15f\n", E2);
    printf("CCSD(TQ-3): %18.15f\n", E3);
    printf("CCSD(TQ-4): %18.15f\n", E4);

    this->put("E(2)", new U(E2));
    this->put("E(3)", new U(E3));
    this->put("E(4)", new U(E4));

    printf("\n");
    printf("T1[2]: %15.12f\n", T_2(1)({1,0},{0,1}).norm(2));
    printf("T2[2]: %15.12f\n", T_2(2)({1,0},{0,1}).norm(2));
    printf("W1[2]: %15.12f\n", WABEF_2({1,0},{1,0}).norm(2));
    printf("W2[2]: %15.12f\n", WAMEI_2({1,0},{0,1}).norm(2));
    printf("W3[2]: %15.12f\n", WAMEI_2({1,0},{1,0}).norm(2));
    printf("W4[2]: %15.12f\n", WMNIJ_2({0,1},{0,1}).norm(2));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSD_TQ_N);
REGISTER_TASK(aquarius::cc::CCSD_TQ_N<double>,"ccsd(tq-n)");
