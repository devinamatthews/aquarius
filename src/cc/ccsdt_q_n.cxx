#include "ccsdt_q_n.hpp"

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
CCSDT_Q_N<U>::CCSDT_Q_N(const string& name, Config& config)
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
bool CCSDT_Q_N<U>::run(task::TaskDAG& dag, const Arena& arena)
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

    SpinorbitalTensor<U> DIJ(FMI);
    SpinorbitalTensor<U> DAB(FAE);
    SpinorbitalTensor<U> GIJKL(WMNIJ);
    SpinorbitalTensor<U> GAIBJ(WAMEI);
    SpinorbitalTensor<U> GABCD(WABEF);
    SpinorbitalTensor<U> GIJAK(WMNEJ);
    SpinorbitalTensor<U> GAIBC(WAMEF);

    ExcitationOperator<U,4> T_1("T^(1)", arena, occ, vrt);
    ExcitationOperator<U,4> T_2("T^(2)", arena, occ, vrt);
    ExcitationOperator<U,4> T_3("T^(3)", arena, occ, vrt);
    ExcitationOperator<U,4> T_4("T^(4)", arena, occ, vrt);
    ExcitationOperator<U,4> T_5("T^(5)", arena, occ, vrt);

    DeexcitationOperator<U,4> L_1("L^(1)", arena, occ, vrt);
    DeexcitationOperator<U,4> L_2("L^(2)", arena, occ, vrt);
    DeexcitationOperator<U,4> L_3("L^(3)", arena, occ, vrt);
    DeexcitationOperator<U,4> L_4("L^(4)", arena, occ, vrt);
    DeexcitationOperator<U,4> L_5("L^(5)", arena, occ, vrt);

    ExcitationOperator<U,4> DT_1("DT^(1)", arena, occ, vrt);
    ExcitationOperator<U,4> DT_2("DT^(2)", arena, occ, vrt);
    ExcitationOperator<U,4> DT_3("DT^(3)", arena, occ, vrt);
    ExcitationOperator<U,4> DT_4("DT^(4)", arena, occ, vrt);
    ExcitationOperator<U,4> DT_5("DT^(5)", arena, occ, vrt);

    DeexcitationOperator<U,4> DL_1("DL^(1)", arena, occ, vrt);
    DeexcitationOperator<U,4> DL_2("DL^(2)", arena, occ, vrt);
    DeexcitationOperator<U,4> DL_3("DL^(3)", arena, occ, vrt);
    DeexcitationOperator<U,4> DL_4("DL^(4)", arena, occ, vrt);
    DeexcitationOperator<U,4> DL_5("DL^(5)", arena, occ, vrt);

    SpinorbitalTensor<U> WABCEJK_1(WABCEJK);
    SpinorbitalTensor<U> WABMIJK_1(WABMIJK);

    SpinorbitalTensor<U> FMI_2(FMI);
    SpinorbitalTensor<U> FAE_2(FAE);
    SpinorbitalTensor<U> WAMEI_2(WAMEI);
    SpinorbitalTensor<U> WMNIJ_2(WMNIJ);
    SpinorbitalTensor<U> WABEF_2(WABEF);
    SpinorbitalTensor<U> WABEJ_2(WABEJ);
    SpinorbitalTensor<U> WAMIJ_2(WAMIJ);
    SpinorbitalTensor<U> WABMEJI_2(WABMEJI);
    SpinorbitalTensor<U> WABCEJK_2(WABCEJK);
    SpinorbitalTensor<U> WABMIJK_2(WABMIJK);
    SpinorbitalTensor<U> WAMNIJK_2(WAMNIJK);

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

    SpinorbitalTensor<U> FME_4(FME);
    SpinorbitalTensor<U> WAMEI_4(WAMEI);
    SpinorbitalTensor<U> WMNIJ_4(WMNIJ);
    SpinorbitalTensor<U> WABEF_4(WABEF);
    SpinorbitalTensor<U> WABEJ_4(WABEJ);
    SpinorbitalTensor<U> WAMIJ_4(WAMIJ);
    SpinorbitalTensor<U> WAMEF_4(WAMEF);
    SpinorbitalTensor<U> WMNEJ_4(WMNEJ);
    SpinorbitalTensor<U> WABMEJI_4(WABMEJI);
    SpinorbitalTensor<U> WABCEJK_4(WABCEJK);
    SpinorbitalTensor<U> WABMIJK_4(WABMIJK);
    SpinorbitalTensor<U> WAMNIJK_4(WAMNIJK);

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

    SpinorbitalTensor<U> DIJ_2(FMI);
    SpinorbitalTensor<U> DAB_2(FAE);
    SpinorbitalTensor<U> DAI_2(T(1));
    SpinorbitalTensor<U> GIJKL_2(WMNIJ);
    SpinorbitalTensor<U> GAIBJ_2(WAMEI);
    SpinorbitalTensor<U> GABCD_2(WABEF);
    SpinorbitalTensor<U> GIJAK_2(WMNEJ);
    SpinorbitalTensor<U> GAIBC_2(WAMEF);
    SpinorbitalTensor<U> GAIJK_2(WAMIJ);
    SpinorbitalTensor<U> GABCI_2(WABEJ);
    SpinorbitalTensor<U> GAIJBCD_2("G(aij,bcd)", arena, group, {vrt,occ}, {1,2}, {3,0});
    SpinorbitalTensor<U> GIJKABL_2("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1});
    SpinorbitalTensor<U> GIJKALM_2("G(ijk,alm)", arena, group, {vrt,occ}, {0,3}, {1,2});

    SpinorbitalTensor<U> DIJ_3(FMI);
    SpinorbitalTensor<U> DAB_3(FAE);
    SpinorbitalTensor<U> DAI_3(T(1));
    SpinorbitalTensor<U> GIJKL_3(WMNIJ);
    SpinorbitalTensor<U> GAIBJ_3(WAMEI);
    SpinorbitalTensor<U> GABCD_3(WABEF);
    SpinorbitalTensor<U> GIJAK_3(WMNEJ);
    SpinorbitalTensor<U> GAIBC_3(WAMEF);
    SpinorbitalTensor<U> GAIJK_3(WAMIJ);
    SpinorbitalTensor<U> GABCI_3(WABEJ);
    SpinorbitalTensor<U> GAIJBCD_3("G(aij,bcd)", arena, group, {vrt,occ}, {1,2}, {3,0});
    SpinorbitalTensor<U> GIJKABL_3("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1});
    SpinorbitalTensor<U> GIJKALM_3("G(ijk,alm)", arena, group, {vrt,occ}, {0,3}, {1,2});

    SpinorbitalTensor<U> GAIJBCD_4("G(aij,bcd)", arena, group, {vrt,occ}, {1,2}, {3,0});
    SpinorbitalTensor<U> GIJKABL_4("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1});

    SpinorbitalTensor<U> WABEJ_5(WABEJ);
    SpinorbitalTensor<U> WAMIJ_5(WAMIJ);
    ExcitationOperator<U,4> Z("Z", arena, occ, vrt);

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

    DT_1(4)["abcdijkl"]  =     WABCEJK["abcejk"]*T(2)[    "edil"];
    DT_1(4)["abcdijkl"] -=     WABMIJK["abmijk"]*T(2)[    "cdml"];
    DT_1(4)["abcdijkl"] +=       WABEJ[  "abej"]*T(3)[  "ecdikl"];
    DT_1(4)["abcdijkl"] -=       WAMIJ[  "amij"]*T(3)[  "bcdmkl"];
    DT_1(4)["abcdijkl"] += 0.5*WAMNIJK["amnijk"]*T(3)[  "bcdmnl"];
    DT_1(4)["abcdijkl"] -=     WABMEJI["abmeji"]*T(3)[  "ecdmkl"];

    T_1 = DT_1;
    T_1.weight(D);

    /***************************************************************************
     *
     * T^(2)
     *
     **************************************************************************/

    WABCEJK_1[  "abcejk"]  =  -0.5*WMNEF[  "mnef"]*T_1(4)["abcfmjkn"];
    WABMIJK_1[  "abmijk"]  =   0.5*WMNEF[  "mnef"]*T_1(4)["abefijkn"];

      DT_2(2)[    "abij"]  =  0.25*WMNEF[  "mnef"]*T_1(4)["abefijmn"];

      DT_2(3)[  "abcijk"]  =         FME[    "me"]*T_1(4)["abceijkm"];
      DT_2(3)[  "abcijk"] +=   0.5*WAMEF[  "amef"]*T_1(4)["efbcimjk"];
      DT_2(3)[  "abcijk"] -=   0.5*WMNEJ[  "mnek"]*T_1(4)["abecijmn"];

      DT_2(4)["abcdijkl"]  =   WABCEJK_1["abcejk"]*  T(2)[    "edil"];
      DT_2(4)["abcdijkl"] -=   WABMIJK_1["abmijk"]*  T(2)[    "cdml"];
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

      WABCEJK[  "abcejk"]  =         WABEF[  "abef"]*  T(2)[    "fcjk"];
      WABCEJK[  "abcejk"] -=         WAMEI[  "amek"]*  T(2)[    "bcjm"];
      WABCEJK[  "abcejk"] -=           FME[    "me"]*  T(3)[  "abcmjk"];
      WABCEJK[  "abcejk"] +=     0.5*WMNEJ[  "mnej"]*  T(3)[  "abcmnk"];
      WABCEJK[  "abcejk"] +=         WAMEF[  "amef"]*  T(3)[  "fbcmjk"];

      WABMIJK[  "abmijk"]  =         WAMEI[  "amek"]*  T(2)[    "ebij"];
      WABMIJK[  "abmijk"] -=         WMNIJ[  "nmjk"]*  T(2)[    "abin"];
      WABMIJK[  "abmijk"] +=           FME[    "me"]*  T(3)[  "abeijk"];
      WABMIJK[  "abmijk"] +=     0.5*WAMEF[  "bmef"]*  T(3)[  "aefijk"];
      WABMIJK[  "abmijk"] +=         WMNEJ[  "nmek"]*  T(3)[  "abeijn"];

      WAMNIJK[  "amnijk"]  =         WMNEJ[  "mnek"]*  T(2)[    "aeij"];
      WAMNIJK[  "amnijk"] +=     0.5*WMNEF[  "mnef"]*  T(3)[  "aefijk"];

      WABMEJI[  "abmeji"]  =         WAMEF[  "amef"]*  T(2)[    "bfji"];
      WABMEJI[  "abmeji"] -=         WMNEJ[  "nmei"]*  T(2)[    "abnj"];
      WABMEJI[  "abmeji"] +=         WMNEF[  "mnef"]*  T(3)[  "abfnji"];

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
    WABCEJK_2[  "abcejk"] -=     0.5*WMNEF[  "mnef"]*T_2(4)["abcfmjkn"];

    WABMIJK_2[  "abmijk"]  =       WAMEI_2[  "amek"]*  T(2)[    "ebij"];
    WABMIJK_2[  "abmijk"] -=   0.5*WMNIJ_2[  "nmjk"]*  T(2)[    "abin"];
    WABMIJK_2[  "abmijk"] +=     0.5*WAMEF[  "bmef"]*T_2(3)[  "aefijk"];
    WABMIJK_2[  "abmijk"] +=     0.5*WMNEF[  "mnef"]*T_2(4)["abefijkn"];

    WAMNIJK_2[  "amnijk"]  =     0.5*WMNEF[  "mnef"]*T_2(3)[  "aefijk"];

      DT_3(1)[      "ai"]  =           FME[    "me"]*T_2(2)[    "aeim"];
      DT_3(1)[      "ai"] +=     0.5*WAMEF[  "amef"]*T_2(2)[    "efim"];
      DT_3(1)[      "ai"] -=     0.5*WMNEJ[  "mnei"]*T_2(2)[    "eamn"];
      DT_3(1)[      "ai"] +=    0.25*WMNEF[  "mnef"]*T_2(3)[  "aefimn"];

      DT_3(2)[    "abij"]  =         FAE_2[    "af"]*  T(2)[    "fbij"];
      DT_3(2)[    "abij"] -=         FMI_2[    "ni"]*  T(2)[    "abnj"];
      DT_3(2)[    "abij"] +=           FAE[    "af"]*T_2(2)[    "fbij"];
      DT_3(2)[    "abij"] -=           FMI[    "ni"]*T_2(2)[    "abnj"];
      DT_3(2)[    "abij"] +=     0.5*WABEF[  "abef"]*T_2(2)[    "efij"];
      DT_3(2)[    "abij"] +=     0.5*WMNIJ[  "mnij"]*T_2(2)[    "abmn"];
      DT_3(2)[    "abij"] +=         WAMEI[  "amei"]*T_2(2)[    "ebjm"];
      DT_3(2)[    "abij"] +=           FME[    "me"]*T_2(3)[  "abeijm"];
      DT_3(2)[    "abij"] +=     0.5*WAMEF[  "bmef"]*T_2(3)[  "aefijm"];
      DT_3(2)[    "abij"] -=     0.5*WMNEJ[  "mnej"]*T_2(3)[  "abeinm"];
      DT_3(2)[    "abij"] +=    0.25*WMNEF[  "mnef"]*T_2(4)["abefijmn"];

      DT_3(3)[  "abcijk"]  =       WABEJ_2[  "bcek"]*  T(2)[    "aeij"];
      DT_3(3)[  "abcijk"] -=       WAMIJ_2[  "bmjk"]*  T(2)[    "acim"];
      DT_3(3)[  "abcijk"] +=         FAE_2[    "ce"]*  T(3)[  "abeijk"];
      DT_3(3)[  "abcijk"] -=         FMI_2[    "mk"]*  T(3)[  "abcijm"];
      DT_3(3)[  "abcijk"] +=   0.5*WABEF_2[  "abef"]*  T(3)[  "efcijk"];
      DT_3(3)[  "abcijk"] +=   0.5*WMNIJ_2[  "mnij"]*  T(3)[  "abcmnk"];
      DT_3(3)[  "abcijk"] +=       WAMEI_2[  "amei"]*  T(3)[  "ebcjmk"];
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
      DT_3(4)["abcdijkl"] +=       WABEJ_2[  "abej"]*  T(3)[  "ecdikl"];
      DT_3(4)["abcdijkl"] -=       WAMIJ_2[  "amij"]*  T(3)[  "bcdmkl"];
      DT_3(4)["abcdijkl"] += 0.5*WAMNIJK_2["amnijk"]*  T(3)[  "bcdmnl"];
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

      WAMIJ_2[    "amij"] +=           FME[    "me"]*T_2(2)[    "aeij"];
      WABEJ_2[    "abej"] -=           FME[    "me"]*T_2(2)[    "abmj"];

    T_3 = DT_3;
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
    WABCEJK_3[  "abcejk"] -=     0.5*WMNEF[  "mnef"]*T_3(4)["abcfmjkn"];

    WABMIJK_3[  "abmijk"]  =       WAMEI_3[  "amek"]*  T(2)[    "ebij"];
    WABMIJK_3[  "abmijk"] -=   0.5*WMNIJ_3[  "nmjk"]*  T(2)[    "abin"];
    WABMIJK_3[  "abmijk"] +=   0.5*WAMEF_3[  "bmef"]*  T(3)[  "aefijk"];
    WABMIJK_3[  "abmijk"] +=     0.5*WAMEF[  "bmef"]*T_3(3)[  "aefijk"];
    WABMIJK_3[  "abmijk"] +=     0.5*WMNEF[  "mnef"]*T_3(4)["abefijkn"];

    WAMNIJK_3[  "amnijk"]  =       WMNEJ_3[  "mnek"]*  T(2)[    "aeij"];
    WAMNIJK_3[  "amnijk"] +=     0.5*WMNEF[  "mnef"]*T_3(3)[  "aefijk"];

    WABMEJI_3[  "abmeji"]  =       WAMEF_3[  "amef"]*  T(2)[    "bfji"];
    WABMEJI_3[  "abmeji"] -=       WMNEJ_3[  "nmei"]*  T(2)[    "abnj"];

      DT_4(1)[      "ai"]  =           FAE[    "ae"]*T_3(1)[      "ei"];
      DT_4(1)[      "ai"] -=           FMI[    "mi"]*T_3(1)[      "am"];
      DT_4(1)[      "ai"] -=         WAMEI[  "amei"]*T_3(1)[      "em"];
      DT_4(1)[      "ai"] +=           FME[    "me"]*T_3(2)[    "aeim"];
      DT_4(1)[      "ai"] +=     0.5*WAMEF[  "amef"]*T_3(2)[    "efim"];
      DT_4(1)[      "ai"] -=     0.5*WMNEJ[  "mnei"]*T_3(2)[    "eamn"];
      DT_4(1)[      "ai"] +=    0.25*WMNEF[  "mnef"]*T_3(3)[  "aefimn"];

      DT_4(2)[    "abij"]  =         FAE_3[    "af"]*  T(2)[    "fbij"];
      DT_4(2)[    "abij"] -=         FMI_3[    "ni"]*  T(2)[    "abnj"];
      DT_4(2)[    "abij"] +=         FME_3[    "me"]*  T(3)[  "abeijm"];
      DT_4(2)[    "abij"] +=         WABEJ[  "abej"]*T_3(1)[      "ei"];
      DT_4(2)[    "abij"] -=         WAMIJ[  "amij"]*T_3(1)[      "bm"];
      DT_4(2)[    "abij"] +=           FAE[    "af"]*T_3(2)[    "fbij"];
      DT_4(2)[    "abij"] -=           FMI[    "ni"]*T_3(2)[    "abnj"];
      DT_4(2)[    "abij"] +=     0.5*WABEF[  "abef"]*T_3(2)[    "efij"];
      DT_4(2)[    "abij"] +=     0.5*WMNIJ[  "mnij"]*T_3(2)[    "abmn"];
      DT_4(2)[    "abij"] +=         WAMEI[  "amei"]*T_3(2)[    "ebjm"];
      DT_4(2)[    "abij"] +=           FME[    "me"]*T_3(3)[  "abeijm"];
      DT_4(2)[    "abij"] +=     0.5*WAMEF[  "bmef"]*T_3(3)[  "aefijm"];
      DT_4(2)[    "abij"] -=     0.5*WMNEJ[  "mnej"]*T_3(3)[  "abeinm"];
      DT_4(2)[    "abij"] +=    0.25*WMNEF[  "mnef"]*T_3(4)["abefijmn"];

        FMI_3[      "mi"] +=           FME[    "me"]*T_3(1)[      "ei"];
        FAE_3[      "ae"] -=           FME[    "me"]*T_3(1)[      "am"];

      DT_4(3)[  "abcijk"]  =       WABEJ_3[  "bcek"]*  T(2)[    "aeij"];
      DT_4(3)[  "abcijk"] -=       WAMIJ_3[  "bmjk"]*  T(2)[    "acim"];
      DT_4(3)[  "abcijk"] +=         FAE_3[    "ce"]*  T(3)[  "abeijk"];
      DT_4(3)[  "abcijk"] -=         FMI_3[    "mk"]*  T(3)[  "abcijm"];
      DT_4(3)[  "abcijk"] +=   0.5*WABEF_3[  "abef"]*  T(3)[  "efcijk"];
      DT_4(3)[  "abcijk"] +=   0.5*WMNIJ_3[  "mnij"]*  T(3)[  "abcmnk"];
      DT_4(3)[  "abcijk"] +=       WAMEI_3[  "amei"]*  T(3)[  "ebcjmk"];
      DT_4(3)[  "abcijk"] +=         WABEJ[  "bcek"]*T_3(2)[    "aeij"];
      DT_4(3)[  "abcijk"] -=         WAMIJ[  "bmjk"]*T_3(2)[    "acim"];
      DT_4(3)[  "abcijk"] +=           FAE[    "ce"]*T_3(3)[  "abeijk"];
      DT_4(3)[  "abcijk"] -=           FMI[    "mk"]*T_3(3)[  "abcijm"];
      DT_4(3)[  "abcijk"] +=     0.5*WABEF[  "abef"]*T_3(3)[  "efcijk"];
      DT_4(3)[  "abcijk"] +=     0.5*WMNIJ[  "mnij"]*T_3(3)[  "abcmnk"];
      DT_4(3)[  "abcijk"] +=         WAMEI[  "amei"]*T_3(3)[  "ebcjmk"];
      DT_4(3)[  "abcijk"] +=           FME[    "me"]*T_3(4)["abceijkm"];
      DT_4(3)[  "abcijk"] +=     0.5*WAMEF[  "amef"]*T_3(4)["efbcimjk"];
      DT_4(3)[  "abcijk"] -=     0.5*WMNEJ[  "mnek"]*T_3(4)["abecijmn"];

      WABEJ_3[    "abej"] -=         FME_3[    "me"]*  T(2)[    "abmj"];

      DT_4(4)["abcdijkl"]  =     WABCEJK_1["abcejk"]*T_2(2)[    "edil"];
      DT_4(4)["abcdijkl"] -=     WABMIJK_1["abmijk"]*T_2(2)[    "cdml"];
      DT_4(4)["abcdijkl"] +=         FAE_2[    "ae"]*T_1(4)["ebcdijkl"];
      DT_4(4)["abcdijkl"] -=         FMI_2[    "mi"]*T_1(4)["abcdmjkl"];
      DT_4(4)["abcdijkl"] +=   0.5*WABEF_2[  "abef"]*T_1(4)["efcdijkl"];
      DT_4(4)["abcdijkl"] +=   0.5*WMNIJ_2[  "mnij"]*T_1(4)["abcdmnkl"];
      DT_4(4)["abcdijkl"] +=       WAMEI_2[  "amei"]*T_1(4)["ebcdjmkl"];
      DT_4(4)["abcdijkl"] +=     WABCEJK_3["abcejk"]*  T(2)[    "edil"];
      DT_4(4)["abcdijkl"] -=     WABMIJK_3["abmijk"]*  T(2)[    "cdml"];
      DT_4(4)["abcdijkl"] +=       WABEJ_3[  "abej"]*  T(3)[  "ecdikl"];
      DT_4(4)["abcdijkl"] -=       WAMIJ_3[  "amij"]*  T(3)[  "bcdmkl"];
      DT_4(4)["abcdijkl"] += 0.5*WAMNIJK_3["amnijk"]*  T(3)[  "bcdmnl"];
      DT_4(4)["abcdijkl"] -=     WABMEJI_3["abmeji"]*  T(3)[  "ecdmkl"];
      DT_4(4)["abcdijkl"] +=       WABCEJK["abcejk"]*T_3(2)[    "edil"];
      DT_4(4)["abcdijkl"] -=       WABMIJK["abmijk"]*T_3(2)[    "cdml"];
      DT_4(4)["abcdijkl"] +=         WABEJ[  "abej"]*T_3(3)[  "ecdikl"];
      DT_4(4)["abcdijkl"] -=         WAMIJ[  "amij"]*T_3(3)[  "bcdmkl"];
      DT_4(4)["abcdijkl"] +=   0.5*WAMNIJK["amnijk"]*T_3(3)[  "bcdmnl"];
      DT_4(4)["abcdijkl"] -=       WABMEJI["abmeji"]*T_3(3)[  "ecdmkl"];
      DT_4(4)["abcdijkl"] +=           FAE[    "ae"]*T_3(4)["ebcdijkl"];
      DT_4(4)["abcdijkl"] -=           FMI[    "mi"]*T_3(4)["abcdmjkl"];
      DT_4(4)["abcdijkl"] +=     0.5*WABEF[  "abef"]*T_3(4)["efcdijkl"];
      DT_4(4)["abcdijkl"] +=     0.5*WMNIJ[  "mnij"]*T_3(4)["abcdmnkl"];
      DT_4(4)["abcdijkl"] +=         WAMEI[  "amei"]*T_3(4)["ebcdjmkl"];

      WAMIJ_3[    "amij"] +=           FME[    "me"]*T_3(2)[    "aeij"];
      WABEJ_3[    "abej"] -=           FME[    "me"]*T_3(2)[    "abmj"];

    T_4 = DT_4;
    T_4.weight(D);

    /***************************************************************************
     *
     * T^(5)
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

    WAMNIJK_2[  "amnijk"]  =         WMNEJ[  "mnek"]*T_2(2)[    "aeij"];
    WAMNIJK_2[  "amnijk"] +=     0.5*WMNEF[  "mnef"]*T_2(3)[  "aefijk"];

    WABMEJI_2[  "abmeji"]  =         WAMEF[  "amef"]*T_2(2)[    "bfji"];
    WABMEJI_2[  "abmeji"] -=         WMNEJ[  "nmei"]*T_2(2)[    "abnj"];
    WABMEJI_2[  "abmeji"] +=     0.5*WMNEF[  "mnef"]*T_2(3)[  "abfnji"];

        FME_4[      "me"]  =         WMNEF[  "mnef"]*T_4(1)[      "fn"];

      WMNEJ_4[    "mnej"]  =         WMNEF[  "mnef"]*T_4(1)[      "fj"];

      WAMEF_4[    "amef"]  =        -WMNEF[  "nmef"]*T_4(1)[      "an"];

      WAMIJ_4[    "amij"]  =         WMNEJ[  "nmej"]*T_4(2)[    "aein"];
      WAMIJ_4[    "amij"] +=     0.5*WAMEF[  "amef"]*T_4(2)[    "efij"];
      WAMIJ_4[    "amij"] +=     0.5*WMNEF[  "mnef"]*T_4(3)[  "aefijn"];
      WAMIJ_4[    "amij"] +=         WAMEI[  "amej"]*T_4(1)[      "ei"];
      WAMIJ_4[    "amij"] -=         WMNIJ[  "nmij"]*T_4(1)[      "an"];
      WAMIJ_4[    "amij"] +=         FME_4[    "me"]*  T(2)[    "aeij"];

      WABEJ_4[    "abej"]  =         WAMEF[  "amef"]*T_4(2)[    "fbmj"];
      WABEJ_4[    "abej"] +=     0.5*WMNEJ[  "mnej"]*T_4(2)[    "abmn"];
      WABEJ_4[    "abej"] -=     0.5*WMNEF[  "mnef"]*T_4(3)[  "afbmnj"];
      WABEJ_4[    "abej"] -=         WAMEI[  "amej"]*T_4(1)[      "bm"];
      WABEJ_4[    "abej"] +=         WABEF[  "abef"]*T_4(1)[      "fj"];
      WABEJ_4[    "abej"] -=         FME_4[    "me"]*  T(2)[    "abmj"];

      WABEF_4[    "abef"]  =     0.5*WMNEF[  "mnef"]*T_4(2)[    "abmn"];
      WABEF_4[    "abef"] -=         WAMEF[  "amef"]*T_4(1)[      "bm"];

      WAMEI_4[    "amei"]  =        -WMNEF[  "mnef"]*T_4(2)[    "afin"];
      WAMEI_4[    "amei"] +=         WAMEF[  "amef"]*T_4(1)[      "fi"];
      WAMEI_4[    "amei"] -=         WMNEJ[  "nmei"]*T_4(1)[      "an"];

      WMNIJ_4[    "mnij"]  =     0.5*WMNEF[  "mnef"]*T_4(2)[    "efij"];
      WMNIJ_4[    "mnij"] +=         WMNEJ[  "mnej"]*T_4(1)[      "ei"];

    WABCEJK_4[  "abcejk"]  =   0.5*WABEF_4[  "abef"]*  T(2)[    "fcjk"];
    WABCEJK_4[  "abcejk"] -=     0.5*WMNEF[  "mnef"]*T_4(4)["abcfmjkn"];

    WABMIJK_4[  "abmijk"]  =       WAMEI_4[  "amek"]*  T(2)[    "ebij"];
    WABMIJK_4[  "abmijk"] -=   0.5*WMNIJ_4[  "nmjk"]*  T(2)[    "abin"];
    WABMIJK_4[  "abmijk"] +=   0.5*WAMEF_4[  "bmef"]*  T(3)[  "aefijk"];
    WABMIJK_4[  "abmijk"] +=     0.5*WAMEF[  "bmef"]*T_4(3)[  "aefijk"];
    WABMIJK_4[  "abmijk"] +=     0.5*WMNEF[  "mnef"]*T_4(4)["abefijkn"];

    WAMNIJK_4[  "amnijk"]  =       WMNEJ_4[  "mnek"]*  T(2)[    "aeij"];
    WAMNIJK_4[  "amnijk"] +=     0.5*WMNEF[  "mnef"]*T_4(3)[  "aefijk"];

    WABMEJI_4[  "abmeji"]  =       WAMEF_4[  "amef"]*  T(2)[    "bfji"];
    WABMEJI_4[  "abmeji"] -=       WMNEJ_4[  "nmei"]*  T(2)[    "abnj"];

      DT_5(4)["abcdijkl"]  =     WABCEJK_1["abcejk"]*T_3(2)[    "edil"];
      DT_5(4)["abcdijkl"] -=     WABMIJK_1["abmijk"]*T_3(2)[    "cdml"];
      DT_5(4)["abcdijkl"] +=     WABCEJK_2["abcejk"]*T_2(2)[    "edil"];
      DT_5(4)["abcdijkl"] -=     WABMIJK_2["abmijk"]*T_2(2)[    "cdml"];
      DT_5(4)["abcdijkl"] +=       WABEJ_2[  "abej"]*T_2(3)[  "ecdikl"];
      DT_5(4)["abcdijkl"] -=       WAMIJ_2[  "amij"]*T_2(3)[  "bcdmkl"];
      DT_5(4)["abcdijkl"] += 0.5*WAMNIJK_2["amnijk"]*T_2(3)[  "bcdmnl"];
      DT_5(4)["abcdijkl"] -=     WABMEJI_2["abmeji"]*T_2(3)[  "ecdmkl"];
      DT_5(4)["abcdijkl"] +=         FAE_2[    "ae"]*T_2(4)["ebcdijkl"];
      DT_5(4)["abcdijkl"] -=         FMI_2[    "mi"]*T_2(4)["abcdmjkl"];
      DT_5(4)["abcdijkl"] +=   0.5*WABEF_2[  "abef"]*T_2(4)["efcdijkl"];
      DT_5(4)["abcdijkl"] +=   0.5*WMNIJ_2[  "mnij"]*T_2(4)["abcdmnkl"];
      DT_5(4)["abcdijkl"] +=       WAMEI_2[  "amei"]*T_2(4)["ebcdjmkl"];
      DT_5(4)["abcdijkl"] +=         FAE_3[    "ae"]*T_1(4)["ebcdijkl"];
      DT_5(4)["abcdijkl"] -=         FMI_3[    "mi"]*T_1(4)["abcdmjkl"];
      DT_5(4)["abcdijkl"] +=   0.5*WABEF_3[  "abef"]*T_1(4)["efcdijkl"];
      DT_5(4)["abcdijkl"] +=   0.5*WMNIJ_3[  "mnij"]*T_1(4)["abcdmnkl"];
      DT_5(4)["abcdijkl"] +=       WAMEI_3[  "amei"]*T_1(4)["ebcdjmkl"];
      DT_5(4)["abcdijkl"] +=     WABCEJK_4["abcejk"]*  T(2)[    "edil"];
      DT_5(4)["abcdijkl"] -=     WABMIJK_4["abmijk"]*  T(2)[    "cdml"];
      DT_5(4)["abcdijkl"] +=       WABEJ_4[  "abej"]*  T(3)[  "ecdikl"];
      DT_5(4)["abcdijkl"] -=       WAMIJ_4[  "amij"]*  T(3)[  "bcdmkl"];
      DT_5(4)["abcdijkl"] += 0.5*WAMNIJK_4["amnijk"]*  T(3)[  "bcdmnl"];
      DT_5(4)["abcdijkl"] -=     WABMEJI_4["abmeji"]*  T(3)[  "ecdmkl"];
      DT_5(4)["abcdijkl"] +=       WABCEJK["abcejk"]*T_4(2)[    "edil"];
      DT_5(4)["abcdijkl"] -=       WABMIJK["abmijk"]*T_4(2)[    "cdml"];
      DT_5(4)["abcdijkl"] +=         WABEJ[  "abej"]*T_4(3)[  "ecdikl"];
      DT_5(4)["abcdijkl"] -=         WAMIJ[  "amij"]*T_4(3)[  "bcdmkl"];
      DT_5(4)["abcdijkl"] +=   0.5*WAMNIJK["amnijk"]*T_4(3)[  "bcdmnl"];
      DT_5(4)["abcdijkl"] -=       WABMEJI["abmeji"]*T_4(3)[  "ecdmkl"];
      DT_5(4)["abcdijkl"] +=           FAE[    "ae"]*T_4(4)["ebcdijkl"];
      DT_5(4)["abcdijkl"] -=           FMI[    "mi"]*T_4(4)["abcdmjkl"];
      DT_5(4)["abcdijkl"] +=     0.5*WABEF[  "abef"]*T_4(4)["efcdijkl"];
      DT_5(4)["abcdijkl"] +=     0.5*WMNIJ[  "mnij"]*T_4(4)["abcdmnkl"];
      DT_5(4)["abcdijkl"] +=         WAMEI[  "amei"]*T_4(4)["ebcdjmkl"];

      WAMIJ_4[    "amij"] +=           FME[    "me"]*T_4(2)[    "aeij"];
      WABEJ_4[    "abej"] -=           FME[    "me"]*T_4(2)[    "abmj"];

    T_5 = DT_5;
    T_5.weight(D);

    /***************************************************************************
     *
     * L^(1)
     *
     **************************************************************************/

    DL_1(4)["ijklabcd"]  = WMNEF["ijab"]*L(2)[  "klcd"];
    DL_1(4)["ijklabcd"] +=   FME[  "ia"]*L(3)["jklbcd"];
    DL_1(4)["ijklabcd"] += WAMEF["elcd"]*L(3)["ijkabe"];
    DL_1(4)["ijklabcd"] -= WMNEJ["ijam"]*L(3)["mklbcd"];

    L_1 = DL_1;
    L_1.weight(D);

    /***************************************************************************
     *
     * L^(2)
     *
     **************************************************************************/

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

      DL_2(1)[      "ia"]  =               WMNEF[    "miea"]*    DAI_1[      "em"];
      DL_2(1)[      "ia"] +=  (1.0/ 2.0)*  WMNEF[    "imef"]*  GABCI_1[    "efam"];
      DL_2(1)[      "ia"] -=  (1.0/ 2.0)*  WMNEF[    "mnea"]*  GAIJK_1[    "eimn"];
      DL_2(1)[      "ia"] -=  (1.0/ 2.0)*  WABEF[    "efga"]*  GAIBC_1[    "gief"];
      DL_2(1)[      "ia"] +=               WAMEI[    "eifm"]*  GAIBC_1[    "fmea"];
      DL_2(1)[      "ia"] -=               WAMEI[    "eman"]*  GIJAK_1[    "inem"];
      DL_2(1)[      "ia"] +=  (1.0/ 2.0)*  WMNIJ[    "imno"]*  GIJAK_1[    "noam"];
      DL_2(1)[      "ia"] -=  (1.0/ 2.0)*  WAMEF[    "gief"]*  GABCD_1[    "efga"];
      DL_2(1)[      "ia"] +=               WAMEF[    "fmea"]*  GAIBJ_1[    "eifm"];
      DL_2(1)[      "ia"] -=               WMNEJ[    "inem"]*  GAIBJ_1[    "eman"];
      DL_2(1)[      "ia"] +=  (1.0/ 2.0)*  WMNEJ[    "noam"]*  GIJKL_1[    "imno"];

      DL_2(2)[    "ijab"]  =              -WAMEF[    "fiae"]*  GAIBC_1[    "ejbf"];
      DL_2(2)[    "ijab"] -=               WMNEJ[    "ijem"]*  GAIBC_1[    "emab"];
      DL_2(2)[    "ijab"] -=               WAMEF[    "emab"]*  GIJAK_1[    "ijem"];
      DL_2(2)[    "ijab"] -=               WMNEJ[    "niam"]*  GIJAK_1[    "mjbn"];
      DL_2(2)[    "ijab"] +=  (1.0/ 2.0)*  WMNEF[    "ijef"]*  GABCD_1[    "efab"];
      DL_2(2)[    "ijab"] +=               WMNEF[    "imea"]*  GAIBJ_1[    "ejbm"];
      DL_2(2)[    "ijab"] +=  (1.0/ 2.0)*  WMNEF[    "mnab"]*  GIJKL_1[    "ijmn"];
      DL_2(2)[    "ijab"] +=  (1.0/12.0)*WABCEJK[  "efgbmn"]*   L_1(4)["ijmnaefg"];
      DL_2(2)[    "ijab"] -=  (1.0/12.0)*WABMIJK[  "efimno"]*   L_1(4)["mnojefab"];

      DL_2(3)[  "ijkabc"]  =               WMNEF[    "ijae"]*  GAIBC_1[    "ekbc"];
      DL_2(3)[  "ijkabc"] -=               WMNEF[    "mkbc"]*  GIJAK_1[    "ijam"];
      DL_2(3)[  "ijkabc"] +=  (1.0/ 2.0)*  WMNEF[    "mnbc"]*GIJKALM_1[  "ijkamn"];
      DL_2(3)[  "ijkabc"] -=               WAMEF[    "embc"]*GIJKABL_1[  "ijkaem"];
      DL_2(3)[  "ijkabc"] +=  (1.0/ 2.0)*  WABEJ[    "efcm"]*   L_1(4)["ijkmabef"];
      DL_2(3)[  "ijkabc"] -=  (1.0/ 2.0)*  WAMIJ[    "eknm"]*   L_1(4)["ijmnabce"];
      DL_2(3)[  "ijkabc"] -=  (1.0/ 4.0)*WABMEJI[  "efkcnm"]*   L_1(4)["ijmnabef"];
      DL_2(3)[  "ijkabc"] +=  (1.0/ 6.0)*WAMNIJK[  "eijmno"]*   L_1(4)["mnokeabc"];

      DL_2(4)["ijklabcd"]  =               WMNEF[    "ijae"]*GAIJBCD_1[  "eklbcd"];
      DL_2(4)["ijklabcd"] -=               WMNEF[    "mlcd"]*GIJKABL_1[  "ijkabm"];
      DL_2(4)["ijklabcd"] +=                 FAE[      "ea"]*   L_1(4)["ijklebcd"];
      DL_2(4)["ijklabcd"] -=                 FMI[      "im"]*   L_1(4)["mjklabcd"];
      DL_2(4)["ijklabcd"] +=  (1.0/ 2.0)*  WABEF[    "efab"]*   L_1(4)["ijklefcd"];
      DL_2(4)["ijklabcd"] +=  (1.0/ 2.0)*  WMNIJ[    "ijmn"]*   L_1(4)["mnklabcd"];
      DL_2(4)["ijklabcd"] +=               WAMEI[    "eiam"]*   L_1(4)["mjklbecd"];

    L_2 = DL_2;
    L_2.weight(D);

    /***************************************************************************
     *
     * L^(3)
     *
     **************************************************************************/

        DIJ_2[      "ij"]  =  (1.0/ 12.0)*     T(3)[  "efgjmn"]*   L_2(3)[  "imnefg"];
        DIJ_2[      "ij"] +=  (1.0/144.0)*   T_1(4)["efghjmno"]*   L_1(4)["imnoefgh"];
        DIJ_2[      "ij"] +=  (1.0/ 12.0)*   T_2(3)[  "efgjmn"]*     L(3)[  "imnefg"];

        DAB_2[      "ab"]  = -(1.0/ 12.0)*     T(3)[  "aefmno"]*   L_2(3)[  "mnobef"];
        DAB_2[      "ab"] -=  (1.0/144.0)*   T_1(4)["aefgmnop"]*   L_1(4)["mnopbefg"];
        DAB_2[      "ab"] -=  (1.0/ 12.0)*   T_2(3)[  "aefmno"]*     L(3)[  "mnobef"];

        GIJAK[    "ijak"]  =  (1.0/  2.0)*     T(2)[    "efkm"]*     L(3)[  "ijmaef"];
        GAIBC[    "aibc"]  = -(1.0/  2.0)*     T(2)[    "aemn"]*     L(3)[  "minbce"];

    GIJKABL_2[  "ijkabl"]  =  (1.0/  2.0)*     T(2)[    "eflm"]*   L_2(4)["ijkmabef"];
    GAIJBCD_2[  "aijbcd"]  = -(1.0/  2.0)*     T(2)[    "aemn"]*   L_2(4)["mijnbcde"];
    GIJKALM_2[  "ijkalm"]  =  (1.0/  6.0)*     T(3)[  "efglmn"]*   L_2(4)["ijknaefg"];

      GIJAK_2[    "ijak"]  =  (1.0/ 12.0)*     T(3)[  "efgkmn"]*   L_2(4)["ijmnaefg"];
      GIJAK_2[    "ijak"] +=  (1.0/  2.0)*     T(2)[    "efkm"]*   L_2(3)[  "ijmaef"];
      GIJAK_2[    "ijak"] +=  (1.0/  2.0)*   T_2(2)[    "efkm"]*     L(3)[  "ijmaef"];

      GAIBC_2[    "aibc"]  = -(1.0/ 12.0)*     T(3)[  "aefmno"]*   L_2(4)["minobcef"];
      GAIBC_2[    "aibc"] -=  (1.0/  2.0)*   T_2(2)[    "aemn"]*     L(3)[  "minbce"];

      GAIJK_2[    "aijk"]  =  (1.0/ 12.0)*   T_2(4)["aefgjkmn"]*     L(3)[  "imnefg"];
      GAIJK_2[    "aijk"] +=  (1.0/  4.0)*     T(3)[  "efamnj"]*GIJKABL_2[  "mniefk"];
      GAIJK_2[    "aijk"] +=  (1.0/  6.0)*     T(3)[  "efgjkm"]*GAIJBCD_2[  "aimefg"];

      GABCI_2[    "abci"]  = -(1.0/ 12.0)*   T_2(4)["abefmino"]*     L(3)[  "mnocef"];
      GABCI_2[    "abci"] +=  (1.0/  4.0)*     T(3)[  "efbmni"]*GAIJBCD_2[  "amncef"];
      GABCI_2[    "abci"] +=  (1.0/  6.0)*     T(3)[  "eabmno"]*GIJKABL_2[  "mnoeci"];

        DAI_2[      "ai"]  =  (1.0/  4.0)*     T(3)[  "aefimn"]*   L_2(2)[    "mnef"];
        DAI_2[      "ai"] +=  (1.0/ 36.0)*   T_2(4)["aefgimno"]*     L(3)[  "mnoefg"];
        DAI_2[      "ai"] +=  (1.0/  4.0)*   T_2(3)[  "aefimn"]*     L(2)[    "mnef"];
        DAI_2[      "ai"] -=  (1.0/  2.0)*     T(2)[    "eamn"]*  GIJAK_2[    "mnei"];
        DAI_2[      "ai"] +=  (1.0/  2.0)*     T(2)[    "efim"]*  GAIBC_2[    "amef"];

      GAIBC_2[    "aibc"] -=  (1.0/  2.0)*     T(2)[    "aemn"]*   L_2(3)[  "minbce"];

      GABCD_2[    "abcd"]  =  (1.0/  6.0)*     T(3)[  "abemno"]*   L_2(3)[  "mnocde"];
      GABCD_2[    "abcd"] +=  (1.0/ 48.0)*   T_1(4)["abefmnop"]*   L_1(4)["mnopcdef"];
      GABCD_2[    "abcd"] +=  (1.0/  6.0)*   T_2(3)[  "abemno"]*     L(3)[  "mnocde"];
      GABCD_2[    "abcd"] -=  (1.0/  4.0)*     T(2)[    "bemn"]*GAIJBCD_2[  "amncde"];

      GAIBJ_2[    "aibj"]  = -(1.0/  4.0)*     T(3)[  "aefjmn"]*   L_2(3)[  "imnbef"];
      GAIBJ_2[    "aibj"] -=  (1.0/ 36.0)*   T_1(4)["aefgjmno"]*   L_1(4)["imnobefg"];
      GAIBJ_2[    "aibj"] -=  (1.0/  4.0)*   T_2(3)[  "aefjmn"]*     L(3)[  "imnbef"];
      GAIBJ_2[    "aibj"] +=  (1.0/  2.0)*     T(2)[    "eamn"]*GIJKABL_2[  "imnbej"];

      GIJKL_2[    "ijkl"]  =  (1.0/  6.0)*     T(3)[  "efgklm"]*   L_2(3)[  "ijmefg"];
      GIJKL_2[    "ijkl"] +=  (1.0/ 48.0)*   T_1(4)["efghklmn"]*   L_1(4)["ijmnefgh"];
      GIJKL_2[    "ijkl"] +=  (1.0/  6.0)*   T_2(3)[  "efgklm"]*     L(3)[  "ijmefg"];
      GIJKL_2[    "ijkl"] +=  (1.0/  4.0)*     T(2)[    "efmk"]*GIJKABL_2[  "mijefl"];

      DL_3(1)[      "ia"]  =                  FAE_2[      "ea"]*     L(1)[      "ie"];
      DL_3(1)[      "ia"] -=                  FMI_2[      "im"]*     L(1)[      "ma"];
      DL_3(1)[      "ia"] -=                WAMEI_2[    "eiam"]*     L(1)[      "me"];
      DL_3(1)[      "ia"] +=  (1.0/  2.0)*  WABEJ_2[    "efam"]*     L(2)[    "imef"];
      DL_3(1)[      "ia"] -=  (1.0/  2.0)*  WAMIJ_2[    "eimn"]*     L(2)[    "mnea"];
      DL_3(1)[      "ia"] -=  (1.0/  2.0)*  WABEF_2[    "efga"]*    GAIBC[    "gief"];
      DL_3(1)[      "ia"] +=                WAMEI_2[    "eifm"]*    GAIBC[    "fmea"];
      DL_3(1)[      "ia"] -=                WAMEI_2[    "eman"]*    GIJAK[    "inem"];
      DL_3(1)[      "ia"] +=  (1.0/  2.0)*  WMNIJ_2[    "imno"]*    GIJAK[    "noam"];
      DL_3(1)[      "ia"] +=                    FAE[      "ea"]*   L_2(1)[      "ie"];
      DL_3(1)[      "ia"] -=                    FMI[      "im"]*   L_2(1)[      "ma"];
      DL_3(1)[      "ia"] -=                  WAMEI[    "eiam"]*   L_2(1)[      "me"];
      DL_3(1)[      "ia"] +=  (1.0/  2.0)*    WABEJ[    "efam"]*   L_2(2)[    "imef"];
      DL_3(1)[      "ia"] -=  (1.0/  2.0)*    WAMIJ[    "eimn"]*   L_2(2)[    "mnea"];
      DL_3(1)[      "ia"] +=                    FME[      "ie"]*    DAB_2[      "ea"];
      DL_3(1)[      "ia"] -=                    FME[      "ma"]*    DIJ_2[      "im"];
      DL_3(1)[      "ia"] +=                  WMNEF[    "miea"]*    DAI_2[      "em"];
      DL_3(1)[      "ia"] -=  (1.0/  2.0)*    WABEF[    "efga"]*  GAIBC_2[    "gief"];
      DL_3(1)[      "ia"] +=                  WAMEI[    "eifm"]*  GAIBC_2[    "fmea"];
      DL_3(1)[      "ia"] -=                  WAMEI[    "eman"]*  GIJAK_2[    "inem"];
      DL_3(1)[      "ia"] +=  (1.0/  2.0)*    WMNIJ[    "imno"]*  GIJAK_2[    "noam"];
      DL_3(1)[      "ia"] -=  (1.0/  2.0)*    WAMEF[    "gief"]*  GABCD_2[    "efga"];
      DL_3(1)[      "ia"] +=                  WAMEF[    "fmea"]*  GAIBJ_2[    "eifm"];
      DL_3(1)[      "ia"] -=                  WMNEJ[    "inem"]*  GAIBJ_2[    "eman"];
      DL_3(1)[      "ia"] +=  (1.0/  2.0)*    WMNEJ[    "noam"]*  GIJKL_2[    "imno"];
      DL_3(1)[      "ia"] +=  (1.0/  2.0)*    WMNEF[    "imef"]*  GABCI_2[    "efam"];
      DL_3(1)[      "ia"] -=  (1.0/  2.0)*    WMNEF[    "mnea"]*  GAIJK_2[    "eimn"];

        DIJ_2[      "ij"] +=  (1.0/  2.0)*     T(2)[    "efjm"]*   L_2(2)[    "imef"];
        DIJ_2[      "ij"] +=  (1.0/  2.0)*   T_2(2)[    "efjm"]*     L(2)[    "imef"];

        DAB_2[      "ab"] -=  (1.0/  2.0)*     T(2)[    "aemn"]*   L_2(2)[    "mnbe"];
        DAB_2[      "ab"] -=  (1.0/  2.0)*   T_2(2)[    "aemn"]*     L(2)[    "mnbe"];

      DL_3(1)[      "ia"] -=                  WMNEJ[    "inam"]*    DIJ_2[      "mn"];
      DL_3(1)[      "ia"] -=                  WAMEF[    "fiea"]*    DAB_2[      "ef"];

      DL_3(2)[    "ijab"]  =                  FAE_2[      "ea"]*     L(2)[    "ijeb"];
      DL_3(2)[    "ijab"] -=                  FMI_2[      "im"]*     L(2)[    "mjab"];
      DL_3(2)[    "ijab"] +=  (1.0/  2.0)*  WABEF_2[    "efab"]*     L(2)[    "ijef"];
      DL_3(2)[    "ijab"] +=  (1.0/  2.0)*  WMNIJ_2[    "ijmn"]*     L(2)[    "mnab"];
      DL_3(2)[    "ijab"] +=                WAMEI_2[    "eiam"]*     L(2)[    "mjbe"];
      DL_3(2)[    "ijab"] +=  (1.0/  2.0)*  WABEJ_2[    "efbm"]*     L(3)[  "ijmaef"];
      DL_3(2)[    "ijab"] -=  (1.0/  2.0)*  WAMIJ_2[    "ejnm"]*     L(3)[  "imnabe"];
      DL_3(2)[    "ijab"] +=  (1.0/ 12.0)*WABCEJK_1[  "efgamn"]*   L_1(4)["ijmnebfg"];
      DL_3(2)[    "ijab"] -=  (1.0/ 12.0)*WABMIJK_1[  "efjmno"]*   L_1(4)["mnioefab"];
      DL_3(2)[    "ijab"] +=                    FME[      "ia"]*   L_2(1)[      "jb"];
      DL_3(2)[    "ijab"] +=                  WAMEF[    "ejab"]*   L_2(1)[      "ie"];
      DL_3(2)[    "ijab"] -=                  WMNEJ[    "ijam"]*   L_2(1)[      "mb"];
      DL_3(2)[    "ijab"] +=                    FAE[      "ea"]*   L_2(2)[    "ijeb"];
      DL_3(2)[    "ijab"] -=                    FMI[      "im"]*   L_2(2)[    "mjab"];
      DL_3(2)[    "ijab"] +=  (1.0/  2.0)*    WABEF[    "efab"]*   L_2(2)[    "ijef"];
      DL_3(2)[    "ijab"] +=  (1.0/  2.0)*    WMNIJ[    "ijmn"]*   L_2(2)[    "mnab"];
      DL_3(2)[    "ijab"] +=                  WAMEI[    "eiam"]*   L_2(2)[    "mjbe"];
      DL_3(2)[    "ijab"] +=  (1.0/  2.0)*    WABEJ[    "efbm"]*   L_2(3)[  "ijmaef"];
      DL_3(2)[    "ijab"] -=  (1.0/  2.0)*    WAMIJ[    "ejnm"]*   L_2(3)[  "imnabe"];
      DL_3(2)[    "ijab"] +=  (1.0/ 12.0)*  WABCEJK[  "efgamn"]*   L_2(4)["ijmnebfg"];
      DL_3(2)[    "ijab"] -=  (1.0/ 12.0)*  WABMIJK[  "efjmno"]*   L_2(4)["mnioefab"];
      DL_3(2)[    "ijab"] -=                  WMNEF[    "mjab"]*    DIJ_2[      "im"];
      DL_3(2)[    "ijab"] +=                  WMNEF[    "ijeb"]*    DAB_2[      "ea"];
      DL_3(2)[    "ijab"] +=  (1.0/  2.0)*    WMNEF[    "ijef"]*  GABCD_2[    "efab"];
      DL_3(2)[    "ijab"] +=                  WMNEF[    "imea"]*  GAIBJ_2[    "ejbm"];
      DL_3(2)[    "ijab"] +=  (1.0/  2.0)*    WMNEF[    "mnab"]*  GIJKL_2[    "ijmn"];
      DL_3(2)[    "ijab"] -=                  WAMEF[    "fiae"]*  GAIBC_2[    "ejbf"];
      DL_3(2)[    "ijab"] -=                  WMNEJ[    "ijem"]*  GAIBC_2[    "emab"];
      DL_3(2)[    "ijab"] -=                  WAMEF[    "emab"]*  GIJAK_2[    "ijem"];
      DL_3(2)[    "ijab"] -=                  WMNEJ[    "niam"]*  GIJAK_2[    "mjbn"];

      DL_3(3)[  "ijkabc"]  =                  FAE_2[      "ea"]*     L(3)[  "ijkebc"];
      DL_3(3)[  "ijkabc"] -=                  FMI_2[      "im"]*     L(3)[  "mjkabc"];
      DL_3(3)[  "ijkabc"] +=  (1.0/  2.0)*  WABEF_2[    "efab"]*     L(3)[  "ijkefc"];
      DL_3(3)[  "ijkabc"] +=  (1.0/  2.0)*  WMNIJ_2[    "ijmn"]*     L(3)[  "mnkabc"];
      DL_3(3)[  "ijkabc"] +=                WAMEI_2[    "eiam"]*     L(3)[  "mjkbec"];
      DL_3(3)[  "ijkabc"] +=                  WMNEF[    "ijab"]*   L_2(1)[      "kc"];
      DL_3(3)[  "ijkabc"] +=                    FME[      "ia"]*   L_2(2)[    "jkbc"];
      DL_3(3)[  "ijkabc"] +=                  WAMEF[    "ekbc"]*   L_2(2)[    "ijae"];
      DL_3(3)[  "ijkabc"] -=                  WMNEJ[    "ijam"]*   L_2(2)[    "mkbc"];
      DL_3(3)[  "ijkabc"] +=                    FAE[      "ea"]*   L_2(3)[  "ijkebc"];
      DL_3(3)[  "ijkabc"] -=                    FMI[      "im"]*   L_2(3)[  "mjkabc"];
      DL_3(3)[  "ijkabc"] +=  (1.0/  2.0)*    WABEF[    "efab"]*   L_2(3)[  "ijkefc"];
      DL_3(3)[  "ijkabc"] +=  (1.0/  2.0)*    WMNIJ[    "ijmn"]*   L_2(3)[  "mnkabc"];
      DL_3(3)[  "ijkabc"] +=                  WAMEI[    "eiam"]*   L_2(3)[  "mjkbec"];
      DL_3(3)[  "ijkabc"] +=  (1.0/  2.0)*    WABEJ[    "efam"]*   L_2(4)["ijkmebcf"];
      DL_3(3)[  "ijkabc"] -=  (1.0/  2.0)*    WAMIJ[    "eknm"]*   L_2(4)["ijmnabce"];
      DL_3(3)[  "ijkabc"] -=  (1.0/  4.0)*  WABMEJI[  "efkcnm"]*   L_2(4)["ijmnabef"];
      DL_3(3)[  "ijkabc"] +=  (1.0/  6.0)*  WAMNIJK[  "eijmno"]*   L_2(4)["mnokeabc"];
      DL_3(3)[  "ijkabc"] +=                  WMNEF[    "ijae"]*  GAIBC_2[    "ekbc"];
      DL_3(3)[  "ijkabc"] -=                  WMNEF[    "mkbc"]*  GIJAK_2[    "ijam"];
      DL_3(3)[  "ijkabc"] -=                  WAMEF[    "embc"]*GIJKABL_2[  "ijkaem"];
      DL_3(3)[  "ijkabc"] +=  (1.0/  2.0)*    WMNEF[    "mnbc"]*GIJKALM_2[  "ijkamn"];

      DL_3(4)["ijklabcd"]  =                  WMNEF[    "ijab"]*   L_2(2)[    "klcd"];
      DL_3(4)["ijklabcd"] +=                    FME[      "ia"]*   L_2(3)[  "jklbcd"];
      DL_3(4)["ijklabcd"] +=                  WAMEF[    "ejab"]*   L_2(3)[  "iklecd"];
      DL_3(4)["ijklabcd"] -=                  WMNEJ[    "ijam"]*   L_2(3)[  "mklbcd"];
      DL_3(4)["ijklabcd"] +=                  WMNEF[    "ijae"]*GAIJBCD_2[  "eklbcd"];
      DL_3(4)["ijklabcd"] -=                  WMNEF[    "mlcd"]*GIJKABL_2[  "ijkabm"];
      DL_3(4)["ijklabcd"] +=                    FAE[      "ea"]*   L_2(4)["ijklebcd"];
      DL_3(4)["ijklabcd"] -=                    FMI[      "im"]*   L_2(4)["mjklabcd"];
      DL_3(4)["ijklabcd"] +=  (1.0/  2.0)*    WABEF[    "efab"]*   L_2(4)["ijklefcd"];
      DL_3(4)["ijklabcd"] +=  (1.0/  2.0)*    WMNIJ[    "ijmn"]*   L_2(4)["mnklabcd"];
      DL_3(4)["ijklabcd"] +=                  WAMEI[    "eiam"]*   L_2(4)["mjklbecd"];

    L_3 = DL_3;
    L_3.weight(D);

    /***************************************************************************
     *
     * L^(4)
     *
     **************************************************************************/

    WABCEJK_2[  "abcejk"]  =               -WAMEI_2[    "amek"]*     T(2)[    "bcjm"];
    WABCEJK_2[  "abcejk"] +=                WABEF_2[    "abef"]*     T(2)[    "fcjk"];
    WABCEJK_2[  "abcejk"] -=                  WAMEI[    "amek"]*   T_2(2)[    "bcjm"];
    WABCEJK_2[  "abcejk"] +=                  WABEF[    "abef"]*   T_2(2)[    "fcjk"];
    WABCEJK_2[  "abcejk"] -=                    FME[      "me"]*   T_2(3)[  "abcmjk"];
    WABCEJK_2[  "abcejk"] +=              0.5*WMNEJ[    "mnej"]*   T_2(3)[  "abcmnk"];
    WABCEJK_2[  "abcejk"] +=                  WAMEF[    "amef"]*   T_2(3)[  "fbcmjk"];
    WABCEJK_2[  "abcejk"] -=  (1.0/  2.0)*    WMNEF[    "mnef"]*   T_2(4)["abcfmjkn"];

    WABMIJK_2[  "abmijk"]  =                WAMEI_2[    "amek"]*     T(2)[    "ebij"];
    WABMIJK_2[  "abmijk"] -=                WMNIJ_2[    "nmjk"]*     T(2)[    "abin"];
    WABMIJK_2[  "abmijk"] +=                  WAMEI[    "amek"]*   T_2(2)[    "ebij"];
    WABMIJK_2[  "abmijk"] -=                  WMNIJ[    "nmjk"]*   T_2(2)[    "abin"];
    WABMIJK_2[  "abmijk"] +=                    FME[      "me"]*   T_2(3)[  "abeijk"];
    WABMIJK_2[  "abmijk"] +=              0.5*WAMEF[    "bmef"]*   T_2(3)[  "aefijk"];
    WABMIJK_2[  "abmijk"] +=                  WMNEJ[    "nmek"]*   T_2(3)[  "abeijn"];
    WABMIJK_2[  "abmijk"] +=  (1.0/  2.0)*    WMNEF[    "mnef"]*   T_2(4)["abefijkn"];

    WAMNIJK_2[  "amnijk"]  =  (1.0/  2.0)*    WMNEF[    "mnef"]*   T_2(3)[  "aefijk"];
    WAMNIJK_2[  "amnijk"] +=                  WMNEJ[    "mnek"]*   T_2(2)[    "aeij"];

    WABMEJI_2[  "abmeji"]  =                  WAMEF[    "amef"]*   T_2(2)[    "bfji"];
    WABMEJI_2[  "abmeji"] -=                  WMNEJ[    "nmei"]*   T_2(2)[    "abnj"];
    WABMEJI_2[  "abmeji"] +=                  WMNEF[    "mnef"]*   T_2(3)[  "abfnji"];

    GIJKABL_3[  "ijkabl"]  =  (1.0/  2.0)*     T(2)[    "eflm"]*   L_3(4)["ijkmabef"];
    GIJKABL_3[  "ijkabl"] +=  (1.0/  2.0)*   T_2(2)[    "eflm"]*   L_1(4)["ijkmabef"];

    GAIJBCD_3[  "aijbcd"]  = -(1.0/  2.0)*     T(2)[    "aemn"]*   L_3(4)["mijnbcde"];
    GAIJBCD_3[  "aijbcd"] -=  (1.0/  2.0)*   T_2(2)[    "aemn"]*   L_1(4)["mijnbcde"];

    GIJKALM_3[  "ijkalm"]  =  (1.0/  6.0)*     T(3)[  "efglmn"]*   L_3(4)["ijknaefg"];
    GIJKALM_3[  "ijkalm"] +=  (1.0/  6.0)*   T_2(3)[  "efglmn"]*   L_1(4)["ijknaefg"];

        DIJ_3[      "ij"]  =  (1.0/ 12.0)*     T(3)[  "efgjmn"]*   L_3(3)[  "imnefg"];
        DIJ_3[      "ij"] +=  (1.0/144.0)*   T_1(4)["efghjmno"]*   L_2(4)["imnoefgh"];
        DIJ_3[      "ij"] +=  (1.0/144.0)*   T_2(4)["efghjmno"]*   L_1(4)["imnoefgh"];
        DIJ_3[      "ij"] +=  (1.0/ 12.0)*   T_3(3)[  "efgjmn"]*     L(3)[  "imnefg"];

        DAB_3[      "ab"]  = -(1.0/ 12.0)*     T(3)[  "aefmno"]*   L_3(3)[  "mnobef"];
        DAB_3[      "ab"] -=  (1.0/144.0)*   T_1(4)["aefgmnop"]*   L_2(4)["mnopbefg"];
        DAB_3[      "ab"] -=  (1.0/144.0)*   T_2(4)["aefgmnop"]*   L_1(4)["mnopbefg"];
        DAB_3[      "ab"] -=  (1.0/ 12.0)*   T_3(3)[  "aefmno"]*     L(3)[  "mnobef"];

      GIJAK_3[    "ijak"]  =  (1.0/ 12.0)*     T(3)[  "efgkmn"]*   L_3(4)["ijmnaefg"];
      GIJAK_3[    "ijak"] +=  (1.0/  2.0)*     T(2)[    "efkm"]*   L_3(3)[  "ijmaef"];
      GIJAK_3[    "ijak"] +=  (1.0/ 12.0)*   T_2(3)[  "efgkmn"]*   L_1(4)["ijmnaefg"];
      GIJAK_3[    "ijak"] +=  (1.0/  2.0)*   T_3(2)[    "efkm"]*     L(3)[  "ijmaef"];

      GAIBC_3[    "aibc"]  = -(1.0/ 12.0)*     T(3)[  "aefmno"]*   L_3(4)["minobcef"];
      GAIBC_3[    "aibc"] -=  (1.0/ 12.0)*   T_2(3)[  "aefmno"]*   L_1(4)["minobcef"];
      GAIBC_3[    "aibc"] -=  (1.0/  2.0)*   T_3(2)[    "aemn"]*     L(3)[  "minbce"];

      GAIJK_3[    "aijk"]  =  (1.0/ 12.0)*   T_1(4)["aefgjkmn"]*   L_2(3)[  "imnefg"];
      GAIJK_3[    "aijk"] +=  (1.0/ 12.0)*   T_3(4)["aefgjkmn"]*     L(3)[  "imnefg"];
      GAIJK_3[    "aijk"] +=  (1.0/  4.0)*     T(3)[  "efamnj"]*GIJKABL_3[  "mniefk"];
      GAIJK_3[    "aijk"] +=  (1.0/  6.0)*     T(3)[  "efgjkm"]*GAIJBCD_3[  "aimefg"];
      GAIJK_3[    "aijk"] +=  (1.0/  4.0)*   T_2(3)[  "efamnj"]*GIJKABL_1[  "mniefk"];
      GAIJK_3[    "aijk"] +=  (1.0/  6.0)*   T_2(3)[  "efgjkm"]*GAIJBCD_1[  "aimefg"];

      GABCI_3[    "abci"]  = -(1.0/ 12.0)*   T_1(4)["abefmino"]*   L_2(3)[  "mnocef"];
      GABCI_3[    "abci"] -=  (1.0/ 12.0)*   T_3(4)["abefmino"]*     L(3)[  "mnocef"];
      GABCI_3[    "abci"] +=  (1.0/  4.0)*     T(3)[  "efbmni"]*GAIJBCD_3[  "amncef"];
      GABCI_3[    "abci"] +=  (1.0/  6.0)*     T(3)[  "eabmno"]*GIJKABL_3[  "mnoeci"];
      GABCI_3[    "abci"] +=  (1.0/  4.0)*   T_2(3)[  "efbmni"]*GAIJBCD_1[  "amncef"];
      GABCI_3[    "abci"] +=  (1.0/  6.0)*   T_2(3)[  "eabmno"]*GIJKABL_1[  "mnoeci"];

        DAI_3[      "ai"]  =  (1.0/  4.0)*     T(3)[  "aefimn"]*   L_3(2)[    "mnef"];
        DAI_3[      "ai"] +=  (1.0/ 36.0)*   T_1(4)["aefgimno"]*   L_2(3)[  "mnoefg"];
        DAI_3[      "ai"] +=  (1.0/ 36.0)*   T_3(4)["aefgimno"]*     L(3)[  "mnoefg"];
        DAI_3[      "ai"] +=  (1.0/  4.0)*   T_3(3)[  "aefimn"]*     L(2)[    "mnef"];
        DAI_3[      "ai"] -=  (1.0/  2.0)*   T_2(2)[    "eamn"]*  GIJAK_1[    "mnei"];
        DAI_3[      "ai"] +=  (1.0/  2.0)*   T_2(2)[    "efim"]*  GAIBC_1[    "amef"];
        DAI_3[      "ai"] -=  (1.0/  2.0)*     T(2)[    "eamn"]*  GIJAK_3[    "mnei"];
        DAI_3[      "ai"] +=  (1.0/  2.0)*     T(2)[    "efim"]*  GAIBC_3[    "amef"];

      GAIBC_3[    "aibc"] -=  (1.0/  2.0)*     T(2)[    "aemn"]*   L_3(3)[  "minbce"];

      GABCD_3[    "abcd"]  =  (1.0/  6.0)*     T(3)[  "abemno"]*   L_3(3)[  "mnocde"];
      GABCD_3[    "abcd"] +=  (1.0/ 48.0)*   T_1(4)["abefmnop"]*   L_2(4)["mnopcdef"];
      GABCD_3[    "abcd"] +=  (1.0/ 48.0)*   T_2(4)["abefmnop"]*   L_1(4)["mnopcdef"];
      GABCD_3[    "abcd"] +=  (1.0/  6.0)*   T_3(3)[  "abemno"]*     L(3)[  "mnocde"];
      GABCD_3[    "abcd"] -=  (1.0/  4.0)*     T(2)[    "bemn"]*GAIJBCD_3[  "amncde"];
      GABCD_3[    "abcd"] -=  (1.0/  4.0)*   T_2(2)[    "bemn"]*GAIJBCD_1[  "amncde"];

      GAIBJ_3[    "aibj"]  = -(1.0/  4.0)*     T(3)[  "aefjmn"]*   L_3(3)[  "imnbef"];
      GAIBJ_3[    "aibj"] -=  (1.0/ 36.0)*   T_1(4)["aefgjmno"]*   L_2(4)["imnobefg"];
      GAIBJ_3[    "aibj"] -=  (1.0/ 36.0)*   T_2(4)["aefgjmno"]*   L_1(4)["imnobefg"];
      GAIBJ_3[    "aibj"] -=  (1.0/  4.0)*   T_3(3)[  "aefjmn"]*     L(3)[  "imnbef"];
      GAIBJ_3[    "aibj"] +=  (1.0/  2.0)*     T(2)[    "eamn"]*GIJKABL_3[  "imnbej"];
      GAIBJ_3[    "aibj"] +=  (1.0/  2.0)*   T_2(2)[    "efjm"]*GAIJBCD_1[  "aimbef"];

      GIJKL_3[    "ijkl"]  =  (1.0/  6.0)*     T(3)[  "efgklm"]*   L_3(3)[  "ijmefg"];
      GIJKL_3[    "ijkl"] +=  (1.0/ 48.0)*   T_1(4)["efghklmn"]*   L_2(4)["ijmnefgh"];
      GIJKL_3[    "ijkl"] +=  (1.0/ 48.0)*   T_2(4)["efghklmn"]*   L_1(4)["ijmnefgh"];
      GIJKL_3[    "ijkl"] +=  (1.0/  6.0)*   T_3(3)[  "efgklm"]*     L(3)[  "ijmefg"];
      GIJKL_3[    "ijkl"] +=  (1.0/  4.0)*     T(2)[    "efmk"]*GIJKABL_3[  "mijefl"];
      GIJKL_3[    "ijkl"] +=  (1.0/  4.0)*   T_2(2)[    "efmk"]*GIJKABL_1[  "mijefl"];

          DIJ[      "ij"]  =  (1.0/ 12.0)*     T(3)[  "efgjmn"]*     L(3)[  "imnefg"];
          DAB[      "ab"]  = -(1.0/ 12.0)*     T(3)[  "aefmno"]*     L(3)[  "mnobef"];

        GABCD[    "abcd"]  =  (1.0/  6.0)*     T(3)[  "abemno"]*     L(3)[  "mnocde"];
        GAIBJ[    "aibj"]  = -(1.0/  4.0)*     T(3)[  "aefjmn"]*     L(3)[  "imnbef"];
        GIJKL[    "ijkl"]  =  (1.0/  6.0)*     T(3)[  "efgklm"]*     L(3)[  "ijmefg"];

      DL_4(1)[      "ia"]  =                  FME_3[      "ia"];
      DL_4(1)[      "ia"] +=                  FAE_3[      "ea"]*     L(1)[      "ie"];
      DL_4(1)[      "ia"] -=                  FMI_3[      "im"]*     L(1)[      "ma"];
      DL_4(1)[      "ia"] -=                WAMEI_3[    "eiam"]*     L(1)[      "me"];
      DL_4(1)[      "ia"] +=  (1.0/  2.0)*  WABEJ_3[    "efam"]*     L(2)[    "imef"];
      DL_4(1)[      "ia"] -=  (1.0/  2.0)*  WAMIJ_3[    "eimn"]*     L(2)[    "mnea"];
      DL_4(1)[      "ia"] +=                  FME_3[      "ie"]*      DAB[      "ea"];
      DL_4(1)[      "ia"] -=                  FME_3[      "ma"]*      DIJ[      "im"];
      DL_4(1)[      "ia"] -=  (1.0/  2.0)*  WABEF_3[    "efga"]*    GAIBC[    "gief"];
      DL_4(1)[      "ia"] +=                WAMEI_3[    "eifm"]*    GAIBC[    "fmea"];
      DL_4(1)[      "ia"] -=                WAMEI_3[    "eman"]*    GIJAK[    "inem"];
      DL_4(1)[      "ia"] +=  (1.0/  2.0)*  WMNIJ_3[    "imno"]*    GIJAK[    "noam"];
      DL_4(1)[      "ia"] -=  (1.0/  2.0)*  WAMEF_3[    "gief"]*    GABCD[    "efga"];
      DL_4(1)[      "ia"] +=                WAMEF_3[    "fmea"]*    GAIBJ[    "eifm"];
      DL_4(1)[      "ia"] -=                WMNEJ_3[    "inem"]*    GAIBJ[    "eman"];
      DL_4(1)[      "ia"] +=  (1.0/  2.0)*  WMNEJ_3[    "noam"]*    GIJKL[    "imno"];
      DL_4(1)[      "ia"] +=                    FAE[      "ea"]*   L_3(1)[      "ie"];
      DL_4(1)[      "ia"] -=                    FMI[      "im"]*   L_3(1)[      "ma"];
      DL_4(1)[      "ia"] -=                  WAMEI[    "eiam"]*   L_3(1)[      "me"];
      DL_4(1)[      "ia"] -=  (1.0/  2.0)*  WABEF_2[    "efga"]*  GAIBC_1[    "gief"];
      DL_4(1)[      "ia"] +=                WAMEI_2[    "eifm"]*  GAIBC_1[    "fmea"];
      DL_4(1)[      "ia"] -=                WAMEI_2[    "eman"]*  GIJAK_1[    "inem"];
      DL_4(1)[      "ia"] +=  (1.0/  2.0)*  WMNIJ_2[    "imno"]*  GIJAK_1[    "noam"];
      DL_4(1)[      "ia"] +=  (1.0/  2.0)*    WABEJ[    "efam"]*   L_3(2)[    "imef"];
      DL_4(1)[      "ia"] -=  (1.0/  2.0)*    WAMIJ[    "eimn"]*   L_3(2)[    "mnea"];
      DL_4(1)[      "ia"] +=                    FME[      "ie"]*    DAB_3[      "ea"];
      DL_4(1)[      "ia"] -=                    FME[      "ma"]*    DIJ_3[      "im"];
      DL_4(1)[      "ia"] +=                  WMNEF[    "miea"]*    DAI_3[      "em"];
      DL_4(1)[      "ia"] -=  (1.0/  2.0)*    WABEF[    "efga"]*  GAIBC_3[    "gief"];
      DL_4(1)[      "ia"] +=                  WAMEI[    "eifm"]*  GAIBC_3[    "fmea"];
      DL_4(1)[      "ia"] -=                  WAMEI[    "eman"]*  GIJAK_3[    "inem"];
      DL_4(1)[      "ia"] +=  (1.0/  2.0)*    WMNIJ[    "imno"]*  GIJAK_3[    "noam"];
      DL_4(1)[      "ia"] -=  (1.0/  2.0)*    WAMEF[    "gief"]*  GABCD_3[    "efga"];
      DL_4(1)[      "ia"] +=                  WAMEF[    "fmea"]*  GAIBJ_3[    "eifm"];
      DL_4(1)[      "ia"] -=                  WMNEJ[    "inem"]*  GAIBJ_3[    "eman"];
      DL_4(1)[      "ia"] +=  (1.0/  2.0)*    WMNEJ[    "noam"]*  GIJKL_3[    "imno"];
      DL_4(1)[      "ia"] +=  (1.0/  2.0)*    WMNEF[    "imef"]*  GABCI_3[    "efam"];
      DL_4(1)[      "ia"] -=  (1.0/  2.0)*    WMNEF[    "mnea"]*  GAIJK_3[    "eimn"];

        DIJ_3[      "ij"] +=  (1.0/  2.0)*     T(2)[    "efjm"]*   L_3(2)[    "imef"];
        DIJ_3[      "ij"] +=  (1.0/  2.0)*   T_3(2)[    "efjm"]*     L(2)[    "imef"];

        DAB_3[      "ab"] -=  (1.0/  2.0)*     T(2)[    "aemn"]*   L_3(2)[    "mnbe"];
        DAB_3[      "ab"] -=  (1.0/  2.0)*   T_3(2)[    "aemn"]*     L(2)[    "mnbe"];

      DL_4(1)[      "ia"] -=                  WMNEJ[    "inam"]*    DIJ_3[      "mn"];
      DL_4(1)[      "ia"] -=                  WAMEF[    "fiea"]*    DAB_3[      "ef"];

      DL_4(2)[    "ijab"]  =                  FME_3[      "ia"]*     L(1)[      "jb"];
      DL_4(2)[    "ijab"] +=                WAMEF_3[    "ejab"]*     L(1)[      "ie"];
      DL_4(2)[    "ijab"] -=                WMNEJ_3[    "ijam"]*     L(1)[      "mb"];
      DL_4(2)[    "ijab"] +=                  FAE_3[      "ea"]*     L(2)[    "ijeb"];
      DL_4(2)[    "ijab"] -=                  FMI_3[      "im"]*     L(2)[    "mjab"];
      DL_4(2)[    "ijab"] +=  (1.0/  2.0)*  WABEF_3[    "efab"]*     L(2)[    "ijef"];
      DL_4(2)[    "ijab"] +=  (1.0/  2.0)*  WMNIJ_3[    "ijmn"]*     L(2)[    "mnab"];
      DL_4(2)[    "ijab"] +=                WAMEI_3[    "eiam"]*     L(2)[    "mjbe"];
      DL_4(2)[    "ijab"] +=  (1.0/  2.0)*  WABEJ_3[    "efbm"]*     L(3)[  "ijmaef"];
      DL_4(2)[    "ijab"] -=  (1.0/  2.0)*  WAMIJ_3[    "ejnm"]*     L(3)[  "imnabe"];
      DL_4(2)[    "ijab"] -=                WAMEF_3[    "fiae"]*    GAIBC[    "ejbf"];
      DL_4(2)[    "ijab"] -=                WMNEJ_3[    "ijem"]*    GAIBC[    "emab"];
      DL_4(2)[    "ijab"] -=                WAMEF_3[    "emab"]*    GIJAK[    "ijem"];
      DL_4(2)[    "ijab"] -=                WMNEJ_3[    "niam"]*    GIJAK[    "mjbn"];
      DL_4(2)[    "ijab"] +=  (1.0/ 12.0)*WABCEJK_2[  "efgamn"]*   L_1(4)["ijmnebfg"];
      DL_4(2)[    "ijab"] -=  (1.0/ 12.0)*WABMIJK_2[  "efjmno"]*   L_1(4)["mnioefab"];
      DL_4(2)[    "ijab"] +=  (1.0/ 12.0)*WABCEJK_1[  "efgamn"]*   L_2(4)["ijmnebfg"];
      DL_4(2)[    "ijab"] -=  (1.0/ 12.0)*WABMIJK_1[  "efjmno"]*   L_2(4)["mnioefab"];
      DL_4(2)[    "ijab"] +=                    FME[      "ia"]*   L_3(1)[      "jb"];
      DL_4(2)[    "ijab"] +=                  WAMEF[    "ejab"]*   L_3(1)[      "ie"];
      DL_4(2)[    "ijab"] -=                  WMNEJ[    "ijam"]*   L_3(1)[      "mb"];
      DL_4(2)[    "ijab"] +=                    FAE[      "ea"]*   L_3(2)[    "ijeb"];
      DL_4(2)[    "ijab"] -=                    FMI[      "im"]*   L_3(2)[    "mjab"];
      DL_4(2)[    "ijab"] +=  (1.0/  2.0)*    WABEF[    "efab"]*   L_3(2)[    "ijef"];
      DL_4(2)[    "ijab"] +=  (1.0/  2.0)*    WMNIJ[    "ijmn"]*   L_3(2)[    "mnab"];
      DL_4(2)[    "ijab"] +=                  WAMEI[    "eiam"]*   L_3(2)[    "mjbe"];
      DL_4(2)[    "ijab"] +=  (1.0/  2.0)*    WABEJ[    "efbm"]*   L_3(3)[  "ijmaef"];
      DL_4(2)[    "ijab"] -=  (1.0/  2.0)*    WAMIJ[    "ejnm"]*   L_3(3)[  "imnabe"];
      DL_4(2)[    "ijab"] +=  (1.0/ 12.0)*  WABCEJK[  "efgamn"]*   L_3(4)["ijmnebfg"];
      DL_4(2)[    "ijab"] -=  (1.0/ 12.0)*  WABMIJK[  "efjmno"]*   L_3(4)["mnioefab"];
      DL_4(2)[    "ijab"] -=                  WMNEF[    "mjab"]*    DIJ_3[      "im"];
      DL_4(2)[    "ijab"] +=                  WMNEF[    "ijeb"]*    DAB_3[      "ea"];
      DL_4(2)[    "ijab"] +=  (1.0/  2.0)*    WMNEF[    "ijef"]*  GABCD_3[    "efab"];
      DL_4(2)[    "ijab"] +=                  WMNEF[    "imea"]*  GAIBJ_3[    "ejbm"];
      DL_4(2)[    "ijab"] +=  (1.0/  2.0)*    WMNEF[    "mnab"]*  GIJKL_3[    "ijmn"];
      DL_4(2)[    "ijab"] -=                  WAMEF[    "fiae"]*  GAIBC_3[    "ejbf"];
      DL_4(2)[    "ijab"] -=                  WMNEJ[    "ijem"]*  GAIBC_3[    "emab"];
      DL_4(2)[    "ijab"] -=                  WAMEF[    "emab"]*  GIJAK_3[    "ijem"];
      DL_4(2)[    "ijab"] -=                  WMNEJ[    "niam"]*  GIJAK_3[    "mjbn"];

      DL_4(3)[  "ijkabc"]  =                  FME_3[      "ia"]*     L(2)[    "jkbc"];
      DL_4(3)[  "ijkabc"] +=                WAMEF_3[    "ekbc"]*     L(2)[    "ijae"];
      DL_4(3)[  "ijkabc"] -=                WMNEJ_3[    "ijam"]*     L(2)[    "mkbc"];
      DL_4(3)[  "ijkabc"] +=                  FAE_3[      "ea"]*     L(3)[  "ijkebc"];
      DL_4(3)[  "ijkabc"] -=                  FMI_3[      "im"]*     L(3)[  "mjkabc"];
      DL_4(3)[  "ijkabc"] +=  (1.0/  2.0)*  WABEF_3[    "efab"]*     L(3)[  "ijkefc"];
      DL_4(3)[  "ijkabc"] +=  (1.0/  2.0)*  WMNIJ_3[    "ijmn"]*     L(3)[  "mnkabc"];
      DL_4(3)[  "ijkabc"] +=                WAMEI_3[    "eiam"]*     L(3)[  "mjkbec"];
      DL_4(3)[  "ijkabc"] +=  (1.0/  2.0)*  WABEJ_2[    "efam"]*   L_1(4)["ijkmebcf"];
      DL_4(3)[  "ijkabc"] -=  (1.0/  2.0)*  WAMIJ_2[    "eknm"]*   L_1(4)["ijmnabce"];
      DL_4(3)[  "ijkabc"] -=  (1.0/  4.0)*WABMEJI_2[  "efkcnm"]*   L_1(4)["ijmnabef"];
      DL_4(3)[  "ijkabc"] +=  (1.0/  6.0)*WAMNIJK_2[  "eijmno"]*   L_1(4)["mnokeabc"];
      DL_4(3)[  "ijkabc"] +=                  WMNEF[    "ijab"]*   L_3(1)[      "kc"];
      DL_4(3)[  "ijkabc"] +=                    FME[      "ia"]*   L_3(2)[    "jkbc"];
      DL_4(3)[  "ijkabc"] +=                  WAMEF[    "ekbc"]*   L_3(2)[    "ijae"];
      DL_4(3)[  "ijkabc"] -=                  WMNEJ[    "ijam"]*   L_3(2)[    "mkbc"];
      DL_4(3)[  "ijkabc"] +=                    FAE[      "ea"]*   L_3(3)[  "ijkebc"];
      DL_4(3)[  "ijkabc"] -=                    FMI[      "im"]*   L_3(3)[  "mjkabc"];
      DL_4(3)[  "ijkabc"] +=  (1.0/  2.0)*    WABEF[    "efab"]*   L_3(3)[  "ijkefc"];
      DL_4(3)[  "ijkabc"] +=  (1.0/  2.0)*    WMNIJ[    "ijmn"]*   L_3(3)[  "mnkabc"];
      DL_4(3)[  "ijkabc"] +=                  WAMEI[    "eiam"]*   L_3(3)[  "mjkbec"];
      DL_4(3)[  "ijkabc"] +=  (1.0/  2.0)*    WABEJ[    "efam"]*   L_3(4)["ijkmebcf"];
      DL_4(3)[  "ijkabc"] -=  (1.0/  2.0)*    WAMIJ[    "eknm"]*   L_3(4)["ijmnabce"];
      DL_4(3)[  "ijkabc"] -=  (1.0/  4.0)*  WABMEJI[  "efkcnm"]*   L_3(4)["ijmnabef"];
      DL_4(3)[  "ijkabc"] +=  (1.0/  6.0)*  WAMNIJK[  "eijmno"]*   L_3(4)["mnokeabc"];
      DL_4(3)[  "ijkabc"] +=                  WMNEF[    "ijae"]*  GAIBC_3[    "ekbc"];
      DL_4(3)[  "ijkabc"] -=                  WMNEF[    "mkbc"]*  GIJAK_3[    "ijam"];
      DL_4(3)[  "ijkabc"] -=                  WAMEF[    "embc"]*GIJKABL_3[  "ijkaem"];
      DL_4(3)[  "ijkabc"] +=  (1.0/  2.0)*    WMNEF[    "mnbc"]*GIJKALM_3[  "ijkamn"];

      DL_4(4)["ijklabcd"]  =                  FME_3[      "ia"]*     L(3)[  "jklbcd"];
      DL_4(4)["ijklabcd"] +=                WAMEF_3[    "ejab"]*     L(3)[  "iklecd"];
      DL_4(4)["ijklabcd"] -=                WMNEJ_3[    "ijam"]*     L(3)[  "mklbcd"];
      DL_4(4)["ijklabcd"] +=                  FAE_2[      "ea"]*   L_1(4)["ijklebcd"];
      DL_4(4)["ijklabcd"] -=                  FMI_2[      "im"]*   L_1(4)["mjklabcd"];
      DL_4(4)["ijklabcd"] +=  (1.0/  2.0)*  WABEF_2[    "efab"]*   L_1(4)["ijklefcd"];
      DL_4(4)["ijklabcd"] +=  (1.0/  2.0)*  WMNIJ_2[    "ijmn"]*   L_1(4)["mnklabcd"];
      DL_4(4)["ijklabcd"] +=                WAMEI_2[    "eiam"]*   L_1(4)["mjklbecd"];
      DL_4(4)["ijklabcd"] +=                  WMNEF[    "ijab"]*   L_3(2)[    "klcd"];
      DL_4(4)["ijklabcd"] +=                    FME[      "ia"]*   L_3(3)[  "jklbcd"];
      DL_4(4)["ijklabcd"] +=                  WAMEF[    "ejab"]*   L_3(3)[  "iklecd"];
      DL_4(4)["ijklabcd"] -=                  WMNEJ[    "ijam"]*   L_3(3)[  "mklbcd"];
      DL_4(4)["ijklabcd"] +=                  WMNEF[    "ijae"]*GAIJBCD_3[  "eklbcd"];
      DL_4(4)["ijklabcd"] -=                  WMNEF[    "mlcd"]*GIJKABL_3[  "ijkabm"];
      DL_4(4)["ijklabcd"] +=                    FAE[      "ea"]*   L_3(4)["ijklebcd"];
      DL_4(4)["ijklabcd"] -=                    FMI[      "im"]*   L_3(4)["mjklabcd"];
      DL_4(4)["ijklabcd"] +=  (1.0/  2.0)*    WABEF[    "efab"]*   L_3(4)["ijklefcd"];
      DL_4(4)["ijklabcd"] +=  (1.0/  2.0)*    WMNIJ[    "ijmn"]*   L_3(4)["mnklabcd"];
      DL_4(4)["ijklabcd"] +=                  WAMEI[    "eiam"]*   L_3(4)["mjklbecd"];

    L_4 = DL_4;
    L_4.weight(D);

    /***************************************************************************
     *
     * L^(5)
     *
     **************************************************************************/

    GIJKABL_4[  "ijkabl"]  =  (1.0/  2.0)*   T_3(2)[    "eflm"]*   L_1(4)["ijkmabef"];
    GIJKABL_4[  "ijkabl"] +=  (1.0/  2.0)*   T_2(2)[    "eflm"]*   L_2(4)["ijkmabef"];
    GIJKABL_4[  "ijkabl"] +=  (1.0/  2.0)*     T(2)[    "eflm"]*   L_4(4)["ijkmabef"];

    GAIJBCD_4[  "aijbcd"]  = -(1.0/  2.0)*   T_3(2)[    "aemn"]*   L_1(4)["mijnbcde"];
    GAIJBCD_4[  "aijbcd"] -=  (1.0/  2.0)*   T_2(2)[    "aemn"]*   L_2(4)["mijnbcde"];
    GAIJBCD_4[  "aijbcd"] -=  (1.0/  2.0)*     T(2)[    "aemn"]*   L_4(4)["mijnbcde"];

      DL_5(4)["ijklabcd"]  =                  FME_4[      "ia"]*     L(3)[  "jklbcd"];
      DL_5(4)["ijklabcd"] +=                WAMEF_4[    "ejab"]*     L(3)[  "iklecd"];
      DL_5(4)["ijklabcd"] -=                WMNEJ_4[    "ijam"]*     L(3)[  "mklbcd"];
      DL_5(4)["ijklabcd"] +=                  FAE_3[      "ea"]*   L_1(4)["ijklebcd"];
      DL_5(4)["ijklabcd"] -=                  FMI_3[      "im"]*   L_1(4)["mjklabcd"];
      DL_5(4)["ijklabcd"] +=  (1.0/  2.0)*  WABEF_3[    "efab"]*   L_1(4)["ijklefcd"];
      DL_5(4)["ijklabcd"] +=  (1.0/  2.0)*  WMNIJ_3[    "ijmn"]*   L_1(4)["mnklabcd"];
      DL_5(4)["ijklabcd"] +=                WAMEI_3[    "eiam"]*   L_1(4)["mjklbecd"];
      DL_5(4)["ijklabcd"] +=                  FAE_2[      "ea"]*   L_2(4)["ijklebcd"];
      DL_5(4)["ijklabcd"] -=                  FMI_2[      "im"]*   L_2(4)["mjklabcd"];
      DL_5(4)["ijklabcd"] +=  (1.0/  2.0)*  WABEF_2[    "efab"]*   L_2(4)["ijklefcd"];
      DL_5(4)["ijklabcd"] +=  (1.0/  2.0)*  WMNIJ_2[    "ijmn"]*   L_2(4)["mnklabcd"];
      DL_5(4)["ijklabcd"] +=                WAMEI_2[    "eiam"]*   L_2(4)["mjklbecd"];
      DL_5(4)["ijklabcd"] +=                  WMNEF[    "ijab"]*   L_4(2)[    "klcd"];
      DL_5(4)["ijklabcd"] +=                    FME[      "ia"]*   L_4(3)[  "jklbcd"];
      DL_5(4)["ijklabcd"] +=                  WAMEF[    "ejab"]*   L_4(3)[  "iklecd"];
      DL_5(4)["ijklabcd"] -=                  WMNEJ[    "ijam"]*   L_4(3)[  "mklbcd"];
      DL_5(4)["ijklabcd"] +=                  WMNEF[    "ijae"]*GAIJBCD_4[  "eklbcd"];
      DL_5(4)["ijklabcd"] -=                  WMNEF[    "mlcd"]*GIJKABL_4[  "ijkabm"];
      DL_5(4)["ijklabcd"] +=                    FAE[      "ea"]*   L_4(4)["ijklebcd"];
      DL_5(4)["ijklabcd"] -=                    FMI[      "im"]*   L_4(4)["mjklabcd"];
      DL_5(4)["ijklabcd"] +=  (1.0/  2.0)*    WABEF[    "efab"]*   L_4(4)["ijklefcd"];
      DL_5(4)["ijklabcd"] +=  (1.0/  2.0)*    WMNIJ[    "ijmn"]*   L_4(4)["mnklabcd"];
      DL_5(4)["ijklabcd"] +=                  WAMEI[    "eiam"]*   L_4(4)["mjklbecd"];

    L_5 = DL_5;
    L_5.weight(D);

    /***************************************************************************
     *           _
     * <0|(1+L)[[H,T^(1)],T^(3)]|0>
     *
     **************************************************************************/

    Z(3)["abcijk"]  =       FME_3[    "me"]*T_1(4)["abceijkm"];
    Z(3)["abcijk"] += 0.5*WAMEF_3[  "amef"]*T_1(4)["efbcimjk"];
    Z(3)["abcijk"] -= 0.5*WMNEJ_3[  "mnek"]*T_1(4)["abecijmn"];

    U E0113 = (1.0/36.0)*scalar(L(3)["mnoefg"]*Z(3)["efgmno"]);

    /***************************************************************************
     *               _
     * 1/2 <0|(1+L)[[H,T^(2)],T^(2)]|0>
     *
     **************************************************************************/

    WABEJ_2[  "abej"] +=         FME[    "me"]*T_2(2)[  "abmj"];

       Z(2)[  "abij"]  =       FAE_2[    "af"]*T_2(2)[  "fbij"];
       Z(2)[  "abij"] -=       FMI_2[    "ni"]*T_2(2)[  "abnj"];
       Z(2)[  "abij"] += 0.5*WMNIJ_2[  "mnij"]*T_2(2)[  "abmn"];
       Z(2)[  "abij"] += 0.5*WAMEI_2[  "amei"]*T_2(2)[  "ebjm"];

       Z(3)["abcijk"]  =     WABEJ_2[  "bcek"]*T_2(2)[  "aeij"];
       Z(3)["abcijk"] -=     WAMIJ_2[  "bmjk"]*T_2(2)[  "acim"];
       Z(3)["abcijk"] +=       FAE_2[    "ce"]*T_2(3)["abeijk"];
       Z(3)["abcijk"] -=       FMI_2[    "mk"]*T_2(3)["abcijm"];
       Z(3)["abcijk"] += 0.5*WABEF_2[  "abef"]*T_2(3)["efcijk"];
       Z(3)["abcijk"] += 0.5*WMNIJ_2[  "mnij"]*T_2(3)["abcmnk"];
       Z(3)["abcijk"] +=     WAMEI_2[  "amei"]*T_2(3)["ebcjmk"];

    WABEJ_2[  "abej"] -=         FME[    "me"]*T_2(2)[  "abmj"];

    U E0122 = (1.0/ 4.0)*scalar(L(2)[  "mnef"]*Z(2)[  "efmn"])+
              (1.0/36.0)*scalar(L(3)["mnoefg"]*Z(3)["efgmno"]);

    /***************************************************************************
     *           _
     * <0|L^(1)[[H,T^(1)],T^(2)]|0>
     *
     **************************************************************************/

    Z(4)["abcdijkl"]  =   WABCEJK_1["abcejk"]*T_2(2)[    "edil"];
    Z(4)["abcdijkl"] -=   WABMIJK_1["abmijk"]*T_2(2)[    "cdml"];
    Z(4)["abcdijkl"] +=       FAE_2[    "ae"]*T_1(4)["ebcdijkl"];
    Z(4)["abcdijkl"] -=       FMI_2[    "mi"]*T_1(4)["abcdmjkl"];
    Z(4)["abcdijkl"] += 0.5*WABEF_2[  "abef"]*T_1(4)["efcdijkl"];
    Z(4)["abcdijkl"] += 0.5*WMNIJ_2[  "mnij"]*T_1(4)["abcdmnkl"];
    Z(4)["abcdijkl"] +=     WAMEI_2[  "amei"]*T_1(4)["ebcdjmkl"];

    U E1112 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*Z(4)["efghmnop"]);

    /***************************************************************************
     *           _
     * <0|(1+L)[[H,T^(1)],T^(4)]|0>
     *
     **************************************************************************/

    Z(3)["abcijk"]  =       FME_4[    "me"]*T_1(4)["abceijkm"];
    Z(3)["abcijk"] += 0.5*WAMEF_4[  "amef"]*T_1(4)["efbcimjk"];
    Z(3)["abcijk"] -= 0.5*WMNEJ_4[  "mnek"]*T_1(4)["abecijmn"];

    U E0114 = (1.0/36.0)*scalar(L(3)["mnoefg"]*Z(3)["efgmno"]);

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

    U E1122 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*Z(4)["efghmnop"]);

    /***************************************************************************
     *           _
     * <0|L^(2)[[H,T^(1)],T^(2)]|0>
     *
     **************************************************************************/

    Z(4)["abcdijkl"]  =   WABCEJK_1["abcejk"]*T_2(2)[    "edil"];
    Z(4)["abcdijkl"] -=   WABMIJK_1["abmijk"]*T_2(2)[    "cdml"];
    Z(4)["abcdijkl"] +=       FAE_2[    "ae"]*T_1(4)["ebcdijkl"];
    Z(4)["abcdijkl"] -=       FMI_2[    "mi"]*T_1(4)["abcdmjkl"];
    Z(4)["abcdijkl"] += 0.5*WABEF_2[  "abef"]*T_1(4)["efcdijkl"];
    Z(4)["abcdijkl"] += 0.5*WMNIJ_2[  "mnij"]*T_1(4)["abcdmnkl"];
    Z(4)["abcdijkl"] +=     WAMEI_2[  "amei"]*T_1(4)["ebcdjmkl"];

    U E2112 = (1.0/576.0)*scalar(L_2(4)["mnopefgh"]*Z(4)["efghmnop"]);

    /***************************************************************************
     *
     * CCSDT(Q-2)
     *
     **************************************************************************/

    U E101 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*DT_1(4)["efghmnop"]);
    U E2 = E101;

    /***************************************************************************
     *
     * CCSDT(Q-3)
     *
     **************************************************************************/

    U E201 = (1.0/576.0)*scalar(L_2(4)["mnopefgh"]*DT_1(4)["efghmnop"]);
    U E102 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*DT_2(4)["efghmnop"]);

    printf("E201: %18.15f\n", E201);
    printf("E102: %18.15f\n", E102);
    printf("\n");

    U E3 = E102;

    /***************************************************************************
     *
     * CCSDT(Q-4)
     *
     **************************************************************************/

    U E301 =  (1.0/576.0)*scalar(L_3(4)["mnopefgh"]*DT_1(4)["efghmnop"]);
    U E202 =  (1.0/  1.0)*scalar(L_2(1)[      "me"]*DT_2(1)[      "em"])
             +(1.0/  4.0)*scalar(L_2(2)[    "mnef"]*DT_2(2)[    "efmn"])
             +(1.0/ 36.0)*scalar(L_2(3)[  "mnoefg"]*DT_2(3)[  "efgmno"])
             +(1.0/576.0)*scalar(L_2(4)["mnopefgh"]*DT_2(4)["efghmnop"]);
    U E103 =  (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*DT_3(4)["efghmnop"]);

    printf("E301: %18.15f\n", E301);
    printf("E202: %18.15f\n", E202);
    printf("E103: %18.15f\n", E103);
    printf("\n");

    U E4 = E103;

    /***************************************************************************
     *
     * CCSDT(Q-5)
     *
     **************************************************************************/

    U E302s = (1.0/  1.0)*scalar(L_3(1)[      "me"]*DT_2(1)[      "em"]);
    U E302d = (1.0/  4.0)*scalar(L_3(2)[    "mnef"]*DT_2(2)[    "efmn"]);
    U E302t = (1.0/ 36.0)*scalar(L_3(3)[  "mnoefg"]*DT_2(3)[  "efgmno"]);
    U E302q = (1.0/576.0)*scalar(L_3(4)["mnopefgh"]*DT_2(4)["efghmnop"]);
    U E203s = (1.0/  1.0)*scalar(L_2(1)[      "me"]*DT_3(1)[      "em"]);
    U E203d = (1.0/  4.0)*scalar(L_2(2)[    "mnef"]*DT_3(2)[    "efmn"]);
    U E203t = (1.0/ 36.0)*scalar(L_2(3)[  "mnoefg"]*DT_3(3)[  "efgmno"]);
    U E203q = (1.0/576.0)*scalar(L_2(4)["mnopefgh"]*DT_3(4)["efghmnop"]);

    U E401 = (1.0/576.0)*scalar(L_4(4)["mnopefgh"]*DT_1(4)["efghmnop"]);
    U E302 = E302s+E302d+E302t+E302q;
    U E203 = E203s+E203d+E203t+E203q;
    U E104 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*DT_4(4)["efghmnop"]);

    printf("E401:  %18.15f\n", E401);
    printf("E302s: %18.15f\n", E302s);
    printf("E302d: %18.15f\n", E302d);
    printf("E302t: %18.15f\n", E302t);
    printf("E302q: %18.15f\n", E302q);
    printf("E302:  %18.15f\n", E302);
    printf("E203s: %18.15f\n", E203s);
    printf("E203d: %18.15f\n", E203d);
    printf("E203t: %18.15f\n", E203t);
    printf("E203q: %18.15f\n", E203q);
    printf("E203:  %18.15f\n", E203);
    printf("E104:  %18.15f\n", E104);
    printf("\n");

    printf("E0113: %18.15f\n", E0113);
    printf("E0122: %18.15f\n", E0122);
    printf("E1112: %18.15f\n", E1112);
    printf("\n");

    printf("E401:             %18.15f\n", E401);
    printf("E302+E0113+E1112: %18.15f\n", E302+E0113+E1112);
    printf("\n");

    printf("E302:               %18.15f\n", E302);
    printf("E203+2*E0122+E1112: %18.15f\n", E203+2*E0122+E1112);
    printf("\n");

    printf("E203:             %18.15f\n", E203);
    printf("E104+E0113-E1112: %18.15f\n", E104+E0113-E1112);
    printf("\n");

    U E5 = E203 + E1112 + E0122;

    /***************************************************************************
     *
     * CCSDT(Q-6)
     *
     **************************************************************************/

    U E402s = (1.0/  1.0)*scalar(L_4(1)[      "me"]*DT_2(1)[      "em"]);
    U E402d = (1.0/  4.0)*scalar(L_4(2)[    "mnef"]*DT_2(2)[    "efmn"]);
    U E402t = (1.0/ 36.0)*scalar(L_4(3)[  "mnoefg"]*DT_2(3)[  "efgmno"]);
    U E402q = (1.0/576.0)*scalar(L_4(4)["mnopefgh"]*DT_2(4)["efghmnop"]);
    U E303s = (1.0/  1.0)*scalar(L_3(1)[      "me"]*DT_3(1)[      "em"]);
    U E303d = (1.0/  4.0)*scalar(L_3(2)[    "mnef"]*DT_3(2)[    "efmn"]);
    U E303t = (1.0/ 36.0)*scalar(L_3(3)[  "mnoefg"]*DT_3(3)[  "efgmno"]);
    U E303q = (1.0/576.0)*scalar(L_3(4)["mnopefgh"]*DT_3(4)["efghmnop"]);
    U E204s = (1.0/  1.0)*scalar(L_2(1)[      "me"]*DT_4(1)[      "em"]);
    U E204d = (1.0/  4.0)*scalar(L_2(2)[    "mnef"]*DT_4(2)[    "efmn"]);
    U E204t = (1.0/ 36.0)*scalar(L_2(3)[  "mnoefg"]*DT_4(3)[  "efgmno"]);
    U E204q = (1.0/576.0)*scalar(L_2(4)["mnopefgh"]*DT_4(4)["efghmnop"]);

    U E501 = (1.0/576.0)*scalar(L_5(4)["mnopefgh"]*DT_1(4)["efghmnop"]);
    U E402 = E402s+E402d+E402t+E402q;
    U E303 = E303s+E303d+E303t+E303q;
    U E204 = E204s+E204d+E204t+E204q;
    U E105 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*DT_5(4)["efghmnop"]);

    printf("E501:  %18.15f\n", E501);
    printf("E402s: %18.15f\n", E402s);
    printf("E402d: %18.15f\n", E402d);
    printf("E402t: %18.15f\n", E402t);
    printf("E402q: %18.15f\n", E402q);
    printf("E402:  %18.15f\n", E402);
    printf("E303s: %18.15f\n", E303s);
    printf("E303d: %18.15f\n", E303d);
    printf("E303t: %18.15f\n", E303t);
    printf("E303q: %18.15f\n", E303q);
    printf("E303:  %18.15f\n", E303);
    printf("E204s: %18.15f\n", E204s);
    printf("E204d: %18.15f\n", E204d);
    printf("E204t: %18.15f\n", E204t);
    printf("E204q: %18.15f\n", E204q);
    printf("E204:  %18.15f\n", E204);
    printf("E105:  %18.15f\n", E105);
    printf("\n");

    printf("E0114: %18.15f\n", E0114);
    printf("E0123: %18.15f\n", E0123);
    printf("E1113: %18.15f\n", E1113);
    printf("E1122: %18.15f\n", E1122);
    printf("E2112: %18.15f\n", E2112);
    printf("\n");

    printf("E501:                   %18.15f\n", E501);
    printf("E402+E0114+E1113+E2112: %18.15f\n", E402+E0114+E1113+E2112);
    printf("\n");

    printf("E402:                     %18.15f\n", E402);
    printf("E303+E0123+2*E1122+E2112: %18.15f\n", E303+E0123+2*E1122+E2112);
    printf("\n");

    printf("E303:                   %18.15f\n", E303);
    printf("E204+E0123+E1113-E2112: %18.15f\n", E204+E0123+E1113-E2112);
    printf("\n");

    printf("E204:                   %18.15f\n", E204);
    printf("E105+E0114-E1113-E1122: %18.15f\n", E105+E0114-E1113-E1122);
    printf("\n");

    U E6 = E204 + E1122 + E1113 + E0123;

    printf("CCSDT(Q-2): %18.15f\n", E2);
    printf("CCSDT(Q-3): %18.15f\n", E3);
    printf("CCSDT(Q-4): %18.15f\n", E4);
    printf("CCSDT(Q-5): %18.15f\n", E5);
    printf("CCSDT(Q-6): %18.15f\n", E6);

    this->put("E(2)", new U(E2));
    this->put("E(3)", new U(E3));
    this->put("E(4)", new U(E4));
    this->put("E(5)", new U(E5));
    this->put("E(6)", new U(E6));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSDT_Q_N);
REGISTER_TASK(aquarius::cc::CCSDT_Q_N<double>,"ccsdt(q-n)");
