#include "ccsd_tq_n_opt.hpp"

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
CCSD_TQ_N_opt<U>::CCSD_TQ_N_opt(const string& name, Config& config)
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
bool CCSD_TQ_N_opt<U>::run(task::TaskDAG& dag, const Arena& arena)
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

    ExcitationOperator<U,3> Z("Z", arena, occ, vrt);

    /***************************************************************************
     *
     * L^(1)
     *
     **************************************************************************/

    L_1(3)[  "ijkabc"]  = WMNEF["ijab"]*L(1)[  "kc"];
    L_1(3)[  "ijkabc"] +=   FME[  "ia"]*L(2)["jkbc"];
    L_1(3)[  "ijkabc"] += WAMEF["ekbc"]*L(2)["ijae"];
    L_1(3)[  "ijkabc"] -= WMNEJ["ijam"]*L(2)["mkbc"];

    L_1(4)["ijklabcd"]  = WMNEF["ijab"]*L(2)["klcd"];

    L_1.weight(D);

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

     T_1(3)[  "abcijk"]  =   WTWABEJ[  "bcek"]*T(2)["aeij"];
     T_1(3)[  "abcijk"] -=     WAMIJ[  "bmjk"]*T(2)["acim"];

    U E2_3 = (1.0/ 36.0)*scalar(L_1(3)[  "mnoefg"]*T_1(3)[  "efgmno"]);

     T_1(4)["abcdijkl"]  =   WABCEJK["abcejk"]*T(2)["edil"];
     T_1(4)["abcdijkl"] -=   WABMIJK["abmijk"]*T(2)["cdml"];

    U E2_4 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_1(4)["efghmnop"]);

    T_1.weight(D);

    /***************************************************************************
     *
     * T^(2)
     *
     **************************************************************************/

      WAMIJ_1[    "amij"]  =   0.5*WMNEF[  "mnef"]*T_1(3)[  "aefijn"];
      WABEJ_1[    "abej"]  =  -0.5*WMNEF[  "mnef"]*T_1(3)[  "afbmnj"];

    WABMIJK_1[  "abmijk"] +=   0.5*WAMEF[  "bmef"]*T_1(3)[  "aefijk"];

       T_2(1)[      "ai"]  =  0.25*WMNEF[  "mnef"]*T_1(3)[  "aefimn"];

       T_2(2)[    "abij"]  =   0.5*WAMEF[  "bmef"]*T_1(3)[  "aefijm"];
       T_2(2)[    "abij"] -=   0.5*WMNEJ[  "mnej"]*T_1(3)[  "abeinm"];
       T_2(2)[    "abij"] +=         FME[    "me"]*T_1(3)[  "abeijm"];
       T_2(2)[    "abij"] +=  0.25*WMNEF[  "mnef"]*T_1(4)["abefijmn"];

       T_2(3)[  "abcijk"]  =     WABEJ_1[  "bcek"]*  T(2)[    "aeij"];
       T_2(3)[  "abcijk"] -=     WAMIJ_1[  "bmjk"]*  T(2)[    "acim"];
       T_2(3)[  "abcijk"] +=         FAE[    "ce"]*T_1(3)[  "abeijk"];
       T_2(3)[  "abcijk"] -=         FMI[    "mk"]*T_1(3)[  "abcijm"];
       T_2(3)[  "abcijk"] +=   0.5*WABEF[  "abef"]*T_1(3)[  "efcijk"];
       T_2(3)[  "abcijk"] +=   0.5*WMNIJ[  "mnij"]*T_1(3)[  "abcmnk"];
       T_2(3)[  "abcijk"] +=       WAMEI[  "amei"]*T_1(3)[  "ebcjmk"];

    U E3_3_3 = (1.0/ 36.0)*scalar(L_1(3)[  "mnoefg"]*T_2(3)[  "efgmno"]);

       T_2(3)[  "abcijk"] +=         FME[    "me"]*T_1(4)["abceijkm"];
       T_2(3)[  "abcijk"] +=   0.5*WAMEF[  "amef"]*T_1(4)["efbcimjk"];
       T_2(3)[  "abcijk"] -=   0.5*WMNEJ[  "mnek"]*T_1(4)["abecijmn"];

    U E3_3_34 = (1.0/ 36.0)*scalar(L_1(3)[  "mnoefg"]*T_2(3)[  "efgmno"]);
    U E3_3_4 = E3_3_34-E3_3_3;

       T_2(4)["abcdijkl"]  =  -WABMIJK_1["abmijk"]*  T(2)[    "cdml"];
       T_2(4)["abcdijkl"] +=       WABEJ[  "abej"]*T_1(3)[  "ecdikl"];
       T_2(4)["abcdijkl"] -=       WAMIJ[  "amij"]*T_1(3)[  "bcdmkl"];
       T_2(4)["abcdijkl"] += 0.5*WAMNIJK["amnijk"]*T_1(3)[  "bcdmnl"];
       T_2(4)["abcdijkl"] -=     WABMEJI["abmeji"]*T_1(3)[  "ecdmkl"];

    U E3_4_3 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_2(4)["efghmnop"]);

    WABCEJK_1[  "abcejk"]  =  -0.5*WMNEF[  "mnef"]*T_1(4)["abcfmjkn"];
    WABMIJK_1[  "abmijk"]  =   0.5*WMNEF[  "mnef"]*T_1(4)["abefijkn"];

       T_2(4)["abcdijkl"] +=   WABCEJK_1["abcejk"]*  T(2)[    "edil"];
       T_2(4)["abcdijkl"] -=   WABMIJK_1["abmijk"]*  T(2)[    "cdml"];
       T_2(4)["abcdijkl"] +=         FAE[    "ae"]*T_1(4)["ebcdijkl"];
       T_2(4)["abcdijkl"] -=         FMI[    "mi"]*T_1(4)["abcdmjkl"];
       T_2(4)["abcdijkl"] +=   0.5*WABEF[  "abef"]*T_1(4)["efcdijkl"];
       T_2(4)["abcdijkl"] +=   0.5*WMNIJ[  "mnij"]*T_1(4)["abcdmnkl"];
       T_2(4)["abcdijkl"] +=       WAMEI[  "amei"]*T_1(4)["ebcdjmkl"];

    U E3_4_34 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_2(4)["efghmnop"]);
    U E3_4_4 = E3_4_34-E3_4_3;

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
      WAMIJ_2[    "amij"] +=         FME_2[    "me"]*  T(2)[    "aeij"];

      WABEJ_2[    "abej"]  =         WAMEF[  "amef"]*T_2(2)[    "fbmj"];
      WABEJ_2[    "abej"] +=     0.5*WMNEJ[  "mnej"]*T_2(2)[    "abmn"];
      WABEJ_2[    "abej"] +=         WABEF[  "abef"]*T_2(1)[      "fj"];
      WABEJ_2[    "abej"] -=         WAMEI[  "amej"]*T_2(1)[      "bm"];

      WABCEJK[  "abcejk"]  =         WABEF[  "abef"]*  T(2)[    "fcjk"];
      WABCEJK[  "abcejk"] -=         WAMEI[  "amek"]*  T(2)[    "bcjm"];

      WABMIJK[  "abmijk"]  =         WAMEI[  "amek"]*  T(2)[    "ebij"];
      WABMIJK[  "abmijk"] -=         WMNIJ[  "nmjk"]*  T(2)[    "abin"];

      WAMNIJK[  "amnijk"]  =         WMNEJ[  "mnek"]*  T(2)[    "aeij"];

      WABMEJI[  "abmeji"]  =         WAMEF[  "amef"]*  T(2)[    "bfji"];
      WABMEJI[  "abmeji"] -=         WMNEJ[  "nmei"]*  T(2)[    "abnj"];

      WAMEI_2[    "amei"]  =         WMNEF[  "mnef"]*T_2(2)[    "afni"];
      WAMEI_2[    "amei"] +=         WAMEF[  "amef"]*T_2(1)[      "fi"];
      WAMEI_2[    "amei"] -=         WMNEJ[  "nmei"]*T_2(1)[      "an"];

      WMNIJ_2[    "mnij"]  =     0.5*WMNEF[  "mnef"]*T_2(2)[    "efij"];
      WMNIJ_2[    "mnij"] +=         WMNEJ[  "mnej"]*T_2(1)[      "ei"];

      WABEF_2[    "abef"]  =     0.5*WMNEF[  "mnef"]*T_2(2)[    "abmn"];
      WABEF_2[    "abef"] -=         WAMEF[  "amef"]*T_2(1)[      "bm"];

    WABCEJK_2[  "abcejk"]  =   0.5*WABEF_2[  "abef"]*  T(2)[    "fcjk"];

    WABMIJK_2[  "abmijk"]  =       WAMEI_2[  "amek"]*  T(2)[    "ebij"];
    WABMIJK_2[  "abmijk"] -=   0.5*WMNIJ_2[  "nmjk"]*  T(2)[    "abin"];

       T_3(3)[  "abcijk"]  =       WABEJ_2[  "bcek"]*  T(2)[    "aeij"];
       T_3(3)[  "abcijk"] -=       WAMIJ_2[  "bmjk"]*  T(2)[    "acim"];
       T_3(3)[  "abcijk"] +=         WABEJ[  "bcek"]*T_2(2)[    "aeij"];
       T_3(3)[  "abcijk"] -=         WAMIJ[  "bmjk"]*T_2(2)[    "acim"];

    U E4_3_12 = (1.0/ 36.0)*scalar(L_1(3)[  "mnoefg"]*T_3(3)[  "efgmno"]);

      WAMIJ_2[    "amij"]  =     0.5*WMNEF[  "mnef"]*T_2(3)[  "aefijn"];
      WABEJ_2[    "abej"]  =    -0.5*WMNEF[  "mnef"]*T_2(3)[  "afbmnj"];

       T_3(3)[  "abcijk"] +=       WABEJ_2[  "bcek"]*  T(2)[    "aeij"];
       T_3(3)[  "abcijk"] -=       WAMIJ_2[  "bmjk"]*  T(2)[    "acim"];
       T_3(3)[  "abcijk"] +=           FAE[    "ce"]*T_2(3)[  "abeijk"];
       T_3(3)[  "abcijk"] -=           FMI[    "mk"]*T_2(3)[  "abcijm"];
       T_3(3)[  "abcijk"] +=     0.5*WABEF[  "abef"]*T_2(3)[  "efcijk"];
       T_3(3)[  "abcijk"] +=     0.5*WMNIJ[  "mnij"]*T_2(3)[  "abcmnk"];
       T_3(3)[  "abcijk"] +=         WAMEI[  "amei"]*T_2(3)[  "ebcjmk"];

    U E4_3_123 = (1.0/ 36.0)*scalar(L_1(3)[  "mnoefg"]*T_3(3)[  "efgmno"]);
    U E4_3_3 = E4_3_123-E4_3_12;

       T_3(3)[  "abcijk"] +=           FME[    "me"]*T_2(4)["abceijkm"];
       T_3(3)[  "abcijk"] +=     0.5*WAMEF[  "amef"]*T_2(4)["efbcimjk"];
       T_3(3)[  "abcijk"] -=     0.5*WMNEJ[  "mnek"]*T_2(4)["abecijmn"];

    U E4_3_1234 = (1.0/ 36.0)*scalar(L_1(3)[  "mnoefg"]*T_3(3)[  "efgmno"]);
    U E4_3_4 = E4_3_1234-E4_3_123;

       T_3(4)["abcdijkl"]  =     WABCEJK_2["abcejk"]*  T(2)[    "edil"];
       T_3(4)["abcdijkl"] -=     WABMIJK_2["abmijk"]*  T(2)[    "cdml"];
       T_3(4)["abcdijkl"] +=       WABCEJK["abcejk"]*T_2(2)[    "edil"];
       T_3(4)["abcdijkl"] -=       WABMIJK["abmijk"]*T_2(2)[    "cdml"];

    U E4_4_12 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_3(4)["efghmnop"]);

    WABMIJK_2[  "abmijk"]  =     0.5*WAMEF[  "bmef"]*T_2(3)[  "aefijk"];

       T_3(4)["abcdijkl"] -=     WABMIJK_2["abmijk"]*  T(2)[    "cdml"];
       T_3(4)["abcdijkl"] +=         WABEJ[  "abej"]*T_2(3)[  "ecdikl"];
       T_3(4)["abcdijkl"] -=         WAMIJ[  "amij"]*T_2(3)[  "bcdmkl"];
       T_3(4)["abcdijkl"] +=   0.5*WAMNIJK["amnijk"]*T_2(3)[  "bcdmnl"];
       T_3(4)["abcdijkl"] -=       WABMEJI["abmeji"]*T_2(3)[  "ecdmkl"];

    U E4_4_123 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_3(4)["efghmnop"]);
    U E4_4_3 = E4_4_123-E4_4_12;

    WABCEJK_2[  "abcejk"]  =    -0.5*WMNEF[  "mnef"]*T_2(4)["abcfmjkn"];
    WABMIJK_2[  "abmijk"]  =     0.5*WMNEF[  "mnef"]*T_2(4)["abefijkn"];

       T_3(4)["abcdijkl"] +=     WABCEJK_2["abcejk"]*  T(2)[    "edil"];
       T_3(4)["abcdijkl"] -=     WABMIJK_2["abmijk"]*  T(2)[    "cdml"];
       T_3(4)["abcdijkl"] +=           FAE[    "ae"]*T_2(4)["ebcdijkl"];
       T_3(4)["abcdijkl"] -=           FMI[    "mi"]*T_2(4)["abcdmjkl"];
       T_3(4)["abcdijkl"] +=     0.5*WABEF[  "abef"]*T_2(4)["efcdijkl"];
       T_3(4)["abcdijkl"] +=     0.5*WMNIJ[  "mnij"]*T_2(4)["abcdmnkl"];
       T_3(4)["abcdijkl"] +=         WAMEI[  "amei"]*T_2(4)["ebcdjmkl"];

    U E4_4_1234 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_3(4)["efghmnop"]);
    U E4_4_4 = E4_4_1234-E4_4_123;

    /***************************************************************************
     *
     * Complete W^(2)
     *
     **************************************************************************/

    WMNEJ_2["mnej"]  =  WMNEF["mnef"]*T_2(1)[  "fj"];
    WAMEF_2["amef"]  = -WMNEF["nmef"]*T_2(1)[  "an"];
    WABEJ_2["abej"] -=  FME_2[  "me"]*  T(2)["abmj"];

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

    WAMNIJK_1["amnijk"]  =            0.5*WMNEF[  "mnef"]*T_1(3)[  "aefijk"];
    WABMEJI_1["abmeji"]  =            0.5*WMNEF[  "mnef"]*T_1(3)[  "abfnji"];

         Z(3)["abcijk"]  = (1.0/ 2.0)*  WABEJ_1[  "efcm"]*L_1(4)["ijkmabef"];
         Z(3)["abcijk"] -= (1.0/ 2.0)*  WAMIJ_1[  "eknm"]*L_1(4)["ijmnabce"];
         Z(3)["abcijk"] -= (1.0/ 4.0)*WABMEJI_1["efkcnm"]*L_1(4)["ijmnabef"];
         Z(3)["abcijk"] += (1.0/ 6.0)*WAMNIJK_1["eijmno"]*L_1(4)["mnokeabc"];

    U E1111 = (1.0/36.0)*scalar(T_1(3)["efgmno"]*Z(3)["efgmno"]);

    U E2 = E2_3+E2_4;

    printf("\n");
    printf("<0|L3^[1]Hbar|0>: %18.15f\n", E2_3);
    printf("<0|L4^[1]Hbar|0>: %18.15f\n", E2_4);
    printf("CCSD(TQ-2) Total: %18.15f\n", E2);

    U E3 = E3_3_3+E3_3_4+E3_4_3+E3_4_4;

    printf("\n");
    printf("<0|L3^[1][Hbar,T3^[1]]|0>: %18.15f\n", E3_3_3);
    printf("<0|L3^[1][Hbar,T4^[1]]|0>: %18.15f\n", E3_3_4);
    printf("<0|L4^[1][Hbar,T3^[1]]|0>: %18.15f\n", E3_4_3);
    printf("<0|L4^[1][Hbar,T4^[1]]|0>: %18.15f\n", E3_4_4);
    printf("CCSD(TQ-3) Total:          %18.15f\n", E3);

    U E4 = E4_3_12+E4_3_3+E4_3_4+E4_4_12+E4_4_3+E4_4_4+E0112+E1111;

    printf("\n");
    printf("<0|L3^[1][Hbar,T1^[2]+T2^[2]]|0>:       %18.15f\n", E4_3_12);
    printf("<0|L3^[1][Hbar,T3^[2]]|0>:              %18.15f\n", E4_3_3);
    printf("<0|L3^[1][Hbar,T4^[2]]|0>:              %18.15f\n", E4_3_4);
    printf("<0|L4^[1][Hbar,T1^[2]+T2^[2]]|0>:       %18.15f\n", E4_4_12);
    printf("<0|L4^[1][Hbar,T3^[2]]|0>:              %18.15f\n", E4_4_3);
    printf("<0|L4^[1][Hbar,T4^[2]]|0>:              %18.15f\n", E4_4_4);
    printf("<0|(1+L)[[Hbar,T3^[1]],T1^[2]]|0>:      %18.15f\n", E0112);
    printf("1/2 <0|L4^[1][[Hbar,T3^[1]],T3^[1]]|0>: %18.15f\n", E1111);
    printf("CCSD(TQ-4) Total:                       %18.15f\n", E4);

    printf("\n");

    this->put("E(2)", new U(E2));
    this->put("E(3)", new U(E3));
    this->put("E(4)", new U(E4));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSD_TQ_N_opt);
REGISTER_TASK(aquarius::cc::CCSD_TQ_N_opt<double>,"ccsd(tq-n)_opt");
