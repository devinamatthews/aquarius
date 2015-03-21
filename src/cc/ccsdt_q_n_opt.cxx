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

    DeexcitationOperator<U,4> L_1("L^(1)", arena, occ, vrt);

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

    SpinorbitalTensor<U> GAIJBCD_2("G(aij,bcd)", arena, group, {vrt,occ}, {1,2}, {3,0});
    SpinorbitalTensor<U> GIJKABL_2("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1});

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

      WABCEJK[  "abcejk"]  =         WABEF[  "abef"]*  T(2)[    "fcjk"];
      WABCEJK[  "abcejk"] -=         WAMEI[  "amek"]*  T(2)[    "bcjm"];
      WABCEJK[  "abcejk"] +=     0.5*WMNEJ[  "mnej"]*  T(3)[  "abcmnk"];
      WABCEJK[  "abcejk"] +=         WAMEF[  "amef"]*  T(3)[  "fbcmjk"];

      WABMIJK[  "abmijk"]  =         WAMEI[  "amek"]*  T(2)[    "ebij"];
      WABMIJK[  "abmijk"] -=         WMNIJ[  "nmjk"]*  T(2)[    "abin"];
      WABMIJK[  "abmijk"] +=     0.5*WAMEF[  "bmef"]*  T(3)[  "aefijk"];
      WABMIJK[  "abmijk"] +=         WMNEJ[  "nmek"]*  T(3)[  "abeijn"];

      WAMNIJK[  "amnijk"]  =         WMNEJ[  "mnek"]*  T(2)[    "aeij"];
      WAMNIJK[  "amnijk"] +=     0.5*WMNEF[  "mnef"]*  T(3)[  "aefijk"];

      WABMEJI[  "abmeji"]  =         WAMEF[  "amef"]*  T(2)[    "bfji"];
      WABMEJI[  "abmeji"] -=         WMNEJ[  "nmei"]*  T(2)[    "abnj"];
      WABMEJI[  "abmeji"] +=         WMNEF[  "mnef"]*  T(3)[  "abfnji"];

      WAMIJ_2[    "amij"]  =         WMNEJ[  "nmej"]*T_2(2)[    "aein"];
      WAMIJ_2[    "amij"] +=           FME[    "me"]*T_2(2)[    "aeij"];
      WAMIJ_2[    "amij"] +=     0.5*WAMEF[  "amef"]*T_2(2)[    "efij"];
      WAMIJ_2[    "amij"] +=     0.5*WMNEF[  "mnef"]*T_2(3)[  "aefijn"];

      WABEJ_2[    "abej"]  =         WAMEF[  "amef"]*T_2(2)[    "fbmj"];
      WABEJ_2[    "abej"] -=           FME[    "me"]*T_2(2)[    "abmj"];
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

    WABCEJK_2[  "abcejk"]  =    -0.5*WMNEF[  "mnef"]*T_2(4)["abcfmjkn"];
    WABMIJK_2[  "abmijk"]  =     0.5*WMNEF[  "mnef"]*T_2(4)["abefijkn"];

       T_3(4)["abcdijkl"] +=     WABCEJK_2["abcejk"]*  T(2)[    "edil"];
       T_3(4)["abcdijkl"] -=     WABMIJK_2["abmijk"]*  T(2)[    "cdml"];
       T_3(4)["abcdijkl"] +=           FAE[    "ae"]*T_2(4)["ebcdijkl"];
       T_3(4)["abcdijkl"] -=           FMI[    "mi"]*T_2(4)["abcdmjkl"];
       T_3(4)["abcdijkl"] +=     0.5*WABEF[  "abef"]*T_2(4)["efcdijkl"];
       T_3(4)["abcdijkl"] +=     0.5*WMNIJ[  "mnij"]*T_2(4)["abcdmnkl"];
       T_3(4)["abcdijkl"] +=         WAMEI[  "amei"]*T_2(4)["ebcdjmkl"];

    U E4_234 = (1.0/576.0)*scalar(L_1(4)["mnopefgh"]*T_3(4)["efghmnop"]);
    U E4_4 = E4_234-E4_23;

    printf("\n");
    printf("CCSDT(Q-2): %18.15f\n", E2);

    printf("\n");
    printf("CCSDT(Q-3): %18.15f\n", E3);

    U E4 = E4_23+E4_4;

    printf("\n");
    printf("<0|L4^[1][Hbar,T2^[2]+T3^[2]]|0>: %18.15f\n", E4_23);
    printf("<0|L4^[1][Hbar,T4^[2]]|0>:        %18.15f\n", E4_4);
    printf("CCSDT(Q-4) Total:                 %18.15f\n", E4);

    printf("\n");

    this->put("E(2)", new U(E2));
    this->put("E(3)", new U(E3));
    this->put("E(4)", new U(E4));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSDT_Q_N_opt);
REGISTER_TASK(aquarius::cc::CCSDT_Q_N_opt<double>,"ccsdt(q-n)_opt");
