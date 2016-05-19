#include "ccsd_t_n.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;

namespace aquarius
{
namespace cc
{

template <typename U>
CCSD_T_N<U>::CCSD_T_N(const string& name, Config& config)
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
    this->addProduct(Product("double", "E(5)", reqs));
    this->addProduct(Product("double", "E(6)", reqs));
}

template <typename U>
bool CCSD_T_N<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U>>("H");
    const STTwoElectronOperator<U>& Hbar = this->template get<STTwoElectronOperator<U>>("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

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

    ExcitationOperator<U,3> T_1("T^(1)", arena, occ, vrt);
    ExcitationOperator<U,3> T_2("T^(2)", arena, occ, vrt);
    ExcitationOperator<U,3> T_3("T^(3)", arena, occ, vrt);
    ExcitationOperator<U,3> T_4("T^(4)", arena, occ, vrt);
    ExcitationOperator<U,3> T_5("T^(5)", arena, occ, vrt);

    DeexcitationOperator<U,3> L_1("L^(1)", arena, occ, vrt);
    DeexcitationOperator<U,3> L_2("L^(2)", arena, occ, vrt);
    DeexcitationOperator<U,3> L_3("L^(3)", arena, occ, vrt);
    DeexcitationOperator<U,3> L_4("L^(4)", arena, occ, vrt);
    DeexcitationOperator<U,3> L_5("L^(5)", arena, occ, vrt);

    ExcitationOperator<U,3> DT_1("DT^(1)", arena, occ, vrt);
    ExcitationOperator<U,3> DT_2("DT^(2)", arena, occ, vrt);
    ExcitationOperator<U,3> DT_3("DT^(3)", arena, occ, vrt);
    ExcitationOperator<U,3> DT_4("DT^(4)", arena, occ, vrt);
    ExcitationOperator<U,3> DT_5("DT^(5)", arena, occ, vrt);

    DeexcitationOperator<U,3> DL_1("DL^(1)", arena, occ, vrt);
    DeexcitationOperator<U,3> DL_2("DL^(2)", arena, occ, vrt);
    DeexcitationOperator<U,3> DL_3("DL^(3)", arena, occ, vrt);
    DeexcitationOperator<U,3> DL_4("DL^(4)", arena, occ, vrt);
    DeexcitationOperator<U,3> DL_5("DL^(5)", arena, occ, vrt);

    SpinorbitalTensor<U> WTWABEJ(WABEJ);
    WTWABEJ["abej"] += FME["me"]*T(2)["abmj"];

    SpinorbitalTensor<U> WABEJ_1(WABEJ);
    SpinorbitalTensor<U> WAMIJ_1(WAMIJ);

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

    SpinorbitalTensor<U> FME_3(FME);
    SpinorbitalTensor<U> FMI_3(FMI);
    SpinorbitalTensor<U> FAE_3(FAE);
    SpinorbitalTensor<U> WAMEI_3(WAMEI);
    SpinorbitalTensor<U> WMNIJ_3(WMNIJ);
    SpinorbitalTensor<U> WABEF_3(WABEF);
    SpinorbitalTensor<U> WABEJ_3(WABEJ);
    SpinorbitalTensor<U> WAMIJ_3(WAMIJ);
    SpinorbitalTensor<U> WMNEJ_3(WMNEJ);
    SpinorbitalTensor<U> WAMEF_3(WAMEF);

    SpinorbitalTensor<U> FME_4(FME);
    SpinorbitalTensor<U> WABEJ_4(WABEJ);
    SpinorbitalTensor<U> WAMIJ_4(WAMIJ);
    SpinorbitalTensor<U> WMNEJ_4(WMNEJ);
    SpinorbitalTensor<U> WAMEF_4(WAMEF);

    SpinorbitalTensor<U> DIJ_1(FMI);
    SpinorbitalTensor<U> DAB_1(FAE);
    SpinorbitalTensor<U> DAI_1(T(1));
    SpinorbitalTensor<U> GIJAK_1(WMNEJ);
    SpinorbitalTensor<U> GAIBC_1(WAMEF);

    SpinorbitalTensor<U> DIJ_2(FMI);
    SpinorbitalTensor<U> DAB_2(FAE);
    SpinorbitalTensor<U> DTWIJ_2(FMI);
    SpinorbitalTensor<U> DTWAB_2(FAE);
    SpinorbitalTensor<U> DAI_2(T(1));
    SpinorbitalTensor<U> GIJKL_2(WMNIJ);
    SpinorbitalTensor<U> GAIBJ_2(WAMEI);
    SpinorbitalTensor<U> GABCD_2(WABEF);
    SpinorbitalTensor<U> GIJAK_2(WMNEJ);
    SpinorbitalTensor<U> GAIBC_2(WAMEF);

    SpinorbitalTensor<U> DIJ_3(FMI);
    SpinorbitalTensor<U> DAB_3(FAE);
    SpinorbitalTensor<U> DTWIJ_3(FMI);
    SpinorbitalTensor<U> DTWAB_3(FAE);
    SpinorbitalTensor<U> DAI_3(T(1));
    SpinorbitalTensor<U> GIJKL_3(WMNIJ);
    SpinorbitalTensor<U> GAIBJ_3(WAMEI);
    SpinorbitalTensor<U> GABCD_3(WABEF);
    SpinorbitalTensor<U> GIJAK_3(WMNEJ);
    SpinorbitalTensor<U> GAIBC_3(WAMEF);

    SpinorbitalTensor<U> GIJAK_4(WMNEJ);
    SpinorbitalTensor<U> GAIBC_4(WAMEF);

    SpinorbitalTensor<U> Tau(T(2));
    SpinorbitalTensor<U> DIJ(FMI);
    SpinorbitalTensor<U> DAB(FAE);
    SpinorbitalTensor<U> FTWMI_2(FMI);
    SpinorbitalTensor<U> FTWAE_2(FAE);
    SpinorbitalTensor<U> FTWMI_3(FMI);
    SpinorbitalTensor<U> FTWAE_3(FAE);
    SpinorbitalTensor<U> WTWAMEI_2(WAMEI);
    SpinorbitalTensor<U> WTWABEJ_2(WABEJ);
    SpinorbitalTensor<U> WTWAMIJ_2(WAMIJ);
    SpinorbitalTensor<U> WTWABEJ_3(WABEJ);
    SpinorbitalTensor<U> WTWAMIJ_3(WAMIJ);
    SpinorbitalTensor<U> WTWABEJ_4(WABEJ);
    SpinorbitalTensor<U> WTWAMIJ_4(WAMIJ);
    ExcitationOperator<U,3> Z("Z", arena, occ, vrt);

    /***************************************************************************
     *
     * T^(1)
     *
     **************************************************************************/

    DT_1(3)["abcijk"]  = WTWABEJ["bcek"]*T(2)["aeij"];
    DT_1(3)["abcijk"] -=   WAMIJ["bmjk"]*T(2)["acim"];

    T_1 = DT_1;
    T_1.weight(D);

    /***************************************************************************
     *
     * T^(2)
     *
     **************************************************************************/

    WAMIJ_1[  "amij"]  =  0.5*WMNEF["mnef"]*T_1(3)["aefijn"];
    WABEJ_1[  "abej"]  = -0.5*WMNEF["mnef"]*T_1(3)["afbmnj"];

    DT_2(1)[    "ai"]  = 0.25*WMNEF["mnef"]*T_1(3)["aefimn"];

    DT_2(2)[  "abij"]  =  0.5*WAMEF["bmef"]*T_1(3)["aefijm"];
    DT_2(2)[  "abij"] -=  0.5*WMNEJ["mnej"]*T_1(3)["abeinm"];
    DT_2(2)[  "abij"] +=        FME[  "me"]*T_1(3)["abeijm"];

    DT_2(3)["abcijk"]  =    WABEJ_1["bcek"]*  T(2)[  "aeij"];
    DT_2(3)["abcijk"] -=    WAMIJ_1["bmjk"]*  T(2)[  "acim"];
    DT_2(3)["abcijk"] +=        FAE[  "ce"]*T_1(3)["abeijk"];
    DT_2(3)["abcijk"] -=        FMI[  "mk"]*T_1(3)["abcijm"];
    DT_2(3)["abcijk"] +=  0.5*WABEF["abef"]*T_1(3)["efcijk"];
    DT_2(3)["abcijk"] +=  0.5*WMNIJ["mnij"]*T_1(3)["abcmnk"];
    DT_2(3)["abcijk"] +=      WAMEI["amei"]*T_1(3)["ebcjmk"];

    T_2 = DT_2;
    T_2.weight(D);

    /***************************************************************************
     *
     * T^(3)
     *
     **************************************************************************/

      FME_2[    "me"]  =       WMNEF["mnef"]*T_2(1)[    "fn"];

      FMI_2[    "mi"]  =   0.5*WMNEF["mnef"]*T_2(2)[  "efin"];
      FMI_2[    "mi"] +=       WMNEJ["nmfi"]*T_2(1)[    "fn"];

      FAE_2[    "ae"]  =  -0.5*WMNEF["mnef"]*T_2(2)[  "afmn"];
      FAE_2[    "ae"] +=       WAMEF["amef"]*T_2(1)[    "fm"];

    WAMIJ_2[  "amij"]  =       WMNEJ["nmej"]*T_2(2)[  "aein"];
    WAMIJ_2[  "amij"] +=   0.5*WAMEF["amef"]*T_2(2)[  "efij"];
    WAMIJ_2[  "amij"] -=       WMNIJ["nmij"]*T_2(1)[    "an"];
    WAMIJ_2[  "amij"] +=       WAMEI["amej"]*T_2(1)[    "ei"];
    WAMIJ_2[  "amij"] +=   0.5*WMNEF["mnef"]*T_2(3)["aefijn"];
    WAMIJ_2[  "amij"] +=       FME_2[  "me"]*  T(2)[  "aeij"];

    WABEJ_2[  "abej"]  =       WAMEF["amef"]*T_2(2)[  "fbmj"];
    WABEJ_2[  "abej"] +=   0.5*WMNEJ["mnej"]*T_2(2)[  "abmn"];
    WABEJ_2[  "abej"] +=       WABEF["abef"]*T_2(1)[    "fj"];
    WABEJ_2[  "abej"] -=       WAMEI["amej"]*T_2(1)[    "bm"];
    WABEJ_2[  "abej"] -=   0.5*WMNEF["mnef"]*T_2(3)["afbmnj"];

    DT_3(1)[    "ai"]  =         FAE[  "ae"]*T_2(1)[    "ei"];
    DT_3(1)[    "ai"] -=         FMI[  "mi"]*T_2(1)[    "am"];
    DT_3(1)[    "ai"] -=       WAMEI["amei"]*T_2(1)[    "em"];
    DT_3(1)[    "ai"] +=         FME[  "me"]*T_2(2)[  "aeim"];
    DT_3(1)[    "ai"] +=   0.5*WAMEF["amef"]*T_2(2)[  "efim"];
    DT_3(1)[    "ai"] -=   0.5*WMNEJ["mnei"]*T_2(2)[  "eamn"];
    DT_3(1)[    "ai"] +=  0.25*WMNEF["mnef"]*T_2(3)["aefimn"];

    DT_3(2)[  "abij"]  =       FAE_2[  "af"]*  T(2)[  "fbij"];
    DT_3(2)[  "abij"] -=       FMI_2[  "ni"]*  T(2)[  "abnj"];
    DT_3(2)[  "abij"] +=       WABEJ["abej"]*T_2(1)[    "ei"];
    DT_3(2)[  "abij"] -=       WAMIJ["amij"]*T_2(1)[    "bm"];
    DT_3(2)[  "abij"] +=         FAE[  "af"]*T_2(2)[  "fbij"];
    DT_3(2)[  "abij"] -=         FMI[  "ni"]*T_2(2)[  "abnj"];
    DT_3(2)[  "abij"] +=   0.5*WABEF["abef"]*T_2(2)[  "efij"];
    DT_3(2)[  "abij"] +=   0.5*WMNIJ["mnij"]*T_2(2)[  "abmn"];
    DT_3(2)[  "abij"] +=       WAMEI["amei"]*T_2(2)[  "ebjm"];
    DT_3(2)[  "abij"] +=   0.5*WAMEF["bmef"]*T_2(3)["aefijm"];
    DT_3(2)[  "abij"] -=   0.5*WMNEJ["mnej"]*T_2(3)["abeinm"];
    DT_3(2)[  "abij"] +=         FME[  "me"]*T_2(3)["abeijm"];

    DT_3(3)["abcijk"]  =     WABEJ_2["bcek"]*  T(2)[  "aeij"];
    DT_3(3)["abcijk"] -=     WAMIJ_2["bmjk"]*  T(2)[  "acim"];
    DT_3(3)["abcijk"] +=       WABEJ["bcek"]*T_2(2)[  "aeij"];
    DT_3(3)["abcijk"] -=       WAMIJ["bmjk"]*T_2(2)[  "acim"];
    DT_3(3)["abcijk"] +=         FAE[  "ce"]*T_2(3)["abeijk"];
    DT_3(3)["abcijk"] -=         FMI[  "mk"]*T_2(3)["abcijm"];
    DT_3(3)["abcijk"] +=   0.5*WABEF["abef"]*T_2(3)["efcijk"];
    DT_3(3)["abcijk"] +=   0.5*WMNIJ["mnij"]*T_2(3)["abcmnk"];
    DT_3(3)["abcijk"] +=       WAMEI["amei"]*T_2(3)["ebcjmk"];

    T_3 = DT_3;
    T_3.weight(D);

    /***************************************************************************
     *
     * T^(4)
     *
     **************************************************************************/

    FTWMI_2 = FMI_2;
    FTWAE_2 = FAE_2;

      FMI_2[    "mi"] +=         FME[  "me"]*T_2(1)[    "ei"];
      FAE_2[    "ae"] -=         FME[  "me"]*T_2(1)[    "am"];

    WAMEI_2[  "amei"]  =       WMNEF["mnef"]*T_2(2)[  "afni"];
    WAMEI_2[  "amei"] +=       WAMEF["amef"]*T_2(1)[    "fi"];
    WAMEI_2[  "amei"] -=       WMNEJ["nmei"]*T_2(1)[    "an"];

    WMNIJ_2[  "mnij"]  =   0.5*WMNEF["mnef"]*T_2(2)[  "efij"];
    WMNIJ_2[  "mnij"] +=       WMNEJ["mnej"]*T_2(1)[    "ei"];

    WABEF_2[  "abef"]  =   0.5*WMNEF["mnef"]*T_2(2)[  "abmn"];
    WABEF_2[  "abef"] -=       WAMEF["amef"]*T_2(1)[    "bm"];

      FME_3[    "me"]  =       WMNEF["mnef"]*T_3(1)[    "fn"];

      FMI_3[    "mi"]  =   0.5*WMNEF["mnef"]*T_3(2)[  "efin"];
      FMI_3[    "mi"] +=       WMNEJ["nmfi"]*T_3(1)[    "fn"];

      FAE_3[    "ae"]  =  -0.5*WMNEF["mnef"]*T_3(2)[  "afmn"];
      FAE_3[    "ae"] +=       WAMEF["amef"]*T_3(1)[    "fm"];

    WAMIJ_3[  "amij"]  =       WMNEJ["nmej"]*T_3(2)[  "aein"];
    WAMIJ_3[  "amij"] +=   0.5*WAMEF["amef"]*T_3(2)[  "efij"];
    WAMIJ_3[  "amij"] -=       WMNIJ["nmij"]*T_3(1)[    "an"];
    WAMIJ_3[  "amij"] +=       WAMEI["amej"]*T_3(1)[    "ei"];
    WAMIJ_3[  "amij"] +=   0.5*WMNEF["mnef"]*T_3(3)["aefijn"];
    WAMIJ_3[  "amij"] +=       FME_3[  "me"]*  T(2)[  "aeij"];

    WABEJ_3[  "abej"]  =       WAMEF["amef"]*T_3(2)[  "fbmj"];
    WABEJ_3[  "abej"] +=   0.5*WMNEJ["mnej"]*T_3(2)[  "abmn"];
    WABEJ_3[  "abej"] +=       WABEF["abef"]*T_3(1)[    "fj"];
    WABEJ_3[  "abej"] -=       WAMEI["amej"]*T_3(1)[    "bm"];
    WABEJ_3[  "abej"] -=   0.5*WMNEF["mnef"]*T_3(3)["afbmnj"];

    DT_4(1)[    "ai"]  =         FAE[  "ae"]*T_3(1)[    "ei"];
    DT_4(1)[    "ai"] -=         FMI[  "mi"]*T_3(1)[    "am"];
    DT_4(1)[    "ai"] -=       WAMEI["amei"]*T_3(1)[    "em"];
    DT_4(1)[    "ai"] +=         FME[  "me"]*T_3(2)[  "aeim"];
    DT_4(1)[    "ai"] +=   0.5*WAMEF["amef"]*T_3(2)[  "efim"];
    DT_4(1)[    "ai"] -=   0.5*WMNEJ["mnei"]*T_3(2)[  "eamn"];
    DT_4(1)[    "ai"] +=  0.25*WMNEF["mnef"]*T_3(3)["aefimn"];

    DT_4(2)[  "abij"]  =       FAE_3[  "af"]*  T(2)[  "fbij"];
    DT_4(2)[  "abij"] -=       FMI_3[  "ni"]*  T(2)[  "abnj"];
    DT_4(2)[  "abij"] +=       FME_2[  "me"]*T_1(3)["abeijm"];
    DT_4(2)[  "abij"] +=     WABEJ_1["abej"]*T_2(1)[    "ei"];
    DT_4(2)[  "abij"] -=     WAMIJ_1["amij"]*T_2(1)[    "bm"];
    DT_4(2)[  "abij"] +=       WABEJ["abej"]*T_3(1)[    "ei"];
    DT_4(2)[  "abij"] -=       WAMIJ["amij"]*T_3(1)[    "bm"];
    DT_4(2)[  "abij"] +=         FAE[  "af"]*T_3(2)[  "fbij"];
    DT_4(2)[  "abij"] -=         FMI[  "ni"]*T_3(2)[  "abnj"];
    DT_4(2)[  "abij"] +=   0.5*WABEF["abef"]*T_3(2)[  "efij"];
    DT_4(2)[  "abij"] +=   0.5*WMNIJ["mnij"]*T_3(2)[  "abmn"];
    DT_4(2)[  "abij"] +=       WAMEI["amei"]*T_3(2)[  "ebjm"];
    DT_4(2)[  "abij"] +=   0.5*WAMEF["bmef"]*T_3(3)["aefijm"];
    DT_4(2)[  "abij"] -=   0.5*WMNEJ["mnej"]*T_3(3)["abeinm"];
    DT_4(2)[  "abij"] +=         FME[  "me"]*T_3(3)["abeijm"];

    DT_4(3)["abcijk"]  =     WABEJ_3["bcek"]*  T(2)[  "aeij"];
    DT_4(3)["abcijk"] -=     WAMIJ_3["bmjk"]*  T(2)[  "acim"];
    DT_4(3)["abcijk"] +=       FAE_2[  "ce"]*T_1(3)["abeijk"];
    DT_4(3)["abcijk"] -=       FMI_2[  "mk"]*T_1(3)["abcijm"];
    DT_4(3)["abcijk"] += 0.5*WABEF_2["abef"]*T_1(3)["efcijk"];
    DT_4(3)["abcijk"] += 0.5*WMNIJ_2["mnij"]*T_1(3)["abcmnk"];
    DT_4(3)["abcijk"] +=     WAMEI_2["amei"]*T_1(3)["ebcjmk"];
    DT_4(3)["abcijk"] +=     WABEJ_1["bcek"]*T_2(2)[  "aeij"];
    DT_4(3)["abcijk"] -=     WAMIJ_1["bmjk"]*T_2(2)[  "acim"];
    DT_4(3)["abcijk"] +=       WABEJ["bcek"]*T_3(2)[  "aeij"];
    DT_4(3)["abcijk"] -=       WAMIJ["bmjk"]*T_3(2)[  "acim"];
    DT_4(3)["abcijk"] +=         FAE[  "ce"]*T_3(3)["abeijk"];
    DT_4(3)["abcijk"] -=         FMI[  "mk"]*T_3(3)["abcijm"];
    DT_4(3)["abcijk"] +=   0.5*WABEF["abef"]*T_3(3)["efcijk"];
    DT_4(3)["abcijk"] +=   0.5*WMNIJ["mnij"]*T_3(3)["abcmnk"];
    DT_4(3)["abcijk"] +=       WAMEI["amei"]*T_3(3)["ebcjmk"];

    T_4 = DT_4;
    T_4.weight(D);

    /***************************************************************************
     *
     * T^(5)
     *
     **************************************************************************/

          Tau[  "abij"]  =  0.5*T_2(1)[  "ai"]*T_2(1)[    "bj"];

      WMNEJ_2[  "mnej"]  =       WMNEF["mnef"]*T_2(1)[    "fj"];
      WAMEF_2[  "amef"]  =      -WMNEF["nmef"]*T_2(1)[    "an"];
      WABEJ_2[  "abej"] -=       FME_2[  "me"]*  T(2)[  "abmj"];

    WTWAMIJ_2 = WAMIJ_2;
    WTWAMIJ_2[  "amij"] +=         FME[  "me"]*T_2(2)[  "aeij"];

    FTWMI_3 = FMI_3;
    FTWAE_3 = FAE_3;

        FMI_3[    "mi"] +=         FME[  "me"]*T_3(1)[    "ei"];
        FAE_3[    "ae"] -=         FME[  "me"]*T_3(1)[    "am"];

      WAMEI_3[  "amei"]  =       WMNEF["mnef"]*T_3(2)[  "afni"];
      WAMEI_3[  "amei"] +=       WAMEF["amef"]*T_3(1)[    "fi"];
      WAMEI_3[  "amei"] -=       WMNEJ["nmei"]*T_3(1)[    "an"];

      WMNIJ_3[  "mnij"]  =   0.5*WMNEF["mnef"]*T_3(2)[  "efij"];
      WMNIJ_3[  "mnij"] +=       WMNEJ["mnej"]*T_3(1)[    "ei"];

      WABEF_3[  "abef"]  =   0.5*WMNEF["mnef"]*T_3(2)[  "abmn"];
      WABEF_3[  "abef"] -=       WAMEF["amef"]*T_3(1)[    "bm"];

        FME_4[    "me"]  =       WMNEF["mnef"]*T_4(1)[    "fn"];

      WAMIJ_4[  "amij"]  =       WMNEJ["nmej"]*T_4(2)[  "aein"];
      WAMIJ_4[  "amij"] +=   0.5*WAMEF["amef"]*T_4(2)[  "efij"];
      WAMIJ_4[  "amij"] -=       WMNIJ["nmij"]*T_4(1)[    "an"];
      WAMIJ_4[  "amij"] +=       WAMEI["amej"]*T_4(1)[    "ei"];
      WAMIJ_4[  "amij"] +=   0.5*WMNEF["mnef"]*T_4(3)["aefijn"];
      WAMIJ_4[  "amij"] +=       FME_4[  "me"]*  T(2)[  "aeij"];
      WAMIJ_4[  "amij"] -=     WMNIJ_2["nmij"]*T_2(1)[    "an"];
      WAMIJ_4[  "amij"] +=     WMNEJ_2["nmej"]*T_2(2)[  "aein"];
      WAMIJ_4[  "amij"] +=   0.5*WAMEF["amef"]*   Tau[  "efij"];

      WABEJ_4[  "abej"]  =       WAMEF["amef"]*T_4(2)[  "fbmj"];
      WABEJ_4[  "abej"] +=   0.5*WMNEJ["mnej"]*T_4(2)[  "abmn"];
      WABEJ_4[  "abej"] +=       WABEF["abef"]*T_4(1)[    "fj"];
      WABEJ_4[  "abej"] -=       WAMEI["amej"]*T_4(1)[    "bm"];
      WABEJ_4[  "abej"] -=   0.5*WMNEF["mnef"]*T_4(3)["afbmnj"];
      WABEJ_4[  "abej"] +=     WABEF_2["abef"]*T_2(1)[    "fj"];
      WABEJ_4[  "abej"] +=     WAMEF_2["amef"]*T_2(2)[  "fbmj"];
      WABEJ_4[  "abej"] +=   0.5*WMNEJ["mnej"]*   Tau[  "abmn"];

      DT_5(3)["abcijk"]  =     WABEJ_4["bcek"]*  T(2)[  "aeij"];
      DT_5(3)["abcijk"] -=     WAMIJ_4["bmjk"]*  T(2)[  "acim"];
      DT_5(3)["abcijk"] +=       FAE_3[  "ce"]*T_1(3)["abeijk"];
      DT_5(3)["abcijk"] -=       FMI_3[  "mk"]*T_1(3)["abcijm"];
      DT_5(3)["abcijk"] += 0.5*WABEF_3["abef"]*T_1(3)["efcijk"];
      DT_5(3)["abcijk"] += 0.5*WMNIJ_3["mnij"]*T_1(3)["abcmnk"];
      DT_5(3)["abcijk"] +=     WAMEI_3["amei"]*T_1(3)["ebcjmk"];
      DT_5(3)["abcijk"] +=     WABEJ_2["bcek"]*T_2(2)[  "aeij"];
      DT_5(3)["abcijk"] -=   WTWAMIJ_2["bmjk"]*T_2(2)[  "acim"];
      DT_5(3)["abcijk"] +=       FAE_2[  "ce"]*T_2(3)["abeijk"];
      DT_5(3)["abcijk"] -=       FMI_2[  "mk"]*T_2(3)["abcijm"];
      DT_5(3)["abcijk"] += 0.5*WABEF_2["abef"]*T_2(3)["efcijk"];
      DT_5(3)["abcijk"] += 0.5*WMNIJ_2["mnij"]*T_2(3)["abcmnk"];
      DT_5(3)["abcijk"] +=     WAMEI_2["amei"]*T_2(3)["ebcjmk"];
      DT_5(3)["abcijk"] +=     WABEJ_1["bcek"]*T_3(2)[  "aeij"];
      DT_5(3)["abcijk"] -=     WAMIJ_1["bmjk"]*T_3(2)[  "acim"];
      DT_5(3)["abcijk"] +=       WABEJ["bcek"]*T_4(2)[  "aeij"];
      DT_5(3)["abcijk"] -=       WAMIJ["bmjk"]*T_4(2)[  "acim"];
      DT_5(3)["abcijk"] +=         FAE[  "ce"]*T_4(3)["abeijk"];
      DT_5(3)["abcijk"] -=         FMI[  "mk"]*T_4(3)["abcijm"];
      DT_5(3)["abcijk"] +=   0.5*WABEF["abef"]*T_4(3)["efcijk"];
      DT_5(3)["abcijk"] +=   0.5*WMNIJ["mnij"]*T_4(3)["abcmnk"];
      DT_5(3)["abcijk"] +=       WAMEI["amei"]*T_4(3)["ebcjmk"];

    T_5 = DT_5;
    T_5.weight(D);

    /***************************************************************************
     *
     * Complete W^(3) and W^(4)
     *
     **************************************************************************/

    WMNEJ_3["mnej"]  =  WMNEF["mnef"]*T_3(1)[  "fj"];
    WAMEF_3["amef"]  = -WMNEF["nmef"]*T_3(1)[  "an"];
    WABEJ_3["abej"] -=  FME_3[  "me"]*  T(2)["abmj"];

    WMNEJ_4["mnej"]  =  WMNEF["mnef"]*T_4(1)[  "fj"];
    WAMEF_4["amef"]  = -WMNEF["nmef"]*T_4(1)[  "an"];
    WABEJ_4["abej"] -=  FME_4[  "me"]*  T(2)["abmj"];

    /***************************************************************************
     *
     * L^(1)
     *
     **************************************************************************/

    DL_1(3)["ijkabc"]  = WMNEF["ijab"]*L(1)[  "kc"];
    DL_1(3)["ijkabc"] +=   FME[  "ia"]*L(2)["jkbc"];
    DL_1(3)["ijkabc"] += WAMEF["ekbc"]*L(2)["ijae"];
    DL_1(3)["ijkabc"] -= WMNEJ["ijam"]*L(2)["mkbc"];

    L_1 = DL_1;
    L_1.weight(D);

    /***************************************************************************
     *
     * L^(2)
     *
     **************************************************************************/

    GIJAK_1[  "ijak"]  =  (1.0/ 2.0)*   T(2)[  "efkm"]* L_1(3)["ijmaef"];
    GAIBC_1[  "aibc"]  = -(1.0/ 2.0)*   T(2)[  "aemn"]* L_1(3)["minbce"];

      DAI_1[    "ai"]  =  (1.0/ 4.0)* T_1(3)["efamni"]*   L(2)[  "mnef"];
      DAI_1[    "ai"] -=  (1.0/ 2.0)*   T(2)[  "eamn"]*GIJAK_1[  "mnei"];

    DL_2(1)[    "ia"]  =  (1.0/ 2.0)*WABEJ_1[  "efam"]*   L(2)[  "imef"];
    DL_2(1)[    "ia"] -=  (1.0/ 2.0)*WAMIJ_1[  "eimn"]*   L(2)[  "mnea"];
    DL_2(1)[    "ia"] +=               WMNEF[  "miea"]*  DAI_1[    "em"];
    DL_2(1)[    "ia"] -=  (1.0/ 2.0)*  WABEF[  "efga"]*GAIBC_1[  "gief"];
    DL_2(1)[    "ia"] +=               WAMEI[  "eifm"]*GAIBC_1[  "fmea"];
    DL_2(1)[    "ia"] -=               WAMEI[  "eman"]*GIJAK_1[  "inem"];
    DL_2(1)[    "ia"] +=  (1.0/ 2.0)*  WMNIJ[  "imno"]*GIJAK_1[  "noam"];

    DL_2(2)[  "ijab"]  = -             WAMEF[  "fiae"]*GAIBC_1[  "ejbf"];
    DL_2(2)[  "ijab"] -=               WMNEJ[  "ijem"]*GAIBC_1[  "emab"];
    DL_2(2)[  "ijab"] -=               WAMEF[  "emab"]*GIJAK_1[  "ijem"];
    DL_2(2)[  "ijab"] -=               WMNEJ[  "niam"]*GIJAK_1[  "mjbn"];
    DL_2(2)[  "ijab"] +=  (1.0/ 2.0)*  WABEJ[  "efbm"]* L_1(3)["ijmaef"];
    DL_2(2)[  "ijab"] -=  (1.0/ 2.0)*  WAMIJ[  "ejnm"]* L_1(3)["imnabe"];

    DL_2(3)["ijkabc"]  =               WMNEF[  "ijae"]*GAIBC_1[  "ekbc"];
    DL_2(3)["ijkabc"] -=               WMNEF[  "mkbc"]*GIJAK_1[  "ijam"];
    DL_2(3)["ijkabc"] +=                 FAE[    "ea"]* L_1(3)["ijkebc"];
    DL_2(3)["ijkabc"] -=                 FMI[    "im"]* L_1(3)["mjkabc"];
    DL_2(3)["ijkabc"] +=  (1.0/ 2.0)*  WABEF[  "efab"]* L_1(3)["ijkefc"];
    DL_2(3)["ijkabc"] +=  (1.0/ 2.0)*  WMNIJ[  "ijmn"]* L_1(3)["mnkabc"];
    DL_2(3)["ijkabc"] +=               WAMEI[  "eiam"]* L_1(3)["mjkbec"];

    L_2 = DL_2;
    L_2.weight(D);

    /***************************************************************************
     *
     * L^(3)
     *
     **************************************************************************/

        DIJ[    "ij"]  =  (1.0/ 2.0)*   T(2)[  "efjm"]*   L(2)[  "imef"];
        DAB[    "ab"]  = -(1.0/ 2.0)*   T(2)[  "aemn"]*   L(2)[  "mnbe"];

    DTWIJ_2[    "ij"]  =  (1.0/12.0)* T_1(3)["efgjmn"]* L_1(3)["imnefg"];
    DTWIJ_2[    "ij"] +=  (1.0/ 2.0)* T_2(2)[  "efjm"]*   L(2)[  "imef"];

    DTWAB_2[    "ab"]  = -(1.0/12.0)* T_1(3)["aefmno"]* L_1(3)["mnobef"];
    DTWAB_2[    "ab"] -=  (1.0/ 2.0)* T_2(2)[  "aemn"]*   L(2)[  "mnbe"];

      DIJ_2[    "ij"]  =             DTWIJ_2[    "ij"];
      DIJ_2[    "ij"] +=  (1.0/ 2.0)*   T(2)[  "efjm"]* L_2(2)[  "imef"];

      DAB_2[    "ab"]  =             DTWAB_2[    "ab"];
      DAB_2[    "ab"] -=  (1.0/ 2.0)*   T(2)[  "aemn"]* L_2(2)[  "mnbe"];

    GABCD_2[  "abcd"]  =  (1.0/ 6.0)* T_1(3)["abemno"]* L_1(3)["mnocde"];
    GAIBJ_2[  "aibj"]  = -(1.0/ 4.0)* T_1(3)["aefjmn"]* L_1(3)["imnbef"];
    GIJKL_2[  "ijkl"]  =  (1.0/ 6.0)* T_1(3)["efgklm"]* L_1(3)["ijmefg"];

    GIJAK_2[  "ijak"]  =  (1.0/ 2.0)*   T(2)[  "efkm"]* L_2(3)["ijmaef"];
    GAIBC_2[  "aibc"]  = -(1.0/ 2.0)*   T(2)[  "aemn"]* L_2(3)["minbce"];

      DAI_2[    "ai"]  =  (1.0/ 4.0)* T_2(3)["efamni"]*   L(2)[  "mnef"];
      DAI_2[    "ai"] +=              T_2(1)[    "ei"]*    DAB[    "ae"];
      DAI_2[    "ai"] -=              T_2(1)[    "am"]*    DIJ[    "mi"];
      DAI_2[    "ai"] -=  (1.0/ 2.0)*   T(2)[  "eamn"]*GIJAK_2[  "mnei"];

    DL_3(1)[    "ia"]  =               FME_2[    "ia"];
    DL_3(1)[    "ia"] +=               FAE_2[    "ea"]*   L(1)[    "ie"];
    DL_3(1)[    "ia"] -=               FMI_2[    "im"]*   L(1)[    "ma"];
    DL_3(1)[    "ia"] -=             WAMEI_2[  "eiam"]*   L(1)[    "me"];
    DL_3(1)[    "ia"] +=  (1.0/ 2.0)*WABEJ_2[  "efam"]*   L(2)[  "imef"];
    DL_3(1)[    "ia"] -=  (1.0/ 2.0)*WAMIJ_2[  "eimn"]*   L(2)[  "mnea"];
    DL_3(1)[    "ia"] +=                 FME[    "ie"]*DTWAB_2[    "ea"];
    DL_3(1)[    "ia"] -=                 FME[    "ma"]*DTWIJ_2[    "im"];
    DL_3(1)[    "ia"] -=               WMNEJ[  "inam"]*  DIJ_2[    "mn"];
    DL_3(1)[    "ia"] -=               WAMEF[  "fiea"]*  DAB_2[    "ef"];
    DL_3(1)[    "ia"] +=               WMNEF[  "miea"]*  DAI_2[    "em"];
    DL_3(1)[    "ia"] -=  (1.0/ 2.0)*  WABEF[  "efga"]*GAIBC_2[  "gief"];
    DL_3(1)[    "ia"] +=               WAMEI[  "eifm"]*GAIBC_2[  "fmea"];
    DL_3(1)[    "ia"] -=               WAMEI[  "eman"]*GIJAK_2[  "inem"];
    DL_3(1)[    "ia"] +=  (1.0/ 2.0)*  WMNIJ[  "imno"]*GIJAK_2[  "noam"];
    DL_3(1)[    "ia"] -=  (1.0/ 2.0)*  WAMEF[  "gief"]*GABCD_2[  "efga"];
    DL_3(1)[    "ia"] +=               WAMEF[  "fmea"]*GAIBJ_2[  "eifm"];
    DL_3(1)[    "ia"] -=               WMNEJ[  "inem"]*GAIBJ_2[  "eman"];
    DL_3(1)[    "ia"] +=  (1.0/ 2.0)*  WMNEJ[  "noam"]*GIJKL_2[  "imno"];
    DL_3(1)[    "ia"] +=                 FAE[    "ea"]* L_2(1)[    "ie"];
    DL_3(1)[    "ia"] -=                 FMI[    "im"]* L_2(1)[    "ma"];
    DL_3(1)[    "ia"] -=               WAMEI[  "eiam"]* L_2(1)[    "me"];
    DL_3(1)[    "ia"] +=  (1.0/ 2.0)*  WABEJ[  "efam"]* L_2(2)[  "imef"];
    DL_3(1)[    "ia"] -=  (1.0/ 2.0)*  WAMIJ[  "einm"]* L_2(2)[  "mnae"];

    DL_3(2)[  "ijab"]  =               FME_2[    "ia"]*   L(1)[    "jb"];
    DL_3(2)[  "ijab"] +=             WAMEF_2[  "ejab"]*   L(1)[    "ie"];
    DL_3(2)[  "ijab"] -=             WMNEJ_2[  "ijam"]*   L(1)[    "mb"];
    DL_3(2)[  "ijab"] +=               FAE_2[    "ea"]*   L(2)[  "ijeb"];
    DL_3(2)[  "ijab"] -=               FMI_2[    "im"]*   L(2)[  "mjab"];
    DL_3(2)[  "ijab"] +=  (1.0/ 2.0)*WABEF_2[  "efab"]*   L(2)[  "ijef"];
    DL_3(2)[  "ijab"] +=  (1.0/ 2.0)*WMNIJ_2[  "ijmn"]*   L(2)[  "mnab"];
    DL_3(2)[  "ijab"] +=             WAMEI_2[  "eiam"]*   L(2)[  "mjbe"];
    DL_3(2)[  "ijab"] +=  (1.0/ 2.0)*WABEJ_1[  "efbm"]* L_1(3)["ijmaef"];
    DL_3(2)[  "ijab"] -=  (1.0/ 2.0)*WAMIJ_1[  "ejnm"]* L_1(3)["imnabe"];
    DL_3(2)[  "ijab"] +=                 FME[    "ia"]* L_2(1)[    "jb"];
    DL_3(2)[  "ijab"] +=               WAMEF[  "ejab"]* L_2(1)[    "ie"];
    DL_3(2)[  "ijab"] -=               WMNEJ[  "ijam"]* L_2(1)[    "mb"];
    DL_3(2)[  "ijab"] -=               WMNEF[  "mjab"]*  DIJ_2[    "im"];
    DL_3(2)[  "ijab"] +=               WMNEF[  "ijeb"]*  DAB_2[    "ea"];
    DL_3(2)[  "ijab"] +=  (1.0/ 2.0)*  WMNEF[  "ijef"]*GABCD_2[  "efab"];
    DL_3(2)[  "ijab"] +=               WMNEF[  "imea"]*GAIBJ_2[  "ejbm"];
    DL_3(2)[  "ijab"] +=  (1.0/ 2.0)*  WMNEF[  "mnab"]*GIJKL_2[  "ijmn"];
    DL_3(2)[  "ijab"] -=               WAMEF[  "fiae"]*GAIBC_2[  "ejbf"];
    DL_3(2)[  "ijab"] -=               WMNEJ[  "ijem"]*GAIBC_2[  "emab"];
    DL_3(2)[  "ijab"] -=               WAMEF[  "emab"]*GIJAK_2[  "ijem"];
    DL_3(2)[  "ijab"] -=               WMNEJ[  "niam"]*GIJAK_2[  "mjbn"];
    DL_3(2)[  "ijab"] +=                 FAE[    "ea"]* L_2(2)[  "ijeb"];
    DL_3(2)[  "ijab"] -=                 FMI[    "im"]* L_2(2)[  "mjab"];
    DL_3(2)[  "ijab"] +=  (1.0/ 2.0)*  WABEF[  "efab"]* L_2(2)[  "ijef"];
    DL_3(2)[  "ijab"] +=  (1.0/ 2.0)*  WMNIJ[  "ijmn"]* L_2(2)[  "mnab"];
    DL_3(2)[  "ijab"] +=               WAMEI[  "eiam"]* L_2(2)[  "mjbe"];
    DL_3(2)[  "ijab"] +=  (1.0/ 2.0)*  WABEJ[  "efbm"]* L_2(3)["ijmaef"];
    DL_3(2)[  "ijab"] -=  (1.0/ 2.0)*  WAMIJ[  "ejnm"]* L_2(3)["imnabe"];

    DL_3(3)["ijkabc"]  =               FME_2[    "ia"]*   L(2)[  "jkbc"];
    DL_3(3)["ijkabc"] +=             WAMEF_2[  "ekbc"]*   L(2)[  "ijae"];
    DL_3(3)["ijkabc"] -=             WMNEJ_2[  "ijam"]*   L(2)[  "mkbc"];
    DL_3(3)["ijkabc"] +=               WMNEF[  "ijab"]* L_2(1)[    "kc"];
    DL_3(3)["ijkabc"] +=                 FME[    "ia"]* L_2(2)[  "jkbc"];
    DL_3(3)["ijkabc"] +=               WAMEF[  "ekbc"]* L_2(2)[  "ijae"];
    DL_3(3)["ijkabc"] -=               WMNEJ[  "ijam"]* L_2(2)[  "mkbc"];
    DL_3(3)["ijkabc"] +=               WMNEF[  "ijae"]*GAIBC_2[  "ekbc"];
    DL_3(3)["ijkabc"] -=               WMNEF[  "mkbc"]*GIJAK_2[  "ijam"];
    DL_3(3)["ijkabc"] +=                 FAE[    "ea"]* L_2(3)["ijkebc"];
    DL_3(3)["ijkabc"] -=                 FMI[    "im"]* L_2(3)["mjkabc"];
    DL_3(3)["ijkabc"] +=  (1.0/ 2.0)*  WABEF[  "efab"]* L_2(3)["ijkefc"];
    DL_3(3)["ijkabc"] +=  (1.0/ 2.0)*  WMNIJ[  "ijmn"]* L_2(3)["mnkabc"];
    DL_3(3)["ijkabc"] +=               WAMEI[  "eiam"]* L_2(3)["mjkbec"];

    L_3 = DL_3;
    L_3.weight(D);

    /***************************************************************************
     *
     * L^(4)
     *
     **************************************************************************/

    GIJAK_3[  "ijak"]  =  (1.0/ 2.0)*     T(2)[  "efkm"]* L_3(3)["ijmaef"];
    GIJAK_3[  "ijak"] +=  (1.0/ 2.0)*   T_2(2)[  "efkm"]* L_1(3)["ijmaef"];

    GAIBC_3[  "aibc"]  = -(1.0/ 2.0)*   T_2(2)[  "aemn"]* L_1(3)["minbce"];

    DTWIJ_3[    "ij"]  =  (1.0/12.0)*   T_1(3)["efgjmn"]* L_2(3)["imnefg"];
    DTWIJ_3[    "ij"] +=  (1.0/12.0)*   T_2(3)["efgjmn"]* L_1(3)["imnefg"];
    DTWIJ_3[    "ij"] +=  (1.0/ 2.0)*   T_3(2)[  "efjm"]*   L(2)[  "imef"];

    DTWAB_3[    "ab"]  = -(1.0/12.0)*   T_1(3)["aefmno"]* L_2(3)["mnobef"];
    DTWAB_3[    "ab"] -=  (1.0/12.0)*   T_2(3)["aefmno"]* L_1(3)["mnobef"];
    DTWAB_3[    "ab"] -=  (1.0/ 2.0)*   T_3(2)[  "aemn"]*   L(2)[  "mnbe"];

      DIJ_3 = DTWIJ_3;
      DIJ_3[    "ij"] +=  (1.0/ 2.0)*     T(2)[  "efjm"]* L_3(2)[  "imef"];

      DAB_3 = DTWAB_3;
      DAB_3[    "ab"] -=  (1.0/ 2.0)*     T(2)[  "aemn"]* L_3(2)[  "mnbe"];

      DAI_3[    "ai"]  =  (1.0/ 4.0)*   T_3(3)["efamni"]*   L(2)[  "mnef"];
      DAI_3[    "ai"] +=  (1.0/ 4.0)*   T_1(3)["efamni"]* L_2(2)[  "mnef"];
      DAI_3[    "ai"] +=                T_3(1)[    "ei"]*    DAB[    "ae"];
      DAI_3[    "ai"] -=                T_3(1)[    "am"]*    DIJ[    "mi"];
      DAI_3[    "ai"] -=  (1.0/ 2.0)*     T(2)[  "eamn"]*GIJAK_3[  "mnei"];
      DAI_3[    "ai"] +=  (1.0/ 2.0)*     T(2)[  "efim"]*GAIBC_3[  "amef"];

    GAIBC_3[  "aibc"] -=  (1.0/ 2.0)*     T(2)[  "aemn"]* L_3(3)["minbce"];

    GABCD_3[  "abcd"]  =  (1.0/ 6.0)*   T_2(3)["abemno"]* L_1(3)["mnocde"];
    GABCD_3[  "abcd"] +=  (1.0/ 6.0)*   T_1(3)["abemno"]* L_2(3)["mnocde"];

    GAIBJ_3[  "aibj"]  = -(1.0/ 4.0)*   T_2(3)["aefjmn"]* L_1(3)["imnbef"];
    GAIBJ_3[  "aibj"] -=  (1.0/ 4.0)*   T_1(3)["aefjmn"]* L_2(3)["imnbef"];

    GIJKL_3[  "ijkl"]  =  (1.0/ 6.0)*   T_2(3)["efgklm"]* L_1(3)["ijmefg"];
    GIJKL_3[  "ijkl"] +=  (1.0/ 6.0)*   T_1(3)["efgklm"]* L_2(3)["ijmefg"];

    DL_4(1)[    "ia"]  =                 FME_3[    "ia"];
    DL_4(1)[    "ia"] +=                 FAE_3[    "ea"]*   L(1)[    "ie"];
    DL_4(1)[    "ia"] -=                 FMI_3[    "im"]*   L(1)[    "ma"];
    DL_4(1)[    "ia"] -=               WAMEI_3[  "eiam"]*   L(1)[    "me"];
    DL_4(1)[    "ia"] +=  (1.0/ 2.0)*  WABEJ_3[  "efam"]*   L(2)[  "imef"];
    DL_4(1)[    "ia"] -=  (1.0/ 2.0)*  WAMIJ_3[  "eimn"]*   L(2)[  "mnea"];
    DL_4(1)[    "ia"] -=  (1.0/ 2.0)*  WABEF_2[  "efga"]*GAIBC_1[  "gief"];
    DL_4(1)[    "ia"] +=               WAMEI_2[  "eifm"]*GAIBC_1[  "fmea"];
    DL_4(1)[    "ia"] -=               WAMEI_2[  "eman"]*GIJAK_1[  "inem"];
    DL_4(1)[    "ia"] +=  (1.0/ 2.0)*  WMNIJ_2[  "imno"]*GIJAK_1[  "noam"];
    DL_4(1)[    "ia"] +=  (1.0/ 2.0)*  WABEJ_1[  "efam"]* L_2(2)[  "imef"];
    DL_4(1)[    "ia"] -=  (1.0/ 2.0)*  WAMIJ_1[  "eimn"]* L_2(2)[  "mnea"];
    DL_4(1)[    "ia"] +=                   FME[    "ie"]*DTWAB_3[    "ea"];
    DL_4(1)[    "ia"] -=                   FME[    "ma"]*DTWIJ_3[    "im"];
    DL_4(1)[    "ia"] -=                 WMNEJ[  "inam"]*  DIJ_3[    "mn"];
    DL_4(1)[    "ia"] -=                 WAMEF[  "fiea"]*  DAB_3[    "ef"];
    DL_4(1)[    "ia"] +=                 WMNEF[  "miea"]*  DAI_3[    "em"];
    DL_4(1)[    "ia"] -=  (1.0/ 2.0)*    WABEF[  "efga"]*GAIBC_3[  "gief"];
    DL_4(1)[    "ia"] +=                 WAMEI[  "eifm"]*GAIBC_3[  "fmea"];
    DL_4(1)[    "ia"] -=                 WAMEI[  "eman"]*GIJAK_3[  "inem"];
    DL_4(1)[    "ia"] +=  (1.0/ 2.0)*    WMNIJ[  "imno"]*GIJAK_3[  "noam"];
    DL_4(1)[    "ia"] -=  (1.0/ 2.0)*    WAMEF[  "gief"]*GABCD_3[  "efga"];
    DL_4(1)[    "ia"] +=                 WAMEF[  "fmea"]*GAIBJ_3[  "eifm"];
    DL_4(1)[    "ia"] -=                 WMNEJ[  "inem"]*GAIBJ_3[  "eman"];
    DL_4(1)[    "ia"] +=  (1.0/ 2.0)*    WMNEJ[  "noam"]*GIJKL_3[  "imno"];
    DL_4(1)[    "ia"] +=                   FAE[    "ea"]* L_3(1)[    "ie"];
    DL_4(1)[    "ia"] -=                   FMI[    "im"]* L_3(1)[    "ma"];
    DL_4(1)[    "ia"] -=                 WAMEI[  "eiam"]* L_3(1)[    "me"];
    DL_4(1)[    "ia"] +=  (1.0/ 2.0)*    WABEJ[  "efam"]* L_3(2)[  "imef"];
    DL_4(1)[    "ia"] -=  (1.0/ 2.0)*    WAMIJ[  "einm"]* L_3(2)[  "mnae"];

    DL_4(2)[  "ijab"]  =                 FME_3[    "ia"]*   L(1)[    "jb"];
    DL_4(2)[  "ijab"] +=               WAMEF_3[  "ejab"]*   L(1)[    "ie"];
    DL_4(2)[  "ijab"] -=               WMNEJ_3[  "ijam"]*   L(1)[    "mb"];
    DL_4(2)[  "ijab"] +=                 FAE_3[    "ea"]*   L(2)[  "ijeb"];
    DL_4(2)[  "ijab"] -=                 FMI_3[    "im"]*   L(2)[  "mjab"];
    DL_4(2)[  "ijab"] +=  (1.0/ 2.0)*  WABEF_3[  "efab"]*   L(2)[  "ijef"];
    DL_4(2)[  "ijab"] +=  (1.0/ 2.0)*  WMNIJ_3[  "ijmn"]*   L(2)[  "mnab"];
    DL_4(2)[  "ijab"] +=               WAMEI_3[  "eiam"]*   L(2)[  "mjbe"];
    DL_4(2)[  "ijab"] -=               WAMEF_2[  "fiae"]*GAIBC_1[  "ejbf"];
    DL_4(2)[  "ijab"] -=               WMNEJ_2[  "ijem"]*GAIBC_1[  "emab"];
    DL_4(2)[  "ijab"] -=               WAMEF_2[  "emab"]*GIJAK_1[  "ijem"];
    DL_4(2)[  "ijab"] -=               WMNEJ_2[  "niam"]*GIJAK_1[  "mjbn"];
    DL_4(2)[  "ijab"] +=  (1.0/ 2.0)*  WABEJ_2[  "efbm"]* L_1(3)["ijmaef"];
    DL_4(2)[  "ijab"] -=  (1.0/ 2.0)*WTWAMIJ_2[  "ejnm"]* L_1(3)["imnabe"];
    DL_4(2)[  "ijab"] +=  (1.0/ 2.0)*  WABEJ_1[  "efbm"]* L_2(3)["ijmaef"];
    DL_4(2)[  "ijab"] -=  (1.0/ 2.0)*  WAMIJ_1[  "ejnm"]* L_2(3)["imnabe"];
    DL_4(2)[  "ijab"] -=                 WMNEF[  "mjab"]*  DIJ_3[    "im"];
    DL_4(2)[  "ijab"] +=                 WMNEF[  "ijeb"]*  DAB_3[    "ea"];
    DL_4(2)[  "ijab"] +=  (1.0/ 2.0)*    WMNEF[  "ijef"]*GABCD_3[  "efab"];
    DL_4(2)[  "ijab"] +=                 WMNEF[  "imea"]*GAIBJ_3[  "ejbm"];
    DL_4(2)[  "ijab"] +=  (1.0/ 2.0)*    WMNEF[  "mnab"]*GIJKL_3[  "ijmn"];
    DL_4(2)[  "ijab"] -=                   FME[    "mb"]*GIJAK_3[  "ijam"];
    DL_4(2)[  "ijab"] -=                 WAMEF[  "fiae"]*GAIBC_3[  "ejbf"];
    DL_4(2)[  "ijab"] -=                 WMNEJ[  "ijem"]*GAIBC_3[  "emab"];
    DL_4(2)[  "ijab"] -=                 WAMEF[  "emab"]*GIJAK_3[  "ijem"];
    DL_4(2)[  "ijab"] -=                 WMNEJ[  "niam"]*GIJAK_3[  "mjbn"];
    DL_4(2)[  "ijab"] +=                   FME[    "ia"]* L_3(1)[    "jb"];
    DL_4(2)[  "ijab"] +=                 WAMEF[  "ejab"]* L_3(1)[    "ie"];
    DL_4(2)[  "ijab"] -=                 WMNEJ[  "ijam"]* L_3(1)[    "mb"];
    DL_4(2)[  "ijab"] +=                   FAE[    "ea"]* L_3(2)[  "ijeb"];
    DL_4(2)[  "ijab"] -=                   FMI[    "im"]* L_3(2)[  "mjab"];
    DL_4(2)[  "ijab"] +=  (1.0/ 2.0)*    WABEF[  "efab"]* L_3(2)[  "ijef"];
    DL_4(2)[  "ijab"] +=  (1.0/ 2.0)*    WMNIJ[  "ijmn"]* L_3(2)[  "mnab"];
    DL_4(2)[  "ijab"] +=                 WAMEI[  "eiam"]* L_3(2)[  "mjbe"];
    DL_4(2)[  "ijab"] +=  (1.0/ 2.0)*  WTWABEJ[  "efbm"]* L_3(3)["ijmaef"];
    DL_4(2)[  "ijab"] -=  (1.0/ 2.0)*    WAMIJ[  "ejnm"]* L_3(3)["imnabe"];

    DL_4(3)["ijkabc"]  =                 FME_3[    "ia"]*   L(2)[  "jkbc"];
    DL_4(3)["ijkabc"] +=               WAMEF_3[  "ekbc"]*   L(2)[  "ijae"];
    DL_4(3)["ijkabc"] -=               WMNEJ_3[  "ijam"]*   L(2)[  "mkbc"];
    DL_4(3)["ijkabc"] +=                 FAE_2[    "ea"]* L_1(3)["ijkebc"];
    DL_4(3)["ijkabc"] -=                 FMI_2[    "im"]* L_1(3)["mjkabc"];
    DL_4(3)["ijkabc"] +=  (1.0/ 2.0)*  WABEF_2[  "efab"]* L_1(3)["ijkefc"];
    DL_4(3)["ijkabc"] +=  (1.0/ 2.0)*  WMNIJ_2[  "ijmn"]* L_1(3)["mnkabc"];
    DL_4(3)["ijkabc"] +=               WAMEI_2[  "eiam"]* L_1(3)["mjkbec"];
    DL_4(3)["ijkabc"] +=                 WMNEF[  "ijab"]* L_3(1)[    "kc"];
    DL_4(3)["ijkabc"] +=                   FME[    "ia"]* L_3(2)[  "jkbc"];
    DL_4(3)["ijkabc"] +=                 WAMEF[  "ekbc"]* L_3(2)[  "ijae"];
    DL_4(3)["ijkabc"] -=                 WMNEJ[  "ijam"]* L_3(2)[  "mkbc"];
    DL_4(3)["ijkabc"] +=                 WMNEF[  "ijae"]*GAIBC_3[  "ekbc"];
    DL_4(3)["ijkabc"] -=                 WMNEF[  "mkbc"]*GIJAK_3[  "ijam"];
    DL_4(3)["ijkabc"] +=                   FAE[    "ea"]* L_3(3)["ijkebc"];
    DL_4(3)["ijkabc"] -=                   FMI[    "im"]* L_3(3)["mjkabc"];
    DL_4(3)["ijkabc"] +=  (1.0/ 2.0)*    WABEF[  "efab"]* L_3(3)["ijkefc"];
    DL_4(3)["ijkabc"] +=  (1.0/ 2.0)*    WMNIJ[  "ijmn"]* L_3(3)["mnkabc"];
    DL_4(3)["ijkabc"] +=                 WAMEI[  "eiam"]* L_3(3)["mjkbec"];

    L_4 = DL_4;
    L_4.weight(D);

    /***************************************************************************
     *
     * L^(4)
     *
     **************************************************************************/

    GIJAK_4[  "ijak"]  =  (1.0/ 2.0)*     T(2)[  "efkm"]* L_4(3)["ijmaef"];
    GIJAK_4[  "ijak"] +=  (1.0/ 2.0)*   T_2(2)[  "efkm"]* L_2(3)["ijmaef"];
    GIJAK_4[  "ijak"] +=  (1.0/ 2.0)*   T_3(2)[  "efkm"]* L_1(3)["ijmaef"];

    GAIBC_4[  "aibc"]  = -(1.0/ 2.0)*   T_3(2)[  "aemn"]* L_1(3)["minbce"];
    GAIBC_4[  "aibc"] -=  (1.0/ 2.0)*   T_2(2)[  "aemn"]* L_2(3)["minbce"];
    GAIBC_4[  "aibc"] -=  (1.0/ 2.0)*     T(2)[  "aemn"]* L_4(3)["minbce"];

    DL_5(3)["ijkabc"]  =                 FME_4[    "ia"]*   L(2)[  "jkbc"];
    DL_5(3)["ijkabc"] +=               WAMEF_4[  "ekbc"]*   L(2)[  "ijae"];
    DL_5(3)["ijkabc"] -=               WMNEJ_4[  "ijam"]*   L(2)[  "mkbc"];
    DL_5(3)["ijkabc"] +=                 FAE_3[    "ea"]* L_1(3)["ijkebc"];
    DL_5(3)["ijkabc"] -=                 FMI_3[    "im"]* L_1(3)["mjkabc"];
    DL_5(3)["ijkabc"] +=  (1.0/ 2.0)*  WABEF_3[  "efab"]* L_1(3)["ijkefc"];
    DL_5(3)["ijkabc"] +=  (1.0/ 2.0)*  WMNIJ_3[  "ijmn"]* L_1(3)["mnkabc"];
    DL_5(3)["ijkabc"] +=               WAMEI_3[  "eiam"]* L_1(3)["mjkbec"];
    DL_5(3)["ijkabc"] +=                 FME_2[    "ia"]* L_2(2)[  "jkbc"];
    DL_5(3)["ijkabc"] +=               WAMEF_2[  "ekbc"]* L_2(2)[  "ijae"];
    DL_5(3)["ijkabc"] -=               WMNEJ_2[  "ijam"]* L_2(2)[  "mkbc"];
    DL_5(3)["ijkabc"] +=                 FAE_2[    "ea"]* L_2(3)["ijkebc"];
    DL_5(3)["ijkabc"] -=                 FMI_2[    "im"]* L_2(3)["mjkabc"];
    DL_5(3)["ijkabc"] +=  (1.0/ 2.0)*  WABEF_2[  "efab"]* L_2(3)["ijkefc"];
    DL_5(3)["ijkabc"] +=  (1.0/ 2.0)*  WMNIJ_2[  "ijmn"]* L_2(3)["mnkabc"];
    DL_5(3)["ijkabc"] +=               WAMEI_2[  "eiam"]* L_2(3)["mjkbec"];
    DL_5(3)["ijkabc"] +=                 WMNEF[  "ijab"]* L_4(1)[    "kc"];
    DL_5(3)["ijkabc"] +=                   FME[    "ia"]* L_4(2)[  "jkbc"];
    DL_5(3)["ijkabc"] +=                 WAMEF[  "ekbc"]* L_4(2)[  "ijae"];
    DL_5(3)["ijkabc"] -=                 WMNEJ[  "ijam"]* L_4(2)[  "mkbc"];
    DL_5(3)["ijkabc"] +=                 WMNEF[  "ijae"]*GAIBC_4[  "ekbc"];
    DL_5(3)["ijkabc"] -=                 WMNEF[  "mkbc"]*GIJAK_4[  "ijam"];
    DL_5(3)["ijkabc"] +=                   FAE[    "ea"]* L_4(3)["ijkebc"];
    DL_5(3)["ijkabc"] -=                   FMI[    "im"]* L_4(3)["mjkabc"];
    DL_5(3)["ijkabc"] +=  (1.0/ 2.0)*    WABEF[  "efab"]* L_4(3)["ijkefc"];
    DL_5(3)["ijkabc"] +=  (1.0/ 2.0)*    WMNIJ[  "ijmn"]* L_4(3)["mnkabc"];
    DL_5(3)["ijkabc"] +=                 WAMEI[  "eiam"]* L_4(3)["mjkbec"];

    L_5 = DL_5;
    L_5.weight(D);

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
     * <0|(1+L)[[H,T^(1)],T^(3)]|0>
     *
     **************************************************************************/

    Z(2)["abij"]  =       FME_3[  "me"]*T_1(3)["eabmij"];
    Z(2)["abij"] += 0.5*WAMEF_3["bmef"]*T_1(3)["aefijm"];
    Z(2)["abij"] -= 0.5*WMNEJ_3["nmej"]*T_1(3)["abeimn"];

    U E0113 = 0.25*scalar(L(2)["mnef"]*Z(2)["efmn"]);

    /***************************************************************************
     *           _
     * <0|(1+L)[[H,T^(1)],T^(4)]|0>
     *
     **************************************************************************/

    Z(2)["abij"]  =       FME_4[  "me"]*T_1(3)["eabmij"];
    Z(2)["abij"] += 0.5*WAMEF_4["bmef"]*T_1(3)["aefijm"];
    Z(2)["abij"] -= 0.5*WMNEJ_4["nmej"]*T_1(3)["abeimn"];

    U E0114 = 0.25*scalar(L(2)["mnef"]*Z(2)["efmn"]);

    /***************************************************************************
     *               _
     * 1/2 <0|(1+L)[[H,T^(2)],T^(2)]|0>
     *
     **************************************************************************/

          Tau["abij"]  =  0.5*T_2(1)[  "ai"]*T_2(1)[    "bj"];

    WTWAMIJ_2["amij"]  =   0.5*WAMEF["amef"]*T_2(2)[  "efij"];
    WTWAMIJ_2["amij"] -=   0.5*WMNIJ["nmij"]*T_2(1)[    "an"];
    WTWAMIJ_2["amij"] +=       WAMEI["amej"]*T_2(1)[    "ei"];
    WTWAMIJ_2["amij"] +=       FME_2[  "me"]*  T(2)[  "aeij"];
    WTWAMIJ_2["amij"] +=   0.5*WMNEF["mnef"]*T_2(3)["aefijn"];

    WTWABEJ_2["abej"]  =   0.5*WABEF["abef"]*T_2(1)[    "fj"];
    WTWABEJ_2["abej"] -=       FME_2[  "me"]*  T(2)[  "abmj"];
    WTWABEJ_2["abej"] -=   0.5*WMNEF["mnef"]*T_2(3)["afbmnj"];

    WTWAMEI_2["amei"]  =   0.5*WMNEF["mnef"]*T_2(2)[  "afni"];
    WTWAMEI_2["amei"] +=       WAMEF["amef"]*T_2(1)[    "fi"];
    WTWAMEI_2["amei"] -=       WMNEJ["nmei"]*T_2(1)[    "an"];

         Z(1)[  "ai"]  =     FTWAE_2[  "ae"]*T_2(1)[    "ei"];
         Z(1)[  "ai"] -=       FMI_2[  "mi"]*T_2(1)[    "am"];
         Z(1)[  "ai"] +=       FME_2[  "me"]*T_2(2)[  "eami"];

         Z(2)["abij"]  =   WTWABEJ_2["abej"]*T_2(1)[    "ei"];
         Z(2)["abij"] -=   WTWAMIJ_2["amij"]*T_2(1)[    "bm"];
         Z(2)["abij"] +=       FAE_2[  "ae"]*T_2(2)[  "ebij"];
         Z(2)["abij"] -=       FMI_2[  "mi"]*T_2(2)[  "abmj"];
         Z(2)["abij"] += 0.5*WMNIJ_2["mnij"]*T_2(2)[  "abmn"];
         Z(2)["abij"] +=   WTWAMEI_2["amei"]*T_2(2)[  "ebjm"];
         Z(2)["abij"] +=       FME_2[  "me"]*T_2(3)["eabmij"];

    U E0122 =      scalar( L(1)[  "me"]*Z(1)[  "em"]) +
              0.25*scalar( L(2)["mnef"]*Z(2)["efmn"]) +
              0.25*scalar(WMNEF["mnef"]* Tau["efmn"]);

    /***************************************************************************
     *           _
     * <0|(1+L)[[H,T^(2)],T^(3)]|0>
     *
     **************************************************************************/

          Tau["abij"]  =      T_2(1)[  "ai"]*T_3(1)[    "bj"];

    WTWAMIJ_3["amij"]  =       FME_3[  "me"]*  T(2)[  "aeij"];
    WTWAMIJ_3["amij"] +=   0.5*WMNEF["mnef"]*T_3(3)["aefijn"];

    WTWABEJ_3["abej"]  =      -FME_3[  "me"]*  T(2)[  "abmj"];
    WTWABEJ_3["abej"] -=   0.5*WMNEF["mnef"]*T_3(3)["afbmnj"];

         Z(1)[  "ai"]  =     FTWAE_3[  "ae"]*T_2(1)[    "ei"];
         Z(1)[  "ai"] -=     FTWMI_3[  "mi"]*T_2(1)[    "am"];
         Z(1)[  "ai"] +=       FME_3[  "me"]*T_2(2)[  "eami"];
         Z(1)[  "ai"] +=       FAE_2[  "ae"]*T_3(1)[    "ei"];
         Z(1)[  "ai"] -=       FMI_2[  "mi"]*T_3(1)[    "am"];
         Z(1)[  "ai"] +=       FME_2[  "me"]*T_3(2)[  "eami"];

         Z(2)["abij"]  =   WTWABEJ_3["abej"]*T_2(1)[    "ei"];
         Z(2)["abij"] -=   WTWAMIJ_3["amij"]*T_2(1)[    "bm"];
         Z(2)["abij"] +=       FAE_3[  "ae"]*T_2(2)[  "ebij"];
         Z(2)["abij"] -=       FMI_3[  "mi"]*T_2(2)[  "abmj"];
         Z(2)["abij"] +=       FME_3[  "me"]*T_2(3)["eabmij"];
         Z(2)["abij"] +=     WABEJ_2["abej"]*T_3(1)[    "ei"];
         Z(2)["abij"] -=     WAMIJ_2["amij"]*T_3(1)[    "bm"];
         Z(2)["abij"] +=       FAE_2[  "ae"]*T_3(2)[  "ebij"];
         Z(2)["abij"] -=       FMI_2[  "mi"]*T_3(2)[  "abmj"];
         Z(2)["abij"] += 0.5*WABEF_2["abef"]*T_3(2)[  "efij"];
         Z(2)["abij"] += 0.5*WMNIJ_2["mnij"]*T_3(2)[  "abmn"];
         Z(2)["abij"] +=     WAMEI_2["amei"]*T_3(2)[  "ebjm"];
         Z(2)["abij"] +=       FME_2[  "me"]*T_3(3)["eabmij"];

    U E0123 =      scalar( L(1)[  "me"]*Z(1)[  "em"]) +
              0.25*scalar( L(2)["mnef"]*Z(2)["efmn"]) +
              0.25*scalar(WMNEF["mnef"]* Tau["efmn"]);

    /***************************************************************************
     *           _
     * <0|L^(1)[[H,T^(1)],T^(2)]|0>
     *
     **************************************************************************/

    Z(3)["abcijk"]  =     WABEJ_1["bcek"]*T_2(2)[  "aeij"];
    Z(3)["abcijk"] -=     WAMIJ_1["bmjk"]*T_2(2)[  "acim"];
    Z(3)["abcijk"] +=       FAE_2[  "ae"]*T_1(3)["ebcijk"];
    Z(3)["abcijk"] -=       FMI_2[  "mi"]*T_1(3)["abcmjk"];
    Z(3)["abcijk"] += 0.5*WABEF_2["abef"]*T_1(3)["efcijk"];
    Z(3)["abcijk"] += 0.5*WMNIJ_2["mnij"]*T_1(3)["abcmnk"];
    Z(3)["abcijk"] +=     WAMEI_2["amei"]*T_1(3)["becmjk"];

    U E1112 = (1.0/36.0)*scalar(L_1(3)["mnoefg"]*Z(3)["efgmno"]);

    /***************************************************************************
     *           _
     * <0|L^(1)[[H,T^(1)],T^(3)]|0>
     *
     **************************************************************************/

    Z(3)["abcijk"]  =     WABEJ_1["bcek"]*T_3(2)[  "aeij"];
    Z(3)["abcijk"] -=     WAMIJ_1["bmjk"]*T_3(2)[  "acim"];
    Z(3)["abcijk"] +=       FAE_3[  "ae"]*T_1(3)["ebcijk"];
    Z(3)["abcijk"] -=       FMI_3[  "mi"]*T_1(3)["abcmjk"];
    Z(3)["abcijk"] += 0.5*WABEF_3["abef"]*T_1(3)["efcijk"];
    Z(3)["abcijk"] += 0.5*WMNIJ_3["mnij"]*T_1(3)["abcmnk"];
    Z(3)["abcijk"] +=     WAMEI_3["amei"]*T_1(3)["becmjk"];

    U E1113 = (1.0/36.0)*scalar(L_1(3)["mnoefg"]*Z(3)["efgmno"]);

    /***************************************************************************
     *               _
     * 1/2 <0|L^(1)[[H,T^(2)],T^(2)]|0>
     *
     **************************************************************************/

          Tau[  "abij"]  =  0.5*T_2(1)[  "ai"]*T_2(1)[    "bj"];

    WTWAMIJ_2 = WAMIJ_2;
    WTWAMIJ_2[  "amij"] +=         FME[  "me"]*T_2(2)[  "aeij"];

    WTWAMIJ_4[  "amij"]  =     WMNEJ_2["nmej"]*T_2(2)[  "aein"];
    WTWAMIJ_4[  "amij"] +=   0.5*WAMEF["amef"]*   Tau[  "efij"];
    WTWAMIJ_4[  "amij"] -=     WMNIJ_2["nmij"]*T_2(1)[    "an"];

    WTWABEJ_4[  "abej"]  =     WAMEF_2["amef"]*T_2(2)[  "fbmj"];
    WTWABEJ_4[  "abej"] +=   0.5*WMNEJ["mnej"]*   Tau[  "abmn"];
    WTWABEJ_4[  "abej"] +=     WABEF_2["abef"]*T_2(1)[    "fj"];

         Z(3)["abcijk"]  =   WTWABEJ_4["bcek"]*  T(2)[  "aeij"];
         Z(3)["abcijk"] -=   WTWAMIJ_4["bmjk"]*  T(2)[  "acim"];
         Z(3)["abcijk"] +=     WABEJ_2["bcek"]*T_2(2)[  "aeij"];
         Z(3)["abcijk"] -=   WTWAMIJ_2["bmjk"]*T_2(2)[  "acim"];
         Z(3)["abcijk"] +=       FAE_2[  "ce"]*T_2(3)["abeijk"];
         Z(3)["abcijk"] -=       FMI_2[  "mk"]*T_2(3)["abcijm"];
         Z(3)["abcijk"] += 0.5*WABEF_2["abef"]*T_2(3)["efcijk"];
         Z(3)["abcijk"] += 0.5*WMNIJ_2["mnij"]*T_2(3)["abcmnk"];
         Z(3)["abcijk"] +=     WAMEI_2["amei"]*T_2(3)["ebcjmk"];

    U E1122 = (1.0/36.0)*scalar(L_1(3)["mnoefg"]*Z(3)["efgmno"]);

    /***************************************************************************
     *           _
     * <0|L^(2)[[H,T^(1)],T^(2)]|0>
     *
     **************************************************************************/

    Z(2)[  "abij"]  =       FME_2[  "me"]*T_1(3)["eabmij"];
    Z(2)[  "abij"] += 0.5*WAMEF_2["bmef"]*T_1(3)["aefijm"];
    Z(2)[  "abij"] -= 0.5*WMNEJ_2["nmej"]*T_1(3)["abeimn"];

    Z(3)["abcijk"]  =     WABEJ_1["bcek"]*T_2(2)[  "aeij"];
    Z(3)["abcijk"] -=     WAMIJ_1["bmjk"]*T_2(2)[  "acim"];
    Z(3)["abcijk"] +=       FAE_2[  "ae"]*T_1(3)["ebcijk"];
    Z(3)["abcijk"] -=       FMI_2[  "mi"]*T_1(3)["abcmjk"];
    Z(3)["abcijk"] += 0.5*WABEF_2["abef"]*T_1(3)["efcijk"];
    Z(3)["abcijk"] += 0.5*WMNIJ_2["mnij"]*T_1(3)["abcmnk"];
    Z(3)["abcijk"] +=     WAMEI_2["amei"]*T_1(3)["becmjk"];

    U E2112 = (1.0/ 4.0)*scalar(L_2(2)[  "mnef"]*Z(2)[  "efmn"])+
              (1.0/36.0)*scalar(L_2(3)["mnoefg"]*Z(3)["efgmno"]);

    /***************************************************************************
     *
     * CCSD(T-2)
     *
     **************************************************************************/

    U E101 = (1.0/36.0)*scalar(L_1(3)["mnoefg"]*DT_1(3)["efgmno"]);
    U E2 = E101;

    /***************************************************************************
     *
     * CCSD(T-3)
     *
     **************************************************************************/

    U E201 = (1.0/36.0)*scalar(L_2(3)["mnoefg"]*DT_1(3)["efgmno"]);
    U E102 = (1.0/36.0)*scalar(L_1(3)["mnoefg"]*DT_2(3)["efgmno"]);

    printf("E201: %18.15f\n", E201);
    printf("E102: %18.15f\n", E102);
    printf("\n");

    U E3 = E102;

    /***************************************************************************
     *
     * CCSD(T-4)
     *
     **************************************************************************/

    U E301 =  (1.0/36.0)*scalar(L_3(3)["mnoefg"]*DT_1(3)["efgmno"]);
    U E202 =  (1.0/ 1.0)*scalar(L_2(1)[    "me"]*DT_2(1)[    "em"])
             +(1.0/ 4.0)*scalar(L_2(2)[  "mnef"]*DT_2(2)[  "efmn"])
             +(1.0/36.0)*scalar(L_2(3)["mnoefg"]*DT_2(3)["efgmno"]);
    U E103 =  (1.0/36.0)*scalar(L_1(3)["mnoefg"]*DT_3(3)["efgmno"]);

    printf("E301: %18.15f\n", E301);
    printf("E202: %18.15f\n", E202);
    printf("E103: %18.15f\n", E103);
    printf("\n");

    printf("E0112: %18.15f\n", E0112);
    printf("\n");

    printf("E301:       %18.15f\n", E301);
    printf("E202+E0112: %18.15f\n", E202+E0112);
    printf("\n");

    printf("E202:       %18.15f\n", E202);
    printf("E103+E0112: %18.15f\n", E103+E0112);
    printf("\n");

    U E4 = E202;

    /***************************************************************************
     *
     * CCSD(T-5)
     *
     **************************************************************************/

    U E401 =  (1.0/36.0)*scalar(L_4(3)["mnoefg"]*DT_1(3)["efgmno"]);
    U E302 =  (1.0/ 1.0)*scalar(L_3(1)[    "me"]*DT_2(1)[    "em"])
             +(1.0/ 4.0)*scalar(L_3(2)[  "mnef"]*DT_2(2)[  "efmn"])
             +(1.0/36.0)*scalar(L_3(3)["mnoefg"]*DT_2(3)["efgmno"]);
    U E203 =  (1.0/ 1.0)*scalar(L_2(1)[    "me"]*DT_3(1)[    "em"])
             +(1.0/ 4.0)*scalar(L_2(2)[  "mnef"]*DT_3(2)[  "efmn"])
             +(1.0/36.0)*scalar(L_2(3)["mnoefg"]*DT_3(3)["efgmno"]);
    U E104 =  (1.0/36.0)*scalar(L_1(3)["mnoefg"]*DT_4(3)["efgmno"]);

    printf("E401: %18.15f\n", E401);
    printf("E302: %18.15f\n", E302);
    printf("E203: %18.15f\n", E203);
    printf("E104: %18.15f\n", E104);
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
     * CCSD(T-6)
     *
     **************************************************************************/

    U E501 =  (1.0/36.0)*scalar(L_5(3)["mnoefg"]*DT_1(3)["efgmno"]);
    U E402 =  (1.0/ 1.0)*scalar(L_4(1)[    "me"]*DT_2(1)[    "em"])
             +(1.0/ 4.0)*scalar(L_4(2)[  "mnef"]*DT_2(2)[  "efmn"])
             +(1.0/36.0)*scalar(L_4(3)["mnoefg"]*DT_2(3)["efgmno"]);
    U E303 =  (1.0/ 1.0)*scalar(L_3(1)[    "me"]*DT_3(1)[    "em"])
             +(1.0/ 4.0)*scalar(L_3(2)[  "mnef"]*DT_3(2)[  "efmn"])
             +(1.0/36.0)*scalar(L_3(3)["mnoefg"]*DT_3(3)["efgmno"]);
    U E204 =  (1.0/ 1.0)*scalar(L_2(1)[    "me"]*DT_4(1)[    "em"])
             +(1.0/ 4.0)*scalar(L_2(2)[  "mnef"]*DT_4(2)[  "efmn"])
             +(1.0/36.0)*scalar(L_2(3)["mnoefg"]*DT_4(3)["efgmno"]);
    U E105 =  (1.0/36.0)*scalar(L_1(3)["mnoefg"]*DT_5(3)["efgmno"]);

    printf("E501: %18.15f\n", E501);
    printf("E402: %18.15f\n", E402);
    printf("E303: %18.15f\n", E303);
    printf("E204: %18.15f\n", E204);
    printf("E105: %18.15f\n", E105);
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

    printf("CCSD(T-2): %18.15f\n", E2);
    printf("CCSD(T-3): %18.15f\n", E3);
    printf("CCSD(T-4): %18.15f\n", E4);
    printf("CCSD(T-5): %18.15f\n", E5);
    printf("CCSD(T-6): %18.15f\n", E6);

    this->put("E(2)", new U(E2));
    this->put("E(3)", new U(E3));
    this->put("E(4)", new U(E4));
    this->put("E(5)", new U(E5));
    this->put("E(6)", new U(E6));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSD_T_N);
REGISTER_TASK(aquarius::cc::CCSD_T_N<double>,"ccsd(t-n)");
