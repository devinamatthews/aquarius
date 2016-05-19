#include "ccsd_t_n_opt.hpp"

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
CCSD_T_N_opt<U>::CCSD_T_N_opt(const string& name, Config& config)
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
bool CCSD_T_N_opt<U>::run(task::TaskDAG& dag, const Arena& arena)
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

    DeexcitationOperator<U,3> L_1("L^(1)", arena, occ, vrt);
    DeexcitationOperator<U,3> L_2("L^(2)", arena, occ, vrt);

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

    SpinorbitalTensor<U> DAI_1(T(1));
    SpinorbitalTensor<U> GIJAK_1(WMNEJ);
    SpinorbitalTensor<U> GAIBC_1(WAMEF);

    SpinorbitalTensor<U> Tau(T(2));
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
     * L^(1)
     *
     **************************************************************************/

    L_1(3)["ijkabc"]  = WMNEF["ijab"]*L(1)[  "kc"];
    L_1(3)["ijkabc"] +=   FME[  "ia"]*L(2)["jkbc"];
    L_1(3)["ijkabc"] += WAMEF["ekbc"]*L(2)["ijae"];
    L_1(3)["ijkabc"] -= WMNEJ["ijam"]*L(2)["mkbc"];

    L_1.weight(D);

    /***************************************************************************
     *
     * T^(1)
     *
     **************************************************************************/

    T_1(3)["abcijk"]  = WTWABEJ["bcek"]*T(2)["aeij"];
    T_1(3)["abcijk"] -=   WAMIJ["bmjk"]*T(2)["acim"];

    U E2 = (1.0/36.0)*scalar(L_1(3)["mnoefg"]*T_1(3)["efgmno"]);

    T_1.weight(D);

    /***************************************************************************
     *
     * L^(2)
     *
     **************************************************************************/

    WAMIJ_1[  "amij"]  =  0.5*WMNEF["mnef"]*T_1(3)["aefijn"];
    WABEJ_1[  "abej"]  = -0.5*WMNEF["mnef"]*T_1(3)["afbmnj"];

    GIJAK_1[  "ijak"]  =  (1.0/ 2.0)*   T(2)[  "efkm"]* L_1(3)["ijmaef"];
    GAIBC_1[  "aibc"]  = -(1.0/ 2.0)*   T(2)[  "aemn"]* L_1(3)["minbce"];

      DAI_1[    "ai"]  =  (1.0/ 4.0)* T_1(3)["efamni"]*   L(2)[  "mnef"];
      DAI_1[    "ai"] -=  (1.0/ 2.0)*   T(2)[  "eamn"]*GIJAK_1[  "mnei"];

     L_2(1)[    "ia"]  =  (1.0/ 2.0)*WABEJ_1[  "efam"]*   L(2)[  "imef"];
     L_2(1)[    "ia"] -=  (1.0/ 2.0)*WAMIJ_1[  "eimn"]*   L(2)[  "mnea"];
     L_2(1)[    "ia"] +=               WMNEF[  "miea"]*  DAI_1[    "em"];
     L_2(1)[    "ia"] -=  (1.0/ 2.0)*  WABEF[  "efga"]*GAIBC_1[  "gief"];
     L_2(1)[    "ia"] +=               WAMEI[  "eifm"]*GAIBC_1[  "fmea"];
     L_2(1)[    "ia"] -=               WAMEI[  "eman"]*GIJAK_1[  "inem"];
     L_2(1)[    "ia"] +=  (1.0/ 2.0)*  WMNIJ[  "imno"]*GIJAK_1[  "noam"];

     L_2(2)[  "ijab"]  = -             WAMEF[  "fiae"]*GAIBC_1[  "ejbf"];
     L_2(2)[  "ijab"] -=               WMNEJ[  "ijem"]*GAIBC_1[  "emab"];
     L_2(2)[  "ijab"] -=               WAMEF[  "emab"]*GIJAK_1[  "ijem"];
     L_2(2)[  "ijab"] -=               WMNEJ[  "niam"]*GIJAK_1[  "mjbn"];
     L_2(2)[  "ijab"] +=  (1.0/ 2.0)*  WABEJ[  "efbm"]* L_1(3)["ijmaef"];
     L_2(2)[  "ijab"] -=  (1.0/ 2.0)*  WAMIJ[  "ejnm"]* L_1(3)["imnabe"];

     L_2(3)["ijkabc"]  =               WMNEF[  "ijae"]*GAIBC_1[  "ekbc"];
     L_2(3)["ijkabc"] -=               WMNEF[  "mkbc"]*GIJAK_1[  "ijam"];
     L_2(3)["ijkabc"] +=                 FAE[    "ea"]* L_1(3)["ijkebc"];
     L_2(3)["ijkabc"] -=                 FMI[    "im"]* L_1(3)["mjkabc"];
     L_2(3)["ijkabc"] +=  (1.0/ 2.0)*  WABEF[  "efab"]* L_1(3)["ijkefc"];
     L_2(3)["ijkabc"] +=  (1.0/ 2.0)*  WMNIJ[  "ijmn"]* L_1(3)["mnkabc"];
     L_2(3)["ijkabc"] +=               WAMEI[  "eiam"]* L_1(3)["mjkbec"];

    L_2.weight(D);

    /***************************************************************************
     *
     * T^(2)
     *
     **************************************************************************/

     T_2(1)[    "ai"]  = 0.25*WMNEF["mnef"]*T_1(3)["aefimn"];

     T_2(2)[  "abij"]  =  0.5*WAMEF["bmef"]*T_1(3)["aefijm"];
     T_2(2)[  "abij"] -=  0.5*WMNEJ["mnej"]*T_1(3)["abeinm"];
     T_2(2)[  "abij"] +=        FME[  "me"]*T_1(3)["abeijm"];

     T_2(3)["abcijk"]  =    WABEJ_1["bcek"]*  T(2)[  "aeij"];
     T_2(3)["abcijk"] -=    WAMIJ_1["bmjk"]*  T(2)[  "acim"];
     T_2(3)["abcijk"] +=        FAE[  "ce"]*T_1(3)["abeijk"];
     T_2(3)["abcijk"] -=        FMI[  "mk"]*T_1(3)["abcijm"];
     T_2(3)["abcijk"] +=  0.5*WABEF["abef"]*T_1(3)["efcijk"];
     T_2(3)["abcijk"] +=  0.5*WMNIJ["mnij"]*T_1(3)["abcmnk"];
     T_2(3)["abcijk"] +=      WAMEI["amei"]*T_1(3)["ebcjmk"];

    U E3 = (1.0/36.0)*scalar(L_1(3)["mnoefg"]*T_2(3)["efgmno"]);

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
    WAMIJ_2[  "amij"] +=       FME_2[  "me"]*  T(2)[  "aeij"];

    WABEJ_2[  "abej"]  =       WAMEF["amef"]*T_2(2)[  "fbmj"];
    WABEJ_2[  "abej"] +=   0.5*WMNEJ["mnej"]*T_2(2)[  "abmn"];
    WABEJ_2[  "abej"] +=       WABEF["abef"]*T_2(1)[    "fj"];
    WABEJ_2[  "abej"] -=       WAMEI["amej"]*T_2(1)[    "bm"];

     T_3(1)[    "ai"]  =         FAE[  "ae"]*T_2(1)[    "ei"];
     T_3(1)[    "ai"] -=         FMI[  "mi"]*T_2(1)[    "am"];
     T_3(1)[    "ai"] -=       WAMEI["amei"]*T_2(1)[    "em"];
     T_3(1)[    "ai"] +=         FME[  "me"]*T_2(2)[  "aeim"];
     T_3(1)[    "ai"] +=   0.5*WAMEF["amef"]*T_2(2)[  "efim"];
     T_3(1)[    "ai"] -=   0.5*WMNEJ["mnei"]*T_2(2)[  "eamn"];

     T_3(2)[  "abij"]  =       FAE_2[  "af"]*  T(2)[  "fbij"];
     T_3(2)[  "abij"] -=       FMI_2[  "ni"]*  T(2)[  "abnj"];
     T_3(2)[  "abij"] +=       WABEJ["abej"]*T_2(1)[    "ei"];
     T_3(2)[  "abij"] -=       WAMIJ["amij"]*T_2(1)[    "bm"];
     T_3(2)[  "abij"] +=         FAE[  "af"]*T_2(2)[  "fbij"];
     T_3(2)[  "abij"] -=         FMI[  "ni"]*T_2(2)[  "abnj"];
     T_3(2)[  "abij"] +=   0.5*WABEF["abef"]*T_2(2)[  "efij"];
     T_3(2)[  "abij"] +=   0.5*WMNIJ["mnij"]*T_2(2)[  "abmn"];
     T_3(2)[  "abij"] +=       WAMEI["amei"]*T_2(2)[  "ebjm"];

    U E5_12_12 = (1.0/ 1.0)*scalar(L_2(1)[    "me"]*T_3(1)[    "em"])+
                 (1.0/ 4.0)*scalar(L_2(2)[  "mnef"]*T_3(2)[  "efmn"]);

     T_3(1)[    "ai"] +=  0.25*WMNEF["mnef"]*T_2(3)["aefimn"];

     T_3(2)[  "abij"] +=   0.5*WAMEF["bmef"]*T_2(3)["aefijm"];
     T_3(2)[  "abij"] -=   0.5*WMNEJ["mnej"]*T_2(3)["abeinm"];
     T_3(2)[  "abij"] +=         FME[  "me"]*T_2(3)["abeijm"];

    U E5_12_123 = (1.0/ 1.0)*scalar(L_2(1)[    "me"]*T_3(1)[    "em"])+
                  (1.0/ 4.0)*scalar(L_2(2)[  "mnef"]*T_3(2)[  "efmn"]);
    U E5_12_3 = E5_12_123-E5_12_12;

     T_3(3)["abcijk"]  =     WABEJ_2["bcek"]*  T(2)[  "aeij"];
     T_3(3)["abcijk"] -=     WAMIJ_2["bmjk"]*  T(2)[  "acim"];
     T_3(3)["abcijk"] +=       WABEJ["bcek"]*T_2(2)[  "aeij"];
     T_3(3)["abcijk"] -=       WAMIJ["bmjk"]*T_2(2)[  "acim"];

    U E4_3_12 = (1.0/36.0)*scalar(L_1(3)["mnoefg"]*T_3(3)["efgmno"]);

    U E5_3_12 = (1.0/36.0)*scalar(L_2(3)["mnoefg"]*T_3(3)["efgmno"]);

    WAMIJ_2[  "amij"]  =   0.5*WMNEF["mnef"]*T_2(3)["aefijn"];
    WABEJ_2[  "abej"]  =  -0.5*WMNEF["mnef"]*T_2(3)["afbmnj"];

     T_3(3)["abcijk"] +=     WABEJ_2["bcek"]*  T(2)[  "aeij"];
     T_3(3)["abcijk"] -=     WAMIJ_2["bmjk"]*  T(2)[  "acim"];
     T_3(3)["abcijk"] +=         FAE[  "ce"]*T_2(3)["abeijk"];
     T_3(3)["abcijk"] -=         FMI[  "mk"]*T_2(3)["abcijm"];
     T_3(3)["abcijk"] +=   0.5*WABEF["abef"]*T_2(3)["efcijk"];
     T_3(3)["abcijk"] +=   0.5*WMNIJ["mnij"]*T_2(3)["abcmnk"];
     T_3(3)["abcijk"] +=       WAMEI["amei"]*T_2(3)["ebcjmk"];

    U E4_3_123 = (1.0/36.0)*scalar(L_1(3)["mnoefg"]*T_3(3)["efgmno"]);
    U E4_3_3 = E4_3_123-E4_3_12;

    U E5_3_123 = (1.0/36.0)*scalar(L_2(3)["mnoefg"]*T_3(3)["efgmno"]);
    U E5_3_3 = E5_3_123-E5_3_12;

    T_3.weight(D);

    /***************************************************************************
     *
     * T^(4)
     *
     **************************************************************************/

      FME_3[    "me"]  =       WMNEF["mnef"]*T_3(1)[    "fn"];

      FMI_3[    "mi"]  =   0.5*WMNEF["mnef"]*T_3(2)[  "efin"];
      FMI_3[    "mi"] +=       WMNEJ["nmfi"]*T_3(1)[    "fn"];

      FAE_3[    "ae"]  =  -0.5*WMNEF["mnef"]*T_3(2)[  "afmn"];
      FAE_3[    "ae"] +=       WAMEF["amef"]*T_3(1)[    "fm"];

    WAMIJ_3[  "amij"]  =       WMNEJ["nmej"]*T_3(2)[  "aein"];
    WAMIJ_3[  "amij"] +=   0.5*WAMEF["amef"]*T_3(2)[  "efij"];
    WAMIJ_3[  "amij"] -=       WMNIJ["nmij"]*T_3(1)[    "an"];
    WAMIJ_3[  "amij"] +=       WAMEI["amej"]*T_3(1)[    "ei"];
    WAMIJ_3[  "amij"] +=       FME_3[  "me"]*  T(2)[  "aeij"];

    WABEJ_3[  "abej"]  =       WAMEF["amef"]*T_3(2)[  "fbmj"];
    WABEJ_3[  "abej"] +=   0.5*WMNEJ["mnej"]*T_3(2)[  "abmn"];
    WABEJ_3[  "abej"] +=       WABEF["abef"]*T_3(1)[    "fj"];
    WABEJ_3[  "abej"] -=       WAMEI["amej"]*T_3(1)[    "bm"];

     T_4(1)[    "ai"]  =         FAE[  "ae"]*T_3(1)[    "ei"];
     T_4(1)[    "ai"] -=         FMI[  "mi"]*T_3(1)[    "am"];
     T_4(1)[    "ai"] -=       WAMEI["amei"]*T_3(1)[    "em"];
     T_4(1)[    "ai"] +=         FME[  "me"]*T_3(2)[  "aeim"];
     T_4(1)[    "ai"] +=   0.5*WAMEF["amef"]*T_3(2)[  "efim"];
     T_4(1)[    "ai"] -=   0.5*WMNEJ["mnei"]*T_3(2)[  "eamn"];

     T_4(2)[  "abij"]  =       FAE_3[  "af"]*  T(2)[  "fbij"];
     T_4(2)[  "abij"] -=       FMI_3[  "ni"]*  T(2)[  "abnj"];
     T_4(2)[  "abij"] +=       WABEJ["abej"]*T_3(1)[    "ei"];
     T_4(2)[  "abij"] -=       WAMIJ["amij"]*T_3(1)[    "bm"];
     T_4(2)[  "abij"] +=         FAE[  "af"]*T_3(2)[  "fbij"];
     T_4(2)[  "abij"] -=         FMI[  "ni"]*T_3(2)[  "abnj"];
     T_4(2)[  "abij"] +=   0.5*WABEF["abef"]*T_3(2)[  "efij"];
     T_4(2)[  "abij"] +=   0.5*WMNIJ["mnij"]*T_3(2)[  "abmn"];
     T_4(2)[  "abij"] +=       WAMEI["amei"]*T_3(2)[  "ebjm"];

    U E6_12_12 = (1.0/ 1.0)*scalar(L_2(1)[    "me"]*T_4(1)[    "em"])+
                 (1.0/ 4.0)*scalar(L_2(2)[  "mnef"]*T_4(2)[  "efmn"]);

     T_4(1)[    "ai"] +=  0.25*WMNEF["mnef"]*T_3(3)["aefimn"];

     T_4(2)[  "abij"] +=   0.5*WAMEF["bmef"]*T_3(3)["aefijm"];
     T_4(2)[  "abij"] -=   0.5*WMNEJ["mnej"]*T_3(3)["abeinm"];
     T_4(2)[  "abij"] +=         FME[  "me"]*T_3(3)["abeijm"];

    U E6_12_123 = (1.0/ 1.0)*scalar(L_2(1)[    "me"]*T_4(1)[    "em"])+
                  (1.0/ 4.0)*scalar(L_2(2)[  "mnef"]*T_4(2)[  "efmn"]);
    U E6_12_3 = E6_12_123-E6_12_12;

     T_4(3)["abcijk"]  =     WABEJ_3["bcek"]*  T(2)[  "aeij"];
     T_4(3)["abcijk"] -=     WAMIJ_3["bmjk"]*  T(2)[  "acim"];
     T_4(3)["abcijk"] +=       WABEJ["bcek"]*T_3(2)[  "aeij"];
     T_4(3)["abcijk"] -=       WAMIJ["bmjk"]*T_3(2)[  "acim"];

    U E6_3_12 = (1.0/36.0)*scalar(L_2(3)["mnoefg"]*T_4(3)["efgmno"]);

    WAMIJ_3[  "amij"]  =   0.5*WMNEF["mnef"]*T_3(3)["aefijn"];
    WABEJ_3[  "abej"]  =  -0.5*WMNEF["mnef"]*T_3(3)["afbmnj"];

     T_4(3)["abcijk"] +=     WABEJ_3["bcek"]*  T(2)[  "aeij"];
     T_4(3)["abcijk"] -=     WAMIJ_3["bmjk"]*  T(2)[  "acim"];
     T_4(3)["abcijk"] +=         FAE[  "ce"]*T_3(3)["abeijk"];
     T_4(3)["abcijk"] -=         FMI[  "mk"]*T_3(3)["abcijm"];
     T_4(3)["abcijk"] +=   0.5*WABEF["abef"]*T_3(3)["efcijk"];
     T_4(3)["abcijk"] +=   0.5*WMNIJ["mnij"]*T_3(3)["abcmnk"];
     T_4(3)["abcijk"] +=       WAMEI["amei"]*T_3(3)["ebcjmk"];

    U E6_3_123 = (1.0/36.0)*scalar(L_2(3)["mnoefg"]*T_4(3)["efgmno"]);
    U E6_3_3 = E6_3_123-E6_3_12;

    /***************************************************************************
     *
     * Complete W^(2) and W^(3)
     *
     **************************************************************************/

    FTWMI_2 = FMI_2;
    FTWAE_2 = FAE_2;

      FMI_2[  "mi"] +=        FME[  "me"]*T_2(1)[    "ei"];
      FAE_2[  "ae"] -=        FME[  "me"]*T_2(1)[    "am"];

    WAMEI_2["amei"]  =      WMNEF["mnef"]*T_2(2)[  "afni"];
    WAMEI_2["amei"] +=      WAMEF["amef"]*T_2(1)[    "fi"];
    WAMEI_2["amei"] -=      WMNEJ["nmei"]*T_2(1)[    "an"];

    WMNIJ_2["mnij"]  =  0.5*WMNEF["mnef"]*T_2(2)[  "efij"];
    WMNIJ_2["mnij"] +=      WMNEJ["mnej"]*T_2(1)[    "ei"];

    WABEF_2["abef"]  =  0.5*WMNEF["mnef"]*T_2(2)[  "abmn"];
    WABEF_2["abef"] -=      WAMEF["amef"]*T_2(1)[    "bm"];

    WMNEJ_2["mnej"]  =      WMNEF["mnef"]*T_2(1)[    "fj"];

    WAMEF_2["amef"]  = -    WMNEF["nmef"]*T_2(1)[    "an"];

    WAMIJ_2["amij"]  =      WMNEJ["nmej"]*T_2(2)[  "aein"];
    WAMIJ_2["amij"] +=  0.5*WAMEF["amef"]*T_2(2)[  "efij"];
    WAMIJ_2["amij"] -=      WMNIJ["nmij"]*T_2(1)[    "an"];
    WAMIJ_2["amij"] +=      WAMEI["amej"]*T_2(1)[    "ei"];
    WAMIJ_2["amij"] +=  0.5*WMNEF["mnef"]*T_2(3)["aefijn"];
    WAMIJ_2["amij"] +=      FME_2[  "me"]*  T(2)[  "aeij"];

    WABEJ_2["abej"]  =      WAMEF["amef"]*T_2(2)[  "fbmj"];
    WABEJ_2["abej"] +=  0.5*WMNEJ["mnej"]*T_2(2)[  "abmn"];
    WABEJ_2["abej"] +=      WABEF["abef"]*T_2(1)[    "fj"];
    WABEJ_2["abej"] -=      WAMEI["amej"]*T_2(1)[    "bm"];
    WABEJ_2["abej"] -=  0.5*WMNEF["mnef"]*T_2(3)["afbmnj"];
    WABEJ_2["abej"] -=      FME_2[  "me"]*  T(2)[  "abmj"];

    FTWMI_3 = FMI_3;
    FTWAE_3 = FAE_3;

      FMI_3[  "mi"] +=        FME[  "me"]*T_3(1)[    "ei"];
      FAE_3[  "ae"] -=        FME[  "me"]*T_3(1)[    "am"];

    WAMEI_3["amei"]  =      WMNEF["mnef"]*T_3(2)[  "afni"];
    WAMEI_3["amei"] +=      WAMEF["amef"]*T_3(1)[    "fi"];
    WAMEI_3["amei"] -=      WMNEJ["nmei"]*T_3(1)[    "an"];

    WMNIJ_3["mnij"]  =  0.5*WMNEF["mnef"]*T_3(2)[  "efij"];
    WMNIJ_3["mnij"] +=      WMNEJ["mnej"]*T_3(1)[    "ei"];

    WABEF_3["abef"]  =  0.5*WMNEF["mnef"]*T_3(2)[  "abmn"];
    WABEF_3["abef"] -=      WAMEF["amef"]*T_3(1)[    "bm"];

    WMNEJ_3["mnej"]  =      WMNEF["mnef"]*T_3(1)[    "fj"];
    WAMEF_3["amef"]  = -    WMNEF["nmef"]*T_3(1)[    "an"];
    WABEJ_3["abej"] -=      FME_3[  "me"]*  T(2)[  "abmj"];

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

    printf("\n");
    printf("CCSD(T-2): %18.15f\n", E2);

    printf("\n");
    printf("CCSD(T-3): %18.15f\n", E3);

    U E4 = E4_3_12+E4_3_3+E0112;

    printf("\n");
    printf("<0|L3^[1][Hbar,T1^[2]+T2^[2]]|0>:  %18.15f\n", E4_3_12);
    printf("<0|L3^[1][Hbar,T3^[2]]|0>:         %18.15f\n", E4_3_3);
    printf("<0|(1+L)[[Hbar,T3^[1]],T1^[2]]|0>: %18.15f\n", E0112);
    printf("CCSD(T-4) Total:                   %18.15f\n", E4);

    U E5 = E5_12_12+E5_12_3+E5_3_12+E5_3_3+E0122+E1112;

    printf("\n");
    printf("<0|(L1^[2]+L2^[2])[Hbar,T1^[2]+T2^[2]]|0>: %18.15f\n", E5_12_12);
    printf("<0|(L1^[2]+L2^[2])[Hbar,T3^[2]]|0>:        %18.15f\n", E5_12_3);
    printf("<0|L3^[2][Hbar,T1^[2]+T2^[2]]|0>:          %18.15f\n", E5_3_12);
    printf("<0|L3^[2][Hbar,T3^[2]]|0>:                 %18.15f\n", E5_3_3);
    printf("1/2 <0|(1+L)[[Hbar,T^[2]],T^[2]]|0>:       %18.15f\n", E0122);
    printf("<0|L3^[1][[Hbar,T3^[1]],T^[2]]|0>:         %18.15f\n", E1112);
    printf("CCSD(T-5) Total:                           %18.15f\n", E5);

    U E6 = E6_12_12+E6_12_3+E6_3_12+E6_3_3+E0123+E1113+E1122+E2112;

    printf("\n");
    printf("<0|(L1^[2]+L2^[2])[Hbar,T1^[3]+T2^[3]]|0>: %18.15f\n", E6_12_12);
    printf("<0|(L1^[2]+L2^[2])[Hbar,T3^[3]]|0>:        %18.15f\n", E6_12_3);
    printf("<0|L3^[2][Hbar,T1^[3]+T2^[3]]|0>:          %18.15f\n", E6_3_12);
    printf("<0|L3^[2][Hbar,T3^[3]]|0>:                 %18.15f\n", E6_3_3);
    printf("<0|(1+L)[[Hbar,T^[2]],T^[3]]|0>:           %18.15f\n", E0123);
    printf("<0|L3^[1][[Hbar,T3^[1]],T^[3]]|0>:         %18.15f\n", E1113);
    printf("1/2 <0|L3^[1][[Hbar,T^[2]],T^[2]]|0>:      %18.15f\n", E1122);
    printf("<0|L^[2][[Hbar,T3^[1]],T^[2]]|0>:          %18.15f\n", E2112);
    printf("CCSD(T-6) Total:                           %18.15f\n", E6);

    printf("\n");

    this->put("E(2)", new U(E2));
    this->put("E(3)", new U(E3));
    this->put("E(4)", new U(E4));
    this->put("E(5)", new U(E5));
    this->put("E(6)", new U(E6));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSD_T_N_opt);
REGISTER_TASK(aquarius::cc::CCSD_T_N_opt<double>,"ccsd(t-n)_opt");
