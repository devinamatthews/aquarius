#include "e_ccsd_t_n.hpp"

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
E_CCSD_T_N<U>::E_CCSD_T_N(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("ccsd.T", "T"));
    this->addProduct(Product("double", "E(3)", reqs));
    this->addProduct(Product("double", "E(4)", reqs));
    this->addProduct(Product("double", "E(5)", reqs));
    this->addProduct(Product("double", "E(6)", reqs));
    this->addProduct(Product("double", "E(7)", reqs));
}

template <typename U>
bool E_CCSD_T_N<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U>>("H");
    const STTwoElectronOperator<U>& Hbar = this->template get<STTwoElectronOperator<U>>("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    Denominator<U> D(H);
    const ExcitationOperator<U,2>& T = this->template get<ExcitationOperator<U,2>>("T");

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
    ExcitationOperator<U,3> T_6("T^(6)", arena, occ, vrt);

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
    SpinorbitalTensor<U> FTWAE_2(FAE);
    SpinorbitalTensor<U> WTWAMEI_2(WAMEI);
    SpinorbitalTensor<U> WTWABEJ_2(WABEJ);
    SpinorbitalTensor<U> WTWAMIJ_2(WAMIJ);

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
    SpinorbitalTensor<U> FTWAE_3(FAE);
    SpinorbitalTensor<U> WTWABEJ_3(WABEJ);
    SpinorbitalTensor<U> WTWAMIJ_3(WAMIJ);

    SpinorbitalTensor<U> Tau_4(T(2));
    SpinorbitalTensor<U> FME_4(FME);
    SpinorbitalTensor<U> FMI_4(FMI);
    SpinorbitalTensor<U> FAE_4(FAE);
    SpinorbitalTensor<U> WAMEI_4(WAMEI);
    SpinorbitalTensor<U> WMNIJ_4(WMNIJ);
    SpinorbitalTensor<U> WABEF_4(WABEF);
    SpinorbitalTensor<U> WABEJ_4(WABEJ);
    SpinorbitalTensor<U> WAMIJ_4(WAMIJ);
    SpinorbitalTensor<U> WMNEJ_4(WMNEJ);
    SpinorbitalTensor<U> WAMEF_4(WAMEF);
    SpinorbitalTensor<U> WTWABEJ_4(WABEJ);
    SpinorbitalTensor<U> WTWAMIJ_4(WAMIJ);

    SpinorbitalTensor<U> FMI_5(FMI);
    SpinorbitalTensor<U> FAE_5(FAE);

    /***************************************************************************
     *
     * T^(1)
     *
     **************************************************************************/

    T_1(3)["abcijk"]  = WTWABEJ["bcek"]*T(2)["aeij"];
    T_1(3)["abcijk"] -=   WAMIJ["bmjk"]*T(2)["acim"];

    T_1.weight(D);

    /***************************************************************************
     *
     * T^(2)
     *
     **************************************************************************/

    WAMIJ_1[  "amij"]  =  0.5*WMNEF["mnef"]*T_1(3)["aefijn"];

    WABEJ_1[  "abej"]  = -0.5*WMNEF["mnef"]*T_1(3)["afbmnj"];

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

    T_2.weight(D);

    U E3 = (1.0/1.0)*scalar(  FME[  "me"]*T_2(1)[  "em"])+
           (1.0/4.0)*scalar(WMNEF["mnef"]*T_2(2)["efmn"]);

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
    WAMIJ_2[  "amij"] +=   0.5*WMNEF["mnef"]*T_2(3)["aefijn"];

    WABEJ_2[  "abej"]  =       WAMEF["amef"]*T_2(2)[  "fbmj"];
    WABEJ_2[  "abej"] +=   0.5*WMNEJ["mnej"]*T_2(2)[  "abmn"];
    WABEJ_2[  "abej"] +=       WABEF["abef"]*T_2(1)[    "fj"];
    WABEJ_2[  "abej"] -=       WAMEI["amej"]*T_2(1)[    "bm"];
    WABEJ_2[  "abej"] -=   0.5*WMNEF["mnef"]*T_2(3)["afbmnj"];

     T_3(1)[    "ai"]  =         FAE[  "ae"]*T_2(1)[    "ei"];
     T_3(1)[    "ai"] -=         FMI[  "mi"]*T_2(1)[    "am"];
     T_3(1)[    "ai"] -=       WAMEI["amei"]*T_2(1)[    "em"];
     T_3(1)[    "ai"] +=         FME[  "me"]*T_2(2)[  "aeim"];
     T_3(1)[    "ai"] +=   0.5*WAMEF["amef"]*T_2(2)[  "efim"];
     T_3(1)[    "ai"] -=   0.5*WMNEJ["mnei"]*T_2(2)[  "eamn"];
     T_3(1)[    "ai"] +=  0.25*WMNEF["mnef"]*T_2(3)["aefimn"];

     T_3(2)[  "abij"]  =       FAE_2[  "af"]*  T(2)[  "fbij"];
     T_3(2)[  "abij"] -=       FMI_2[  "ni"]*  T(2)[  "abnj"];
     T_3(2)[  "abij"] +=       WABEJ["abej"]*T_2(1)[    "ei"];
     T_3(2)[  "abij"] -=       WAMIJ["amij"]*T_2(1)[    "bm"];
     T_3(2)[  "abij"] +=         FAE[  "af"]*T_2(2)[  "fbij"];
     T_3(2)[  "abij"] -=         FMI[  "ni"]*T_2(2)[  "abnj"];
     T_3(2)[  "abij"] +=   0.5*WABEF["abef"]*T_2(2)[  "efij"];
     T_3(2)[  "abij"] +=   0.5*WMNIJ["mnij"]*T_2(2)[  "abmn"];
     T_3(2)[  "abij"] +=       WAMEI["amei"]*T_2(2)[  "ebjm"];
     T_3(2)[  "abij"] +=   0.5*WAMEF["bmef"]*T_2(3)["aefijm"];
     T_3(2)[  "abij"] -=   0.5*WMNEJ["mnej"]*T_2(3)["abeinm"];
     T_3(2)[  "abij"] +=         FME[  "me"]*T_2(3)["abeijm"];

     T_3(3)["abcijk"]  =     WABEJ_2["bcek"]*  T(2)[  "aeij"];
     T_3(3)["abcijk"] -=     WAMIJ_2["bmjk"]*  T(2)[  "acim"];
     T_3(3)["abcijk"] +=       WABEJ["bcek"]*T_2(2)[  "aeij"];
     T_3(3)["abcijk"] -=       WAMIJ["bmjk"]*T_2(2)[  "acim"];
     T_3(3)["abcijk"] +=         FAE[  "ce"]*T_2(3)["abeijk"];
     T_3(3)["abcijk"] -=         FMI[  "mk"]*T_2(3)["abcijm"];
     T_3(3)["abcijk"] +=   0.5*WABEF["abef"]*T_2(3)["efcijk"];
     T_3(3)["abcijk"] +=   0.5*WMNIJ["mnij"]*T_2(3)["abcmnk"];
     T_3(3)["abcijk"] +=       WAMEI["amei"]*T_2(3)["ebcjmk"];

    T_3.weight(D);

    U E4 = (1.0/1.0)*scalar(  FME[  "me"]*T_3(1)[  "em"])+
           (1.0/4.0)*scalar(WMNEF["mnef"]*T_3(2)["efmn"]);

    /***************************************************************************
     *
     * T^(4)
     *
     **************************************************************************/

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
    WAMIJ_3[  "amij"] +=       FME_3[  "me"]*  T(2)[  "aeij"];
    WAMIJ_3[  "amij"] +=   0.5*WMNEF["mnef"]*T_3(3)["aefijn"];

    WABEJ_3[  "abej"]  =       WAMEF["amef"]*T_3(2)[  "fbmj"];
    WABEJ_3[  "abej"] +=   0.5*WMNEJ["mnej"]*T_3(2)[  "abmn"];
    WABEJ_3[  "abej"] +=       WABEF["abef"]*T_3(1)[    "fj"];
    WABEJ_3[  "abej"] -=       WAMEI["amej"]*T_3(1)[    "bm"];
    WABEJ_3[  "abej"] -=   0.5*WMNEF["mnef"]*T_3(3)["afbmnj"];

     T_4(1)[    "ai"]  =         FAE[  "ae"]*T_3(1)[    "ei"];
     T_4(1)[    "ai"] -=         FMI[  "mi"]*T_3(1)[    "am"];
     T_4(1)[    "ai"] -=       WAMEI["amei"]*T_3(1)[    "em"];
     T_4(1)[    "ai"] +=         FME[  "me"]*T_3(2)[  "aeim"];
     T_4(1)[    "ai"] +=   0.5*WAMEF["amef"]*T_3(2)[  "efim"];
     T_4(1)[    "ai"] -=   0.5*WMNEJ["mnei"]*T_3(2)[  "eamn"];
     T_4(1)[    "ai"] +=  0.25*WMNEF["mnef"]*T_3(3)["aefimn"];

     T_4(2)[  "abij"]  =       FAE_3[  "af"]*  T(2)[  "fbij"];
     T_4(2)[  "abij"] -=       FMI_3[  "ni"]*  T(2)[  "abnj"];
     T_4(2)[  "abij"] +=       FME_2[  "me"]*T_1(3)["abeijm"];
     T_4(2)[  "abij"] +=     WABEJ_1["abej"]*T_2(1)[    "ei"];
     T_4(2)[  "abij"] -=     WAMIJ_1["amij"]*T_2(1)[    "bm"];
     T_4(2)[  "abij"] +=       WABEJ["abej"]*T_3(1)[    "ei"];
     T_4(2)[  "abij"] -=       WAMIJ["amij"]*T_3(1)[    "bm"];
     T_4(2)[  "abij"] +=         FAE[  "af"]*T_3(2)[  "fbij"];
     T_4(2)[  "abij"] -=         FMI[  "ni"]*T_3(2)[  "abnj"];
     T_4(2)[  "abij"] +=   0.5*WABEF["abef"]*T_3(2)[  "efij"];
     T_4(2)[  "abij"] +=   0.5*WMNIJ["mnij"]*T_3(2)[  "abmn"];
     T_4(2)[  "abij"] +=       WAMEI["amei"]*T_3(2)[  "ebjm"];
     T_4(2)[  "abij"] +=   0.5*WAMEF["bmef"]*T_3(3)["aefijm"];
     T_4(2)[  "abij"] -=   0.5*WMNEJ["mnej"]*T_3(3)["abeinm"];
     T_4(2)[  "abij"] +=         FME[  "me"]*T_3(3)["abeijm"];

     T_4(3)["abcijk"]  =     WABEJ_3["bcek"]*  T(2)[  "aeij"];
     T_4(3)["abcijk"] -=     WAMIJ_3["bmjk"]*  T(2)[  "acim"];
     T_4(3)["abcijk"] +=       FAE_2[  "ae"]*T_1(3)["ebcijk"];
     T_4(3)["abcijk"] -=       FMI_2[  "mi"]*T_1(3)["abcmjk"];
     T_4(3)["abcijk"] += 0.5*WABEF_2["abef"]*T_1(3)["efcijk"];
     T_4(3)["abcijk"] += 0.5*WMNIJ_2["mnij"]*T_1(3)["abcmnk"];
     T_4(3)["abcijk"] +=     WAMEI_2["amei"]*T_1(3)["becmjk"];
     T_4(3)["abcijk"] +=     WABEJ_1["bcek"]*T_2(2)[  "aeij"];
     T_4(3)["abcijk"] -=     WAMIJ_1["bmjk"]*T_2(2)[  "acim"];
     T_4(3)["abcijk"] +=       WABEJ["bcek"]*T_3(2)[  "aeij"];
     T_4(3)["abcijk"] -=       WAMIJ["bmjk"]*T_3(2)[  "acim"];
     T_4(3)["abcijk"] +=         FAE[  "ce"]*T_3(3)["abeijk"];
     T_4(3)["abcijk"] -=         FMI[  "mk"]*T_3(3)["abcijm"];
     T_4(3)["abcijk"] +=   0.5*WABEF["abef"]*T_3(3)["efcijk"];
     T_4(3)["abcijk"] +=   0.5*WMNIJ["mnij"]*T_3(3)["abcmnk"];
     T_4(3)["abcijk"] +=       WAMEI["amei"]*T_3(3)["ebcjmk"];

    T_4.weight(D);

    U E5 = (1.0/2.0)*scalar(FME_2[  "me"]*T_2(1)[  "em"])+
           (1.0/1.0)*scalar(  FME[  "me"]*T_4(1)[  "em"])+
           (1.0/4.0)*scalar(WMNEF["mnef"]*T_4(2)["efmn"]);

    /***************************************************************************
     *
     * T^(5)
     *
     **************************************************************************/

        Tau_4[  "abij"]  =  0.5*T_2(1)[  "ai"]*T_2(1)[    "bj"];

      FTWAE_2[    "ae"]  =  -0.5*WMNEF["mnef"]*T_2(2)[  "afmn"];
      FTWAE_2[    "ae"] +=       WAMEF["amef"]*T_2(1)[    "fm"];

    WTWAMIJ_2[  "amij"]  =   0.5*WAMEF["amef"]*T_2(2)[  "efij"];
    WTWAMIJ_2[  "amij"] -=   0.5*WMNIJ["nmij"]*T_2(1)[    "an"];
    WTWAMIJ_2[  "amij"] +=       WAMEI["amej"]*T_2(1)[    "ei"];
    WTWAMIJ_2[  "amij"] +=       FME_2[  "me"]*  T(2)[  "aeij"];
    WTWAMIJ_2[  "amij"] +=   0.5*WMNEF["mnef"]*T_2(3)["aefijn"];

    WTWABEJ_2[  "abej"]  =   0.5*WABEF["abef"]*T_2(1)[    "fj"];
    WTWABEJ_2[  "abej"] -=       FME_2[  "me"]*  T(2)[  "abmj"];
    WTWABEJ_2[  "abej"] -=   0.5*WMNEF["mnef"]*T_2(3)["afbmnj"];

    WTWAMEI_2[  "amei"]  =   0.5*WMNEF["mnef"]*T_2(2)[  "afni"];
    WTWAMEI_2[  "amei"] +=       WAMEF["amef"]*T_2(1)[    "fi"];
    WTWAMEI_2[  "amei"] -=       WMNEJ["nmei"]*T_2(1)[    "an"];

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

        FMI_4[    "mi"]  =   0.5*WMNEF["mnef"]*T_4(2)[  "efin"];
        FMI_4[    "mi"] +=       WMNEJ["nmfi"]*T_4(1)[    "fn"];

        FAE_4[    "ae"]  =  -0.5*WMNEF["mnef"]*T_4(2)[  "afmn"];
        FAE_4[    "ae"] +=       WAMEF["amef"]*T_4(1)[    "fm"];

      WAMIJ_4[  "amij"]  =       WMNEJ["nmej"]*T_4(2)[  "aein"];
      WAMIJ_4[  "amij"] +=   0.5*WAMEF["amef"]*T_4(2)[  "efij"];
      WAMIJ_4[  "amij"] -=       WMNIJ["nmij"]*T_4(1)[    "an"];
      WAMIJ_4[  "amij"] +=       WAMEI["amej"]*T_4(1)[    "ei"];
      WAMIJ_4[  "amij"] +=       FME_4[  "me"]*  T(2)[  "aeij"];
      WAMIJ_4[  "amij"] +=   0.5*WMNEF["mnef"]*T_4(3)["aefijn"];
      WAMIJ_4[  "amij"] +=     WMNEJ_2["nmej"]*T_2(2)[  "aein"];
      WAMIJ_4[  "amij"] +=   0.5*WAMEF["amef"]* Tau_4[  "efij"];
      WAMIJ_4[  "amij"] -=     WMNIJ_2["nmij"]*T_2(1)[    "an"];

      WABEJ_4[  "abej"]  =       WAMEF["amef"]*T_4(2)[  "fbmj"];
      WABEJ_4[  "abej"] +=   0.5*WMNEJ["mnej"]*T_4(2)[  "abmn"];
      WABEJ_4[  "abej"] +=       WABEF["abef"]*T_4(1)[    "fj"];
      WABEJ_4[  "abej"] -=       WAMEI["amej"]*T_4(1)[    "bm"];
      WABEJ_4[  "abej"] -=   0.5*WMNEF["mnef"]*T_4(3)["afbmnj"];
      WABEJ_4[  "abej"] +=     WAMEF_2["amef"]*T_2(2)[  "fbmj"];
      WABEJ_4[  "abej"] +=   0.5*WMNEJ["mnej"]* Tau_4[  "abmn"];
      WABEJ_4[  "abej"] +=     WABEF_2["abef"]*T_2(1)[    "fj"];

       T_5(1)[    "ai"]  =     FTWAE_2[  "ae"]*T_2(1)[    "ei"];
       T_5(1)[    "ai"] -=       FMI_2[  "mi"]*T_2(1)[    "am"];
       T_5(1)[    "ai"] +=       FME_2[  "me"]*T_2(2)[  "eami"];
       T_5(1)[    "ai"] +=         FAE[  "ae"]*T_4(1)[    "ei"];
       T_5(1)[    "ai"] -=         FMI[  "mi"]*T_4(1)[    "am"];
       T_5(1)[    "ai"] -=       WAMEI["amei"]*T_4(1)[    "em"];
       T_5(1)[    "ai"] +=         FME[  "me"]*T_4(2)[  "aeim"];
       T_5(1)[    "ai"] +=   0.5*WAMEF["amef"]*T_4(2)[  "efim"];
       T_5(1)[    "ai"] -=   0.5*WMNEJ["mnei"]*T_4(2)[  "eamn"];
       T_5(1)[    "ai"] +=  0.25*WMNEF["mnef"]*T_4(3)["aefimn"];

       T_5(2)[  "abij"]  =       FAE_4[  "af"]*  T(2)[  "fbij"];
       T_5(2)[  "abij"] -=       FMI_4[  "ni"]*  T(2)[  "abnj"];
       T_5(2)[  "abij"] +=       FME_3[  "me"]*T_1(3)["abeijm"];
       T_5(2)[  "abij"] +=   WTWABEJ_2["abej"]*T_2(1)[    "ei"];
       T_5(2)[  "abij"] -=   WTWAMIJ_2["amij"]*T_2(1)[    "bm"];
       T_5(2)[  "abij"] +=       FAE_2[  "ae"]*T_2(2)[  "ebij"];
       T_5(2)[  "abij"] -=       FMI_2[  "mi"]*T_2(2)[  "abmj"];
       T_5(2)[  "abij"] += 0.5*WMNIJ_2["mnij"]*T_2(2)[  "abmn"];
       T_5(2)[  "abij"] +=   WTWAMEI_2["amei"]*T_2(2)[  "ebjm"];
       T_5(2)[  "abij"] +=       FME_2[  "me"]*T_2(3)["eabmij"];
       T_5(2)[  "abij"] +=     WABEJ_1["abej"]*T_3(1)[    "ei"];
       T_5(2)[  "abij"] -=     WAMIJ_1["amij"]*T_3(1)[    "bm"];
       T_5(2)[  "abij"] +=       WABEJ["abej"]*T_4(1)[    "ei"];
       T_5(2)[  "abij"] -=       WAMIJ["amij"]*T_4(1)[    "bm"];
       T_5(2)[  "abij"] +=         FAE[  "af"]*T_4(2)[  "fbij"];
       T_5(2)[  "abij"] -=         FMI[  "ni"]*T_4(2)[  "abnj"];
       T_5(2)[  "abij"] +=   0.5*WABEF["abef"]*T_4(2)[  "efij"];
       T_5(2)[  "abij"] +=   0.5*WMNIJ["mnij"]*T_4(2)[  "abmn"];
       T_5(2)[  "abij"] +=       WAMEI["amei"]*T_4(2)[  "ebjm"];
       T_5(2)[  "abij"] +=   0.5*WAMEF["bmef"]*T_4(3)["aefijm"];
       T_5(2)[  "abij"] -=   0.5*WMNEJ["mnej"]*T_4(3)["abeinm"];
       T_5(2)[  "abij"] +=         FME[  "me"]*T_4(3)["abeijm"];

    WTWAMIJ_2[  "amij"]  =       WMNEJ["nmej"]*T_2(2)[  "aein"];
    WTWAMIJ_2[  "amij"] +=   0.5*WAMEF["amef"]*T_2(2)[  "efij"];
    WTWAMIJ_2[  "amij"] -=       WMNIJ["nmij"]*T_2(1)[    "an"];
    WTWAMIJ_2[  "amij"] +=       WAMEI["amej"]*T_2(1)[    "ei"];
    WTWAMIJ_2[  "amij"] +=   0.5*WMNEF["mnef"]*T_2(3)["aefijn"];
    WTWAMIJ_2[  "amij"] +=       FME_2[  "me"]*  T(2)[  "aeij"];
    WTWAMIJ_2[  "amij"] +=         FME[  "me"]*T_2(2)[  "aeij"];

       T_5(3)["abcijk"]  =     WABEJ_4["bcek"]*  T(2)[  "aeij"];
       T_5(3)["abcijk"] -=     WAMIJ_4["bmjk"]*  T(2)[  "acim"];
       T_5(3)["abcijk"] +=       FAE_3[  "ae"]*T_1(3)["ebcijk"];
       T_5(3)["abcijk"] -=       FMI_3[  "mi"]*T_1(3)["abcmjk"];
       T_5(3)["abcijk"] += 0.5*WABEF_3["abef"]*T_1(3)["efcijk"];
       T_5(3)["abcijk"] += 0.5*WMNIJ_3["mnij"]*T_1(3)["abcmnk"];
       T_5(3)["abcijk"] +=     WAMEI_3["amei"]*T_1(3)["becmjk"];
       T_5(3)["abcijk"] +=     WABEJ_2["bcek"]*T_2(2)[  "aeij"];
       T_5(3)["abcijk"] -=   WTWAMIJ_2["bmjk"]*T_2(2)[  "acim"];
       T_5(3)["abcijk"] +=       FAE_2[  "ce"]*T_2(3)["abeijk"];
       T_5(3)["abcijk"] -=       FMI_2[  "mk"]*T_2(3)["abcijm"];
       T_5(3)["abcijk"] += 0.5*WABEF_2["abef"]*T_2(3)["efcijk"];
       T_5(3)["abcijk"] += 0.5*WMNIJ_2["mnij"]*T_2(3)["abcmnk"];
       T_5(3)["abcijk"] +=     WAMEI_2["amei"]*T_2(3)["ebcjmk"];
       T_5(3)["abcijk"] +=     WABEJ_1["bcek"]*T_3(2)[  "aeij"];
       T_5(3)["abcijk"] -=     WAMIJ_1["bmjk"]*T_3(2)[  "acim"];
       T_5(3)["abcijk"] +=       WABEJ["bcek"]*T_4(2)[  "aeij"];
       T_5(3)["abcijk"] -=       WAMIJ["bmjk"]*T_4(2)[  "acim"];
       T_5(3)["abcijk"] +=         FAE[  "ce"]*T_4(3)["abeijk"];
       T_5(3)["abcijk"] -=         FMI[  "mk"]*T_4(3)["abcijm"];
       T_5(3)["abcijk"] +=   0.5*WABEF["abef"]*T_4(3)["efcijk"];
       T_5(3)["abcijk"] +=   0.5*WMNIJ["mnij"]*T_4(3)["abcmnk"];
       T_5(3)["abcijk"] +=       WAMEI["amei"]*T_4(3)["ebcjmk"];

    T_5.weight(D);

    U E6 = (1.0/1.0)*scalar(FME_2[  "me"]*T_3(1)[  "em"])+
           (1.0/1.0)*scalar(  FME[  "me"]*T_5(1)[  "em"])+
           (1.0/4.0)*scalar(WMNEF["mnef"]*T_5(2)["efmn"]);

    /***************************************************************************
     *
     * T^(6)
     *
     **************************************************************************/

    WTWAMIJ_2[  "amij"]  =       WMNEJ["nmej"]*T_2(2)[  "aein"];
    WTWAMIJ_2[  "amij"] +=   0.5*WAMEF["amef"]*T_2(2)[  "efij"];
    WTWAMIJ_2[  "amij"] -=       WMNIJ["nmij"]*T_2(1)[    "an"];
    WTWAMIJ_2[  "amij"] +=       WAMEI["amej"]*T_2(1)[    "ei"];
    WTWAMIJ_2[  "amij"] +=       FME_2[  "me"]*  T(2)[  "aeij"];
    WTWAMIJ_2[  "amij"] +=   0.5*WMNEF["mnef"]*T_2(3)["aefijn"];

    WTWABEJ_2[  "abej"]  =       WAMEF["amef"]*T_2(2)[  "fbmj"];
    WTWABEJ_2[  "abej"] +=   0.5*WMNEJ["mnej"]*T_2(2)[  "abmn"];
    WTWABEJ_2[  "abej"] +=       WABEF["abef"]*T_2(1)[    "fj"];
    WTWABEJ_2[  "abej"] -=       WAMEI["amej"]*T_2(1)[    "bm"];
    WTWABEJ_2[  "abej"] -=       FME_2[  "me"]*  T(2)[  "abmj"];
    WTWABEJ_2[  "abej"] -=   0.5*WMNEF["mnef"]*T_2(3)["afbmnj"];

      FTWAE_3[    "ae"]  =  -0.5*WMNEF["mnef"]*T_3(2)[  "afmn"];
      FTWAE_3[    "ae"] +=       WAMEF["amef"]*T_3(1)[    "fm"];

    WTWAMIJ_3[  "amij"]  =       FME_3[  "me"]*  T(2)[  "aeij"];
    WTWAMIJ_3[  "amij"] +=   0.5*WMNEF["mnef"]*T_3(3)["aefijn"];

    WTWABEJ_3[  "abej"]  = -     FME_3[  "me"]*  T(2)[  "abmj"];
    WTWABEJ_3[  "abej"] -=   0.5*WMNEF["mnef"]*T_3(3)["afbmnj"];

        FMI_4[    "mi"] +=         FME[  "me"]*T_4(1)[    "ei"];
        FAE_4[    "ae"] -=         FME[  "me"]*T_4(1)[    "am"];

      WAMEI_4[  "amei"]  =       WMNEF["mnef"]*T_4(2)[  "afni"];
      WAMEI_4[  "amei"] +=       WAMEF["amef"]*T_4(1)[    "fi"];
      WAMEI_4[  "amei"] -=       WMNEJ["nmei"]*T_4(1)[    "an"];

      WMNIJ_4[  "mnij"]  =   0.5*WMNEF["mnef"]*T_4(2)[  "efij"];
      WMNIJ_4[  "mnij"] +=       WMNEJ["mnej"]*T_4(1)[    "ei"];

      WABEF_4[  "abef"]  =   0.5*WMNEF["mnef"]*T_4(2)[  "abmn"];
      WABEF_4[  "abef"] -=       WAMEF["amef"]*T_4(1)[    "bm"];

        FMI_5[    "mi"]  =   0.5*WMNEF["mnef"]*T_5(2)[  "efin"];
        FMI_5[    "mi"] +=       WMNEJ["nmfi"]*T_5(1)[    "fn"];

        FAE_5[    "ae"]  =  -0.5*WMNEF["mnef"]*T_5(2)[  "afmn"];
        FAE_5[    "ae"] +=       WAMEF["amef"]*T_5(1)[    "fm"];

       T_6(1)[    "ai"]  =       FAE_3[  "ae"]*T_2(1)[    "ei"];
       T_6(1)[    "ai"] -=       FMI_3[  "mi"]*T_2(1)[    "am"];
       T_6(1)[    "ai"] +=       FME_3[  "me"]*T_2(2)[  "aeim"];
       T_6(1)[    "ai"] +=       FAE_2[  "ae"]*T_3(1)[    "ei"];
       T_6(1)[    "ai"] -=       FMI_2[  "mi"]*T_3(1)[    "am"];
       T_6(1)[    "ai"] +=       FME_2[  "me"]*T_3(2)[  "aeim"];
       T_6(1)[    "ai"] +=         FAE[  "ae"]*T_5(1)[    "ei"];
       T_6(1)[    "ai"] -=         FMI[  "mi"]*T_5(1)[    "am"];
       T_6(1)[    "ai"] -=       WAMEI["amei"]*T_5(1)[    "em"];
       T_6(1)[    "ai"] +=         FME[  "me"]*T_5(2)[  "aeim"];
       T_6(1)[    "ai"] +=   0.5*WAMEF["amef"]*T_5(2)[  "efim"];
       T_6(1)[    "ai"] -=   0.5*WMNEJ["mnei"]*T_5(2)[  "eamn"];
       T_6(1)[    "ai"] +=  0.25*WMNEF["mnef"]*T_5(3)["aefimn"];

       T_6(2)[  "abij"]  =       FAE_5[  "af"]*  T(2)[  "fbij"];
       T_6(2)[  "abij"] -=       FMI_5[  "ni"]*  T(2)[  "abnj"];
       T_6(2)[  "abij"] +=       FME_4[  "me"]*T_1(3)["abeijm"];
       T_6(2)[  "abij"] +=   WTWABEJ_3["abej"]*T_2(1)[    "ei"];
       T_6(2)[  "abij"] -=   WTWAMIJ_3["amij"]*T_2(1)[    "bm"];
       T_6(2)[  "abij"] +=     FTWAE_3[  "af"]*T_2(2)[  "fbij"];
       T_6(2)[  "abij"] -=       FMI_3[  "ni"]*T_2(2)[  "abnj"];
       T_6(2)[  "abij"] +=       FME_3[  "me"]*T_2(3)["abeijm"];
       T_6(2)[  "abij"] +=   WTWABEJ_2["abej"]*T_3(1)[    "ei"];
       T_6(2)[  "abij"] -=   WTWAMIJ_2["amij"]*T_3(1)[    "bm"];
       T_6(2)[  "abij"] +=     FTWAE_2[  "af"]*T_3(2)[  "fbij"];
       T_6(2)[  "abij"] -=       FMI_2[  "ni"]*T_3(2)[  "abnj"];
       T_6(2)[  "abij"] += 0.5*WABEF_2["abef"]*T_3(2)[  "efij"];
       T_6(2)[  "abij"] += 0.5*WMNIJ_2["mnij"]*T_3(2)[  "abmn"];
       T_6(2)[  "abij"] +=     WAMEI_2["amei"]*T_3(2)[  "ebjm"];
       T_6(2)[  "abij"] +=       FME_2[  "me"]*T_3(3)["abeijm"];
       T_6(2)[  "abij"] +=     WABEJ_1["abej"]*T_4(1)[    "ei"];
       T_6(2)[  "abij"] -=     WAMIJ_1["amij"]*T_4(1)[    "bm"];
       T_6(2)[  "abij"] +=       WABEJ["abej"]*T_5(1)[    "ei"];
       T_6(2)[  "abij"] -=       WAMIJ["amij"]*T_5(1)[    "bm"];
       T_6(2)[  "abij"] +=         FAE[  "af"]*T_5(2)[  "fbij"];
       T_6(2)[  "abij"] -=         FMI[  "ni"]*T_5(2)[  "abnj"];
       T_6(2)[  "abij"] +=   0.5*WABEF["abef"]*T_5(2)[  "efij"];
       T_6(2)[  "abij"] +=   0.5*WMNIJ["mnij"]*T_5(2)[  "abmn"];
       T_6(2)[  "abij"] +=       WAMEI["amei"]*T_5(2)[  "ebjm"];
       T_6(2)[  "abij"] +=   0.5*WAMEF["bmef"]*T_5(3)["aefijm"];
       T_6(2)[  "abij"] -=   0.5*WMNEJ["mnej"]*T_5(3)["abeinm"];
       T_6(2)[  "abij"] +=         FME[  "me"]*T_5(3)["abeijm"];

    T_6.weight(D);

    U E7 = (1.0/1.0)*scalar(FME_2[  "me"]*T_4(1)[  "em"])+
           (1.0/2.0)*scalar(FME_3[  "me"]*T_3(1)[  "em"])+
           (1.0/1.0)*scalar(  FME[  "me"]*T_6(1)[  "em"])+
           (1.0/4.0)*scalar(WMNEF["mnef"]*T_6(2)["efmn"]);

    printf("\n");
    printf("E-CCSD(T-3): %18.15f\n", E3);
    printf("\n");
    printf("E-CCSD(T-4): %18.15f\n", E4);
    printf("\n");
    printf("E-CCSD(T-5): %18.15f\n", E5);
    printf("\n");
    printf("E-CCSD(T-6): %18.15f\n", E6);
    printf("\n");
    printf("E-CCSD(T-7): %18.15f\n", E7);
    printf("\n");

    this->put("E(3)", new U(E3));
    this->put("E(4)", new U(E4));
    this->put("E(5)", new U(E5));
    this->put("E(6)", new U(E6));
    this->put("E(7)", new U(E7));

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::E_CCSD_T_N);
REGISTER_TASK(aquarius::cc::E_CCSD_T_N<double>,"e-ccsd(t-n)");
