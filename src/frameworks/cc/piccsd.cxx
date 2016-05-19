#include "piccsd.hpp"

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
PiCCSD<U>::PiCCSD(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd(t)_l.S", "S"));
    reqs.push_back(Requirement("ccsd(t)_l.Upsilon", "Upsilon"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("ccsd.T", "T"));
    reqs.push_back(Requirement("ccsd.L", "L"));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("ccsd(t)_l.Pi", "Pi", reqs));
}

template <typename U>
bool PiCCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    this->put   (   "Pi", new DeexcitationOperator<U,2>("Pi", arena, occ, vrt));
    this->puttmp(    "D", new Denominator         <U  >(H));
    this->puttmp(    "Z", new DeexcitationOperator<U,2>("Z", arena, occ, vrt));
    this->puttmp("Omega", new DeexcitationOperator<U,2>("Omega", arena, occ, vrt));
    this->puttmp(  "GIM", new SpinorbitalTensor   <U  >("G(im)", H.getIJ()));
    this->puttmp(  "GEA", new SpinorbitalTensor   <U  >("G(ea)", H.getAB()));
    this->puttmp(  "XME", new SpinorbitalTensor   <U  >("X(me)", H.getIA()));
    this->puttmp(  "XEA", new SpinorbitalTensor   <U  >("X(ea)", H.getAB()));
    this->puttmp(  "XIM", new SpinorbitalTensor   <U  >("X(im)", H.getIJ()));
    this->puttmp("XEIAM", new SpinorbitalTensor   <U  >("X(ei,am)", H.getAIBJ()));
    this->puttmp("XIJMN", new SpinorbitalTensor   <U  >("X(ij,mn)", H.getIJKL()));
    this->puttmp("XEIMN", new SpinorbitalTensor   <U  >("X(ei,mn)", H.getAIJK()));
    this->puttmp("XEFAM", new SpinorbitalTensor   <U  >("X(ef,am)", H.getABCI()));
    this->puttmp(  "DIA", new SpinorbitalTensor   <U  >("D(ia)", H.getIA()));
    this->puttmp(  "DAI", new SpinorbitalTensor   <U  >("D(ai)", H.getAI()));
    this->puttmp(  "DAB", new SpinorbitalTensor   <U  >("D(ab)", H.getAB()));
    this->puttmp(  "DIJ", new SpinorbitalTensor   <U  >("D(ij)", H.getIJ()));
    this->puttmp("GIJKL", new SpinorbitalTensor   <U  >("G(ij,kl)", H.getIJKL()));
    this->puttmp("GIJAK", new SpinorbitalTensor   <U  >("G(ij,ak)", H.getIJAK()));
    this->puttmp("GAIBC", new SpinorbitalTensor   <U  >("G(ai,bc)", H.getAIBC()));

    auto& T   = this->template get   <ExcitationOperator  <U,2>>("T");
    auto& Ups = this->template get   <ExcitationOperator  <U,2>>("Upsilon");
    auto& L   = this->template get   <DeexcitationOperator<U,2>>("L");
    auto& S   = this->template get   <DeexcitationOperator<U,2>>("S");
    auto& P   = this->template get   <DeexcitationOperator<U,2>>("Pi");
    auto& O   = this->template gettmp<DeexcitationOperator<U,2>>("Omega");
    auto& Z   = this->template gettmp<DeexcitationOperator<U,2>>("Z");
    auto& D   = this->template gettmp<Denominator         <U  >>("D");

    auto&   FME =   H.getIA();
    auto&   FAE =   H.getAB();
    auto&   FMI =   H.getIJ();
    auto& WMNEF = H.getIJAB();
    auto& WAMEF = H.getAIBC();
    auto& WABEJ = H.getABCI();
    auto& WABEF = H.getABCD();
    auto& WMNIJ = H.getIJKL();
    auto& WMNEJ = H.getIJAK();
    auto& WAMIJ = H.getAIJK();
    auto& WAMEI = H.getAIBJ();

    auto&   XME = this->template gettmp<SpinorbitalTensor<U>>(  "XME");
    auto&   XEA = this->template gettmp<SpinorbitalTensor<U>>(  "XEA");
    auto&   XIM = this->template gettmp<SpinorbitalTensor<U>>(  "XIM");
    auto& XEIAM = this->template gettmp<SpinorbitalTensor<U>>("XEIAM");
    auto& XIJMN = this->template gettmp<SpinorbitalTensor<U>>("XIJMN");
    auto& XEIMN = this->template gettmp<SpinorbitalTensor<U>>("XEIMN");
    auto& XEFAM = this->template gettmp<SpinorbitalTensor<U>>("XEFAM");
    auto&   DIA = this->template gettmp<SpinorbitalTensor<U>>(  "DIA");
    auto&   DAI = this->template gettmp<SpinorbitalTensor<U>>(  "DAI");
    auto&   DAB = this->template gettmp<SpinorbitalTensor<U>>(  "DAB");
    auto&   DIJ = this->template gettmp<SpinorbitalTensor<U>>(  "DIJ");
    auto& GIJKL = this->template gettmp<SpinorbitalTensor<U>>("GIJKL");
    auto& GIJAK = this->template gettmp<SpinorbitalTensor<U>>("GIJAK");
    auto& GAIBC = this->template gettmp<SpinorbitalTensor<U>>("GAIBC");

    Z(0) = 0;
    P(0) = 0;
    P(1) = 0;
    P(2) = 0;

    U LUps =      scalar(L(1)[  "ia"]*Ups(1)[  "ai"]) +
             0.25*scalar(L(2)["ijab"]*Ups(2)["abij"]);

    print("LU: %.15f\n", real(LUps));

      XME[  "me"]  =       WMNEF["mnef"]*Ups(1)[  "fn"];

      XEA[  "ea"]  = -       FME[  "ma"]*Ups(1)[  "em"];
      XEA[  "ea"] +=       WAMEF["emaf"]*Ups(1)[  "fm"];
      XEA[  "ea"] -=   0.5*WMNEF["mnaf"]*Ups(2)["efmn"];

      XIM[  "im"]  =         FME[  "ie"]*Ups(1)[  "em"];
      XIM[  "im"] +=       WMNEJ["niem"]*Ups(1)[  "en"];
      XIM[  "im"] +=   0.5*WMNEF["inef"]*Ups(2)["efmn"];

    XEIAM["eiam"]  = -     WMNEF["inaf"]*Ups(2)["efmn"];

    XIJMN["ijmn"]  =   0.5*WMNEF["ijef"]*Ups(2)["efmn"];

    XEIMN["eimn"]  =       WMNEJ["oifn"]*Ups(2)["efmo"];
    XEIMN["eimn"] +=   0.5*WAMEF["eifg"]*Ups(2)["fgmn"];
    XEIMN["eimn"] +=         XME[  "if"]*  T(2)["efmn"];

    XEFAM["efam"]  =       WAMEF["enag"]*Ups(2)["gfnm"];
    XEFAM["efam"] +=   0.5*WMNEJ["noam"]*Ups(2)["efno"];
    XEFAM["efam"] -=         XME[  "na"]*  T(2)["efnm"];

      DIA[  "me"]  =        L(2)["mnef"]*Ups(1)[  "fn"];

      DAB[  "ea"]  = - 0.5* L(2)["mnaf"]*Ups(2)["efmn"];
      DIJ[  "im"]  =   0.5* L(2)["inef"]*Ups(2)["efmn"];

    GIJKL["ijmn"]  =   0.5* L(2)["ijef"]*Ups(2)["efmn"];

    GIJAK["ijak"]  =        L(2)["ijae"]*Ups(1)[  "ek"];

    GAIBC["aibc"]  = -      L(2)["mibc"]*Ups(1)[  "am"];

      DAI[  "ai"]  = - 0.5* T(2)["eamn"]* GIJAK["mnei"];
      DAI[  "ai"] +=   0.5* T(2)["efim"]* GAIBC["amef"];

     O(1)[  "ia"]  =         XME[  "ia"];
     O(1)[  "ia"] +=  LUps*  FME[  "ia"];
     O(1)[  "ia"] +=       WMNEF["miea"]*   DAI[  "em"];
     O(1)[  "ia"] +=         FAE[  "ea"]*   DIA[  "ie"];
     O(1)[  "ia"] -=         FMI[  "im"]*   DIA[  "ma"];
     O(1)[  "ia"] +=         XEA[  "ea"]*  L(1)[  "ie"];
     O(1)[  "ia"] -=         XIM[  "im"]*  L(1)[  "ma"];
     O(1)[  "ia"] +=         FME[  "ie"]*   DAB[  "ea"];
     O(1)[  "ia"] -=         FME[  "ma"]*   DIJ[  "im"];
     O(1)[  "ia"] -=       XEIAM["eiam"]*  L(1)[  "me"];
     O(1)[  "ia"] -=       WAMEI["eiam"]*   DIA[  "me"];
     O(1)[  "ia"] +=   0.5*XEFAM["efam"]*  L(2)["imef"];
     O(1)[  "ia"] -=   0.5*XEIMN["eimn"]*  L(2)["mnea"];
     O(1)[  "ia"] -=   0.5*WABEF["efga"]* GAIBC["gief"];
     O(1)[  "ia"] +=       WAMEI["eifm"]* GAIBC["fmea"];
     O(1)[  "ia"] -=       WAMEI["eman"]* GIJAK["inem"];
     O(1)[  "ia"] +=   0.5*WMNIJ["imno"]* GIJAK["noam"];

      DAB[  "ea"] -=        L(1)[  "ma"]*Ups(1)[  "em"];
      DIJ[  "im"] +=        L(1)[  "ie"]*Ups(1)[  "em"];

     O(1)[  "ia"] -=       WAMEF["fiea"]*   DAB[  "ef"];
     O(1)[  "ia"] -=       WMNEJ["inam"]*   DIJ[  "mn"];

     O(2)["ijab"]  =        S(2)["ijab"];
     O(2)["ijab"] +=  LUps*WMNEF["ijab"];
     O(2)["ijab"] +=         XME[  "ia"]*  L(1)[  "jb"];
     O(2)["ijab"] +=         FME[  "ia"]*   DIA[  "jb"];
     O(2)["ijab"] +=         XEA[  "ea"]*  L(2)["ijeb"];
     O(2)["ijab"] -=         XIM[  "im"]*  L(2)["mjab"];
     O(2)["ijab"] +=       WAMEF["ejab"]*   DIA[  "ie"];
     O(2)["ijab"] -=       WMNEJ["ijam"]*   DIA[  "mb"];
     O(2)["ijab"] +=       WMNEF["ijeb"]*   DAB[  "ea"];
     O(2)["ijab"] -=       WMNEF["mjab"]*   DIJ[  "im"];
     O(2)["ijab"] -=       XEIAM["eiam"]*  L(2)["mjeb"];
     O(2)["ijab"] +=   0.5*XIJMN["ijmn"]*  L(2)["mnab"];
     O(2)["ijab"] +=   0.5*WMNEF["mnab"]* GIJKL["ijmn"];
     O(2)["ijab"] -=       WAMEF["fiae"]* GAIBC["ejbf"];
     O(2)["ijab"] -=       WMNEJ["ijem"]* GAIBC["emab"];
     O(2)["ijab"] -=       WAMEF["emab"]* GIJAK["ijem"];
     O(2)["ijab"] -=       WMNEJ["niam"]* GIJAK["mjbn"];

    print("DIA: %.15f\n", DIA({0,0},{0,0}).norm(2));
    print("FAE: %.15f\n", FAE({0,0},{0,0}).norm(2));
    print("FMI: %.15f\n", FMI({0,0},{0,0}).norm(2));
    print("WMAEI: %.15f\n", WAMEI({1,0},{0,1}).norm(2));

    print("Omega1: %.15f\n", O(1)({0,0},{0,0}).norm(2));
    print("Omega2: %.15f\n", O(2)({0,1},{1,0}).norm(2));

    Iterative<U>::run(dag, arena);

    print("Pi1: %.15f\n", P(1)({0,0},{0,0}).norm(2));
    print("Pi2: %.15f\n", P(2)({0,1},{1,0}).norm(2));

    this->put("convergence", new U(this->conv()));

    return true;
}

template <typename U>
void PiCCSD<U>::iterate(const Arena& arena)
{
    const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    auto&   FME =   H.getIA();
    auto&   FAE =   H.getAB();
    auto&   FMI =   H.getIJ();
    auto& WMNEF = H.getIJAB();
    auto& WAMEF = H.getAIBC();
    auto& WABEJ = H.getABCI();
    auto& WABEF = H.getABCD();
    auto& WMNIJ = H.getIJKL();
    auto& WMNEJ = H.getIJAK();
    auto& WAMIJ = H.getAIJK();
    auto& WAMEI = H.getAIBJ();

    auto& T = this->template get   <ExcitationOperator  <U,2>>("T");
    auto& P = this->template get   <DeexcitationOperator<U,2>>("Pi");
    auto& O = this->template gettmp<DeexcitationOperator<U,2>>("Omega");
    auto& D = this->template gettmp<Denominator         <U  >>("D");
    auto& Z = this->template gettmp<DeexcitationOperator<U,2>>("Z");

    auto& GIM = this->template gettmp<SpinorbitalTensor<U>>("GIM");
    auto& GEA = this->template gettmp<SpinorbitalTensor<U>>("GEA");

    /***************************************************************************
     *
     * Intermediates for Pi-CCSD
     */
    GIM["mn"] =  0.5*T(2)["efno"]*P(2)["moef"];
    GEA["ef"] = -0.5*T(2)["egmn"]*P(2)["mnfg"];
    /*
     **************************************************************************/

    /***************************************************************************
     *
     * Pi-CCSD iteration
     */
    Z(1)[  "ia"]  =      O(1)[  "ia"];
    Z(1)[  "ia"] +=       FAE[  "ea"]*P(1)[  "ie"];
    Z(1)[  "ia"] -=       FMI[  "im"]*P(1)[  "ma"];
    Z(1)[  "ia"] -=     WAMEI["eiam"]*P(1)[  "me"];
    Z(1)[  "ia"] += 0.5*WABEJ["efam"]*P(2)["imef"];
    Z(1)[  "ia"] -= 0.5*WAMIJ["eimn"]*P(2)["mnea"];
    Z(1)[  "ia"] -=     WMNEJ["inam"]* GIM[  "mn"];
    Z(1)[  "ia"] -=     WAMEF["fiea"]* GEA[  "ef"];

    Z(2)["ijab"]  =      O(2)["ijab"];
    Z(2)["ijab"] +=       FME[  "ia"]*P(1)[  "jb"];
    Z(2)["ijab"] +=     WAMEF["ejab"]*P(1)[  "ie"];
    Z(2)["ijab"] -=     WMNEJ["ijam"]*P(1)[  "mb"];
    Z(2)["ijab"] +=       FAE[  "ea"]*P(2)["ijeb"];
    Z(2)["ijab"] -=       FMI[  "im"]*P(2)["mjab"];
    Z(2)["ijab"] += 0.5*WABEF["efab"]*P(2)["ijef"];
    Z(2)["ijab"] += 0.5*WMNIJ["ijmn"]*P(2)["mnab"];
    Z(2)["ijab"] +=     WAMEI["eiam"]*P(2)["mjbe"];
    Z(2)["ijab"] -=     WMNEF["mjab"]* GIM[  "im"];
    Z(2)["ijab"] +=     WMNEF["ijeb"]* GEA[  "ea"];
    /*
     **************************************************************************/

    Z.weight(D);
    P += Z;

    this->conv() = Z.norm(00);

    diis.extrapolate(P, Z);
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::PiCCSD);
REGISTER_TASK(aquarius::cc::PiCCSD<double>,"piccsd",spec);
