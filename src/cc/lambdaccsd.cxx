#include "lambdaccsd.hpp"

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
LambdaCCSD<U>::LambdaCCSD(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("ccsd.T", "T"));
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("ccsd.L", "L", reqs));
}

template <typename U>
bool LambdaCCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const STTwoElectronOperator<U>& H = this->template get<STTwoElectronOperator<U> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    this->put   (  "L", new DeexcitationOperator<U,2>("L", arena, occ, vrt));
    this->puttmp(  "D", new Denominator         <U  >(H));
    this->puttmp(  "Z", new DeexcitationOperator<U,2>("Z", arena, occ, vrt));
    this->puttmp("GIM", new SpinorbitalTensor   <U  >("G(im)", H.getIJ()));
    this->puttmp("GEA", new SpinorbitalTensor   <U  >("G(ea)", H.getAB()));

    ExcitationOperator  <U,2>& T = this->template get   <ExcitationOperator  <U,2>>("T");
    DeexcitationOperator<U,2>& L = this->template get   <DeexcitationOperator<U,2>>("L");
    DeexcitationOperator<U,2>& Z = this->template gettmp<DeexcitationOperator<U,2>>("Z");
    Denominator         <U  >& D = this->template gettmp<Denominator         <U  >>("D");

    Z(0) = 0;
    L(0) = 1;
    L(1)[  "ia"] = T(1)[  "ai"];
    L(2)["ijab"] = T(2)["abij"];

    Iterative<U>::run(dag, arena);

    this->put("energy", new U(this->energy()));
    this->put("convergence", new U(this->conv()));

    return true;
}

template <typename U>
void LambdaCCSD<U>::iterate(const Arena& arena)
{
    const STTwoElectronOperator<U>& H = this->template get<STTwoElectronOperator<U> >("Hbar");

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

    ExcitationOperator  <U,2>& T = this->template get   <ExcitationOperator  <U,2>>("T");
    DeexcitationOperator<U,2>& L = this->template get   <DeexcitationOperator<U,2>>("L");
    Denominator         <U  >& D = this->template gettmp<Denominator         <U  >>("D");
    DeexcitationOperator<U,2>& Z = this->template gettmp<DeexcitationOperator<U,2>>("Z");

    SpinorbitalTensor<U>& GIM = this->template gettmp<SpinorbitalTensor<U>>("GIM");
    SpinorbitalTensor<U>& GEA = this->template gettmp<SpinorbitalTensor<U>>("GEA");

    /***************************************************************************
     *
     * Intermediates for Lambda-CCSD
     */
    GIM["mn"]  =  0.5*T(2)["efno"]*L(2)["moef"];
    GEA["ef"]  = -0.5*T(2)["egmn"]*L(2)["mnfg"];
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
    Z(1)[  "ia"] -=     WMNEJ["inam"]* GIM[  "mn"];
    Z(1)[  "ia"] -=     WAMEF["fiea"]* GEA[  "ef"];

    Z(2)["ijab"]  =     WMNEF["ijab"];
    Z(2)["ijab"] +=       FME[  "ia"]*L(1)[  "jb"];
    Z(2)["ijab"] +=     WAMEF["ejab"]*L(1)[  "ie"];
    Z(2)["ijab"] -=     WMNEJ["ijam"]*L(1)[  "mb"];
    Z(2)["ijab"] +=       FAE[  "ea"]*L(2)["ijeb"];
    Z(2)["ijab"] -=       FMI[  "im"]*L(2)["mjab"];
    Z(2)["ijab"] += 0.5*WABEF["efab"]*L(2)["ijef"];
    Z(2)["ijab"] += 0.5*WMNIJ["ijmn"]*L(2)["mnab"];
    Z(2)["ijab"] +=     WAMEI["eiam"]*L(2)["mjbe"];
    Z(2)["ijab"] -=     WMNEF["mjab"]* GIM[  "im"];
    Z(2)["ijab"] +=     WMNEF["ijeb"]* GEA[  "ea"];
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::LambdaCCSD);
REGISTER_TASK(aquarius::cc::LambdaCCSD<double>,"lambdaccsd",spec);
