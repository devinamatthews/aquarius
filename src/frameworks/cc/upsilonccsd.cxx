#include "upsilonccsd.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;

namespace aquarius
{
namespace cc
{

template <typename Type>
UpsilonCCSD<Type>::UpsilonCCSD(const string& name, Config& config)
: Iterative<Type>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd(t)_l.U", "U"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("ccsd.T", "T"));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("ccsd(t)_l.Upsilon", "Upsilon", reqs));
}

template <typename Type>
bool UpsilonCCSD<Type>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& H = this->template get<TwoElectronOperator<Type>>("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    auto& Ups = this->put   ("Upsilon", new ExcitationOperator<Type,2>("Upsilon", arena, occ, vrt));
    auto& Z   = this->puttmp(      "Z", new ExcitationOperator<Type,2>("Z", arena, occ, vrt));
    auto& D   = this->puttmp(      "D", new Denominator       <Type  >(H));

    this->puttmp("XAE", new SpinorbitalTensor<Type>("X(ae)", H.getAB()));
    this->puttmp("XMI", new SpinorbitalTensor<Type>("X(mi)", H.getIJ()));

      Z(0) = 0;
    Ups(0) = 0;
    Ups(1) = 0;
    Ups(2) = 0;

    Iterative<Type>::run(dag, arena);

    print("Upsilon1: %.15f\n", Ups(1)({0,0},{0,0}).norm(2));
    print("Upsilon2: %.15f\n", Ups(2)({1,0},{0,1}).norm(2));

    this->put("convergence", new Type(this->conv()));

    return true;
}

template <typename Type>
void UpsilonCCSD<Type>::iterate(const Arena& arena)
{
    const auto& H = this->template get<TwoElectronOperator<Type>>("Hbar");

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

    auto& T   = this->template get   <ExcitationOperator<Type,2>>("T");
    auto& Ups = this->template get   <ExcitationOperator<Type,2>>("Upsilon");
    auto& U   = this->template get   <ExcitationOperator<Type,2>>("U");
    auto& D   = this->template gettmp<Denominator       <Type  >>("D");
    auto& Z   = this->template gettmp<ExcitationOperator<Type,2>>("Z");

    auto& XAE = this->template gettmp<SpinorbitalTensor<Type>>("XAE");
    auto& XMI = this->template gettmp<SpinorbitalTensor<Type>>("XMI");

    /**************************************************************************
     *
     * Intermediates for Upsilon-CCSD
     */
     XMI[  "mi"]  =  0.5*WMNEF["mnef"]*Ups(2)["efin"];
     XMI[  "mi"] +=      WMNEJ["nmfi"]*Ups(1)[  "fn"];

     XAE[  "ae"]  = -0.5*WMNEF["mnef"]*Ups(2)["afmn"];
     XAE[  "ae"] +=      WAMEF["amef"]*Ups(1)[  "fm"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * Upsilon-CCSD Iteration
     */
    Z(1)[  "ai"]  =       U(1)[  "ai"];
    Z(1)[  "ai"] +=        FAE[  "ae"]*Ups(1)[  "ei"];
    Z(1)[  "ai"] -=        FMI[  "mi"]*Ups(1)[  "am"];
    Z(1)[  "ai"] -=      WAMEI["amei"]*Ups(1)[  "em"];
    Z(1)[  "ai"] +=        FME[  "me"]*Ups(2)["aeim"];
    Z(1)[  "ai"] +=  0.5*WAMEF["amef"]*Ups(2)["efim"];
    Z(1)[  "ai"] -=  0.5*WMNEJ["mnei"]*Ups(2)["eamn"];

    Z(2)["abij"]  =       U(2)["abij"];
    Z(2)["abij"] +=      WABEJ["abej"]*Ups(1)[  "ei"];
    Z(2)["abij"] -=      WAMIJ["amij"]*Ups(1)[  "bm"];
    Z(2)["abij"] +=        FAE[  "ae"]*Ups(2)["ebij"];
    Z(2)["abij"] -=        FMI[  "mi"]*Ups(2)["abmj"];
    Z(2)["abij"] +=        XAE[  "ae"]*  T(2)["ebij"];
    Z(2)["abij"] -=        XMI[  "mi"]*  T(2)["abmj"];
    Z(2)["abij"] +=  0.5*WABEF["abef"]*Ups(2)["efij"];
    Z(2)["abij"] +=  0.5*WMNIJ["mnij"]*Ups(2)["abmn"];
    Z(2)["abij"] +=      WAMEI["amei"]*Ups(2)["ebjm"];
    /*
     *************************************************************************/

    Z.weight(D);
    Ups += Z;

    this->conv() = Z.norm(00);

    diis.extrapolate(Ups, Z);
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::UpsilonCCSD);
REGISTER_TASK(aquarius::cc::UpsilonCCSD<double>,"upsilonccsd",spec);
