#include "perturbedccsd.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;

namespace aquarius
{
namespace cc
{

template <typename U>
PerturbedCCSD<U>::PerturbedCCSD(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.T", "T"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("1epert", "A"));
    reqs.push_back(Requirement("double", "omega"));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("ccsd.TA", "TA", reqs));
}

template <typename U>
bool PerturbedCCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& A = this->template get<OneElectronOperator  <U>>("A");
    const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    auto& T = this->template get<ExcitationOperator<U,2>>("T");

    auto& TA = this->put   ("TA", new ExcitationOperator  <U,2>("T^A", arena, occ, vrt));
    //auto& X  = this->puttmp( "X", new STExcitationOperator<U,2>("X", A, T));
    auto& Z  = this->puttmp( "Z", new ExcitationOperator  <U,2>("Z", arena, occ, vrt));
    auto& D  = this->puttmp( "D", new Denominator         <U  >(H));

    omega = this->template get<U>("omega");

    //X(0) = 0;
    TA.weight(D, omega);

    return true;
}

template <typename U>
void PerturbedCCSD<U>::iterate(const Arena& arena)
{
    const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    auto& TA = this->template get   <ExcitationOperator<U,2>>("TA");
    auto& D  = this->template gettmp<Denominator       <U  >>( "D");
    auto& X  = this->template gettmp<ExcitationOperator<U,2>>( "X");
    auto& Z  = this->template gettmp<ExcitationOperator<U,2>>( "Z");

    Z = X;
    //TODO: H.contract(TA, Z);

    Z.weight(D, omega);
    TA += Z;

    this->conv() = Z.norm(00);

    diis.extrapolate(TA, Z);
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::PerturbedCCSD);
REGISTER_TASK(aquarius::cc::PerturbedCCSD<double>, "perturbedccsd");
