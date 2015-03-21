#include "perturbedlambdaccsd.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;

namespace aquarius
{
namespace cc
{

template <typename U>
PerturbedLambdaCCSD<U>::PerturbedLambdaCCSD(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.L", "L"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("1epert", "A"));
    reqs.push_back(Requirement("double", "omega"));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("ccsd.LA", "LA", reqs));
}

template <typename U>
bool PerturbedLambdaCCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const PerturbedSTTwoElectronOperator<U>& A = this->template get<PerturbedSTTwoElectronOperator<U> >("A");
    const STTwoElectronOperator<U>& H = this->template get<STTwoElectronOperator<U> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    this->put   ("LA", new DeexcitationOperator<U,2>("L^A", arena, occ, vrt));
    this->puttmp( "D", new DeexcitationOperator<U,2>("D", arena, occ, vrt));
    this->puttmp( "N", new DeexcitationOperator<U,2>("N", arena, occ, vrt));
    this->puttmp( "Z", new DeexcitationOperator<U,2>("Z", arena, occ, vrt));

    U omega = this->template get<U>("omega");

    DeexcitationOperator<U,2>&  L = this->template get   <DeexcitationOperator<U,2> >( "L");
    DeexcitationOperator<U,2>& LA = this->template get   <DeexcitationOperator<U,2> >("LA");
    DeexcitationOperator<U,2>&  D = this->template gettmp<DeexcitationOperator<U,2> >( "D");
    DeexcitationOperator<U,2>&  N = this->template gettmp<DeexcitationOperator<U,2> >( "N");

    D(0) = (U)1.0;
    D(1)["ia"]  = H.getIJ()["ii"];
    D(1)["ia"] -= H.getAB()["aa"];
    D(2)["ijab"]  = H.getIJ()["ii"];
    D(2)["ijab"] += H.getIJ()["jj"];
    D(2)["ijab"] -= H.getAB()["aa"];
    D(2)["ijab"] -= H.getAB()["bb"];

    D -= omega;
    D = 1/D;

    N(1) = (U)0.0;
    N(2) = (U)0.0;
    //TODO: A.contract(L, N);
    N(0) = (U)0.0;

    LA = N*D;

    Iterative<U>::run(dag, arena);

    this->put("convergence", new U(this->conv()));

    return true;
}

template <typename U>
void PerturbedLambdaCCSD<U>::iterate(const Arena& arena)
{
    const STTwoElectronOperator<U>& H = this->template get<STTwoElectronOperator<U> >("Hbar");

    DeexcitationOperator<U,2>& LA = this->template get   <DeexcitationOperator<U,2> >("LA");
    DeexcitationOperator<U,2>&  D = this->template gettmp<DeexcitationOperator<U,2> >( "D");
    DeexcitationOperator<U,2>&  N = this->template gettmp<DeexcitationOperator<U,2> >( "N");
    DeexcitationOperator<U,2>&  Z = this->template gettmp<DeexcitationOperator<U,2> >( "Z");

    Z = N;
    //TODO: H.contract(LA, Z);

     Z *= D;
    LA += Z;

    this->conv() = Z.norm(00);

    diis.extrapolate(LA, Z);
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::PerturbedLambdaCCSD);
REGISTER_TASK(aquarius::cc::PerturbedLambdaCCSD<double>, "perturbedlambdaccsd");
