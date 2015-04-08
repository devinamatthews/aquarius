#include "lccd.hpp"

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
LCCD<U>::LCCD(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("double", "S2", reqs));
    this->addProduct(Product("double", "multiplicity", reqs));
    this->addProduct(Product("lccd.T", "T", reqs));
    this->addProduct(Product("lccd.Hbar", "Hbar", reqs));
}

template <typename U>
bool LCCD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& H = this->template get<TwoElectronOperator<U>>("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    auto& T = this->put   (  "T", new ExcitationOperator<U,2>("T", arena, occ, vrt));
    auto& Z = this->puttmp(  "Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));
    auto& D = this->puttmp(  "D", new Denominator       <U  >(H));

    Z(0) = 0;
    T(0) = 0;
    T(1) = 0;
    T(2) = H.getABIJ();

    T.weight(D);

    CTF_Timer_epoch ep(this->name.c_str());
    ep.begin();
    Iterative<U>::run(dag, arena);
    ep.end();

    this->put("energy", new U(this->energy()));
    this->put("convergence", new U(this->conv()));

    /*
    if (isUsed("S2") || isUsed("multiplicity"))
    {
        double s2 = this->template getProjectedS2(occ, vrt, T(1), T(2));
        double mult = sqrt(4*s2+1);

        this->put("S2", new Scalar(arena, s2));
        this->put("multiplicity", new Scalar(arena, mult));
    }
    */

    if (this->isUsed("Hbar"))
    {
        this->put("Hbar", new STTwoElectronOperator<U>("Hbar", H, T, true));
    }

    return true;
}

template <typename U>
void LCCD<U>::iterate(const Arena& arena)
{
    const auto& H = this->template get<TwoElectronOperator<U>>("H");

    const SpinorbitalTensor<U>&   fAE =   H.getAB();
    const SpinorbitalTensor<U>&   fMI =   H.getIJ();
    const SpinorbitalTensor<U>& VABIJ = H.getABIJ();
    const SpinorbitalTensor<U>& VABEF = H.getABCD();
    const SpinorbitalTensor<U>& VMNIJ = H.getIJKL();
    const SpinorbitalTensor<U>& VAMEI = H.getAIBJ();

    auto& T = this->template get   <ExcitationOperator<U,2>>("T");
    auto& D = this->template gettmp<Denominator       <U  >>("D");
    auto& Z = this->template gettmp<ExcitationOperator<U,2>>("Z");

    /**************************************************************************
     *
     * LCCD Iteration
     */
    Z(2)["abij"]  =     VABIJ["abij"];
    Z(2)["abij"] +=       fAE[  "af"]*T(2)["fbij"];
    Z(2)["abij"] -=       fMI[  "ni"]*T(2)["abnj"];
    Z(2)["abij"] += 0.5*VABEF["abef"]*T(2)["efij"];
    Z(2)["abij"] += 0.5*VMNIJ["mnij"]*T(2)["abmn"];
    Z(2)["abij"] +=     VAMEI["amei"]*T(2)["ebjm"];
    /*
     *************************************************************************/

    Z.weight(D);
    T += Z;

    this->energy() = 0.25*real(scalar(H.getABIJ()*T(2)));
    this->conv() = Z.norm(00);

    diis.extrapolate(T, Z);
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::LCCD);
REGISTER_TASK(aquarius::cc::LCCD<double>,"lccd",spec);
