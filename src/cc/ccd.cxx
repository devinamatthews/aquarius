#include "ccd.hpp"

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
CCD<U>::CCD(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("double", "S2", reqs));
    this->addProduct(Product("double", "multiplicity", reqs));
    this->addProduct(Product("ccd.T", "T", reqs));
    this->addProduct(Product("ccd.Hbar", "Hbar", reqs));
}

template <typename U>
bool CCD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U> >("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    this->put   (  "T", new ExcitationOperator<U,2>("T", arena, occ, vrt));
    this->puttmp(  "Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));
    this->puttmp(  "D", new Denominator<U>(H));

    this->puttmp(  "FAE", new SpinorbitalTensor<U>(    "F(ae)",   H.getAB()));
    this->puttmp(  "FMI", new SpinorbitalTensor<U>(    "F(mi)",   H.getIJ()));
    this->puttmp("WMNIJ", new SpinorbitalTensor<U>( "W(mn,ij)", H.getIJKL()));
    this->puttmp("WAMEI", new SpinorbitalTensor<U>("W~(am,ei)", H.getAIBJ()));

    ExcitationOperator<U,2>& T = this->template get   <ExcitationOperator<U,2> >(  "T");
    Denominator<U>&          D = this->template gettmp<Denominator<U> >         (  "D");
    ExcitationOperator<U,2>& Z = this->template gettmp<ExcitationOperator<U,2> >(  "Z");

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
void CCD<U>::iterate(const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U> >("H");

    const SpinorbitalTensor<U>&   fAE =   H.getAB();
    const SpinorbitalTensor<U>&   fMI =   H.getIJ();
    const SpinorbitalTensor<U>& VABIJ = H.getABIJ();
    const SpinorbitalTensor<U>& VMNEF = H.getIJAB();
    const SpinorbitalTensor<U>& VABEF = H.getABCD();
    const SpinorbitalTensor<U>& VMNIJ = H.getIJKL();
    const SpinorbitalTensor<U>& VAMEI = H.getAIBJ();

    ExcitationOperator<U,2>& T = this->template get   <ExcitationOperator<U,2> >(  "T");
    Denominator<U>&          D = this->template gettmp<Denominator<U>          >(  "D");
    ExcitationOperator<U,2>& Z = this->template gettmp<ExcitationOperator<U,2> >(  "Z");

    SpinorbitalTensor<U>&   FAE = this->template gettmp<SpinorbitalTensor<U> >(  "FAE");
    SpinorbitalTensor<U>&   FMI = this->template gettmp<SpinorbitalTensor<U> >(  "FMI");
    SpinorbitalTensor<U>& WMNIJ = this->template gettmp<SpinorbitalTensor<U> >("WMNIJ");
    SpinorbitalTensor<U>& WAMEI = this->template gettmp<SpinorbitalTensor<U> >("WAMEI");

    /**************************************************************************
     *
     * Intermediates for CCD
     */
      FMI[  "mi"]  =       fMI[  "mi"];
      FMI[  "mi"] += 0.5*VMNEF["mnef"]*T(2)["efin"];

      FAE[  "ae"]  =       fAE[  "ae"];
      FAE[  "ae"] -= 0.5*VMNEF["mnef"]*T(2)["afmn"];

    WMNIJ["mnij"]  =     VMNIJ["mnij"];
    WMNIJ["mnij"] += 0.5*VMNEF["mnef"]*T(2)["efij"];

    WAMEI["amei"]  =     VAMEI["amei"];
    WAMEI["amei"] += 0.5*VMNEF["mnef"]*T(2)["afni"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * CCD Iteration
     */
    Z(2)["abij"]  =     VABIJ["abij"];
    Z(2)["abij"] +=       FAE[  "af"]*T(2)["fbij"];
    Z(2)["abij"] -=       FMI[  "ni"]*T(2)["abnj"];
    Z(2)["abij"] += 0.5*VABEF["abef"]*T(2)["efij"];
    Z(2)["abij"] += 0.5*WMNIJ["mnij"]*T(2)["abmn"];
    Z(2)["abij"] +=     WAMEI["amei"]*T(2)["ebjm"];
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCD);
REGISTER_TASK(aquarius::cc::CCD<double>,"ccd",spec);
