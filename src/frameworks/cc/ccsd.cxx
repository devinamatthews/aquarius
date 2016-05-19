#include "ccsd.hpp"

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
CCSD<U>::CCSD(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    this->addProduct(Product("double", "mp2", reqs));
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("double", "S2", reqs));
    this->addProduct(Product("double", "multiplicity", reqs));
    this->addProduct(Product("ccsd.T", "T", reqs));
    this->addProduct(Product("ccsd.Hbar", "Hbar", reqs));
}

template <typename U>
bool CCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& H = this->template get<TwoElectronOperator<U>>("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    auto& T   = this->put   (  "T", new ExcitationOperator<U,2>("T", arena, occ, vrt));
    auto& Z   = this->puttmp(  "Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));
    auto& Tau = this->puttmp("Tau", new SpinorbitalTensor <U  >("Tau", H.getABIJ()));
    auto& D   = this->puttmp(  "D", new Denominator       <U  >(H));

    this->puttmp(  "FAE", new SpinorbitalTensor<U>(    "F(ae)",   H.getAB()));
    this->puttmp(  "FMI", new SpinorbitalTensor<U>(    "F(mi)",   H.getIJ()));
    this->puttmp(  "FME", new SpinorbitalTensor<U>(    "F(me)",   H.getIA()));
    this->puttmp("WMNIJ", new SpinorbitalTensor<U>( "W(mn,ij)", H.getIJKL()));
    this->puttmp("WMNEJ", new SpinorbitalTensor<U>( "W(mn,ej)", H.getIJAK()));
    this->puttmp("WAMIJ", new SpinorbitalTensor<U>("W~(am,ij)", H.getAIJK()));
    this->puttmp("WAMEI", new SpinorbitalTensor<U>("W~(am,ei)", H.getAIBJ()));

    Z(0) = (U)0.0;
    T(0) = (U)0.0;
    T(1) = H.getAI();
    T(2) = H.getABIJ();

    T.weight(D);

    Tau["abij"]  = T(2)["abij"];
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    double mp2 = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));
    Logger::log(arena) << "MP2 energy = " << setprecision(15) << mp2 << endl;
    this->put("mp2", new U(mp2));

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
void CCSD<U>::iterate(const Arena& arena)
{
    const auto& H = this->template get<TwoElectronOperator<U>>("H");

    const SpinorbitalTensor<U>&   fAI =   H.getAI();
    const SpinorbitalTensor<U>&   fME =   H.getIA();
    const SpinorbitalTensor<U>&   fAE =   H.getAB();
    const SpinorbitalTensor<U>&   fMI =   H.getIJ();
    const SpinorbitalTensor<U>& VABIJ = H.getABIJ();
    const SpinorbitalTensor<U>& VMNEF = H.getIJAB();
    const SpinorbitalTensor<U>& VAMEF = H.getAIBC();
    const SpinorbitalTensor<U>& VABEJ = H.getABCI();
    const SpinorbitalTensor<U>& VABEF = H.getABCD();
    const SpinorbitalTensor<U>& VMNIJ = H.getIJKL();
    const SpinorbitalTensor<U>& VMNEJ = H.getIJAK();
    const SpinorbitalTensor<U>& VAMIJ = H.getAIJK();
    const SpinorbitalTensor<U>& VAMEI = H.getAIBJ();

    auto& T   = this->template get   <ExcitationOperator<U,2>>(  "T");
    auto& D   = this->template gettmp<Denominator       <U  >>(  "D");
    auto& Z   = this->template gettmp<ExcitationOperator<U,2>>(  "Z");
    auto& Tau = this->template gettmp<SpinorbitalTensor <U  >>("Tau");

    auto&   FME = this->template gettmp<SpinorbitalTensor<U>>(  "FME");
    auto&   FAE = this->template gettmp<SpinorbitalTensor<U>>(  "FAE");
    auto&   FMI = this->template gettmp<SpinorbitalTensor<U>>(  "FMI");
    auto& WMNIJ = this->template gettmp<SpinorbitalTensor<U>>("WMNIJ");
    auto& WMNEJ = this->template gettmp<SpinorbitalTensor<U>>("WMNEJ");
    auto& WAMIJ = this->template gettmp<SpinorbitalTensor<U>>("WAMIJ");
    auto& WAMEI = this->template gettmp<SpinorbitalTensor<U>>("WAMEI");

    Tau["abij"]  = T(2)["abij"];
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    /**************************************************************************
     *
     * Intermediates for CCSD
     */
      FME[  "me"]  =       fME[  "me"];
      FME[  "me"] +=     VMNEF["mnef"]*T(1)[  "fn"];

      FMI[  "mi"]  =       fMI[  "mi"];
      FMI[  "mi"] += 0.5*VMNEF["mnef"]*T(2)["efin"];
      FMI[  "mi"] +=       FME[  "me"]*T(1)[  "ei"];
      FMI[  "mi"] +=     VMNEJ["nmfi"]*T(1)[  "fn"];

      FAE[  "ae"]  =       fAE[  "ae"];
      FAE[  "ae"] -= 0.5*VMNEF["mnef"]*T(2)["afmn"];
      FAE[  "ae"] -=       FME[  "me"]*T(1)[  "am"];
      FAE[  "ae"] +=     VAMEF["amef"]*T(1)[  "fm"];

    WMNIJ["mnij"]  =     VMNIJ["mnij"];
    WMNIJ["mnij"] += 0.5*VMNEF["mnef"]* Tau["efij"];
    WMNIJ["mnij"] +=     VMNEJ["mnej"]*T(1)[  "ei"];

    WMNEJ["mnej"]  =     VMNEJ["mnej"];
    WMNEJ["mnej"] +=     VMNEF["mnef"]*T(1)[  "fj"];

    WAMIJ["amij"]  =     VAMIJ["amij"];
    WAMIJ["amij"] += 0.5*VAMEF["amef"]* Tau["efij"];
    WAMIJ["amij"] +=     VAMEI["amej"]*T(1)[  "ei"];

    WAMEI["amei"]  =     VAMEI["amei"];
    WAMEI["amei"] += 0.5*VMNEF["mnef"]*T(2)["afni"];
    WAMEI["amei"] +=     VAMEF["amef"]*T(1)[  "fi"];
    WAMEI["amei"] -=     WMNEJ["nmei"]*T(1)[  "an"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * CCSD Iteration
     */
    Z(1)[  "ai"]  =       fAI[  "ai"];
    Z(1)[  "ai"] +=       fAE[  "ae"]*T(1)[  "ei"];
    Z(1)[  "ai"] -=       FMI[  "mi"]*T(1)[  "am"];
    Z(1)[  "ai"] -=     VAMEI["amei"]*T(1)[  "em"];
    Z(1)[  "ai"] +=       FME[  "me"]*T(2)["aeim"];
    Z(1)[  "ai"] += 0.5*VAMEF["amef"]* Tau["efim"];
    Z(1)[  "ai"] -= 0.5*WMNEJ["mnei"]*T(2)["eamn"];

    Z(2)["abij"]  =     VABIJ["abij"];
    Z(2)["abij"] +=     VABEJ["abej"]*T(1)[  "ei"];
    Z(2)["abij"] -=     WAMIJ["amij"]*T(1)[  "bm"];
    Z(2)["abij"] +=       FAE[  "ae"]*T(2)["ebij"];
    Z(2)["abij"] -=       FMI[  "mi"]*T(2)["abmj"];
    Z(2)["abij"] += 0.5*VABEF["abef"]* Tau["efij"];
    Z(2)["abij"] += 0.5*WMNIJ["mnij"]* Tau["abmn"];
    Z(2)["abij"] +=     WAMEI["amei"]*T(2)["ebjm"];
    /*
     *************************************************************************/

    Z.weight(D);
    T += Z;

    Tau["abij"]  = T(2)["abij"];
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];
    this->energy() = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));
    this->conv() = Z.norm(00);

    diis.extrapolate(T, Z);
}

/*
template <typename U>
double CCSD<U>::getProjectedS2(const MOSpace<U>& occ, const MOSpace<U>& vrt,
                               const SpinorbitalTensor<U>& T1,
                               const SpinorbitalTensor<U>& T2)
{
    const Arena& arena = occ.arena;

    int N = occ.nao;
    int nI = occ.nalpha;
    int ni = occ.nbeta;
    int nA = vrt.nalpha;
    int na = vrt.nbeta;

    vector<int> shapeNN = {NS,NS};
    vector<int> shapeNNNN = {NS,NS,NS,NS};
    vector<int> sizeAI = {nA,nI};
    vector<int> sizeAi = {nA,ni};
    vector<int> sizeaI = {na,nI};
    vector<int> sizeai = {na,ni};
    vector<int> sizeIi = {nI,ni};
    vector<int> sizeIn = {nI,N};
    vector<int> sizein = {ni,N};
    vector<int> sizeAaIi = {nA,na,nI,ni};

    const CTFTensor<U>& CA = vrt.Calpha;
    const CTFTensor<U>& Ca = vrt.Cbeta;
    const CTFTensor<U>& CI = occ.Calpha;
    const CTFTensor<U>& Ci = occ.Cbeta;

    //TODO
    CTFTensor<U> S(arena, 2, {N,N}, shapeNN, true);

    CTFTensor<U> DAI(arena, 2, sizeAI, shapeNN, false);
    CTFTensor<U> DAi(arena, 2, sizeAi, shapeNN, false);
    CTFTensor<U> DaI(arena, 2, sizeaI, shapeNN, false);
    CTFTensor<U> Dai(arena, 2, sizeai, shapeNN, false);
    CTFTensor<U> DIj(arena, 2, sizeIi, shapeNN, false);
    CTFTensor<U> DAbIj(arena, 4, sizeAaIi, shapeNNNN, false);
    CTFTensor<U> tmp1(arena, 2, sizeIn, shapeNN, false);
    CTFTensor<U> tmp2(arena, 2, sizein, shapeNN, false);

    tmp1["Iq"] = CI["pI"]*S["pq"];
    DIj["Ij"] = tmp1["Iq"]*Ci["qj"];
    DaI["aI"] = tmp1["Iq"]*Ca["qa"];

    tmp2["iq"] = Ci["pi"]*S["pq"];
    DAi["Ai"] = tmp2["iq"]*CA["qA"];

    DAI["AI"] = DAi["Aj"]*DIj["Ij"];
    Dai["ai"] = DaI["aJ"]*DIj["Ji"];
    DAbIj["AbIj"] = DAi["Aj"]*DaI["bI"];

    const CTFTensor<U>& T1A = T1(1,0,0,1);
    const CTFTensor<U>& T1B = T1(0,0,0,0);
    CTFTensor<U> TauAB(T2(1,0,0,1));

    TauAB["AbIj"] += T1A["AI"]*T1B["bj"];

    U S2 = (U)0;

    U S2T11 = -scalar(DAI*T1A);
    U S2T12 = -scalar(Dai*T1B);
    U S2T2 = -scalar(DAbIj*TauAB);

    return aquarius::abs(S2+S2T11+S2T12+S2T2);
}
*/

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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSD);
REGISTER_TASK(aquarius::cc::CCSD<double>,"ccsd",spec);
