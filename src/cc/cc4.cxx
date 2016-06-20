#include "cc4.hpp"

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
CC4<U>::CC4(const string& name, Config& config)
: Iterative<U>(name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    this->addProduct(Product("double", "mp2", reqs));
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("double", "S2", reqs));
    this->addProduct(Product("double", "multiplicity", reqs));
    this->addProduct(Product("cc4.T", "T", reqs));
    this->addProduct(Product("cc4.Hbar", "Hbar", reqs));
}

template <typename U>
bool CC4<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const auto& H = this->template get<TwoElectronOperator<U>>("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    auto& T   = this->put   (   "T", new ExcitationOperator<U,3>("T", arena, occ, vrt));
    auto& Z   = this->puttmp(   "Z", new ExcitationOperator<U,3>("Z", arena, occ, vrt));
    auto& Tau = this->puttmp( "Tau", new SpinorbitalTensor <U  >("Tau", H.getABIJ()));
    auto& T12 = this->puttmp("T1^2", new SpinorbitalTensor <U  >("T1^2", H.getABIJ()));
    auto& D   = this->puttmp(   "D", new Denominator       <U  >(H));

    this->puttmp(  "FAE", new SpinorbitalTensor<U>(   "F(ae)",   H.getAB()));
    this->puttmp(  "FMI", new SpinorbitalTensor<U>(   "F(mi)",   H.getIJ()));
    this->puttmp(  "FME", new SpinorbitalTensor<U>(   "F(me)",   H.getIA()));
    this->puttmp("WMNIJ", new SpinorbitalTensor<U>("W(mn,ij)", H.getIJKL()));
    this->puttmp("WMNEJ", new SpinorbitalTensor<U>("W(mn,ej)", H.getIJAK()));
    this->puttmp("WAMIJ", new SpinorbitalTensor<U>("W(am,ij)", H.getAIJK()));
    this->puttmp("WAMEI", new SpinorbitalTensor<U>("W(am,ei)", H.getAIBJ()));
    this->puttmp("WABEF", new SpinorbitalTensor<U>("W(ab,ef)", H.getABCD()));
    this->puttmp("WABEJ", new SpinorbitalTensor<U>("W(ab,ej)", H.getABCI()));
    this->puttmp("WAMEF", new SpinorbitalTensor<U>("W(am,ef)", H.getAIBC()));

    this->puttmp("T4", new SpinorbitalTensor<U>("T(abcd,ijkl)", arena,
                                               H.getABIJ().getGroup(),
                                               {vrt, occ}, {4, 0},
                                               {0, 4}));

    this->puttmp("WABCEJK", new SpinorbitalTensor<U>("W~(abc,ejk)", arena,
                                               H.getABIJ().getGroup(),
                                               {vrt, occ}, {3, 0},
                                               {1, 2}));
    this->puttmp("WABMIJK", new SpinorbitalTensor<U>("W~(abm,ijk)", arena,
                                               H.getABIJ().getGroup(),
                                               {vrt, occ}, {2, 1},
                                               {0, 3}));

    Z(0) = (U)0.0;
    T(0) = (U)0.0;
    T(1) = H.getAI();
    T(2) = H.getABIJ();
    T(3) = (U)0.0;

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

    if (this->isUsed("Hbar"))
    {
        this->put("Hbar", new STTwoElectronOperator<U>("Hbar", H, T, true));
    }

    return true;
}

template <typename U>
void CC4<U>::iterate(const Arena& arena)
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

    auto& T   = this->template get   <ExcitationOperator<U,3>>(   "T");
    auto& D   = this->template gettmp<Denominator       <U  >>(   "D");
    auto& Z   = this->template gettmp<ExcitationOperator<U,3>>(   "Z");
    auto& Tau = this->template gettmp<SpinorbitalTensor <U  >>( "Tau");
    auto& T12 = this->template gettmp<SpinorbitalTensor <U  >>("T1^2");

    auto& T4 = this->template gettmp<SpinorbitalTensor<U>>("T4");

    auto&     FME = this->template gettmp<SpinorbitalTensor<U>>(    "FME");
    auto&     FAE = this->template gettmp<SpinorbitalTensor<U>>(    "FAE");
    auto&     FMI = this->template gettmp<SpinorbitalTensor<U>>(    "FMI");
    auto&   WMNIJ = this->template gettmp<SpinorbitalTensor<U>>(  "WMNIJ");
    auto&   WMNEJ = this->template gettmp<SpinorbitalTensor<U>>(  "WMNEJ");
    auto&   WAMIJ = this->template gettmp<SpinorbitalTensor<U>>(  "WAMIJ");
    auto&   WAMEI = this->template gettmp<SpinorbitalTensor<U>>(  "WAMEI");
    auto&   WABEF = this->template gettmp<SpinorbitalTensor<U>>(  "WABEF");
    auto&   WABEJ = this->template gettmp<SpinorbitalTensor<U>>(  "WABEJ");
    auto&   WAMEF = this->template gettmp<SpinorbitalTensor<U>>(  "WAMEF");
    auto& WABCEJK = this->template gettmp<SpinorbitalTensor<U>>("WABCEJK");
    auto& WABMIJK = this->template gettmp<SpinorbitalTensor<U>>("WABMIJK");

    Tau["abij"]  = T(2)["abij"];
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    T12["abij"]  = 0.5*T(1)["ai"]*T(1)["bj"];

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
    Z(2)["abij"] +=       FAE[  "af"]*T(2)["fbij"];
    Z(2)["abij"] -=       FMI[  "ni"]*T(2)["abnj"];
    Z(2)["abij"] += 0.5*VABEF["abef"]* Tau["efij"];
    Z(2)["abij"] += 0.5*WMNIJ["mnij"]* Tau["abmn"];
    Z(2)["abij"] +=     WAMEI["amei"]*T(2)["ebjm"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * Intermediates for CCSDT
     */
    WAMIJ["amij"] +=     WMNEJ["nmej"]*T(2)[  "aein"];
    WAMIJ["amij"] -=     WMNIJ["nmij"]*T(1)[    "an"];
    WAMIJ["amij"] +=       FME[  "me"]*T(2)[  "aeij"];
    WAMIJ["amij"] += 0.5*VMNEF["mnef"]*T(3)["aefijn"];

    WAMEI["amei"] += 0.5*VMNEF["mnef"]*T(2)[  "afni"];
    WAMEI["amei"] += 0.5*WMNEJ["nmei"]*T(1)[    "an"];

    WABEJ["abej"]  =     VABEJ["abej"];
    WABEJ["abej"] +=     VAMEF["amef"]*T(2)[  "fbmj"];
    WABEJ["abej"] += 0.5*WMNEJ["mnej"]*T(2)[  "abmn"];
    WABEJ["abej"] +=     VABEF["abef"]*T(1)[    "fj"];
    WABEJ["abej"] -=     WAMEI["amej"]*T(1)[    "bm"];
    WABEJ["abej"] -= 0.5*VMNEF["mnef"]*T(3)["afbmnj"];

    WAMEI["amei"] -= 0.5*WMNEJ["nmei"]*T(1)[    "an"];

    WABEF["abef"]  =     VABEF["abef"];
    WABEF["abef"] -=     VAMEF["amef"]*T(1)[    "bm"];
    WABEF["abef"] += 0.5*VMNEF["mnef"]* Tau[  "abmn"];

    WAMEF["amef"]  =     VAMEF["amef"];
    WAMEF["amef"] -=     VMNEF["nmef"]*T(1)[    "an"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * CCSDT Iteration
     */
    Z(1)[    "ai"] += 0.25*VMNEF["mnef"]*T(3)["aefimn"];

    Z(2)[  "abij"] +=  0.5*WAMEF["bmef"]*T(3)["aefijm"];
    Z(2)[  "abij"] -=  0.5*WMNEJ["mnej"]*T(3)["abeinm"];
    Z(2)[  "abij"] +=        FME[  "me"]*T(3)["abeijm"];

    Z(3)["abcijk"]  =      WABEJ["bcek"]*T(2)[  "aeij"];
    Z(3)["abcijk"] -=      WAMIJ["bmjk"]*T(2)[  "acim"];
    Z(3)["abcijk"] +=        FAE[  "ce"]*T(3)["abeijk"];
    Z(3)["abcijk"] -=        FMI[  "mk"]*T(3)["abcijm"];
    Z(3)["abcijk"] +=  0.5*WABEF["abef"]*T(3)["efcijk"];
    Z(3)["abcijk"] +=  0.5*WMNIJ["mnij"]*T(3)["abcmnk"];
    Z(3)["abcijk"] +=      WAMEI["amei"]*T(3)["ebcjmk"];
    /*
     **************************************************************************/

    /**************************************************************************
     *
     * Intermediates for CC4
     */
      WAMEI[  "amei"]  =     VAMEI["amei"];
      WAMEI[  "amei"] +=     VAMEF["amef"]*T(1)[  "fi"];
      WAMEI[  "amei"] -=     WMNEJ["nmei"]*T(1)[  "an"];

      WMNIJ[  "mnij"]  =     VMNIJ["mnij"];
      WMNIJ[  "mnij"] +=     VMNEJ["mnej"]*T(1)[  "ei"];
      WMNIJ[  "mnij"] += 0.5*VMNEF["mnef"]* T12["efij"];

      WABEF[  "abef"]  =     VABEF["abef"];
      WABEF[  "abef"] -=     VAMEF["amef"]*T(1)[  "bm"];
      WABEF[  "abef"] += 0.5*VMNEF["mnef"]* T12["abmn"];

      WAMIJ[  "amij"]  =     VAMIJ["amij"];
      WAMIJ[  "amij"] += 0.5*VAMEF["amef"]* T12["efij"];
      WAMIJ[  "amij"] +=     VAMEI["amej"]*T(1)[  "ei"];
      WAMIJ[  "amij"] -=     WMNIJ["nmij"]*T(1)[  "an"];

      WABEJ[  "abej"]  =     VABEJ["abej"];
      WABEJ[  "abej"] += 0.5*VMNEJ["mnej"]* T12["abmn"];
      WABEJ[  "abej"] -=     VAMEI["amej"]*T(1)[  "bm"];
      WABEJ[  "abej"] +=     WABEF["abef"]*T(1)[  "fj"];

    WABCEJK["abcejk"]  = 0.5*WABEF["abef"]*T(2)["fcjk"];

    WABMIJK["abmijk"]  =     WAMEI["amek"]*T(2)["ebij"];
    WABMIJK["abmijk"] -= 0.5*WMNIJ["mnkj"]*T(2)["abin"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * CC4 Iteration
     */
      T4["abcdijkl"]  =     WABCEJK["abcejk"]*T(2)[    "edil"];
      T4["abcdijkl"] -=     WABMIJK["abmijk"]*T(2)[    "cdml"];
      T4["abcdijkl"] +=       WABEJ[  "abej"]*T(3)[  "ecdikl"];
      T4["abcdijkl"] -=       WAMIJ[  "amij"]*T(3)[  "bcdmkl"];

    T4.weight({&D.getDA(), &D.getDI()}, {&D.getDa(), &D.getDi()});

    Z(2)[    "abij"] +=  0.25*VMNEF[  "mnef"]*  T4["abefijmn"];

    Z(3)[  "abcijk"] +=         FME[    "me"]*  T4["abceijkm"];
    Z(3)[  "abcijk"] +=   0.5*WAMEF[  "amef"]*  T4["efbcimjk"];
    Z(3)[  "abcijk"] -=   0.5*WMNEJ[  "mnek"]*  T4["abecijmn"];
    /*
     **************************************************************************/

    Z.weight(D);
    T += Z;

    Tau["abij"]  = T(2)["abij"];
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];
    this->energy() = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));
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
sub_iterations?
    int 0,
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CC4);
REGISTER_TASK(aquarius::cc::CC4<double>,"cc4",spec);
