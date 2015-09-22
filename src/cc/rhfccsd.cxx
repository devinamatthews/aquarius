#include "rhfccsd.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;
using namespace aquarius::convergence;

namespace aquarius
{
namespace cc
{

template <typename U>
RHFCCSD<U>::RHFCCSD(const string& name, Config& config)
: Iterative<U>(name, config), diis_config(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement( "mofock",     "f"));
    reqs.push_back(Requirement("<Ab|Cd>", "VABCD"));
    reqs.push_back(Requirement("<Ab|Ci>", "VABCI"));
    reqs.push_back(Requirement("<Ab|Ij>", "VABIJ"));
    reqs.push_back(Requirement("<Ai|Bj>", "VAIBJ"));
    reqs.push_back(Requirement("<Ai|Jb>", "VAIJB"));
    reqs.push_back(Requirement("<Ai|Jk>", "VAIJK"));
    reqs.push_back(Requirement("<Ij|Kl>", "VIJKL"));
    this->addProduct(Product("double", "mp2", reqs));
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("rhfccsd.T1", "T1", reqs));
    this->addProduct(Product("rhfccsd.T2", "T2", reqs));
    this->addProduct(Product("ccsd.F(AE)", "FAE", reqs));
    this->addProduct(Product("ccsd.F(MI)", "FMI", reqs));
    this->addProduct(Product("ccsd.F(ME)", "FME", reqs));
    this->addProduct(Product("ccsd.W(Am,Ei)", "WAMEI", reqs));
    this->addProduct(Product("ccsd.W(Am,Ie)", "WAMIE", reqs));
    this->addProduct(Product("ccsd.W(Mn,Ej)", "WMNEJ", reqs));
    this->addProduct(Product("ccsd.W(Ab,Ej)", "WABEJ", reqs));
    this->addProduct(Product("ccsd.W(Am,Ij)", "WAMIJ", reqs));
    this->addProduct(Product("ccsd.W(Mn,Ij)", "WMNIJ", reqs));
}

template <typename U>
bool RHFCCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& f = this->template get<OneElectronOperator<U>>("f");

    const Space& occ = f.occ;
    const Space& vrt = f.vrt;
    const vector<int>& nI = occ.nalpha;
    const vector<int>& nA = vrt.nalpha;

    const auto& fAI = f.getAI()({0,0},{0,0});
    const auto& VABIJ = this->template get<SymmetryBlockedTensor<U>>("VABIJ");
    const auto& VAIBJ = this->template get<SymmetryBlockedTensor<U>>("VAIBJ");
    const auto& VAIJB = this->template get<SymmetryBlockedTensor<U>>("VAIJB");

    auto&      T1 = this->put   (     "T1", new SymmetryBlockedTensor<U>(     "T1", arena, occ.group, 2,       {nA,nI},       {NS,NS}, false));
    auto&      T2 = this->put   (     "T2", new SymmetryBlockedTensor<U>(     "T2", arena, occ.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));
    auto&      Z1 = this->puttmp(     "Z1", new SymmetryBlockedTensor<U>(     "Z2", arena, occ.group, 2,       {nA,nI},       {NS,NS}, false));
    auto&      Z2 = this->puttmp(     "Z2", new SymmetryBlockedTensor<U>(     "Z2", arena, occ.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));
    auto&     Tau = this->puttmp(    "Tau", new SymmetryBlockedTensor<U>(    "Tau", arena, occ.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));
    auto&    T2SA = this->puttmp(   "T2SA", new SymmetryBlockedTensor<U>(   "T2SA", arena, occ.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));
    auto&   TauSA = this->puttmp(  "TauSA", new SymmetryBlockedTensor<U>(  "TauSA", arena, occ.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));
    auto& VABIJSA = this->puttmp("VABIJSA", new SymmetryBlockedTensor<U>("VABIJSA", arena, occ.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));
    auto& VAIJBSA = this->puttmp("VAIJBSA", new SymmetryBlockedTensor<U>("VAIJBSA", arena, occ.group, 4, {nA,nI,nI,nA}, {NS,NS,NS,NS}, false));
    auto&       D = this->puttmp(      "D", new Denominator          <U>(f));

    auto&   FAE = this->put(  "FAE", new SymmetryBlockedTensor<U>(   "F(AE)", arena, occ.group, 2,       {nA,nA},       {NS,NS}, false));
    auto&   FMI = this->put(  "FMI", new SymmetryBlockedTensor<U>(   "F(MI)", arena, occ.group, 2,       {nI,nI},       {NS,NS}, false));
    auto&   FME = this->put(  "FME", new SymmetryBlockedTensor<U>(   "F(ME)", arena, occ.group, 2,       {nI,nA},       {NS,NS}, false));
    auto& WMNIJ = this->put("WMNIJ", new SymmetryBlockedTensor<U>("W(Mn,Ij)", arena, occ.group, 4, {nI,nI,nI,nI}, {NS,NS,NS,NS}, false));
    auto& WMNEJ = this->put("WMNEJ", new SymmetryBlockedTensor<U>("W(Mn,Ej)", arena, occ.group, 4, {nI,nI,nA,nI}, {NS,NS,NS,NS}, false));
    auto& WAMIJ = this->put("WAMIJ", new SymmetryBlockedTensor<U>("W(Am,Ij)", arena, occ.group, 4, {nA,nI,nI,nI}, {NS,NS,NS,NS}, false));
    auto& WAMEI = this->put("WAMEI", new SymmetryBlockedTensor<U>("W(Am,Ei)", arena, occ.group, 4, {nA,nI,nA,nI}, {NS,NS,NS,NS}, false));
    auto& WAMIE = this->put("WAMIE", new SymmetryBlockedTensor<U>("W(Am,Ie)", arena, occ.group, 4, {nA,nI,nI,nA}, {NS,NS,NS,NS}, false));

    T1 = fAI;
    T1.weight({&D.getDA(), &D.getDI()});

    T2 = VABIJ;
    T2.weight({&D.getDA(), &D.getDA(), &D.getDI(), &D.getDI()});

    Tau["abij"]  = T2["abij"];
    Tau["abij"] += T1["ai"]*T1["bj"];

    VABIJSA["abij"]  = 2*VABIJ["abij"];
    VABIJSA["abij"] -=   VABIJ["abji"];

    VAIJBSA["aijb"]  = 2*VAIJB["aijb"];
    VAIJBSA["aijb"] -=   VAIBJ["aibj"];

    double mp2 = 2*real(scalar(fAI*T1)) + real(scalar(VABIJSA*Tau));
    Logger::log(arena) << "MP2 energy = " << setprecision(15) << mp2 << endl;
    this->put("mp2", new U(mp2));

    this->puttmp("DIIS", new DIIS<SymmetryBlockedTensor<U>>(diis_config, 2, 2));

    CTF_Timer_epoch ep(this->name.c_str());
    ep.begin();
    Iterative<U>::run(dag, arena);
    ep.end();

    this->put("energy", new U(this->energy()));
    this->put("convergence", new U(this->conv()));

    /*
     * Form Hbar
     */
    {
        const auto& fAB = f.getAB()({0,0},{0,0});
        const auto& fIJ = f.getIJ()({0,0},{0,0});

        const auto& VABCD = this->template get<SymmetryBlockedTensor<U>>("VABCD");
        const auto& VABCI = this->template get<SymmetryBlockedTensor<U>>("VABCI");
        const auto& VAIJK = this->template get<SymmetryBlockedTensor<U>>("VAIJK");
        const auto& VIJKL = this->template get<SymmetryBlockedTensor<U>>("VIJKL");

          Tau["abij"]  =          T2["abij"];
          Tau["abij"] +=          T1[  "ai"]  *T1[  "bj"];

         T2SA["abij"]  =        2*T2["abij"];
         T2SA["abij"] -=          T2["abji"];

          FME[  "me"]  =         fAI[  "em"];
          FME[  "me"] +=     VABIJSA["efmn"]*  T1[  "fn"];

          FMI[  "mi"]  =         fIJ[  "mi"];
          FMI[  "mi"] +=     VABIJSA["efmn"]*  T2["efin"];
          FMI[  "mi"] +=         FME[  "me"]*  T1[  "ei"];
          FMI[  "mi"] +=   2*  VAIJK["finm"]*  T1[  "fn"];
          FMI[  "mi"] -=       VAIJK["fimn"]*  T1[  "fn"];

          FAE[  "ae"]  =         fAB[  "ae"];
          FAE[  "ae"] -=     VABIJSA["efmn"]*  T2["afmn"];
          FAE[  "ae"] -=         FME[  "me"]*  T1[  "am"];
          FAE[  "ae"] +=   2*  VABCI["efam"]*  T1[  "fm"];
          FAE[  "ae"] -=       VABCI["feam"]*  T1[  "fm"];

        WMNIJ["mnij"]  = 0.5*  VIJKL["mnij"];
        WMNIJ["mnij"] += 0.5*  VABIJ["efmn"]* Tau["efij"];
        WMNIJ["mnij"] +=       VAIJK["ejmn"]*  T1[  "ei"];
        WMNIJ["mnij"] +=       WMNIJ["nmji"];

        WMNEJ["mnej"]  =       VAIJK["ejmn"];
        WMNEJ["mnej"] +=       VABIJ["efmn"]*  T1[  "fj"];

        WAMIJ["amij"]  =       VAIJK["amij"];
        WAMIJ["amij"] +=       VABCI["efam"]* Tau["efij"];
        WAMIJ["amij"] +=       VAIBJ["amej"]*  T1[  "ei"];
        WAMIJ["amij"] +=       VAIJB["amie"]*  T1[  "ej"];
        WAMIJ["amij"] -=       WMNIJ["nmij"]*  T1[  "an"];
        WAMIJ["amij"] +=       WMNEJ["nmej"]*T2SA["aein"];
        WAMIJ["amij"] -=       WMNEJ["mnej"]*  T2["aein"];
        WAMIJ["amij"] -=       WMNEJ["mnei"]*  T2["aenj"];
        WAMIJ["amij"] +=         FME[  "me"]*  T2["aeij"];

        WAMIE["amie"]  =     VAIJBSA["amie"];
        WAMIE["amie"] +=     VABIJSA["efmn"]*T2SA["afin"];
        WAMIE["amie"] +=   2*  VABCI["feam"]*  T1[  "fi"];
        WAMIE["amie"] -=       VABCI["efam"]*  T1[  "fi"];

        WAMEI["amei"]  =       VAIBJ["amei"];
        WAMEI["amei"] -=       VABIJ["efnm"]*  T2["afni"];
        WAMEI["amei"] +=       VABCI["efam"]*  T1[  "fi"];

        if (this->isUsed("WABEJ"))
        {
            auto& WABEJ = this->put("WABEJ", new SymmetryBlockedTensor<U>("W(Ab,Ej)", arena, occ.group, 4, {nA,nA,nA,nI}, {NS,NS,NS,NS}, false));

            WABEJ["abej"]  =     VABCI["abej"];
            WABEJ["abej"] +=     WMNEJ["mnej"]* Tau["abmn"];
            WABEJ["abej"] -= 0.5*WAMIE["bmje"]*  T1[  "am"];
            WABEJ["abej"] -= 0.5*WAMEI["bmej"]*  T1[  "am"];
            WABEJ["abej"] -=     WAMEI["amej"]*  T1[  "bm"];
            WABEJ["abej"] +=     VABCD["abef"]*  T1[  "fj"];
            WABEJ["abej"] +=     VABCI["efam"]*T2SA["fbmj"];
            WABEJ["abej"] -=     VABCI["feam"]*  T2["fbmj"];
            WABEJ["abej"] -=     VABCI["febm"]*  T2["afmj"];
            WABEJ["abej"] -=       FME[  "me"]*  T2["abmj"];
        }

        WAMIE["amie"] -=   2*  WMNEJ["mnei"]*  T1[  "an"];
        WAMIE["amie"] +=       WMNEJ["nmei"]*  T1[  "an"];

        WAMEI["amei"] -=       WMNEJ["nmei"]*  T1[  "an"];
    }

    return true;
}

template <typename U>
void RHFCCSD<U>::iterate(const Arena& arena)
{
    const auto& f = this->template get<OneElectronOperator<U>>("f");
    const auto& fAI = f.getAI()({0,0},{0,0});
    const auto& fAB = f.getAB()({0,0},{0,0});
    const auto& fIJ = f.getIJ()({0,0},{0,0});

    const auto& VABCD = this->template get<SymmetryBlockedTensor<U>>("VABCD");
    const auto& VABCI = this->template get<SymmetryBlockedTensor<U>>("VABCI");
    const auto& VABIJ = this->template get<SymmetryBlockedTensor<U>>("VABIJ");
    const auto& VAIBJ = this->template get<SymmetryBlockedTensor<U>>("VAIBJ");
    const auto& VAIJB = this->template get<SymmetryBlockedTensor<U>>("VAIJB");
    const auto& VAIJK = this->template get<SymmetryBlockedTensor<U>>("VAIJK");
    const auto& VIJKL = this->template get<SymmetryBlockedTensor<U>>("VIJKL");

    auto&      T1 = this->template get   <SymmetryBlockedTensor<U>>(     "T1");
    auto&      T2 = this->template get   <SymmetryBlockedTensor<U>>(     "T2");
    auto&      Z1 = this->template gettmp<SymmetryBlockedTensor<U>>(     "Z1");
    auto&      Z2 = this->template gettmp<SymmetryBlockedTensor<U>>(     "Z2");
    auto&     Tau = this->template gettmp<SymmetryBlockedTensor<U>>(    "Tau");
    auto&    T2SA = this->template gettmp<SymmetryBlockedTensor<U>>(   "T2SA");
    auto&   TauSA = this->template gettmp<SymmetryBlockedTensor<U>>(  "TauSA");
    auto& VABIJSA = this->template gettmp<SymmetryBlockedTensor<U>>("VABIJSA");
    auto& VAIJBSA = this->template gettmp<SymmetryBlockedTensor<U>>("VAIJBSA");
    auto&       D = this->template gettmp<Denominator          <U>>(      "D");

    auto&   FAE = this->template get<SymmetryBlockedTensor<U>>(  "FAE");
    auto&   FMI = this->template get<SymmetryBlockedTensor<U>>(  "FMI");
    auto&   FME = this->template get<SymmetryBlockedTensor<U>>(  "FME");
    auto& WMNIJ = this->template get<SymmetryBlockedTensor<U>>("WMNIJ");
    auto& WMNEJ = this->template get<SymmetryBlockedTensor<U>>("WMNEJ");
    auto& WAMIJ = this->template get<SymmetryBlockedTensor<U>>("WAMIJ");
    auto& WAMEI = this->template get<SymmetryBlockedTensor<U>>("WAMEI");
    auto& WAMIE = this->template get<SymmetryBlockedTensor<U>>("WAMIE");

    auto& diis = this->template gettmp<DIIS<SymmetryBlockedTensor<U>>>("DIIS");

    Tau["abij"]  = T2["abij"];
    Tau["abij"] += T1["ai"]*T1["bj"];

    TauSA["abij"]  = 2*Tau["abij"];
    TauSA["abij"] -=   Tau["abji"];

    T2SA["abij"]  = 2*T2["abij"];
    T2SA["abij"] -=   T2["abji"];

    /**************************************************************************
     *
     * Intermediates for CCSD
     */
      FME[  "me"]  =         fAI[  "em"];
      FME[  "me"] +=     VABIJSA["efmn"]*   T1[  "fn"];

      FMI[  "mi"]  =         fIJ[  "mi"];
      FMI[  "mi"] +=     VABIJSA["efmn"]*   T2["efin"];
      FMI[  "mi"] +=         FME[  "me"]*   T1[  "ei"];
      FMI[  "mi"] +=   2*  VAIJK["finm"]*   T1[  "fn"];
      FMI[  "mi"] -=       VAIJK["fimn"]*   T1[  "fn"];

      FAE[  "ae"]  =         fAB[  "ae"];
      FAE[  "ae"] -=     VABIJSA["efmn"]*   T2["afmn"];
      FAE[  "ae"] -=         FME[  "me"]*   T1[  "am"];
      FAE[  "ae"] +=   2*  VABCI["efam"]*   T1[  "fm"];
      FAE[  "ae"] -=       VABCI["feam"]*   T1[  "fm"];

    WMNIJ["mnij"]  = 0.5*  VIJKL["mnij"];
    WMNIJ["mnij"] += 0.5*  VABIJ["efmn"]*  Tau["efij"];
    WMNIJ["mnij"] +=       VAIJK["ejmn"]*   T1[  "ei"];
    WMNIJ["mnij"] +=       WMNIJ["nmji"];

    WMNEJ["mnej"]  =       VAIJK["ejmn"];
    WMNEJ["mnej"] +=       VABIJ["efmn"]*   T1[  "fj"];

    WAMIJ["amij"]  =       VAIJK["amij"];
    WAMIJ["amij"] +=       VABCI["efam"]*  Tau["efij"];
    WAMIJ["amij"] +=       VAIBJ["amej"]*   T1[  "ei"];
    WAMIJ["amij"] +=       VAIJB["amie"]*   T1[  "ej"];

    WAMIE["amie"]  =     VAIJBSA["amie"];
    WAMIE["amie"] += 0.5*VABIJSA["efmn"]* T2SA["afin"];
    WAMIE["amie"] +=   2*  VABCI["feam"]*   T1[  "fi"];
    WAMIE["amie"] -=       VABCI["efam"]*   T1[  "fi"];
    WAMIE["amie"] -=   2*  WMNEJ["mnei"]*   T1[  "an"];
    WAMIE["amie"] +=       WMNEJ["nmei"]*   T1[  "an"];

    WAMEI["amei"]  =       VAIBJ["amei"];
    WAMEI["amei"] -= 0.5*  VABIJ["efnm"]*   T2["afni"];
    WAMEI["amei"] +=       VABCI["efam"]*   T1[  "fi"];
    WAMEI["amei"] -=       WMNEJ["nmei"]*   T1[  "an"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * CCSD Iteration
     */
    Z1[  "ai"]  =       fAI[  "ai"];
    Z1[  "ai"] +=       fAB[  "ae"]*   T1[  "ei"];
    Z1[  "ai"] -=       FMI[  "mi"]*   T1[  "am"];
    Z1[  "ai"] +=   VAIJBSA["amie"]*   T1[  "em"];
    Z1[  "ai"] +=       FME[  "me"]* T2SA["aeim"];
    Z1[  "ai"] +=     VABCI["efam"]*TauSA["efim"];
    Z1[  "ai"] -=     WMNEJ["mnei"]* T2SA["eamn"];

    Z2["abij"]  = 0.5*VABIJ["abij"];
    Z2["abij"] +=     VABCI["abej"]*   T1[  "ei"];
    Z2["abij"] -=     WAMIJ["amij"]*   T1[  "bm"];
    Z2["abij"] +=       FAE[  "ae"]*   T2["ebij"];
    Z2["abij"] -=       FMI[  "mi"]*   T2["abmj"];
    Z2["abij"] += 0.5*VABCD["abef"]*  Tau["efij"];
    Z2["abij"] += 0.5*WMNIJ["mnij"]*  Tau["abmn"];
    Z2["abij"] += 0.5*WAMIE["amie"]* T2SA["ebmj"];
    Z2["abij"] -= 0.5*WAMEI["amei"]*   T2["ebjm"];
    Z2["abij"] -=     WAMEI["amej"]*   T2["ebim"];

    Z2["abij"] +=        Z2["baji"];
    /*
     *************************************************************************/

    Z1.weight({&D.getDA(), &D.getDI()});
    Z2.weight({&D.getDA(), &D.getDA(), &D.getDI(), &D.getDI()});

    T1 += Z1;
    T2 += Z2;

    Tau["abij"]  = T2["abij"];
    Tau["abij"] += T1["ai"]*T1["bj"];

    this->energy() = 2*real(scalar(fAI*T1)) + real(scalar(VABIJSA*Tau));
    this->conv() = max(Z1.norm(00), Z2.norm(00));

    diis.extrapolate(ptr_vector<SymmetryBlockedTensor<U>>{&T1, &T2},
                     ptr_vector<SymmetryBlockedTensor<U>>{&Z1, &Z2});
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::RHFCCSD);
REGISTER_TASK(aquarius::cc::RHFCCSD<double>,"rhfccsd",spec);
