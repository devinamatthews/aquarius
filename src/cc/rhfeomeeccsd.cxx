#include "rhfeomeeccsd.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::convergence;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace cc
{

template <typename U>
struct RHFInnerProd
{
    bool triplet;

    template <typename a_container, typename b_container>
    U operator()(const a_container& a, const b_container& b) const
    {
        if (triplet)
        {
            return scalar(a[0][  "ai"]*b[0][  "ai"]) +
                   scalar(a[1]["abij"]*b[1]["abij"]);
        }
        else
        {
            return 2*scalar(a[0][  "ai"]*b[0][  "ai"]) +
                   2*scalar(a[1]["abij"]*b[1]["abij"]) -
                     scalar(a[1]["abij"]*b[1]["abji"]);
        }
    }
};

template <typename U>
struct RHFWeight
{
    template <typename a_container>
    void operator()(a_container& a, const op::Denominator<U>& D, U omega) const
    {
        a[0].weight({&D.getDA(), &D.getDI()}, omega);
        a[1].weight({&D.getDA(), &D.getDA(), &D.getDI(), &D.getDI()}, omega);
    }
};

template <typename U> using RHFDavidson = Davidson<SymmetryBlockedTensor<U>, RHFInnerProd<U>, RHFWeight<U>>;

template <typename U>
RHFEOMEECCSD<U>::RHFEOMEECCSD(const string& name, Config& config)
: Iterative<U>(name, config),
  davidson_config(config.get("davidson")),
  nroot(config.get<int>("nroot")),
  nsinglet(config.get<int>("nsinglet")),
  ntriplet(config.get<int>("ntriplet"))
{
    vector<Requirement> reqs;
    reqs.emplace_back("molecule", "molecule");
    reqs.emplace_back("mofock", "f");
    reqs.emplace_back("<Ab|Cd>", "VABCD");
    reqs.emplace_back("<Ab|Ci>", "VABCI");
    reqs.emplace_back("<Ab|Ij>", "VABIJ");
    reqs.emplace_back("rhfccsd.T1", "T1");
    reqs.emplace_back("rhfccsd.T2", "T2");
    reqs.emplace_back("ccsd.F(AE)", "FAE");
    reqs.emplace_back("ccsd.F(MI)", "FMI");
    reqs.emplace_back("ccsd.F(ME)", "FME");
    reqs.emplace_back("ccsd.W(Am,Ei)", "WAMEI");
    reqs.emplace_back("ccsd.W(Am,Ie)", "WAMIESA");
    reqs.emplace_back("ccsd.W(Mn,Ej)", "WMNEJ");
    reqs.emplace_back("ccsd.W(Ab,Ej)", "WABEJ");
    reqs.emplace_back("ccsd.W(Am,Ij)", "WAMIJ");
    reqs.emplace_back("ccsd.W(Mn,Ij)", "WMNIJ");
    reqs.emplace_back("rhftda.evals_sing", "evals_sing");
    reqs.emplace_back("rhftda.evals_trip", "evals_trip");
    reqs.emplace_back("rhftda.evecs_sing", "evecs_sing");
    reqs.emplace_back("rhftda.evecs_trip", "evecs_trip");
    this->addProduct("rhfeomeeccsd.energy", "energy", reqs);
    this->addProduct("rhfeomeeccsd.convergence", "convergence", reqs);
    this->addProduct("rhfeomeeccsd.R1", "R1s", reqs);
    this->addProduct("rhfeomeeccsd.R2", "R2s", reqs);
}

template <typename U>
bool RHFEOMEECCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = this->template get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();
    int nirrep = group.getNumIrreps();

    const auto& f = this->template get<OneElectronOperator<U>>("f");
    const Space& occ = f.occ;
    const Space& vrt = f.vrt;
    const vector<int>& nI = occ.nalpha;
    const vector<int>& nA = vrt.nalpha;

    this->puttmp("D", new Denominator<U>(f));
    const auto& evecs_sing = this->template get<vector<unique_vector<SymmetryBlockedTensor<U>>>>("evecs_sing");
    const auto& evecs_trip = this->template get<vector<unique_vector<SymmetryBlockedTensor<U>>>>("evecs_trip");
    const auto& evals_sing = this->template get<vector<vector<U>>>("evals_sing");
    const auto& evals_trip = this->template get<vector<vector<U>>>("evals_trip");

    const auto&   VABIJ = this->template get<SymmetryBlockedTensor<U>>(  "VABIJ");
    const auto&      T1 = this->template get<SymmetryBlockedTensor<U>>(     "T1");
    const auto&      T2 = this->template get<SymmetryBlockedTensor<U>>(     "T2");
    const auto& WAMIESA = this->template get<SymmetryBlockedTensor<U>>("WAMIESA");
    const auto&   WAMEI = this->template get<SymmetryBlockedTensor<U>>(  "WAMEI");

    auto&   WAMIE = this->puttmp(  "WAMIE", new SymmetryBlockedTensor<U>(  "WAMIE", arena, occ.group, 4, {nA,nI,nI,nA}, {NS,NS,NS,NS}, false));
    auto& VABIJSA = this->puttmp("VABIJSA", new SymmetryBlockedTensor<U>("VABIJSA", arena, occ.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));
    auto&     Tau = this->puttmp(    "Tau", new SymmetryBlockedTensor<U>(    "Tau", arena, occ.group, 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));

    VABIJSA["abij"]  = 2*VABIJ["abij"];
    VABIJSA["abij"] -=   VABIJ["abji"];

    Tau["abij"]  = T2["abij"];
    Tau["abij"] += T1["ai"]*T1["bj"];

    WAMIE["amie"]  = 0.5*WAMIESA["amie"];
    WAMIE["amie"] += 0.5*  WAMEI["amei"];

    vector<tuple<U,int,int,int>> tda_sorted;

    int nsing_tot = 0;
    int ntrip_tot = 0;
    for (int i = 0;i < nirrep;i++)
    {
        nsing_tot += evecs_sing[i].size();
        ntrip_tot += evecs_trip[i].size();
    }

    assert(nroot > 0);
    assert(nsinglet >= 0);
    assert(ntriplet >= 0);
    assert(nsinglet+ntriplet <= nroot);
    assert(nsinglet <= nsing_tot);
    assert(ntriplet <= ntrip_tot);
    assert(nroot <= ntrip_tot+nsing_tot);

    for (int i = 0;i < nirrep;i++)
    {
        tda_sorted += zip(evals_sing[i],
                          vector<int>(evals_sing[i].size(), 0),
                          vector<int>(evals_sing[i].size(), i),
                          vector<int>(range<int>(evals_sing[i].size())));
        tda_sorted += zip(evals_trip[i],
                          vector<int>(evals_trip[i].size(), 1),
                          vector<int>(evals_trip[i].size(), i),
                          vector<int>(range<int>(evals_trip[i].size())));
    }

    sort(tda_sorted);

    vector<tuple<U,int,int,int>> roots;

    for (auto& root : tda_sorted)
    {
        if (nsinglet == 0) break;
        if (get<1>(root) == 0)
        {
            nsinglet--;
            nroot--;
            roots.push_back(root);
            get<1>(root) = -1;
        }
    }

    for (auto& root : tda_sorted)
    {
        if (ntriplet == 0) break;
        if (get<1>(root) == 1)
        {
            ntriplet--;
            nroot--;
            roots.push_back(root);
            get<1>(root) = -1;
        }
    }

    for (auto& root : tda_sorted)
    {
        if (nroot == 0) break;
        if (get<1>(root) != -1)
        {
            nroot--;
            roots.push_back(root);
            get<1>(root) = -1;
        }
    }

    auto& R1s = this->put("R1s", new vector<vector<unique_vector<SymmetryBlockedTensor<U>>>>(nirrep));
    auto& R2s = this->put("R2s", new vector<vector<unique_vector<SymmetryBlockedTensor<U>>>>(nirrep));
    auto& Es = this->put("energy", new vector<vector<vector<U>>>(nirrep));
    auto& convs = this->put("convergence", new vector<vector<vector<U>>>(nirrep));

    int idx = 1;
    for (int irrep = 0;irrep < nirrep;irrep++)
    {
        R1s[irrep].resize(2);
        R2s[irrep].resize(2);
        Es[irrep].resize(2);
        convs[irrep].resize(2);

        for (int spin : {0,1})
        {
            quintet = false;
            triplet = spin == 1;

            auto& davidson = this->puttmp("Davidson", new RHFDavidson<U>(davidson_config, 1, 2, RHFInnerProd<U>{triplet}));

            auto& R1 = this->puttmp("R1", new SymmetryBlockedTensor<U>("R1", arena, group, group.getIrrep(irrep), 2,       {nA,nI},       {NS,NS}, false));
            auto& R2 = this->puttmp("R2", new SymmetryBlockedTensor<U>("R2", arena, group, group.getIrrep(irrep), 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));

            this->puttmp(  "XMI", new SymmetryBlockedTensor<U>(   "X(MI)", arena, group, group.getIrrep(irrep), 2,       {nI,nI},       {NS,NS}, false));
            this->puttmp(  "XAE", new SymmetryBlockedTensor<U>(   "X(AE)", arena, group, group.getIrrep(irrep), 2,       {nA,nA},       {NS,NS}, false));
            this->puttmp(  "XME", new SymmetryBlockedTensor<U>(   "X(ME)", arena, group, group.getIrrep(irrep), 2,       {nI,nA},       {NS,NS}, false));
            this->puttmp("XAMIJ", new SymmetryBlockedTensor<U>("X(Am,Ij)", arena, group, group.getIrrep(irrep), 4, {nA,nI,nI,nI}, {NS,NS,NS,NS}, false));
            this->puttmp("XBMJI", new SymmetryBlockedTensor<U>("X(bM,jI)", arena, group, group.getIrrep(irrep), 4, {nA,nI,nI,nI}, {NS,NS,NS,NS}, false));
            this->puttmp("XMNIJ", new SymmetryBlockedTensor<U>("X(Mn,Ij)", arena, group, group.getIrrep(irrep), 4, {nI,nI,nI,nI}, {NS,NS,NS,NS}, false));
            this->puttmp( "R2SA", new SymmetryBlockedTensor<U>(    "R2SA", arena, group, group.getIrrep(irrep), 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));
            this->puttmp(   "Z1", new SymmetryBlockedTensor<U>(      "Z1", arena, group, group.getIrrep(irrep), 2,       {nA,nI},       {NS,NS}, false));
            this->puttmp(   "Z2", new SymmetryBlockedTensor<U>(      "Z2", arena, group, group.getIrrep(irrep), 4, {nA,nA,nI,nI}, {NS,NS,NS,NS}, false));

            for (auto& root : roots)
            {
                if (spin != get<1>(root)) continue;
                if (irrep != get<2>(root)) continue;

                int which = get<3>(root);

                Logger::log(arena) << "Starting root number " << idx << endl;
                Logger::log(arena) << "Guess energy: " << fixed << setprecision(12) << get<0>(root) << endl;

                if (quintet)
                {
                    Logger::log(arena) << "Quintet initial guess" << endl;
                }
                else if (triplet)
                {
                    Logger::log(arena) << "Triplet initial guess" << endl;
                }
                else
                {
                    Logger::log(arena)<< "Singlet initial guess" << endl;
                }

                if (triplet) R1 = evecs_trip[irrep][which];
                else         R1 = evecs_sing[irrep][which];
                R2 = 0;

                Iterative<U>::run(dag, arena);

                if (!this->isConverged())
                {
                    this->error(arena) << "Root " << idx << " did not converge." << endl;
                }

                davidson.getSolution(0, ptr_vector<SymmetryBlockedTensor<U>>{&R1, &R2});
                davidson.nextRoot(1, 2, RHFInnerProd<U>{triplet});
                idx++;

                Es[irrep][spin].push_back(this->energy());
                convs[irrep][spin].push_back(this->conv());
                R1s[irrep][spin].emplace_back(R1);
                R2s[irrep][spin].emplace_back(R2);
            }
        }
    }

    return true;
}

template <typename U>
void RHFEOMEECCSD<U>::iterate(const Arena& arena)
{
    const auto&   VABCD = this->template get<SymmetryBlockedTensor<U>>(  "VABCD");
    const auto&   VABCI = this->template get<SymmetryBlockedTensor<U>>(  "VABCI");
    const auto&   VABIJ = this->template get<SymmetryBlockedTensor<U>>(  "VABIJ");
    const auto&     FAE = this->template get<SymmetryBlockedTensor<U>>(    "FAE");
    const auto&     FMI = this->template get<SymmetryBlockedTensor<U>>(    "FMI");
    const auto&     FME = this->template get<SymmetryBlockedTensor<U>>(    "FME");
    const auto&   WMNIJ = this->template get<SymmetryBlockedTensor<U>>(  "WMNIJ");
    const auto&   WMNEJ = this->template get<SymmetryBlockedTensor<U>>(  "WMNEJ");
    const auto&   WAMIJ = this->template get<SymmetryBlockedTensor<U>>(  "WAMIJ");
    const auto& WAMIESA = this->template get<SymmetryBlockedTensor<U>>("WAMIESA");
    const auto&   WAMEI = this->template get<SymmetryBlockedTensor<U>>(  "WAMEI");
    const auto&   WABEJ = this->template get<SymmetryBlockedTensor<U>>(  "WABEJ");

    auto&    T1 = this->template get   <SymmetryBlockedTensor<U>>(   "T1");
    auto&    T2 = this->template get   <SymmetryBlockedTensor<U>>(   "T2");
    auto&   Tau = this->template gettmp<SymmetryBlockedTensor<U>>(  "Tau");
    auto&    R1 = this->template gettmp<SymmetryBlockedTensor<U>>(   "R1");
    auto&    R2 = this->template gettmp<SymmetryBlockedTensor<U>>(   "R2");
    auto&    Z1 = this->template gettmp<SymmetryBlockedTensor<U>>(   "Z1");
    auto&    Z2 = this->template gettmp<SymmetryBlockedTensor<U>>(   "Z2");
    auto&  R2SA = this->template gettmp<SymmetryBlockedTensor<U>>( "R2SA");
    auto& WAMIE = this->template gettmp<SymmetryBlockedTensor<U>>("WAMIE");

    auto&   XMI = this->template gettmp<SymmetryBlockedTensor<U>>(  "XMI");
    auto&   XAE = this->template gettmp<SymmetryBlockedTensor<U>>(  "XAE");
    auto&   XME = this->template gettmp<SymmetryBlockedTensor<U>>(  "XME");
    auto& XAMIJ = this->template gettmp<SymmetryBlockedTensor<U>>("XAMIJ");
    auto& XBMJI = this->template gettmp<SymmetryBlockedTensor<U>>("XBMJI");
    auto& XMNIJ = this->template gettmp<SymmetryBlockedTensor<U>>("XMNIJ");

    auto& D = this->template gettmp<Denominator<U>>("D");
    auto& davidson = this->template gettmp<RHFDavidson<U>>("Davidson");

    if (quintet)
    {
        XAMIJ["amij"]  =      VABCI["efam"]*  R2["efij"];
        XMNIJ["mnij"]  =      VABIJ["efmn"]*  R2["efij"];

           Z2["abij"]  = -0.5*XAMIJ["amij"]*  T1[  "bm"];
           Z2["abij"] +=  0.5*  FAE[  "ae"]*  R2["ebij"];
           Z2["abij"] -=  0.5*  FMI[  "mi"]*  R2["abmj"];
           Z2["abij"] += 0.25*WMNIJ["mnij"]*  R2["abmn"];
           Z2["abij"] += 0.25*XMNIJ["mnij"]* Tau["abmn"];
           Z2["abij"] += 0.25*VABCD["abef"]*  R2["efij"];
           Z2["abij"] -=      WAMEI["amei"]*  R2["ebmj"];

           Z2["abij"] -= Z2["abji"];
           Z2["abij"] -= Z2["baij"];
    }
    else if (triplet)
    {
          XMI[  "mi"]  =     VABIJ["efmn"]*  R2["efin"];

           Z1[  "ai"]  =       FAE[  "ae"]*  R1[  "ei"];
           Z1[  "ai"] -=       FMI[  "mi"]*  R1[  "am"];
           Z1[  "ai"] -=       XMI[  "mi"]*  T1[  "am"];
           Z1[  "ai"] -=     WAMEI["amei"]*  R1[  "em"];
           Z1[  "ai"] +=       FME[  "me"]*  R2["aeim"];
           Z1[  "ai"] +=     VABCI["efam"]*  R2["efim"];
           Z1[  "ai"] -=     WMNEJ["mnei"]*  R2["aenm"];

          XMI[  "mi"] -=     WMNEJ["mnei"]*  R1[  "en"];

          XME[  "me"]  =  -  VABIJ["efnm"]*  R1[  "fn"];

          XAE[  "ae"]  =  -  VABCI["feam"]*  R1[  "fm"];
          XAE[  "ae"] -=       XME[  "me"]*  T1[  "am"];
          XAE[  "ae"] -=     VABIJ["efmn"]*  R2["afmn"];

        XAMIJ["amij"]  =     VABCI["efam"]*  R2["efij"];
        XBMJI["bmji"]  =     VABCI["febm"]*  R2["efij"];

        XMNIJ["mnij"]  =     VABIJ["efmn"]*  R2["efij"];

           Z2["abij"]  = 2*  WABEJ["abej"]*  R1[  "ei"];
           Z2["abij"] -=     WABEJ["baej"]*  R1[  "ei"];
           Z2["abij"] -=     WABEJ["abei"]*  R1[  "ej"];
           Z2["abij"] -= 2*  WAMIJ["bmji"]*  R1[  "am"];
           Z2["abij"] +=     WAMIJ["bmij"]*  R1[  "am"];
           Z2["abij"] +=     WAMIJ["amji"]*  R1[  "bm"];
           Z2["abij"] -=     XBMJI["bmji"]*  T1[  "am"];
           Z2["abij"] -=     XAMIJ["amij"]*  T1[  "bm"];
           Z2["abij"] +=       FAE[  "ae"]*  R2["ebij"];
           Z2["abij"] +=       FAE[  "be"]*  R2["aeij"];
           Z2["abij"] -=       FMI[  "mi"]*  R2["abmj"];
           Z2["abij"] -=       FMI[  "mj"]*  R2["abim"];
           Z2["abij"] += 2*    XAE[  "ae"]*  T2["ebij"];
           Z2["abij"] -=       XAE[  "ae"]*  T2["beij"];
           Z2["abij"] -=       XAE[  "be"]*  T2["eaij"];
           Z2["abij"] -= 2*    XMI[  "mi"]*  T2["abmj"];
           Z2["abij"] +=       XMI[  "mi"]*  T2["abjm"];
           Z2["abij"] +=       XMI[  "mj"]*  T2["abmi"];
           Z2["abij"] +=     WMNIJ["mnij"]*  R2["abmn"];
           Z2["abij"] +=     XMNIJ["mnij"]* Tau["abmn"];
           Z2["abij"] +=     VABCD["abef"]*  R2["efij"];
           Z2["abij"] +=   WAMIESA["bmje"]*  R2["aeim"];
           Z2["abij"] -=     WAMIE["amje"]*  R2["beim"];
           Z2["abij"] -=     WAMIE["bmie"]*  R2["aejm"];
           Z2["abij"] -=     WAMEI["amei"]*  R2["ebmj"];
           Z2["abij"] -=     WAMEI["amej"]*  R2["ebim"];
           Z2["abij"] -=     WAMEI["bmei"]*  R2["aemj"];
    }
    else
    {
         R2SA["abij"]  =    2*    R2["abij"];
         R2SA["abij"] -=          R2["abji"];

          XMI[  "mi"]  =       VABIJ["efmn"]*R2SA["efin"];

           Z1[  "ai"]  =         FAE[  "ae"]*  R1[  "ei"];
           Z1[  "ai"] -=         FMI[  "mi"]*  R1[  "am"];
           Z1[  "ai"] -=         XMI[  "mi"]*  T1[  "am"];
           Z1[  "ai"] +=     WAMIESA["amie"]*  R1[  "em"];
           Z1[  "ai"] +=         FME[  "me"]*R2SA["aeim"];
           Z1[  "ai"] +=       VABCI["efam"]*R2SA["efim"];
           Z1[  "ai"] -=       WMNEJ["mnei"]*R2SA["eamn"];

          XMI[  "mi"] +=   2*  WMNEJ["nmei"]*  R1[  "en"];
          XMI[  "mi"] -=       WMNEJ["mnei"]*  R1[  "en"];

          XME[  "me"]  =   2*  VABIJ["efmn"]*  R1[  "fn"];
          XME[  "me"] -=       VABIJ["efnm"]*  R1[  "fn"];

          XAE[  "ae"]  =   2*  VABCI["efam"]*  R1[  "fm"];
          XAE[  "ae"] -=       VABCI["feam"]*  R1[  "fm"];
          XAE[  "ae"] -=         XME[  "me"]*  T1[  "am"];
          XAE[  "ae"] -=       VABIJ["efmn"]*R2SA["afmn"];

        XAMIJ["amij"]  =       VABCI["efam"]*  R2["efij"];

        XMNIJ["mnij"]  =       VABIJ["efmn"]*  R2["efij"];

           Z2["abij"]  =       WABEJ["abej"]*  R1[  "ei"];
           Z2["abij"] -=       WAMIJ["amij"]*  R1[  "bm"];
           Z2["abij"] -=       XAMIJ["amij"]*  T1[  "bm"];
           Z2["abij"] +=         FAE[  "ae"]*  R2["ebij"];
           Z2["abij"] -=         FMI[  "mi"]*  R2["abmj"];
           Z2["abij"] +=         XAE[  "ae"]*  T2["ebij"];
           Z2["abij"] -=         XMI[  "mi"]*  T2["abmj"];
           Z2["abij"] += 0.5*  WMNIJ["mnij"]*  R2["abmn"];
           Z2["abij"] += 0.5*  XMNIJ["mnij"]* Tau["abmn"];
           Z2["abij"] += 0.5*  VABCD["abef"]*  R2["efij"];
           Z2["abij"] += 0.5*WAMIESA["amie"]*R2SA["ebmj"];
           Z2["abij"] -= 0.5*  WAMEI["amei"]*  R2["ebjm"];
           Z2["abij"] -=       WAMEI["amej"]*  R2["ebim"];

           Z2["abij"] +=          Z2["baji"];
    }

    this->energy() = davidson.extrapolate(ptr_vector<SymmetryBlockedTensor<U>>{&R1, &R2},
                                          ptr_vector<SymmetryBlockedTensor<U>>{&Z1, &Z2}, D)[0];
    this->conv() = max(Z1.norm(00), Z2.norm(00));
}

}
}

static const char* spec = R"(

nroot?
    int 1,
ntriplet?
    int 0,
nsinglet?
    int 0,
convergence?
    double 1e-9,
max_iterations?
    int 150,
sub_iterations?
    int 0,
conv_type?
    enum { MAXE, RMSE, MAE },
davidson?
{
    order?
            int 10,
    num_reduce?
            int 3,
    compaction?
            enum { discrete, continuous },
}

)";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::RHFEOMEECCSD);
REGISTER_TASK(aquarius::cc::RHFEOMEECCSD<double>, "rhfeomeeccsd",spec);
