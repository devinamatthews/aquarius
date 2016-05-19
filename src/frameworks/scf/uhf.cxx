#include "uhf.hpp"

using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;
using namespace aquarius::time;
using namespace aquarius::op;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace scf
{

template <typename T>
UHF<T>::UHF(const string& name, Config& config)
: Iterative<T>(name, config), frozen_core(config.get<bool>("frozen_core")),
  diis(config.get("diis"), 2)
{
    vector<Requirement> reqs;
    reqs += Requirement("molecule", "molecule");
    reqs += Requirement("ovi", "S");
    reqs += Requirement("1ehamiltonian", "H");
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("double", "S2", reqs));
    this->addProduct(Product("double", "multiplicity", reqs));
    this->addProduct(Product("occspace", "occ", reqs));
    this->addProduct(Product("vrtspace", "vrt", reqs));
    this->addProduct(Product("Ea", "Ea", reqs));
    this->addProduct(Product("Eb", "Eb", reqs));
    this->addProduct(Product("Fa", "Fa", reqs));
    this->addProduct(Product("Fb", "Fb", reqs));
    this->addProduct(Product("Da", "Da", reqs));
    this->addProduct(Product("Db", "Db", reqs));
}

template <typename T>
bool UHF<T>::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = this->template get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();

    vector<int> shapeNN = {NS,NS};
    vector<vector<int>> sizenn = {norb,norb};

    this->put("Fa", new SymmetryBlockedTensor<T>("Fa", arena, group, 2, sizenn, shapeNN, false));
    this->put("Fb", new SymmetryBlockedTensor<T>("Fb", arena, group, 2, sizenn, shapeNN, false));
    this->put("Da", new SymmetryBlockedTensor<T>("Da", arena, group, 2, sizenn, shapeNN, true));
    this->put("Db", new SymmetryBlockedTensor<T>("Db", arena, group, 2, sizenn, shapeNN, true));

    this->puttmp("dF",     new SymmetryBlockedTensor<T>("dF",     arena, group, 2, sizenn, shapeNN, false));
    this->puttmp("Ca",     new SymmetryBlockedTensor<T>("Ca",     arena, group, 2, sizenn, shapeNN, false));
    this->puttmp("Cb",     new SymmetryBlockedTensor<T>("Cb",     arena, group, 2, sizenn, shapeNN, false));
    this->puttmp("dDa",    new SymmetryBlockedTensor<T>("dDa",    arena, group, 2, sizenn, shapeNN, false));
    this->puttmp("dDb",    new SymmetryBlockedTensor<T>("dDb",    arena, group, 2, sizenn, shapeNN, false));
    this->puttmp("S^-1/2", new SymmetryBlockedTensor<T>("S^-1/2", arena, group, 2, sizenn, shapeNN, false));

    occ_alpha.resize(group.getNumIrreps());
    occ_beta.resize(group.getNumIrreps());

    E_alpha.resize(group.getNumIrreps());
    E_beta.resize(group.getNumIrreps());

    for (int i = 0;i < group.getNumIrreps();i++)
    {
        E_alpha[i].resize(norb[i]);
        E_beta[i].resize(norb[i]);
    }

    calcSMinusHalf();

    CTF_Timer_epoch ep(this->name.c_str());
    ep.begin();
    Iterative<T>::run(dag, arena);
    ep.end();

    if (this->isUsed("S2") || this->isUsed("multiplicity"))
    {
        calcS2();
    }

    this->put("energy", new T(this->energy()));
    this->put("convergence", new T(this->conv()));

    int nfrozen = 0;
    if (frozen_core)
    {
        for (vector<Atom>::const_iterator a = molecule.getAtomsBegin();a != molecule.getAtomsEnd();++a)
        {
            int Z = a->getCenter().getElement().getAtomicNumber();
            if      (Z > 86) nfrozen += 31;
            else if (Z > 54) nfrozen += 22;
            else if (Z > 36) nfrozen += 13;
            else if (Z > 18) nfrozen += 9;
            else if (Z > 10) nfrozen += 5;
            else if (Z >  2) nfrozen += 1;
        }
    }

    if (nfrozen > nalpha || nfrozen > nbeta)
        Logger::error(arena) << "There are not enough valence electrons for this multiplicity" << endl;

    vector<pair<real_type_t<T>,int>> E_alpha_occ;
    vector<pair<real_type_t<T>,int>> E_beta_occ;
    for (int i = 0;i < group.getNumIrreps();i++)
    {
        for (int j = 0;j < occ_alpha[i];j++)
        {
            E_alpha_occ.push_back(make_pair(E_alpha[i][j],i));
        }
        for (int j = 0;j < occ_beta[i];j++)
        {
            E_beta_occ.push_back(make_pair(E_beta[i][j],i));
        }
    }
    assert(E_alpha_occ.size() == nalpha);
    assert(E_beta_occ.size() == nbeta);

    sort(E_alpha_occ.begin(), E_alpha_occ.end());
    sort(E_beta_occ.begin(), E_beta_occ.end());

    vector<int> nfrozen_alpha(group.getNumIrreps());
    vector<int> nfrozen_beta(group.getNumIrreps());
    for (int i = 0;i < nfrozen;i++)
    {
        nfrozen_alpha[E_alpha_occ[i].second]++;
        nfrozen_beta[E_beta_occ[i].second]++;
    }

    Logger::log(arena) << "Dropping MOs: " << nfrozen_alpha << ", " << nfrozen_beta << endl;

    vector<int> vrt_alpha(group.getNumIrreps());
    vector<int> vrt_beta(group.getNumIrreps());
    vector<int> vrt0_alpha(group.getNumIrreps());
    vector<int> vrt0_beta(group.getNumIrreps());
    for (int i = 0;i < group.getNumIrreps();i++)
    {
        vrt_alpha[i] = norb[i]-occ_alpha[i];
        vrt_beta[i] = norb[i]-occ_beta[i];
        vrt0_alpha[i] = occ_alpha[i];
        vrt0_beta[i] = occ_beta[i];
        occ_alpha[i] -= nfrozen_alpha[i];
        occ_beta[i] -= nfrozen_beta[i];
    }

    vector<int> zero(norb.size(), 0);
    this->put("occ", new MOSpace<T>(SymmetryBlockedTensor<T>("CI", this->template gettmp<SymmetryBlockedTensor<T>>("Ca"),
                                                             {zero,nfrozen_alpha},
                                                             {norb,occ_alpha}),
                                    SymmetryBlockedTensor<T>("Ci", this->template gettmp<SymmetryBlockedTensor<T>>("Cb"),
                                                             {zero,nfrozen_beta},
                                                             {norb,occ_beta})));

   	this->put("vrt", new MOSpace<T>(SymmetryBlockedTensor<T>("CA", this->template gettmp<SymmetryBlockedTensor<T>>("Ca"),
                                                             {zero,vrt0_alpha},
                                                             {norb,vrt_alpha}),
   	                                SymmetryBlockedTensor<T>("Ca", this->template gettmp<SymmetryBlockedTensor<T>>("Cb"),
                                                             {zero,vrt0_beta},
                                                             {norb,vrt_beta})));

    vector<int> shapeN{NS};
    vector<vector<int>> sizena{norb};
    vector<vector<int>> sizenb{norb};
    for (int i = 0;i < group.getNumIrreps();i++) sizena[0][i] -= nfrozen_alpha[i];
    for (int i = 0;i < group.getNumIrreps();i++) sizenb[0][i] -= nfrozen_beta[i];

    auto& Ea = this->put("Ea", new vector<vector<real_type_t<T>>>(group.getNumIrreps()));
    auto& Eb = this->put("Eb", new vector<vector<real_type_t<T>>>(group.getNumIrreps()));

    for (int i = 0;i < group.getNumIrreps();i++)
    {
        sort(E_alpha[i].begin(), E_alpha[i].end());
        Ea[i].assign(E_alpha[i].begin()+nfrozen_alpha[i], E_alpha[i].end());
    }

    for (int i = 0;i < group.getNumIrreps();i++)
    {
        sort(E_beta[i].begin(), E_beta[i].end());
        Eb[i].assign(E_beta[i].begin()+nfrozen_beta[i], E_beta[i].end());
    }

    return true;
}

template <typename T>
void UHF<T>::iterate(const Arena& arena)
{
    const Molecule& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();
    int norbtot = sum(norb);
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();

    buildFock();
    DIISExtrap();
    calcEnergy();
    diagonalizeFock();

    vector<pair<real_type_t<T>,int>> E_alpha_sorted;
    vector<pair<real_type_t<T>,int>> E_beta_sorted;
    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        occ_alpha[i] = 0;
        occ_beta[i] = 0;
        for (int j = 0;j < norb[i];j++)
        {
            E_alpha_sorted.push_back(make_pair(E_alpha[i][j],i));
        }
        for (int j = 0;j < norb[i];j++)
        {
            E_beta_sorted.push_back(make_pair(E_beta[i][j],i));
        }
    }

    sort(E_alpha_sorted.begin(), E_alpha_sorted.end());
    sort(E_beta_sorted.begin(), E_beta_sorted.end());

    for (int i = 0;i < nalpha;i++)
    {
        occ_alpha[E_alpha_sorted[i].second]++;
    }
    for (int i = 0;i < nbeta;i++)
    {
        occ_beta[E_beta_sorted[i].second]++;
    }

    Logger::log(arena) << "Iteration " << this->iter() << " occupation = " << occ_alpha << ", " << occ_beta << endl;

    calcDensity();

    auto& dDa = this->template gettmp<SymmetryBlockedTensor<T>>("dDa");
    auto& dDb = this->template gettmp<SymmetryBlockedTensor<T>>("dDb");

    switch (this->convtype)
    {
        case Iterative<T>::MAX_ABS:
            this->conv() = max(dDa.norm(00), dDb.norm(00));
            break;
        case Iterative<T>::RMSD:
            this->conv() = (dDa.norm(2)+dDb.norm(2))/sqrt(2*norbtot*norbtot);
            break;
        case Iterative<T>::MAD:
            this->conv() = (dDa.norm(1)+dDb.norm(1))/(2*norbtot*norbtot);
            break;
    }
}

template <typename T>
void UHF<T>::calcS2()
{
    const auto& molecule = this->template get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();

    auto& S = this->template get<SymmetryBlockedTensor<T>>("S");

    vector<int> zero(norb.size(), 0);
    SymmetryBlockedTensor<T> Ca_occ("CI", this->template gettmp<SymmetryBlockedTensor<T>>("Ca"),
                                    {zero,zero}, {norb,occ_alpha});
    SymmetryBlockedTensor<T> Cb_occ("Ci", this->template gettmp<SymmetryBlockedTensor<T>>("Cb"),
                                    {zero,zero}, {norb,occ_beta});

    SymmetryBlockedTensor<T> Delta("Delta", S.arena, group, 2, {{nalpha},{nbeta}}, {NS,NS}, false);
    SymmetryBlockedTensor<T> tmp("tmp", S.arena, group, 2, {{nalpha},norb}, {NS,NS}, false);

    int ndiff = abs(nalpha-nbeta);
    int nmin = min(nalpha, nbeta);

    double S2 = ((ndiff/2)*(ndiff/2+1) + nmin);

    tmp["ai"] = Ca_occ["ja"]*S["ij"];
    Delta["ab"] = tmp["ai"]*Cb_occ["ib"];

    S2 -= aquarius::abs(scalar(Delta*conj(Delta)));

    this->put("S2", new T(S2));
    this->put("multiplicity", new T(sqrt(4*S2+1)));
}

template <typename T>
void UHF<T>::calcEnergy()
{
    const Molecule& molecule = this->template get<Molecule>("molecule");

    auto& H  = this->template get<SymmetryBlockedTensor<T>>("H");
    auto& Fa = this->template get<SymmetryBlockedTensor<T>>("Fa");
    auto& Fb = this->template get<SymmetryBlockedTensor<T>>("Fb");
    auto& Da = this->template get<SymmetryBlockedTensor<T>>("Da");
    auto& Db = this->template get<SymmetryBlockedTensor<T>>("Db");

    /*
     * E = (1/2)Tr[D(F+H)]
     *
     *   = (1/2)Tr[Da*(Fa+H) + Db*(Fb+H)]
     */
    Fa["ab"] += H["ab"];
    Fb["ab"] += H["ab"];
    this->energy()  = molecule.getNuclearRepulsion();
    this->energy() += 0.5*scalar(Da["ab"]*Fa["ab"]);
    this->energy() += 0.5*scalar(Db["ab"]*Fb["ab"]);
    Fa["ab"] -= H["ab"];
    Fb["ab"] -= H["ab"];
}

template <typename T>
void UHF<T>::calcDensity()
{
    const Molecule& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();

    auto& dDa = this->template gettmp<SymmetryBlockedTensor<T>>("dDa");
    auto& dDb = this->template gettmp<SymmetryBlockedTensor<T>>("dDb");
    auto& Da  = this->template get   <SymmetryBlockedTensor<T>>("Da");
    auto& Db  = this->template get   <SymmetryBlockedTensor<T>>("Db");

    vector<int> zero(norb.size(), 0);
    SymmetryBlockedTensor<T> Ca_occ("CI", this->template gettmp<SymmetryBlockedTensor<T>>("Ca"),
                                    {zero,zero}, {norb,occ_alpha});
    SymmetryBlockedTensor<T> Cb_occ("Ci", this->template gettmp<SymmetryBlockedTensor<T>>("Cb"),
                                    {zero,zero}, {norb,occ_beta});

    /*
     * D[ab] = C[ai]*C[bi]
     */
    dDa["ab"]  = Da["ab"];
    dDb["ab"]  = Db["ab"];
     Da["ab"]  = Ca_occ["ai"]*Ca_occ["bi"];
     Db["ab"]  = Cb_occ["ai"]*Cb_occ["bi"];
    dDa["ab"] -= Da["ab"];
    dDb["ab"] -= Db["ab"];
}

template <typename T>
void UHF<T>::DIISExtrap()
{
    auto& S      = this->template get   <SymmetryBlockedTensor<T>>("S");
    auto& Smhalf = this->template gettmp<SymmetryBlockedTensor<T>>("S^-1/2");
    auto& dF     = this->template gettmp<SymmetryBlockedTensor<T>>("dF");
    auto& Fa     = this->template get   <SymmetryBlockedTensor<T>>("Fa");
    auto& Fb     = this->template get   <SymmetryBlockedTensor<T>>("Fb");
    auto& Da     = this->template get   <SymmetryBlockedTensor<T>>("Da");
    auto& Db     = this->template get   <SymmetryBlockedTensor<T>>("Db");

    /*
     * Generate the residual:
     *
     *   R = FDS - SDF
     *
     * Then, convert to the orthonormal basis:
     *
     *   ~    -1/2    -1/2
     *   R = S     R S.
     *
     * In this basis we have
     *
     *   ~    -1/2    -1/2  ~    1/2
     *   F = S     F S    , C = S    C, and
     *
     *   ~   ~ ~T    1/2    T  1/2    1/2    1/2
     *   D = C C  = S    C C  S    = S    D S.
     *
     * And so,
     *
     *   ~    -1/2    -1/2  1/2    1/2    1/2    1/2  -1/2    -1/2
     *   R = S     F S     S    D S    - S    D S    S     F S
     *
     *        ~ ~
     *     = [F,D] = 0 at convergence.
     */
    {
        SymmetryBlockedTensor<T> tmp1("tmp", Fa);
        SymmetryBlockedTensor<T> tmp2("tmp", Fa);

        tmp1["ab"]  =     Fa["ac"]*    Da["cb"];
        tmp2["ab"]  =   tmp1["ac"]*     S["cb"];
        tmp1["ab"]  =      S["ac"]*    Da["cb"];
        tmp2["ab"] -=   tmp1["ac"]*    Fa["cb"];
        tmp1["ab"]  = Smhalf["ac"]*  tmp2["cb"];
          dF["ab"]  =   tmp1["ac"]*Smhalf["cb"];

        tmp1["ab"]  =     Fb["ac"]*    Db["cb"];
        tmp2["ab"]  =   tmp1["ac"]*     S["cb"];
        tmp1["ab"]  =      S["ac"]*    Db["cb"];
        tmp2["ab"] -=   tmp1["ac"]*    Fb["cb"];
        tmp1["ab"]  = Smhalf["ac"]*  tmp2["cb"];
          dF["ab"] +=   tmp1["ac"]*Smhalf["cb"];
    }

    diis.extrapolate(ptr_vector<SymmetryBlockedTensor<T>>{&Fa, &Fb},
                     ptr_vector<SymmetryBlockedTensor<T>>{&dF});
}

INSTANTIATE_SPECIALIZATIONS(UHF);

}
}
