#include "uhf.hpp"

using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;
using namespace aquarius::time;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace scf
{

UHF::UHF(const string& name, Config& config)
: Iterative(name, config), frozen_core(config.get<bool>("frozen_core")),
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

bool UHF::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = this->template get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();

    auto init = TensorInitializer<>("Ca", Field::DOUBLE) <<
                TensorInitializer<PGSYMMETRIC|BOUNDED>(group, {norb,norb}) <<
                TensorInitializer<DISTRIBUTED>(arena);
    Tensor<BOUNDED|PGSYMMETRIC> Ca = this->puttmp<Tensor<>>("Ca", Tensor<PGSYMMETRIC|BOUNDED|DISTRIBUTED>::construct(init));
    Tensor<BOUNDED|PGSYMMETRIC> Cb = this->puttmp<Tensor<>>("Cb", Ca.construct("Cb"));

    this->puttmp<Tensor<>>("dF",     Ca.construct("dF"));
    this->puttmp<Tensor<>>("dDa",    Ca.construct("dDa"));
    this->puttmp<Tensor<>>("dDb",    Ca.construct("dDb"));
    this->puttmp<Tensor<>>("S^-1/2", Ca.construct("S^-1/2"));

    this->put<Tensor<>>("Fa", Ca.construct("Fa"));
    this->put<Tensor<>>("Fb", Ca.construct("Fb"));
    this->put<Tensor<>>("Da", Ca.construct("Da"));
    this->put<Tensor<>>("Db", Ca.construct("Db"));

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
    Iterative::run(dag, arena);
    ep.end();

    this->put("energy", this->energy());
    this->put("convergence", this->conv());

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

    vector<pair<Scalar,int>> E_alpha_occ;
    vector<pair<Scalar,int>> E_beta_occ;
    for (int i = 0;i < group.getNumIrreps();i++)
    {
        for (int j = 0;j < occ_alpha[i];j++)
        {
            E_alpha_occ.emplace_back(E_alpha[i][j], i);
        }
        for (int j = 0;j < occ_beta[i];j++)
        {
            E_beta_occ.emplace_back(E_beta[i][j], i);
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

    Tensor<BOUNDED|PGSYMMETRIC> Ca_occ = this->put<Tensor<>>("CI", Ca.construct("CI", TensorInitializer<PGSYMMETRIC|BOUNDED>(group, {norb,occ_alpha})));
    Tensor<BOUNDED|PGSYMMETRIC> Cb_occ = this->put<Tensor<>>("Ci", Ca.construct("Ci", TensorInitializer<PGSYMMETRIC|BOUNDED>(group, {norb,occ_beta})));
    Tensor<BOUNDED|PGSYMMETRIC> Ca_vrt = this->put<Tensor<>>("CA", Ca.construct("CA", TensorInitializer<PGSYMMETRIC|BOUNDED>(group, {norb,vrt_alpha})));
    Tensor<BOUNDED|PGSYMMETRIC> Cb_vrt = this->put<Tensor<>>("Ca", Ca.construct("Ca", TensorInitializer<PGSYMMETRIC|BOUNDED>(group, {norb,vrt_beta})));

    vector<int> zero(group.getNumIrreps());
    Ca_occ.slice({zero,zero}, Ca);
    Cb_occ.slice({zero,zero}, Cb);
    Ca_vrt.slice({zero,vrt0_alpha}, Ca);
    Cb_vrt.slice({zero,vrt0_beta}, Cb);

    auto& Ea = this->put("Ea", new vector<vector<Scalar>>(group.getNumIrreps()));
    auto& Eb = this->put("Eb", new vector<vector<Scalar>>(group.getNumIrreps()));

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

    if (this->isUsed("S2") || this->isUsed("multiplicity"))
    {
        calcS2();
    }

    return true;
}

void UHF::iterate(const Arena& arena)
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

    vector<pair<Scalar,int>> E_alpha_sorted;
    vector<pair<Scalar,int>> E_beta_sorted;
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

    Tensor<BOUNDED> dDa = this->template gettmp<Tensor<>>("dDa");
    Tensor<BOUNDED> dDb = this->template gettmp<Tensor<>>("dDb");

    switch (this->convtype)
    {
        case Iterative::MAX_ABS:
            this->conv() = max(dDa.norm(00), dDb.norm(00)).to<double>();
            break;
        case Iterative::RMSD:
            this->conv() = (dDa.norm(2)+dDb.norm(2)).to<double>()/sqrt(2*norbtot*norbtot);
            break;
        case Iterative::MAD:
            this->conv() = (dDa.norm(1)+dDb.norm(1)).to<double>()/(2*norbtot*norbtot);
            break;
    }
}

void UHF::calcS2()
{
    const auto& molecule = this->template get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();

    const vector<int>& norb = molecule.getNumOrbitals();

    Tensor<BOUNDED> S = this->template get<Tensor<>>("S");

    Tensor<BOUNDED|PGSYMMETRIC> Ca_occ = this->get<Tensor<>>("CI");
    Tensor<BOUNDED|PGSYMMETRIC> Cb_occ = this->get<Tensor<>>("Ci");

    vector<int> nalpha = Ca_occ.getLengthsPerIrrep()[1];
    vector<int> nbeta = Cb_occ.getLengthsPerIrrep()[1];

    Tensor<BOUNDED> Delta = S.construct("Delta", TensorInitializer<BOUNDED|PGSYMMETRIC>(group, {nalpha,nbeta}));
    Tensor<BOUNDED> tmp = Cb_occ.construct("tmp");

    int ndiff = abs(sum(nalpha) - sum(nbeta));
    int nmin = min(sum(nalpha), sum(nbeta));

    Scalar S2 = double((ndiff/2)*(ndiff/2+1) + nmin);

    tmp["ib"] = S["ij"]*Cb_occ["jb"];
    Delta["ab"] = tmp["ib"]*Ca_occ["ia"];

    S2 -= abs(scalar(Delta*conj(Delta)));

    this->put("S^2", S2);
    this->put("multiplicity", sqrt(4*S2+1));
}

void UHF::calcEnergy()
{
    const Molecule& molecule = this->template get<Molecule>("molecule");

    Tensor<BOUNDED> H  = this->template get<Tensor<>>("H");
    Tensor<BOUNDED> Fa = this->template get<Tensor<>>("Fa");
    Tensor<BOUNDED> Fb = this->template get<Tensor<>>("Fb");
    Tensor<BOUNDED> Da = this->template get<Tensor<>>("Da");
    Tensor<BOUNDED> Db = this->template get<Tensor<>>("Db");

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

void UHF::calcDensity()
{
    const Molecule& molecule = this->template get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();
    const vector<int>& norb = molecule.getNumOrbitals();

    Tensor<BOUNDED> dDa = this->template gettmp<Tensor<>>("dDa");
    Tensor<BOUNDED> dDb = this->template gettmp<Tensor<>>("dDb");
    Tensor<BOUNDED> Da  = this->template get   <Tensor<>>("Da");
    Tensor<BOUNDED> Db  = this->template get   <Tensor<>>("Db");

    Tensor<BOUNDED|PGSYMMETRIC> Ca  = this->template gettmp<Tensor<>>("Ca");
    Tensor<BOUNDED|PGSYMMETRIC> Cb  = this->template gettmp<Tensor<>>("Cb");
    Tensor<BOUNDED|PGSYMMETRIC> Ca_occ = Ca.construct("CI", TensorInitializer<PGSYMMETRIC|BOUNDED>(group, {norb,occ_alpha}));
    Tensor<BOUNDED|PGSYMMETRIC> Cb_occ = Ca.construct("Ci", TensorInitializer<PGSYMMETRIC|BOUNDED>(group, {norb,occ_beta}));

    vector<int> zero(group.getNumIrreps());
    Ca_occ.slice({zero,zero}, Ca);
    Cb_occ.slice({zero,zero}, Cb);

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

void UHF::DIISExtrap()
{
    Tensor<BOUNDED> S      = this->template get   <Tensor<>>("S");
    Tensor<BOUNDED> Smhalf = this->template gettmp<Tensor<>>("S^-1/2");
    Tensor<BOUNDED> dF     = this->template gettmp<Tensor<>>("dF");
    Tensor<BOUNDED> Fa     = this->template get   <Tensor<>>("Fa");
    Tensor<BOUNDED> Fb     = this->template get   <Tensor<>>("Fb");
    Tensor<BOUNDED> Da     = this->template get   <Tensor<>>("Da");
    Tensor<BOUNDED> Db     = this->template get   <Tensor<>>("Db");

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
        auto tmp1 = Fa.construct("tmp");
        auto tmp2 = Fa.construct("tmp");

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

    diis.extrapolate({Fa, Fb}, {dF});
}

}
}
