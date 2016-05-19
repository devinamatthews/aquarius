#include "../../frameworks/scf/cfourscf.hpp"

using namespace aquarius::task;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::op;

namespace aquarius
{
namespace scf
{

template <typename T>
CFOURSCF<T>::CFOURSCF(const string& name, Config& config)
: Task(name, config),
  frozen_core(config.get<bool>("frozen_core")),
  semicanonical(config.get<bool>("semicanonical"))
{
    vector<Requirement> reqs;
    reqs += Requirement("molecule", "molecule");
    this->addProduct(Product("double", "energy", reqs));
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
bool CFOURSCF<T>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& molecule = get<Molecule>("molecule");
    const auto& group = molecule.getGroup();

    vector<int> norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();
    int nirrep = norb.size();
    bool rhf = nalpha == nbeta;

    vector<int> occ_alpha(nirrep);
    vector<int> occ_beta(nirrep);

    vector<vector<double>> E_alpha(nirrep);
    vector<vector<double>> E_beta(nirrep);

    for (int i = 0;i < nirrep;i++)
    {
        E_alpha[i].resize(norb[i]);
        E_beta[i].resize(norb[i]);
    }

    vector<int> shapeNN = {NS,NS};
    vector<vector<int>> sizenn = {norb,norb};

    auto& Fa = this->put   ("Fa", new SymmetryBlockedTensor<T>("Fa", arena, group, 2, sizenn, shapeNN, true));
    auto& Fb = this->put   ("Fb", new SymmetryBlockedTensor<T>("Fb", arena, group, 2, sizenn, shapeNN, true));
    auto& Da = this->put   ("Da", new SymmetryBlockedTensor<T>("Da", arena, group, 2, sizenn, shapeNN, true));
    auto& Db = this->put   ("Db", new SymmetryBlockedTensor<T>("Db", arena, group, 2, sizenn, shapeNN, true));
    auto& Ca = this->puttmp("Ca", new SymmetryBlockedTensor<T>("Ca", arena, group, 2, sizenn, shapeNN, true));
    auto& Cb = this->puttmp("Cb", new SymmetryBlockedTensor<T>("Cb", arena, group, 2, sizenn, shapeNN, true));

    ifstream ifsfock("NEWFOCK");
    ifstream ifsmo("NEWMOS");

    double energy;
    ifsfock >> energy;
    this->put("energy", new T(energy));

    for (int i = 0;i < nirrep;i++) ifsfock >> occ_alpha[i];
    for (int i = 0;i < nirrep;i++) ifsfock >> occ_beta[i];

    for (int spin = 0;spin < (rhf ? 1 : 2);spin++)
    {
        auto& F = (spin == 0 ? Fa : Fb);
        auto& D = (spin == 0 ? Da : Db);
        auto& C = (spin == 0 ? Ca : Cb);
        auto& E = (spin == 0 ? E_alpha : E_beta);
        auto& occ = (spin == 0 ? occ_alpha : occ_beta);

        vector<int> vrt(norb);

        for (int i = 0;i < nirrep;i++)
        {
            vrt[i] -= occ[i];

            matrix<double> fock(norb[i], norb[i]);
            matrix<double>  fmo(norb[i], norb[i]);
            matrix<double>   mo(norb[i], norb[i]);
            matrix<double> dens(norb[i], norb[i]);
            matrix<double>  tmp(norb[i], norb[i]);

            for (int r0 = 0;r0 < norb[i];r0 += 4)
            {
                for (int c = 0;c < norb[i];c++)
                {
                    for (int r = r0;r < min(r0+4,norb[i]);r++)
                    {
                        ifsfock >> fock[c][r];
                        ifsmo >> mo[c][r];
                    }
                }
            }

            gemm('T', 'N', norb[i], norb[i], occ[i],
                 1.0,   mo.data(), norb[i],
                        mo.data(), norb[i],
                 0.0, dens.data(), norb[i]);

            gemm('N', 'T', norb[i], norb[i], norb[i],
                 1.0, fock.data(), norb[i],
                        mo.data(), norb[i],
                 0.0,  tmp.data(), norb[i]);
            gemm('N', 'N', norb[i], norb[i], norb[i],
                 1.0,   mo.data(), norb[i],
                       tmp.data(), norb[i],
                 0.0,  fmo.data(), norb[i]);

            copy(norb[i], fmo.data(), norb[i]+1, E[i].data(), 1);

            if (semicanonical)
            {
                tmp = fmo;

                vector<double> E_tmp(norb[i]);
                vector<int> order = range(norb[i]);

                heev('V', 'U', occ[i], &tmp[     0][     0], norb[i], &E_tmp[     0]);
                heev('V', 'U', vrt[i], &tmp[occ[i]][occ[i]], norb[i], &E_tmp[occ[i]]);

                matrix<double> tmp2(norb[i], norb[i]);

                gemm('T', 'N', occ[i], norb[i], occ[i],
                     1.0,  &tmp[     0][     0], norb[i],
                            &mo[     0][     0], norb[i],
                     0.0, &tmp2[     0][     0], norb[i]);
                gemm('T', 'N', vrt[i], norb[i], vrt[i],
                     1.0,  &tmp[occ[i]][occ[i]], norb[i],
                            &mo[     0][occ[i]], norb[i],
                     0.0, &tmp2[     0][occ[i]], norb[i]);

                cosort(E_tmp.begin(), E_tmp.begin()+occ[i],
                       order.begin(), order.begin()+occ[i]);
                cosort(E_tmp.begin()+occ[i], E_tmp.end(),
                       order.begin()+occ[i], order.end());

                for (int j = 0;j < norb[i];j++)
                {
                    E[i][j] = E_tmp[order[j]];
                    copy(norb[i], &tmp2[0][j], norb[i], &mo[0][j], norb[i]);
                }
            }

            vector<tkv_pair<T>> fpairs;
            vector<tkv_pair<T>> cpairs;
            vector<tkv_pair<T>> dpairs;

            for (int r = 0;r < norb[i];r++)
            {
                for (int c = 0;c < norb[i];c++)
                {
                    fpairs.emplace_back(r+c*norb[i], fock[r][c]);
                    cpairs.emplace_back(r+c*norb[i], mo[r][c]);
                    dpairs.emplace_back(r+c*norb[i], dens[r][c]);
                }
            }

            if (arena.rank == 0)
            {
                F({i,i}).writeRemoteData(fpairs);
                C({i,i}).writeRemoteData(cpairs);
                D({i,i}).writeRemoteData(dpairs);
            }
            else
            {
                F({i,i}).writeRemoteData();
                C({i,i}).writeRemoteData();
                D({i,i}).writeRemoteData();
            }
        }
    }

    if (rhf)
    {
        Fb = Fa;
        Cb = Ca;
        Db = Da;
        E_beta = E_alpha;
        occ_beta = occ_alpha;
    }

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
            E_alpha_occ.emplace_back(E_alpha[i][j], i);
        }
        for (int j = 0;j < occ_beta[i];j++)
        {
            E_beta_occ.emplace_back(E_beta[i][j], i);
        }
    }
    assert(E_alpha_occ.size() == nalpha);
    assert(E_beta_occ.size() == nbeta);

    sort(E_alpha_occ);
    sort(E_beta_occ);

    vector<int> nfrozen_alpha(nirrep);
    vector<int> nfrozen_beta(nirrep);
    for (int i = 0;i < nfrozen;i++)
    {
        nfrozen_alpha[E_alpha_occ[i].second]++;
        nfrozen_beta[E_beta_occ[i].second]++;
    }

    Logger::log(arena) << "Dropping MOs: " << nfrozen_alpha << ", " << nfrozen_beta << endl;

    vector<int> vrt_alpha(nirrep);
    vector<int> vrt_beta(nirrep);
    vector<int> vrt0_alpha(nirrep);
    vector<int> vrt0_beta(nirrep);
    for (int i = 0;i < nirrep;i++)
    {
        vrt_alpha[i] = norb[i]-occ_alpha[i];
        vrt_beta[i] = norb[i]-occ_beta[i];
        vrt0_alpha[i] = occ_alpha[i];
        vrt0_beta[i] = occ_beta[i];
        occ_alpha[i] -= nfrozen_alpha[i];
        occ_beta[i] -= nfrozen_beta[i];
    }

    vector<int> zero(norb.size(), 0);
    auto& occ =
        this->put("occ", new MOSpace<T>(SymmetryBlockedTensor<T>("CI", this->template gettmp<SymmetryBlockedTensor<T>>("Ca"),
                                                                 {zero,nfrozen_alpha},
                                                                 {norb,occ_alpha}),
                                        SymmetryBlockedTensor<T>("Ci", this->template gettmp<SymmetryBlockedTensor<T>>("Cb"),
                                                                 {zero,nfrozen_beta},
                                                                 {norb,occ_beta})));

    auto& vrt =
        this->put("vrt", new MOSpace<T>(SymmetryBlockedTensor<T>("CA", this->template gettmp<SymmetryBlockedTensor<T>>("Ca"),
                                                                 {zero,vrt0_alpha},
                                                                 {norb,vrt_alpha}),
                                        SymmetryBlockedTensor<T>("Ca", this->template gettmp<SymmetryBlockedTensor<T>>("Cb"),
                                                                 {zero,vrt0_beta},
                                                                 {norb,vrt_beta})));

    auto& Ea = this->put("Ea", new vector<vector<real_type_t<T>>>(nirrep));
    auto& Eb = this->put("Eb", new vector<vector<real_type_t<T>>>(nirrep));

    for (int i = 0;i < nirrep;i++)
    {
        sort(E_alpha[i].begin(), E_alpha[i].end());
        Ea[i].assign(E_alpha[i].begin()+nfrozen_alpha[i], E_alpha[i].end());
    }

    for (int i = 0;i < nirrep;i++)
    {
        sort(E_beta[i].begin(), E_beta[i].end());
        Eb[i].assign(E_beta[i].begin()+nfrozen_beta[i], E_beta[i].end());
    }

    return true;
}

}
}

static const char* spec = R"(

    frozen_core?
        bool false,
    semicanonical?
        bool true

)";

INSTANTIATE_SPECIALIZATIONS(aquarius::scf::CFOURSCF);
REGISTER_TASK(CONCAT(aquarius::scf::CFOURSCF<double>), "cfourscf",spec);
