#include "scf.hpp"

using namespace aquarius::task;
using namespace aquarius::input;
using namespace aquarius::tensor;

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
                        ifsfock >> fock[r][c];
                        ifsmo >> mo[r][c];
                    }
                }
            }

            gemm('N', 'T', norb[i], norb[i], occ[i],
                 1.0,   mo.data(), norb[i],
                        mo.data(), norb[i],
                 0.0, dens.data(), norb[i]);

            gemm('N', 'N', norb[i], norb[i], norb[i],
                 1.0, fock.data(), norb[i],
                        mo.data(), norb[i],
                 0.0,  tmp.data(), norb[i]);
            gemm('T', 'N', norb[i], norb[i], norb[i],
                 1.0,   mo.data(), norb[i],
                       tmp.data(), norb[i],
                 0.0,  fmo.data(), norb[i]);

            copy(norb[i], fmo.data(), norb[i]+1, E[i].data(), 1);

            if (semicanonical)
            {
                tmp = fmo;

                heev('V', 'U', occ[i], &tmp[     0][     0], norb[i], &E[i][     0]);
                heev('V', 'U', vrt[i], &tmp[occ[i]][occ[i]], norb[i], &E[i][occ[i]]);

                matrix<double> tmp2 = copy(mo);

                gemm('N', 'N', norb[i], occ[i], occ[i],
                     1.0, &tmp2[     0][     0], norb[i],
                           &tmp[     0][     0], norb[i],
                     0.0,   &mo[     0][     0], norb[i]);
                gemm('N', 'N', norb[i], vrt[i], vrt[i],
                     1.0, &tmp2[occ[i]][     0], norb[i],
                           &tmp[occ[i]][occ[i]], norb[i],
                     0.0,   &mo[occ[i]][     0], norb[i]);
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

    //TODO frozen core
    //TODO vrt, occ
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
