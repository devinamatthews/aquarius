#include "choleskyuhf.hpp"

using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace scf
{

template <typename T, template <typename T_> class WhichUHF>
CholeskyUHF<T,WhichUHF>::CholeskyUHF(const string& name, Config& config)
: WhichUHF<T>(name, config)
{
    for (Product& p : products)
    {
        p.addRequirement(Requirement("cholesky", "cholesky"));
    }
}

template <typename T, template <typename T_> class WhichUHF>
bool CholeskyUHF<T,WhichUHF>::run(TaskDAG& dag, const Arena& arena)
{
    const auto& molecule = this->template get<Molecule>("molecule");
    const auto& chol = this->template get<CholeskyIntegrals<T>>("cholesky");

    const PointGroup& group = molecule.getGroup();

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumAlphaElectrons();

    vector<int> shapeN{NS};
    vector<int> shapeNNN{NS,NS,NS};
    vector<vector<int>> sizer{{chol.getRank()}};
    vector<vector<int>> sizenOr{norb,{nalpha},{chol.getRank()}};
    vector<vector<int>> sizenor{norb,{nbeta},{chol.getRank()}};

    this->puttmp("J",       new SymmetryBlockedTensor<T>("J",    arena, group, 1, sizer,   shapeN,   false));
    this->puttmp("JD",      new SymmetryBlockedTensor<T>("JD",   arena, group, 1, sizer,   shapeN,   false));
    this->puttmp("La_occ",  new SymmetryBlockedTensor<T>("LpI",  arena, group, 3, sizenOr, shapeNNN, false));
    this->puttmp("Lb_occ",  new SymmetryBlockedTensor<T>("Lpi",  arena, group, 3, sizenor, shapeNNN, false));
    this->puttmp("LDa_occ", new SymmetryBlockedTensor<T>("LDpI", arena, group, 3, sizenOr, shapeNNN, false));
    this->puttmp("LDb_occ", new SymmetryBlockedTensor<T>("LDpi", arena, group, 3, sizenor, shapeNNN, false));

    return UHF<T>::run(dag, arena);
}

template <typename T, template <typename T_> class WhichUHF>
void CholeskyUHF<T,WhichUHF>::buildFock()
{
    const auto& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumAlphaElectrons();

    const auto& chol = this->template get<CholeskyIntegrals<T>>("cholesky");

    auto& H  = this->template get<SymmetryBlockedTensor<T>>("H");
    auto& Da = this->template get<SymmetryBlockedTensor<T>>("Da");
    auto& Db = this->template get<SymmetryBlockedTensor<T>>("Db");
    auto& Fa = this->template get<SymmetryBlockedTensor<T>>("Fa");
    auto& Fb = this->template get<SymmetryBlockedTensor<T>>("Fb");

    SymmetryBlockedTensor<T> Ca_occ("CI", this->template gettmp<SymmetryBlockedTensor<T>>("Ca"),
                                    {{0},{0}}, {norb,{nalpha}});
    SymmetryBlockedTensor<T> Cb_occ("Ci", this->template gettmp<SymmetryBlockedTensor<T>>("Cb"),
                                    {{0},{0}}, {norb,{nbeta}});

    auto& J       = this->template gettmp<SymmetryBlockedTensor<T>>("J");
    auto& JD      = this->template gettmp<SymmetryBlockedTensor<T>>("JD");
    auto& La_occ  = this->template gettmp<SymmetryBlockedTensor<T>>("La_occ");
    auto& Lb_occ  = this->template gettmp<SymmetryBlockedTensor<T>>("Lb_occ");
    auto& LDa_occ = this->template gettmp<SymmetryBlockedTensor<T>>("LDa_occ");
    auto& LDb_occ = this->template gettmp<SymmetryBlockedTensor<T>>("LDb_occ");

    /*
     * Coulomb contribution:
     *
     * F[ab] = (Da[cd]+Db[cd])*(ab|cd)
     *
     *       = (Da[cd]+Db[cd])*L[abJ]*D[J]*L[cdJ]
     *
     *       = L(abJ)*{D[J]*{(Da[cd]+Db[cd])*L[cdJ]}}
     *
     *       = L(abJ)*{D[J]*J[J]}
     */
    Da += Db;
    J["J"] = chol.getL()["cdJ"]*Da["cd"];
    JD["J"] = chol.getD()["J"]*J["J"];
    Da -= Db;
    Fa["ab"] = JD["J"]*chol.getL()["abJ"];

    /*
     * Core contribution:
     *
     * F += H
     *
     * Up though this point, Fa = Fb
     */
    Fa += H;
    Fb  = Fa;

    /*
     * Exchange contribution:
     *
     * Fa[ab] -= Da[cd]*(ac|bd)
     *
     *         = C[ci]*C[di]*L[acJ]*D[J]*L[bdJ]
     *
     *         = {C[ci]*L[acJ]}*{D[J]*{C[di]*L[bdJ]}}
     *
     *         = L[aiJ]*{D[J]*L[biJ]}
     */
    La_occ["aiJ"] = chol.getL()["acJ"]*Ca_occ["ci"];
    LDa_occ["aiJ"] = chol.getD()["J"]*La_occ["aiJ"];
    Fa["ab"] -= LDa_occ["aiJ"]*La_occ["biJ"];

    Lb_occ["aiJ"] = chol.getL()["acJ"]*Cb_occ["ci"];
    LDb_occ["aiJ"] = chol.getD()["J"]*Lb_occ["aiJ"];
    Fb["ab"] -= LDb_occ["aiJ"]*Lb_occ["biJ"];
}

}
}

INSTANTIATE_SPECIALIZATIONS_2(aquarius::scf::CholeskyUHF, aquarius::scf::LocalUHF);
INSTANTIATE_SPECIALIZATIONS_2(aquarius::scf::CholeskyUHF, aquarius::scf::ElementalUHF);
REGISTER_TASK(aquarius::scf::CholeskyUHF<double,LocalUHF>, "localcholeskyuhf");
REGISTER_TASK(aquarius::scf::CholeskyUHF<double,ElementalUHF>, "elementalcholeskyuhf");
