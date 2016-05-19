#include "2edensity.hpp"

using namespace aquarius::op;
using namespace aquarius::tensor;

namespace aquarius
{
namespace cc
{

template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const string& name, const MOSpace<U>& occ, const MOSpace<U>& vrt,
                                          const SymmetryBlockedTensor<U>& Da,
                                          const SymmetryBlockedTensor<U>& Db)
: TwoElectronOperator<U>(name, OneElectronDensity<U>(name, occ, vrt, Da, Db))
{
    this->abcd["abcd"]  =     this->ab["ab"]*this->ab["cd"];
    this->abcd["abcd"] += 0.5*this->ab["ac"]*this->ab["bd"];

    this->abci["abci"]  =     this->ab["ab"]*this->ai["ci"];
    this->abci["abci"] += 0.5*this->ab["ac"]*this->ai["bi"];

    this->abij["abij"]  =     this->ab["ab"]*this->ij["ij"];
    this->abij["abij"] += 0.5*this->ai["ai"]*this->ai["bj"];

    this->aibj["aibj"]  =     this->ai["ai"]*this->ai["bj"];
    this->aibj["aibj"] += 0.5*this->ab["ab"]*this->ij["ij"];

    this->ijak["ijak"]  =     this->ij["ik"]*this->ia["ja"];
    this->ijak["ijak"] += 0.5*this->ij["ij"]*this->ia["ka"];

    this->ijkl["ijkl"]  =     this->ij["ij"]*this->ij["kl"];
    this->ijkl["ijkl"] += 0.5*this->ij["ik"]*this->ij["jl"];

    this->aibc["aibc"] = this->abci["bcai"];
    this->aijk["aijk"] = this->ijak["jkai"];
    this->ijab["ijab"] = this->abij["abij"];
}

/*
 * Form the unrelaxed CCSD Density
 */
template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const string& name, const ExcitationOperator<U,2>& T)
: TwoElectronOperator<U>(name, OneElectronDensity<U>(name, T))
{
    this->abij["abij"]  = T(2)["abij"];
    this->abij["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];
}

/*
 * Form the mixed perturbed CCSD Density
 */
template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const string& name, const DeexcitationOperator<U,2>& L,
                                          const ExcitationOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: TwoElectronOperator<U>(name, OneElectronDensity<U>(name, L, T, TA))
{
    TwoElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

    SpinorbitalTensor<U>& IIJ = I.getIJ();
    SpinorbitalTensor<U>& IAB = I.getAB();
    SpinorbitalTensor<U>& IAI = I.getAI();

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    SpinorbitalTensor<U> TauA(TA(2));
    TauA["abij"] += T(1)["ai"]*TA(1)["bj"];

    this->ijak["ijak"] = -L(2)["ijae"]*TA(1)["ek"];

    this->aibc["aibc"] = L(2)["mibc"]*TA(1)["am"];

    this->ijkl["ijkl"] = 0.5*L(2)["ijef"]*TauA["efkl"];

    this->abcd["abcd"] = 0.5*L(2)["mncd"]*TauA["abmn"];

    this->aibj["aibj"] -= L(1)["ib"]*TA(1)["aj"];
    this->aibj["aibj"] -= L(2)["imbe"]*TA(2)["aejm"];

    IAB["ab"]  = 0.5*T(2)["aemn"]*L(2)["mnbe"];
    IIJ["ij"]  = 0.5*T(2)["efim"]*L(2)["jmef"];
    IAI["ai"]  = IAB["ae"]*TA(1)["ei"];
    IAI["ai"] += IIJ["mi"]*TA(1)["am"];
    IAI["ai"] += L(1)["me"]*TA(2)["aeim"];

    this->abci["abci"] += L(1)["mc"]*TA(2)["abmi"];
    this->abci["abci"] += this->ab["ac"]*T(1)["bi"];
    this->abci["abci"] += IAB["ac"]*TA(1)["bi"];
    this->abci["abci"] += this->aibc["amce"]*T(2)["beim"];
    this->abci["abci"] -= this->aibj["amci"]*T(1)["bm"];

    this->aijk["aijk"] -= L(1)["ie"]*TA(2)["aejk"];
    this->aijk["aijk"] += this->ij["ik"]*T(1)["aj"];
    this->aijk["aijk"] -= IIJ["ik"]*TA(1)["aj"];
    this->aijk["aijk"] += this->ijak["miek"]*T(2)["aejm"];

    this->abij["abij"] += TauA["abij"];
    this->abij["abij"] += this->ijkl["mnij"]*Tau["abmn"];
    this->abij["abij"] += this->abcd["abef"]*Tau["efij"];
    this->abij["abij"] += this->ij["mi"]*T(2)["abmj"];
    this->abij["abij"] -= this->ab["ae"]*T(2)["ebij"];
    this->abij["abij"] -= IIJ["mi"]*TA(2)["abmj"];
    this->abij["abij"] -= IAB["ae"]*TA(2)["ebij"];
    this->abij["abij"] -= this->aibj["amej"]*T(2)["ebmj"];
    this->abij["abij"] += this->aijk["bmji"]*T(1)["am"];
    this->abij["abij"] -= this->abci["abej"]*T(1)["ei"];
    this->abij["abij"] += IAI["ai"]*T(1)["bj"];

    this->abci["abci"] -= 0.5*this->ijak["mnci"]*Tau["abmn"];
    this->abci["abci"] -= this->abcd["abce"]*T(1)["ei"];

    this->aijk["aijk"] += 0.5*this->aibc["aief"]*Tau["efjk"];
    this->aijk["aijk"] += this->aibj["aiek"]*T(1)["ej"];
    this->aijk["aijk"] += this->ijkl["mijk"]*T(1)["am"];

    this->aibj["aibj"] += this->ijak["mibj"]*T(1)["am"];
    this->aibj["aibj"] -= this->aibc["aibe"]*T(1)["ej"];
}

/*
 * Form the relaxed CCSD Density
 */
template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const string& name, const DeexcitationOperator<U,2>& L,
                                          const ExcitationOperator<U,2>& T)
: TwoElectronOperator<U>(name, OneElectronDensity<U>(name, L, T))
{
    SpinorbitalTensor<U> aitmp("I", this->ai);

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    this->ijab["ijab"] =  L(2)["ijab"];

    this->ijak["ijak"] = -L(2)["ijae"]*T(1)["ek"];

    this->aibc["aibc"] =  L(2)["mibc"]*T(1)["am"];

    this->ijkl["ijkl"] =  0.5*L(2)["ijef"]*Tau["efkl"];

    this->abcd["abcd"] =  0.5*L(2)["mncd"]*Tau["abmn"];

    this->aibj["aibj"]  = -L(1)["ib"]*T(1)["aj"];
    this->aibj["aibj"] -=  L(2)["imbe"]*T(2)["aejm"];
    this->aibj["aibj"] +=  this->ijak["mibj"]*T(1)["am"];

    this->aijk["aijk"]  = this->ij["ik"]*T(1)["aj"];

    this->abci["abci"]  = L(1)["mc"]*T(2)["abmi"];
    this->abci["abci"] += this->ab["ac"]*T(1)["bi"];
    this->abci["abci"] += this->aibc["amce"]*T(2)["beim"];
    this->abci["abci"] -= 0.5*this->ijak["mnci"]*Tau["abmn"];

    this->abij["abij"]  = Tau["abij"];
    this->abij["abij"] += this->ijkl["mnij"]*Tau["abmn"];
    this->abij["abij"] += this->ij["mi"]*T(2)["abmj"];
    this->abij["abij"] -= this->ab["ae"]*T(2)["ebij"];
    this->abij["abij"] -= this->aibj["amej"]*T(2)["ebmj"];
    this->abij["abij"] += this->aijk["bmji"]*T(1)["am"];
    aitmp["ai"] = -this->aibc["amef"]*T(2)["efim"];
    this->abij["abij"] += aitmp["ai"]*T(1)["bj"];

    this->aijk["aijk"] -= L(1)["ie"]*T(2)["aejk"];
    this->aijk["aijk"] += this->ijak["miek"]*T(2)["aejm"];
    this->aijk["aijk"] += 0.5*this->aibc["aief"]*Tau["efjk"];
}

/*
 * Form the relaxed perturbed CCSD Density
 */
template <typename U>
TwoElectronDensity<U>::TwoElectronDensity(const string& name, const DeexcitationOperator<U,2>& L,
                                          const DeexcitationOperator<U,2>& LA,
                                          const ExcitationOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: TwoElectronOperator<U>(name, OneElectronDensity<U>(name, L, LA, T, TA))
{
    //TODO
    assert(0);
    return;
}

template <typename U>
double TwoElectronDensity<U>::getS2() const
{
    //TODO
    assert(0);
    return 0.0;
}

template <typename U>
double TwoElectronDensity<U>::getProjectedMultiplicity() const
{
    return sqrt(1+4*getS2());
}

INSTANTIATE_SPECIALIZATIONS(TwoElectronDensity);

}
}
