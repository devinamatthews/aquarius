#include "1edensity.hpp"

using namespace aquarius::op;
using namespace aquarius::tensor;

namespace aquarius
{
namespace cc
{

template <typename U>
OneElectronDensity<U>::OneElectronDensity(const string& name, const MOSpace<U>& occ, const MOSpace<U>& vrt,
                                          const SymmetryBlockedTensor<U>& Da,
                                          const SymmetryBlockedTensor<U>& Db)
: OneElectronOperator<U>(name, occ, vrt, Da, Db) {}

/*
 * Form the unrelaxed CCSD density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const string& name, const ExcitationOperator<U,2>& T)
: OneElectronOperator<U>(name, T.arena, T.occ, T.vrt)
{
    this->ai["ai"] = T(1)["ai"];
}

/*
 * Form the partial perturbed CCSD Density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const string& name, const DeexcitationOperator<U,2>& L,
                                          const ExcitationOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: OneElectronOperator<U>(name, T.arena, T.occ, T.vrt)
{
    OneElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

    SpinorbitalTensor<U>& IIJ = I.getIJ();
    SpinorbitalTensor<U>& IAB = I.getAB();

    IAB["ab"] = 0.5*T(2)["aemn"]*L(2)["mnbe"];
    IIJ["ij"] = 0.5*T(2)["efim"]*L(2)["jmef"];

    this->ab["ab"] += TA(1)["am"]*L(1)["mb"];
    this->ab["ab"] += 0.5*TA(2)["aemn"]*L(2)["mnbe"];

    this->ij["ij"] -= TA(1)["ei"]*L(1)["je"];
    this->ij["ij"] -= 0.5*TA(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] += TA(1)["ai"];
    this->ai["ai"] += TA(2)["aeim"]*L(1)["me"];
    this->ai["ai"] += this->ij["mi"]*T(1)["am"];
    this->ai["ai"] -= this->ab["ae"]*T(1)["ei"];
    this->ai["ai"] -= IIJ["mi"]*TA(1)["am"];
    this->ai["ai"] -= IAB["ae"]*TA(1)["ei"];
}

/*
 * Form the relaxed CCSD density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const string& name, const DeexcitationOperator<U,2>& L,
                                          const ExcitationOperator<U,2>& T)
: OneElectronOperator<U>(name, T.arena, T.occ, T.vrt)
{
    this->ia["ia"] += L(1)["ia"];

    this->ab["ab"] += 0.5*T(2)["aemn"]*L(2)["mnbe"];

    this->ij["ij"] -= T(1)["ei"]*L(1)["je"];
    this->ij["ij"] -= 0.5*T(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] += T(1)["ai"];
    this->ai["ai"] += T(2)["aeim"]*L(1)["me"];
    this->ai["ai"] += this->ij["mi"]*T(1)["am"];
    this->ai["ai"] -= this->ab["ae"]*T(1)["ei"];

    this->ab["ab"] += T(1)["am"]*L(1)["mb"];
}

/*
 * Form the relaxed perturbed CCSD Density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const string& name, const DeexcitationOperator<U,2>& L,
                                          const DeexcitationOperator<U,2>& LA,
                                          const ExcitationOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: OneElectronOperator<U>(name, T.arena, T.occ, T.vrt)
{
    OneElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

    SpinorbitalTensor<U>& IIJ = I.getIJ();
    SpinorbitalTensor<U>& IAB = I.getAB();

    this->ia["ia"] += LA(1)["ia"];

    this->ab["ab"] += 0.5*T(2)["aemn"]*LA(2)["mnbe"];

    this->ij["ij"] -= T(1)["ei"]*LA(1)["je"];
    this->ij["ij"] -= 0.5*T(2)["efim"]*LA(2)["jmef"];

    this->ai["ai"] += T(1)["ai"];
    this->ai["ai"] += T(2)["aeim"]*LA(1)["me"];
    this->ai["ai"] += this->ij["mi"]*T(1)["am"];
    this->ai["ai"] -= this->ab["ae"]*T(1)["ei"];

    this->ab["ab"] += T(1)["am"]*LA(1)["mb"];

    IAB["ab"]  = TA(1)["am"]*L(1)["mb"];
    IAB["ab"] += 0.5*TA(2)["aemn"]*L(2)["mnbe"];

    this->ab["ab"] += IAB["ab"];
    this->ij["ij"] -= IIJ["ij"];

    IIJ["ij"]  = TA(1)["ei"]*L(1)["je"];
    IIJ["ij"] += 0.5*TA(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] += TA(1)["ai"];
    this->ai["ai"] += TA(2)["aeim"]*L(1)["me"];
    this->ai["ai"] -= IIJ["mi"]*T(1)["am"];
    this->ai["ai"] -= IAB["ae"]*T(1)["ei"];

    IAB["ab"] = 0.5*T(2)["aemn"]*L(2)["mnbe"];
    IIJ["ij"] = 0.5*T(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] -= IIJ["mi"]*TA(1)["am"];
    this->ai["ai"] -= IAB["ae"]*TA(1)["ei"];
}

INSTANTIATE_SPECIALIZATIONS(OneElectronDensity);

}
}
