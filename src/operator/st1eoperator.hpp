#ifndef _AQUARIUS_OPERATOR_ST1EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_ST1EOPERATOR_HPP_

#include "util/global.hpp"

#include "1eoperator.hpp"
#include "excitationoperator.hpp"

namespace aquarius
{
namespace op
{

/*
 *      _    -T   T       T
 * Form X = e  X e  = (X e ) , up to one-electron terms
 *                          c
 */
template <typename U>
class STOneElectronOperator : public OneElectronOperator<U>
{
    public:
        template <int N>
        STOneElectronOperator(const string& name, const OneElectronOperator<U>& X, const ExcitationOperator<U,N>& T);

        template <int N>
        STOneElectronOperator(const string& name, const TwoElectronOperator<U>& X, const ExcitationOperator<U,N>& T);
};

template <typename U> template <int N>
STOneElectronOperator<U>::STOneElectronOperator(const string& name, const OneElectronOperator<U>& X,
                                                const ExcitationOperator<U,N>& T)
: OneElectronOperator<U>(name, X)
{
    assert(N >= 2 && N <= 4);

    this->ij["mi"] += this->ia["me"]*T(1)["ei"];

    this->ai["ai"] += T(2)["aeim"]*this->ia["me"];
    this->ai["ai"] += T(1)["ei"]*this->ab["ae"];
    this->ai["ai"] -= T(1)["am"]*this->ij["mi"];

    this->ab["ae"] -= this->ia["me"]*T(1)["am"];
}

template <typename U> template <int N>
STOneElectronOperator<U>::STOneElectronOperator(const string& name, const TwoElectronOperator<U>& X,
                                                const ExcitationOperator<U,N>& T)
: OneElectronOperator<U>(name, X)
{
    assert(N >= 2 && N <= 4);

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    this->ia["me"] = this->ijab["mnef"]*T(1)["fn"];

    this->ij["mi"] += 0.5*this->ijab["nmef"]*T(2)["efni"];
    this->ij["mi"] += this->ia["me"]*T(1)["ei"];
    this->ij["mi"] += this->ijak["nmfi"]*T(1)["fn"];

    this->ai["ai"] -= T(1)["em"]*this->aibj["amei"];
    this->ai["ai"] += 0.5*this->aibc["amef"]*Tau["efim"];
    this->ai["ai"] -= 0.5*this->ijak["mnei"]*T(2)["eamn"];
    this->ai["ai"] += T(2)["aeim"]*this->ia["me"];
    this->ai["ai"] += T(1)["ei"]*this->ab["ae"];
    this->ai["ai"] -= T(1)["am"]*this->ij["mi"];

    this->ab["ae"] -= 0.5*this->ijab["mnfe"]*T(2)["famn"];
    this->ab["ae"] -= this->ia["me"]*T(1)["am"];
    this->ab["ae"] += this->aibc["anef"]*T(1)["fn"];

    if (N > 2)
    {
        this->ai["ai"] += 0.25*this->ijab["mnef"]*T(3)["aefimn"];
    }
}

}
}

#endif
