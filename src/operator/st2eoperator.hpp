#ifndef _AQUARIUS_OPERATOR_ST2EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_ST2EOPERATOR_HPP_

#include "util/global.hpp"

#include "tensor/spinorbital_tensor.hpp"

#include "2eoperator.hpp"
#include "deexcitationoperator.hpp"
#include "excitationoperator.hpp"

namespace aquarius
{
namespace op
{

/*
 *      _    -T   T       T
 * Form X = e  X e  = (X e ) , up to two-electron terms
 *                          c
 */
template <typename U>
class STTwoElectronOperator : public TwoElectronOperator<U>
{
    public:
        template <int N>
        STTwoElectronOperator(const string& name, const OneElectronOperator<U>& X, const ExcitationOperator<U,N>& T)
        : TwoElectronOperator<U>(name, X)
        {
            assert(N >= 2 && N <= 4);

            this->ij["mi"] += this->ia["me"]*T(1)["ei"];

            this->ai["ai"] += T(2)["aeim"]*this->ia["me"];
            this->ai["ai"] += T(1)["ei"]*this->ab["ae"];
            this->ai["ai"] -= T(1)["am"]*this->ij["mi"];

            this->ab["ae"] -= this->ia["me"]*T(1)["am"];

            this->abij["abij"] += this->ab["af"]*T(2)["fbij"];
            this->abij["abij"] -= this->ij["ni"]*T(2)["abnj"];

            this->aijk["amij"] += this->ia["me"]*T(2)["aeij"];

            this->abci["abej"] -= this->ia["me"]*T(2)["abmj"];

            if (N > 2)
            {
                this->abij["abij"] += T(3)["abeijm"]*this->ia["me"];
            }
        }

        template <int N>
        STTwoElectronOperator(const string& name, const TwoElectronOperator<U>& X, const ExcitationOperator<U,N>& T, bool isHbar=false)
        : TwoElectronOperator<U>(name, X)
        {
            assert(N >= 2 && N <= 4);

            tensor::SpinorbitalTensor<U> Tau(T(2));
            Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

            this->ia["me"] = this->ijab["mnef"]*T(1)["fn"];

            this->ij["mi"] += 0.5*this->ijab["nmef"]*T(2)["efni"];
            this->ij["mi"] += this->ia["me"]*T(1)["ei"];
            this->ij["mi"] += this->ijak["nmfi"]*T(1)["fn"];

            this->ijkl["mnij"] += 0.5*this->ijab["mnef"]*Tau["efij"];
            this->ijkl["mnij"] += this->ijak["nmei"]*T(1)["ej"];

            this->ijak["mnej"] += this->ijab["mnef"]*T(1)["fj"];

            if (!isHbar)
            {
                this->ai["ai"] -= T(1)["em"]*this->aibj["amei"];
                this->ai["ai"] += 0.5*this->aibc["amef"]*Tau["efim"];
                this->ai["ai"] -= 0.5*this->ijak["mnei"]*T(2)["eamn"];
                this->ai["ai"] += T(2)["aeim"]*this->ia["me"];
                this->ai["ai"] += T(1)["ei"]*this->ab["ae"];
                this->ai["ai"] -= T(1)["am"]*this->ij["mi"];
            }
            else
            {
                this->ai = 0;
            }

            this->ab["ae"] -= 0.5*this->ijab["mnfe"]*T(2)["famn"];
            this->ab["ae"] -= this->ia["me"]*T(1)["am"];
            this->ab["ae"] += this->aibc["anef"]*T(1)["fn"];

            this->aijk["amij"] += 0.5*this->aibc["amef"]*Tau["efij"];
            this->aijk["amij"] += this->aibj["amej"]*T(1)["ei"];

            this->aibj["amei"] -= 0.5*this->ijab["mnef"]*T(2)["afin"];
            this->aibj["amei"] -= this->aibc["amfe"]*T(1)["fi"];
            this->aibj["amei"] -= this->ijak["nmei"]*T(1)["an"];

            if (!isHbar)
            {
                this->abij["abij"] += this->ab["af"]*T(2)["fbij"];
                this->abij["abij"] -= this->ij["ni"]*T(2)["abnj"];
                this->abij["abij"] += this->abci["abej"]*T(1)["ei"];
                this->abij["abij"] -= this->aijk["amij"]*T(1)["bm"];
                this->abij["abij"] += 0.5*this->abcd["abef"]*Tau["efij"];
                this->abij["abij"] += 0.5*this->ijkl["mnij"]*Tau["abmn"];
                this->abij["abij"] -= this->aibj["amei"]*T(2)["ebmj"];
            }
            else
            {
                this->abij = 0;
            }

            this->aijk["amij"] += this->ijak["nmej"]*T(2)["aein"];
            this->aijk["amij"] -= this->ijkl["nmij"]*T(1)["an"];
            this->aijk["amij"] += this->ia["me"]*T(2)["aeij"];

            this->aibj["amei"] -= 0.5*this->ijab["mnef"]*T(2)["afin"];
            this->aibj["amei"] += 0.5*this->ijak["nmei"]*T(1)["an"];

            this->abci["abej"] += 0.5*this->ijak["mnej"]*T(2)["abmn"];
            this->abci["abej"] -= this->aibj["amej"]*T(1)["bm"];
            this->abci["abej"] += this->aibc["amef"]*T(2)["fbmj"];
            this->abci["abej"] += this->abcd["abef"]*T(1)["fj"];
            this->abci["abej"] -= this->ia["me"]*T(2)["abmj"];

            this->aibj["amei"] -= 0.5*this->ijak["nmei"]*T(1)["an"];

            this->abcd["abef"] += 0.5*this->ijab["mnef"]*Tau["abmn"];
            this->abcd["abef"] -= this->aibc["amef"]*T(1)["bm"];

            this->aibc["amef"] -= this->ijab["nmef"]*T(1)["an"];

            if (N > 2)
            {
                if (!isHbar)
                {
                    this->ai["ai"] += 0.25*this->ijab["mnef"]*T(3)["aefimn"];

                    this->abij["abij"] += this->ia["me"]*T(3)["abeijm"];
                    this->abij["abij"] -= 0.5*this->ijak["mnej"]*T(3)["aebimn"];
                    this->abij["abij"] += 0.5*this->aibc["amef"]*T(3)["efbimj"];
                }

                this->abci["abej"] -= 0.5*this->ijab["mnef"]*T(3)["abfmjn"];

                this->aijk["amij"] += 0.5*this->ijab["mnef"]*T(3)["aefijn"];
            }

            if (N > 3 && !isHbar)
            {
                this->abij["abij"] += 0.25*this->ijab["mnef"]*T(4)["abefijmn"];
            }
        }
};

}
}

#endif
