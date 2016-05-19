#ifndef _AQUARIUS_OPERATOR_PERTURBEDST2EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_PERTURBEDST2EOPERATOR_HPP_

#include "util/global.hpp"

#include "2eoperator.hpp"
#include "st2eoperator.hpp"

namespace aquarius
{
namespace op
{

/*
 *      _A    -T  A  T    _   A      A  T     _  A
 * Form X  = e   X  e  + [X, T ] = (X  e ) + (X T ) , up to two-electron terms
 *                                        c        c
 */
template <typename U>
class PerturbedSTTwoElectronOperator : public STTwoElectronOperator<U>
{
    protected:
        const STTwoElectronOperator<U>& X;
        const ExcitationOperator<U,2>& TA;

    public:
        template <int N>
        PerturbedSTTwoElectronOperator(const string& name, const STTwoElectronOperator<U>& X, const OneElectronOperator<U>& XA,
                                       const ExcitationOperator<U,N>& T, const ExcitationOperator<U,N>& TA)
        : STTwoElectronOperator<U>(name, XA, T), X(X), TA(TA)
        {
            OneElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

            tensor::SpinorbitalTensor<U>& IMI = I.getIJ();
            tensor::SpinorbitalTensor<U>& IAE = I.getAB();
            tensor::SpinorbitalTensor<U>& IME = I.getIA();

            IME["me"]  = X.getIJAB()["mnef"]*TA(1)["fn"];

            IMI["mi"]  = X.getIJAK()["nmei"]*TA(1)["en"];
            IMI["mi"] += 0.5*X.getIJAB()["mnef"]*TA(2)["efin"];

            IAE["ae"]  = X.getAIBC()["amef"]*TA(1)["fm"];
            IAE["ae"] -= 0.5*X.getIJAB()["mnef"]*TA(2)["afmn"];

            this->ia["ia"] += IME["ia"];

            this->ij["ij"] += X.getIA()["ie"]*TA(1)["ej"];
            this->ij["ij"] += IMI["ij"];

            this->ab["ab"] -= X.getIA()["mb"]*TA(1)["am"];
            this->ab["ab"] += IAE["ab"];

            this->ai["ai"] += X.getAB()["ae"]*TA(1)["ei"];
            this->ai["ai"] -= X.getIJ()["mi"]*TA(1)["am"];
            this->ai["ai"] += X.getIA()["me"]*TA(2)["aeim"];
            this->ai["ai"] -= X.getAIBJ()["amei"]*TA(1)["em"];
            this->ai["ai"] -= 0.5*X.getIJAK()["nmei"]*TA(2)["aemn"];
            this->ai["ai"] += 0.5*X.getAIBC()["amef"]*TA(2)["efim"];

            this->getIJAK()["ijak"] += X.getIJAB()["ijae"]*TA(1)["ek"];

            this->getAIBC()["aibc"] -= X.getIJAB()["mibc"]*TA(1)["am"];

            this->getIJKL()["ijkl"] += X.getIJAK()["jiek"]*TA(1)["el"];
            this->getIJKL()["ijkl"] += 0.5*X.getIJAB()["ijef"]*TA(2)["efkl"];

            this->getABCD()["abcd"] -= X.getAIBC()["amcd"]*TA(1)["bm"];
            this->getABCD()["abcd"] += 0.5*X.getIJAB()["mncd"]*TA(2)["abmn"];

            this->getAIBJ()["aibj"] += X.getAIBC()["aibe"]*TA(1)["ej"];
            this->getAIBJ()["aibj"] -= X.getIJAK()["mibj"]*TA(1)["am"];
            this->getAIBJ()["aibj"] -= X.getIJAB()["mibe"]*TA(2)["aemj"];

            this->getAIJK()["aijk"] += IME["ie"]*T(2)["aejk"];
            this->getAIJK()["aijk"] += X.getAIBJ()["aiek"]*TA(1)["ej"];
            this->getAIJK()["aijk"] -= X.getIJKL()["mijk"]*TA(1)["am"];
            this->getAIJK()["aijk"] += X.getIJAK()["miek"]*TA(2)["aejm"];
            this->getAIJK()["aijk"] += 0.5*X.getAIBC()["aief"]*TA(2)["efjk"];

            this->getABCI()["abci"] -= IME["mc"]*T(2)["abmi"];
            this->getABCI()["abci"] -= X.getAIBJ()["amci"]*TA(1)["bm"];
            this->getABCI()["abci"] += X.getABCD()["abce"]*TA(1)["ei"];
            this->getABCI()["abci"] += X.getAIBC()["amce"]*TA(2)["beim"];
            this->getABCI()["abci"] += 0.5*X.getIJAK()["mnci"]*TA(2)["abmn"];

            this->abij["abij"] += IAE["ae"]*T(2)["ebij"];
            this->abij["abij"] -= IMI["mi"]*T(2)["abmj"];
            this->abij["abij"] += X.getABCI()["abej"]*TA(1)["ei"];
            this->abij["abij"] -= X.getAIJK()["bmji"]*TA(1)["am"];
            this->abij["abij"] += 0.5*X.getABCD()["abef"]*TA(2)["efij"];
            this->abij["abij"] += 0.5*X.getIJKL()["mnij"]*TA(2)["abmn"];
            this->abij["abij"] -= X.getAIBJ()["amei"]*TA(2)["ebmj"];
        }

        template <int N>
        PerturbedSTTwoElectronOperator(const string& name, const STTwoElectronOperator<U>& X, const TwoElectronOperator<U>& XA,
                                       const ExcitationOperator<U,N>& T, const ExcitationOperator<U,N>& TA)
        : STTwoElectronOperator<U>(name, XA, T), X(X), TA(TA)
        {
            OneElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

            tensor::SpinorbitalTensor<U>& IMI = I.getIJ();
            tensor::SpinorbitalTensor<U>& IAE = I.getAB();
            tensor::SpinorbitalTensor<U>& IME = I.getIA();

            IME["me"]  = X.getIJAB()["mnef"]*TA(1)["fn"];

            IMI["mi"]  = X.getIJAK()["nmei"]*TA(1)["en"];
            IMI["mi"] += 0.5*X.getIJAB()["mnef"]*TA(2)["efin"];

            IAE["ae"]  = X.getAIBC()["amef"]*TA(1)["fm"];
            IAE["ae"] -= 0.5*X.getIJAB()["mnef"]*TA(2)["afmn"];

            this->ia["ia"] += IME["ia"];

            this->ij["ij"] += X.getIA()["ie"]*TA(1)["ej"];
            this->ij["ij"] += IMI["ij"];

            this->ab["ab"] -= X.getIA()["mb"]*TA(1)["am"];
            this->ab["ab"] += IAE["ab"];

            this->ai["ai"] += X.getAB()["ae"]*TA(1)["ei"];
            this->ai["ai"] -= X.getIJ()["mi"]*TA(1)["am"];
            this->ai["ai"] += X.getIA()["me"]*TA(2)["aeim"];
            this->ai["ai"] -= X.getAIBJ()["amei"]*TA(1)["em"];
            this->ai["ai"] -= 0.5*X.getIJAK()["nmei"]*TA(2)["aemn"];
            this->ai["ai"] += 0.5*X.getAIBC()["amef"]*TA(2)["efim"];

            this->getIJAK()["ijak"] += X.getIJAB()["ijae"]*TA(1)["ek"];

            this->getAIBC()["aibc"] -= X.getIJAB()["mibc"]*TA(1)["am"];

            this->getIJKL()["ijkl"] += X.getIJAK()["jiek"]*TA(1)["el"];
            this->getIJKL()["ijkl"] += 0.5*X.getIJAB()["ijef"]*TA(2)["efkl"];

            this->getABCD()["abcd"] -= X.getAIBC()["amcd"]*TA(1)["bm"];
            this->getABCD()["abcd"] += 0.5*X.getIJAB()["mncd"]*TA(2)["abmn"];

            this->getAIBJ()["aibj"] += X.getAIBC()["aibe"]*TA(1)["ej"];
            this->getAIBJ()["aibj"] -= X.getIJAK()["mibj"]*TA(1)["am"];
            this->getAIBJ()["aibj"] -= X.getIJAB()["mibe"]*TA(2)["aemj"];

            this->getAIJK()["aijk"] += IME["ie"]*T(2)["aejk"];
            this->getAIJK()["aijk"] += X.getAIBJ()["aiek"]*TA(1)["ej"];
            this->getAIJK()["aijk"] -= X.getIJKL()["mijk"]*TA(1)["am"];
            this->getAIJK()["aijk"] += X.getIJAK()["miek"]*TA(2)["aejm"];
            this->getAIJK()["aijk"] += 0.5*X.getAIBC()["aief"]*TA(2)["efjk"];

            this->getABCI()["abci"] -= IME["mc"]*T(2)["abmi"];
            this->getABCI()["abci"] -= X.getAIBJ()["amci"]*TA(1)["bm"];
            this->getABCI()["abci"] += X.getABCD()["abce"]*TA(1)["ei"];
            this->getABCI()["abci"] += X.getAIBC()["amce"]*TA(2)["beim"];
            this->getABCI()["abci"] += 0.5*X.getIJAK()["mnci"]*TA(2)["abmn"];

            this->abij["abij"] += IAE["ae"]*T(2)["ebij"];
            this->abij["abij"] -= IMI["mi"]*T(2)["abmj"];
            this->abij["abij"] += X.getABCI()["abej"]*TA(1)["ei"];
            this->abij["abij"] -= X.getAIJK()["bmji"]*TA(1)["am"];
            this->abij["abij"] += 0.5*X.getABCD()["abef"]*TA(2)["efij"];
            this->abij["abij"] += 0.5*X.getIJKL()["mnij"]*TA(2)["abmn"];
            this->abij["abij"] -= X.getAIBJ()["amei"]*TA(2)["ebmj"];
        }
};

}
}

#endif
