/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_OPERATOR_STEXCITATIONOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_STEXCITATIONOPERATOR_HPP_

#include "excitationoperator.hpp"
#include "exponentialoperator.hpp"

namespace aquarius
{
namespace op
{

/*
 *                                  _    -T   T       T
 * Form the pure excitation part of X = e  X e  = (X e )
 *                                                      c
 */
template <typename U, int nex> class STExcitationOperator;

template <typename U>
class STExcitationOperator<U,2> : public ExcitationOperator<U,2>
{
    public:
        STExcitationOperator(const OneElectronOperator<U>& X, const ExponentialOperator<U,2>& T)
        : ExcitationOperator<U,2>(X.uhf)
        {
            transform(X, T, *this);
        }

        STExcitationOperator(const TwoElectronOperator<U>& X, const ExponentialOperator<U,2>& T)
        : ExcitationOperator<U,2>(X.uhf)
        {
            transform(X, T, *this);
        }

        static void transform(const OneElectronOperator<U>& X, const ExponentialOperator<U,2>& T,
                              ExcitationOperator<U,2>& Z)
        {
            OneElectronOperator<U> W(X, OneElectronOperator<U>::AB|
                                        OneElectronOperator<U>::IJ);

            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FAI = W.getAI();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FME = W.getIA();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FAE = W.getAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FMI = W.getIJ();

            Z(0) = FME["me"]*T(1)["em"];

            FMI["mi"] += FME["me"]*T(1)["ei"];

            Z(1)["ai"]  = FAI["ai"];
            Z(1)["ai"] += T(2)["aeim"]*FME["me"];
            Z(1)["ai"] += T(1)["ei"]*FAE["ae"];
            Z(1)["ai"] -= T(1)["am"]*FMI["mi"];

            FAE["ae"] -= FME["me"]*T(1)["am"];

            Z(2)["abij"]  = FAE["af"]*T(2)["fbij"];
            Z(2)["abij"] -= FMI["ni"]*T(2)["abnj"];
        }

        static void transform(const TwoElectronOperator<U>& X, const ExponentialOperator<U,2>& T,
                              ExcitationOperator<U,2>& Z)
        {
            TwoElectronOperator<U> W(const_cast<TwoElectronOperator<U>&>(X), TwoElectronOperator<U>::AB|
                                                                             TwoElectronOperator<U>::IJ|
                                                                             TwoElectronOperator<U>::IA|
                                                                             TwoElectronOperator<U>::IJKL|
                                                                             TwoElectronOperator<U>::IJKA|
                                                                             TwoElectronOperator<U>::IAJK|
                                                                             TwoElectronOperator<U>::AIBJ);

            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FAI = W.getAI();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FME = W.getIA();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FAE = W.getAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FMI = W.getIJ();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABIJ = W.getABIJ();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNEF = W.getIJAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WAMEF = W.getAIBC();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABEJ = W.getABCI();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABEF = W.getABCD();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNIJ = W.getIJKL();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNIE = W.getIJKA();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMBIJ = W.getIAJK();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WAMEI = W.getAIBJ();

            tensor::SpinorbitalTensor< tensor::DistTensor<U> > Tau(T(2));
            Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

            Z(0) = 0;
            //Z(0)  = FME["me"]*T(1)["em"];
            //Z(0) += 0.25*WMNEF["mnef"]*Tau["efmn"];

            FME["me"] += WMNEF["mnef"]*T(1)["fn"];

            FMI["mi"] += 0.5*WMNEF["mnef"]*T(2)["efin"];
            FMI["mi"] += FME["me"]*T(1)["ei"];
            FMI["mi"] += WMNIE["mnif"]*T(1)["fn"];

            WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*Tau["efij"];
            WMNIJ["mnij"] += WMNIE["mnie"]*T(1)["ej"];

            WMNIE["mnie"] += WMNEF["mnfe"]*T(1)["fi"];

            Z(1)["ai"]  = FAI["ai"];
            Z(1)["ai"] -= T(1)["em"]*WAMEI["amei"];
            Z(1)["ai"] += 0.5*WAMEF["amef"]*Tau["efim"];
            Z(1)["ai"] -= 0.5*WMNIE["mnie"]*T(2)["aemn"];
            Z(1)["ai"] += T(2)["aeim"]*FME["me"];
            Z(1)["ai"] += T(1)["ei"]*FAE["ae"];
            Z(1)["ai"] -= T(1)["am"]*FMI["mi"];

            FAE["ae"] -= 0.5*WMNEF["mnef"]*T(2)["afmn"];
            FAE["ae"] -= FME["me"]*T(1)["am"];
            FAE["ae"] += WAMEF["amef"]*T(1)["fm"];

            WMBIJ["mbij"] += 0.5*WAMEF["bmfe"]*Tau["efij"];
            WMBIJ["mbij"] -= WAMEI["bmej"]*T(1)["ei"];

            WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
            WAMEI["amei"] -= WAMEF["amfe"]*T(1)["fi"];
            WAMEI["amei"] += WMNIE["nmie"]*T(1)["an"];

            Z(2)["abij"]  = WABIJ["abij"];
            Z(2)["abij"] += FAE["af"]*T(2)["fbij"];
            Z(2)["abij"] -= FMI["ni"]*T(2)["abnj"];
            Z(2)["abij"] += WABEJ["abej"]*T(1)["ei"];
            Z(2)["abij"] -= WMBIJ["mbij"]*T(1)["am"];
            Z(2)["abij"] += 0.5*WABEF["abef"]*Tau["efij"];
            Z(2)["abij"] += 0.5*WMNIJ["mnij"]*Tau["abmn"];
            Z(2)["abij"] -= WAMEI["amei"]*T(2)["ebmj"];
        }
};

}
}

#endif
