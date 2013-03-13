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

#ifndef _AQUARIUS_CC_EOMCCSD_HPP_
#define _AQUARIUS_CC_EOMCCSD_HPP_

#include "convergence/davidson.hpp"
#include "operator/2eoperator.hpp"

#include "ccsd.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class EOMEECCSD : public Iterative, public op::ExcitationOperator<U,2>
{
    protected:
        op::ExcitationOperator<U,2>& R;
        op::ExcitationOperator<U,2> Z, D;
        const op::TwoElectronOperator<T>& H;
        const op::ExponentialOperator<U,2>& T;
        convergence::Davidson< op::ExcitationOperator<U,2> > davidson;

    public:
        EOMEECCSD(const input::Config& config, const op::TwoElectronOperator<T>& H,
                  const op::ExponentialOperator<U,2>& T)
        : tensor::Tensor<op::ExcitationOperator<U,2>,U>(*this),
          Iterative(config), op::ExcitationOperator<U,2>(H.getSCF()), R(*this),
          Z(this->uhf), D(this->uhf), H(H), T(T), davidson(config.get("davidson"))
        {
            D(0) = 1;
            D(1)["ai"]  = H.getIJ()["ii"];
            D(1)["ai"] -= H.getAB()["aa"];
            D(1)["ai"] -= H.getAIBJ()["aiai"];
            D(2)["abij"]  = H.getIJ()["ii"];
            D(2)["abij"] += H.getIJ()["jj"];
            D(2)["abij"] -= H.getAB()["aa"];
            D(2)["abij"] -= H.getAB()["bb"];
            D(2)["abij"] -= H.getAIBJ()["aiai"];
            D(2)["abij"] -= H.getAIBJ()["bjbj"];
            D(2)["abij"] += H.getABCD()["abab"];
            D(2)["abij"] += H.getIJKL()["ijij"];
            D(2)["abij"] -= H.getIJAB()["imab"]*T(2)["abim"];
            D(2)["abij"] += H.getIJAB()["ijae"]*T(2)["aeij"];

            //TODO: guess
        }

        void _iterate()
        {
            op::OneElectronOperator<U> X(this->uhf);

            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FME = H.getIA();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FAE = H.getAB();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FMI = H.getIJ();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNEF = H.getIJAB();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WAMEF = H.getAIBC();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABEJ = H.getABCI();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABEF = H.getABCD();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNIJ = H.getIJKL();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNIE = H.getIJKA();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMBIJ = H.getIAJK();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WAMEI = H.getAIBJ();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XMI = X.getIJ();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XAE = X.getAB();

            /**************************************************************************
             *
             * X intermediates
             */
            XMI["mi"]  = WMNIE["mnie"]*R(1)["en"];
            XMI["mi"] += 0.5*WMNEF["mnef"]*R(2)["efin"];

            XAE["ae"]  = WAMEF["amef"]*R(1)["fm"];
            XMI["ae"] -= 0.5*WMNEF["mnef"]*R(2)["afmn"];
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * R(1)->R(1) and R(2)->R(1)
             */
            Z(1)["ai"]  = R(1)["ei"]*FAE["ae"];

            Z(1)["ai"] -= R(1)["am"]*FMI["mi"];

            Z(1)["ai"] -= R(1)["em"]*WAMEI["amei"];

            Z(1)["ai"] += R(2)["aeim"]*FME["me"];

            Z(1)["ai"] += 0.5*WAMEF["amef"]*R(2)["efim"];

            Z(1)["ai"] -= 0.5*WMNIE["mnie"]*R(2)["aemn"];
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * R(1)->R(2) and R(2)->R(2)
             */
            Z(2)["abij"]  = FAE["ae"]*R(2)["ebij"];

            Z(2)["abij"] -= FMI["mi"]*R(2)["abmj"];

            Z(2)["abij"] += XAE["ae"]*T(2)["ebij"];

            Z(2)["abij"] -= XMI["mi"]*T(2)["abmj"];

            Z(2)["abij"] += WABEJ["abej"]*R(1)["ei"];

            Z(2)["abij"] -= WMBIJ["mbij"]*R(1)["am"];

            Z(2)["abij"] += 0.5*WABEF["abef"]*R(2)["efij"];

            Z(2)["abij"] += 0.5*WMNIJ["mnij"]*R(2)["abmn"];

            Z(2)["abij"] -= WAMEI["amei"]*R(2)["ebmj"];
            /*
             *************************************************************************/

            energy = davidson.extrapolate(R, Z, D);

            conv =               Z(1)(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,Z(1)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(2).reduce(CTF_OP_MAXABS));
        }
};

}
}

#endif
