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

#ifndef _AQUARIUS_CC_LAMBDACCSD_HPP_
#define _AQUARIUS_CC_LAMBDACCSD_HPP_

#include "operator/2eoperator.hpp"
#include "operator/deexcitationoperator.hpp"

#include "ccsd.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class LambdaCCSD : public Iterative, public op::DeexcitationOperator<U,2>
{
    protected:
        op::DeexcitationOperator<U,2>& L;
        op::DeexcitationOperator<U,2> Z, D, Dreal;
        op::TwoElectronOperator<U> H;
        const CCSD<U>& T;
        convergence::DIIS< op::DeexcitationOperator<U,2> > diis;

    public:
        LambdaCCSD(const input::Config& config, const op::TwoElectronOperator<U>& H,
                   const CCSD<U>& T)
        : tensor::Tensor<op::DeexcitationOperator<U,2>,U>(*this),
          Iterative(config), op::DeexcitationOperator<U,2>(H.getSCF()), L(*this),
          Z(this->uhf), D(this->uhf), Dreal(this->uhf),
          H(const_cast<op::TwoElectronOperator<U>&>(H),
            op::TwoElectronOperator<U>::AB|op::TwoElectronOperator<U>::IJ),
          T(T), diis(config.get("diis"))
        {
            D(0) = 1;
            D(1)["ia"]  = H.getIJ()["ii"];
            D(1)["ia"] -= H.getAB()["aa"];
            D(2)["ijab"]  = H.getIJ()["ii"];
            D(2)["ijab"] += H.getIJ()["jj"];
            D(2)["ijab"] -= H.getAB()["aa"];
            D(2)["ijab"] -= H.getAB()["bb"];

            Dreal = D;
            D = 1/D;

            L(0) = 1;
            L(1) = -H.getIA()*D(1);
            L(2) = -H.getIJAB()*D(2);
        }

        void _iterate()
        {
            op::OneElectronOperator<U> X(this->uhf);

            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FME = H.getIA();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FAE = H.getAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FMI = H.getIJ();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNEF = H.getIJAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WAMEF = H.getAIBC();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABEJ = H.getABCI();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABEF = H.getABCD();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNIJ = H.getIJKL();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNIE = H.getIJKA();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMBIJ = H.getIAJK();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WAMEI = H.getAIBJ();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XMN = X.getIJ();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XEF = X.getAB();

            FAE["aa"] = 0.0;
            FMI["ii"] = 0.0;

            Z(0) = 1;

            tensor::SpinorbitalTensor< tensor::DistTensor<U> > Tau(T(2));
            Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

            /**************************************************************************
             *
             * "Density-like" intermediates
             */
            XMN["mn"] = 0.5*T(2)["efno"]*L(2)["moef"];

            XEF["ef"] = -0.5*L(2)["mnfg"]*T(2)["egmn"];
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * L(1)->L(1) and L(2)->L(1)
             */
            Z(1)["ia"] = -FME["ia"];

            Z(1)["ia"] -= L(1)["me"]*WAMEI["eiam"];

            Z(1)["ia"] += 0.5*WABEJ["efam"]*L(2)["imef"];

            Z(1)["ia"] -= 0.5*WMBIJ["iemn"]*L(2)["mnae"];

            Z(1)["ia"] += L(1)["ie"]*FAE["ea"];

            Z(1)["ia"] -= L(1)["ma"]*FMI["im"];

            Z(1)["ia"] -= XMN["mn"]*WMNIE["nima"];

            Z(1)["ia"] -= XEF["ef"]*WAMEF["fiea"];
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * L(1)->L(2) and L(2)->L(2)
             */
            Z(2)["ijab"] = -WMNEF["ijab"];

            Z(2)["ijab"] += FAE["fa"]*L(2)["ijfb"];

            Z(2)["ijab"] -= FMI["in"]*L(2)["njab"];

            Z(2)["ijab"] += FME["ia"]*L(1)["jb"];

            Z(2)["ijab"] += WAMEF["ejab"]*L(1)["ie"];

            Z(2)["ijab"] -= WMNIE["ijmb"]*L(1)["ma"];

            Z(2)["ijab"] += 0.5*WABEF["efab"]*L(2)["ijef"];

            Z(2)["ijab"] += 0.5*WMNIJ["ijmn"]*L(2)["mnab"];

            Z(2)["ijab"] -= WAMEI["eiam"]*L(2)["mjeb"];

            Z(2)["ijab"] -= XMN["im"]*WMNEF["mjab"];

            Z(2)["ijab"] += XEF["ae"]*WMNEF["ijeb"];
            /*
             *************************************************************************/

            op::DeexcitationOperator<U,2> LH(Z);

            LH += L*T.getEnergy();
            LH -= L*Dreal;

            energy = tensor::scalar(LH*conj(L))/tensor::scalar(L*conj(L));

            Z *= D;
            Z -= L;
            L += Z;

            conv =               Z(1)(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,Z(1)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(2).reduce(CTF_OP_MAXABS));

            diis.extrapolate(L, Z);
        }
};

}
}

#endif
