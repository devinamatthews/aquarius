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

#ifndef _AQUARIUS_CC_PERTURBEDCCSD_HPP_
#define _AQUARIUS_CC_PERTURBEDCCSD_HPP_

#include "operator/2eoperator.hpp"

#include "ccsd.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class PerturbedCCSD : public Iterative, public op::ExcitationOperator<U,2>
{
    protected:
        const op::TwoElectronOperator<U>& H;
        U omega;
        op::ExcitationOperator<U,2>& TX;
        op::ExcitationOperator<U,2> D, Z, X;
        convergence::DIIS< op::ExcitationOperator<U,2> > diis;

    public:
        PerturbedCCSD(const input::Config& config, const op::TwoElectronOperator& H,
                      const CCSD<U>& T, const OneElectronOperator<U>& x, const U omega=0)
        : tensor::Tensor<op::ExcitationOperator<U,2>,U>(*this),
          Iterative(config), op::ExcitationOperator<U,2>(ccsd.uhf),
          H(H), T(T), omega(omega), D(this->uhf), Z(this->uhf), X(this->uhf),
          diis(config.get("diis"))
        {
            D(0) = 1;
            D(1)["ia"]  = moints.getIJ()["ii"];
            D(1)["ia"] -= moints.getAB()["aa"];
            D(2)["ijab"]  = moints.getIJ()["ii"];
            D(2)["ijab"] += moints.getIJ()["jj"];
            D(2)["ijab"] -= moints.getAB()["aa"];
            D(2)["ijab"] -= moints.getAB()["bb"];

            D = 1/D;

            op::OneElectronOperator<U> xt(const_cast<OneElectronOperator<U>&>(x), op::OneElectronOperator<U>::ALL);

            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XTAB = xt.getAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XTIJ = xt.getIJ();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XAB = x.getAB();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XAI = x.getAI();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XIA = x.getIA();
            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XIJ = x.getIJ();

            TX(0) = 0;
            X(0) = 0;

            XTIJ["mi"] += XIA["me"]*T(1)["ei"];
            XTAB["ae"] -= XIA["me"]*T(1)["am"];

            X(1)["ai"]  = XAI["ai"];
            X(1)["ai"] += XAB["ae"]*T(1)["ei"];
            X(1)["ai"] -= XTIJ["mi"]*T(1)["am"];
            X(1)["ai"] += XIA["me"]*T(2)["aeim"];

            X(2)["abij"]  = XTAB["ae"]*T(2)["ebij"];
            X(2)["abij"] -= XTIJ["mi"]*T(2)["abmj"];
        }

        PerturbedCCSD(const input::Config& config, const op::TwoElectronOperator& H,
                      const CCSD<U>& T, const TwoElectronOperator<U>& x, const U omega=0)
        : Iterative(config), op::ExcitationOperator<U,2>(ccsd.uhf),
          H(H), T(T), omega(omega), D(this->uhf), Z(this->uhf), X(this->uhf),
          diis(config.get("diis"))
        {
            D(0) = 1;
            D(1)["ia"]  = moints.getIJ()["ii"];
            D(1)["ia"] -= moints.getAB()["aa"];
            D(2)["ijab"]  = moints.getIJ()["ii"];
            D(2)["ijab"] += moints.getIJ()["jj"];
            D(2)["ijab"] -= moints.getAB()["aa"];
            D(2)["ijab"] -= moints.getAB()["bb"];

            int64_t size;
            U * data;
            data = D(1)(0).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./(data[i]+omega);
            }
            data = D(1)(1).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./(data[i]+omega);
            }
            data = D(2)(0).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./(data[i]+omega);
            }
            data = D(2)(1).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./(data[i]+omega);
            }
            data = D(2)(2).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./(data[i]+omega);
            }

            op::TwoElectronOperator<U> xt(const_cast<op::TwoElectronOperator<U>&>(x),
                                          op::TwoElectronOperator<U>::AB|
                                          op::TwoElectronOperator<U>::IJ|
                                          op::TwoElectronOperator<U>::IA|
                                          op::TwoElectronOperator<U>::IJKL|
                                          op::TwoElectronOperator<U>::IJKA|
                                          op::TwoElectronOperator<U>::IAJK|
                                          op::TwoElectronOperator<U>::AIBJ);

            const tensor::SpinorbitalTensor< tensor::DistTensor<U> >& XAI = x.getAI();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FME = xt.getIA();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FAE = xt.getAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& FMI = xt.getIJ();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNEF = xt.getIJAB();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WAMEF = xt.getAIBC();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABEJ = xt.getABCI();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WABEF = xt.getABCD();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNIJ = xt.getIJKL();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMNIE = xt.getIJKA();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WMBIJ = xt.getIAJK();
            tensor::SpinorbitalTensor< tensor::DistTensor<U> >& WAMEI = xt.getAIBJ();

            TX(0) = 0;
            X(0) = 0;

            tensor::SpinorbitalTensor< tensor::DistTensor<U> > Tau(T(2));
            Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

            FME["me"] += WMNEF["mnef"]*T(1)["fn"];

            FMI["mi"] += 0.5*WMNEF["mnef"]*T(2)["efin"];
            FMI["mi"] += FME["me"]*T(1)["ei"];
            FMI["mi"] += WMNIE["mnif"]*T(1)["fn"];

            WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*Tau["efij"];
            WMNIJ["mnij"] += WMNIE["mnie"]*T(1)["ej"];

            WMNIE["mnie"] += WMNEF["mnfe"]*T(1)["fi"];

            X(1)["ai"]  = XAI["ai"];
            X(1)["ai"] -= T(1)["em"]*WAMEI["amei"];
            X(1)["ai"] += 0.5*WAMEF["amef"]*Tau["efim"];
            X(1)["ai"] -= 0.5*WMNIE["mnie"]*T(2)["aemn"];
            X(1)["ai"] += T(2)["aeim"]*FME["me"];
            X(1)["ai"] += T(1)["ei"]*FAE["ae"];
            X(1)["ai"] -= T(1)["am"]*FMI["mi"];

            FAE["ae"] -= 0.5*WMNEF["mnef"]*T(2)["afmn"];
            FAE["ae"] -= FME["me"]*T(1)["am"];
            FAE["ae"] += WAMEF["amef"]*T(1)["fm"];

            WMBIJ["mbij"] += 0.5*WAMEF["bmfe"]*Tau["efij"];
            WMBIJ["mbij"] -= WAMEI["bmej"]*T(1)["ei"];

            WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
            WAMEI["amei"] -= WAMEF["amfe"]*T(1)["fi"];
            WAMEI["amei"] += WMNIE["nmie"]*T(1)["an"];

            X(2)["abij"]  = WMNEF["ijab"];
            X(2)["abij"] += FAE["af"]*T(2)["fbij"];
            X(2)["abij"] -= FMI["ni"]*T(2)["abnj"];
            X(2)["abij"] += WABEJ["abej"]*T(1)["ei"];
            X(2)["abij"] -= WMBIJ["mbij"]*T(1)["am"];
            X(2)["abij"] += 0.5*WABEF["abef"]*Tau["efij"];
            X(2)["abij"] += 0.5*WMNIJ["mnij"]*Tau["abmn"];
            X(2)["abij"] -= WAMEI["amei"]*T(2)["ebmj"];
        }

        void _iterate()
        {
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

            /**************************************************************************
             *
             * TX(1)->TX(1) and TX(2)->TX(1)
             */
            TX(1)["ai"]  = X(1)["ai"];
            TX(1)["ai"] += FAE["ae"]*TX(1)["ei"];
            TX(1)["ai"] -= FMI["mi"]*TX(1)["am"];
            TX(1)["ai"] -= WAMEI["amei"]*TX(1)["em"];
            TX(1)["ai"] += FME["me"]*TX(2)["aeim"];
            TX(1)["ai"] += 0.5*WAMEF["amef"]*TX(2)["efim"];
            TX(1)["ai"] -= 0.5*WMNIE["mnie"]*TX(2)["aemn"];
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * TX(1)->TX(2) and TX(2)->TX(2)
             */
            TX(2)["abij"]  = X(2)["abij"];
            TX(2)["abij"] += FAE["ae"]*TX(2)["ebij"];
            TX(2)["abij"] -= FMI["mi"]*TX(2)["abmj"];
            TX(2)["abij"] += WABEJ["abej"]*TX(1)["ei"];
            TX(2)["abij"] -= WMBIJ["mbij"]*TX(1)["am"];
            TX(2)["abij"] += 0.5*WABEF["abef"]*TX(2)["efij"];
            TX(2)["abij"] += 0.5*WMNIJ["mnij"]*TX(2)["abmn"];
            TX(2)["abij"] -= WAMEI["amei"]*TX(2)["ebmj"];
            /*
             *************************************************************************/

             Z *=  D;
             Z -= TX;
            TX +=  Z;

            conv =               Z(1)(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,Z(1)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(2).reduce(CTF_OP_MAXABS));

            diis.extrapolate(TX, Z);
        }
};

}
}

#endif
