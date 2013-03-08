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

#ifndef _AQUARIUS_CC_CCSDT_HPP_
#define _AQUARIUS_CC_CCSDT_HPP_

#include "mpi.h"

#include "time/time.hpp"
#include "tensor/spinorbital.hpp"
#include "tensor/dist_tensor.hpp"
#include "scf/moints.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/exponentialoperator.hpp"
#include "convergence/diis.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSDT : public Iterative, public op::ExponentialOperator<U,3>
{
    protected:
        op::ExcitationOperator<U,3>& T;
        op::ExcitationOperator<U,3> D, Z;
        op::TwoElectronOperator<U>& moints;
        convergence::DIIS< op::ExcitationOperator<U,3> > diis;

    public:
        CCSDT(const input::Config& config, op::TwoElectronOperator<U>& moints)
        : Iterative(config), ExcitationOperator<U,3>(moints.getSCF(), false),
          T(*this), D(moints.getSCF(), false), Z(moints.getSCF(), false),
          moints(moints), diis(config.get("diis"))
        {
            D[0] = 1;
            D[1]["ai"]  = moints.getIJ()["ii"];
            D[1]["ai"] -= moints.getAB()["aa"];
            D[2]["abij"]  = moints.getIJ()["ii"];
            D[2]["abij"] += moints.getIJ()["jj"];
            D[2]["abij"] -= moints.getAB()["aa"];
            D[2]["abij"] -= moints.getAB()["bb"];
            D[3]["abcijk"]  = moints.getIJ()["ii"];
            D[3]["abcijk"] -= moints.getIJ()["jj"];
            D[3]["abcijk"] -= moints.getIJ()["kk"];
            D[3]["abcijk"] += moints.getAB()["aa"];
            D[3]["abcijk"] += moints.getAB()["bb"];
            D[3]["abcijk"] += moints.getAB()["cc"];

            int64_t size;
            T * data;
            data = D[1].getSpinCase(0).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./data[i];
            }
            data = D[1].getSpinCase(1).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./data[i];
            }
            data = D[2].getSpinCase(0).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./data[i];
            }
            data = D[2].getSpinCase(1).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./data[i];
            }
            data = D[2].getSpinCase(2).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./data[i];
            }
            data = D[3].getSpinCase(0).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./data[i];
            }
            data = D[3].getSpinCase(1).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./data[i];
            }
            data = D[3].getSpinCase(2).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./data[i];
            }
            data = D[3].getSpinCase(3).getRawData(size);
            for (int i=0; i<size; i++){
              if (fabs(data[i]) > DBL_MIN)
                data[i] = 1./data[i];
            }

            Z[0] = 0;
            T[0] = 0;
            T[1] = moints.getAI()*D[1];
            T[2] = moints.getABIJ()*D[2];
            T[3] = 0;

            tensor::SpinorbitalTensor< tensor::DistTensor<U> > Tau(T[2]);
            Tau["abij"] += 0.5*T[1]["ai"]*T[1]["bj"];

            energy = scalar(moints.getAI()*T[1]) + 0.25*scalar(moints.getABIJ()*Tau);

            conv =          conv,T[1].getSpinCase(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,T[1].getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,T[2].getSpinCase(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,T[2].getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,T[2].getSpinCase(2).reduce(CTF_OP_MAXABS));
        }

        const op::TwoElectronOperator<U>& getMOIntegrals() const
        {
            return moints;
        }

        void _iterate()
        {
            op::TwoElectronOperator H(moints, op::TwoElectronOperator::AB|
                                              op::TwoElectronOperator::IJ|
                                              op::TwoElectronOperator::IA|
                                              op::TwoElectronOperator::AIBC|
                                              op::TwoElectronOperator::ABCI|
                                              op::TwoElectronOperator::ABCD|
                                              op::TwoElectronOperator::IJKL|
                                              op::TwoElectronOperator::IJKA|
                                              op::TwoElectronOperator::IAJK|
                                              op::TwoElectronOperator::AIBJ);

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

            FAE["aa"] = 0.0;
            FMI["ii"] = 0.0;

            tensor::SpinorbitalTensor< tensor::DistTensor<U> > Tau(T[2]);
            Tau["abij"] += 0.5*T[1]["ai"]*T[1]["bj"];

            /**************************************************************************
             *
             * Intermediates for T[1]->T[1] and T[2]->T[1]
             */
            PROFILE_SECTION(calc_FEM)
            FME["me"] = WMNEF["mnef"]*T[1]["fn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_FMI)
            FMI["mi"] = 0.5*WMNEF["nmef"]*T[2]["efni"];
            FMI["mi"] += FME["me"]*T[1]["ei"];
            FMI["mi"] += WMNIE["mnif"]*T[1]["fn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WMNIE)
            WMNIE["mnie"] += WMNEF["mnfe"]*T[1]["fi"];
            PROFILE_STOP
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * T[1]->T[1] and T[2]->T[1]
             */
            PROFILE_SECTION(calc_FAI)
            Z[1]["ai"] = moints.getAI()["ai"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[1]_IN_T[1]_RING)
            Z[1]["ai"] -= T[1]["em"]*WAMEI["amei"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[2]_IN_T[1]_ABCI)
            Z[1]["ai"] += 0.5*WAMEF["amef"]*Tau["efim"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[2]_IN_T[1]_IJKA)
            Z[1]["ai"] -= 0.5*WMNIE["mnie"]*T[2]["aemn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[2]_IN_T[1]_FME)
            Z[1]["ai"] += T[2]["aeim"]*FME["me"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[1]_IN_T[1]_FAE)
            Z[1]["ai"] += T[1]["ei"]*FAE["ae"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[1]_IN_T[1]_FMI)
            Z[1]["ai"] -= T[1]["am"]*FMI["mi"];
            PROFILE_STOP
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * Intermediates for T[1]->T[2] and T[2]->T[2]
             */
            PROFILE_SECTION(calc_FAE)
            FAE["ae"] = -0.5*WMNEF["mnfe"]*T[2]["famn"];
            FAE["ae"] -= FME["me"]*T[1]["am"];
            FAE["ae"] += WABEJ["efan"]*T[1]["fn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WIJKL)
            WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*Tau["efij"];
            WMNIJ["mnij"] += WMNIE["mnie"]*T[1]["ej"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WAIJK)
            WMBIJ["mbij"] += 0.5*WABEJ["febm"]*Tau["efij"];
            WMBIJ["mbij"] -= WAMEI["amej"]*T[1]["ei"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WMBEJ)
            WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T[2]["afin"];
            WAMEI["amei"] -= WABEJ["feam"]*T[1]["fi"];
            WAMEI["amei"] += WMNIE["nmie"]*T[1]["an"];
            PROFILE_STOP
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * T[1]->T[2] and T[2]->T[2]
             */
            PROFILE_SECTION(calc_WMNEF)
            Z[2]["abij"] = WMNEF["ijab"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[2]_IN_T[2]_FAE)
            Z[2]["abij"] += FAE["af"]*T[2]["fbij"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[2]_IN_T[2]_FMI)
            Z[2]["abij"] -= FMI["ni"]*T[2]["abnj"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[1]_IN_T[2]_ABCI)
            Z[2]["abij"] += WABEJ["abej"]*T[1]["ei"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[1]_IN_T[2]_IJKA)
            Z[2]["abij"] -= WMBIJ["mbij"]*T[1]["am"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[2]_IN_T[2]_ABCD)
            Z[2]["abij"] += 0.5*WABEF["abef"]*Tau["efij"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[2]_IN_T[2]_IJKL)
            Z[2]["abij"] += 0.5*WMNIJ["mnij"]*Tau["abmn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[2]_IN_T[2]_RING)
            Z[2]["abij"] -= WAMEI["amei"]*T[2]["ebmj"];
            PROFILE_STOP
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * Intermediates for CCSDT
             */
            PROFILE_SECTION(calc_WAIJK2)
            WMBIJ["mbij"] += WMNIE["mnie"]*T[2]["bejn"];
            WMBIJ["mbij"] -= WMNIJ["mnij"]*T[1]["bn"];
            WMBIJ["mbij"] += FME["me"]*T[2]["ebij"];
            //WMBIJ["mbij"] += 0.5*WMNEF["mnef"]*T[3]["efbinj"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WMBEJ2)
            WAMEI["amei"] *= 2;
            WAMEI["amei"] -= moints.getAIBJ()["amei"];
            WAMEI["amei"] += WAMEF["amfe"]*T[1]["fi"];
            WAMEI["amei"] -= 1.5*WMNIE["nmie"]*T[1]["an"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WABCI)
            WABEJ["abej"] += WAMEF["amef"]*T[2]["fbmj"];
            WABEJ["abej"] += 0.5*WMNIE["nmje"]*T[2]["abmn"];
            WABEJ["abej"] += WABEF["abef"]*T[1]["fj"];
            WABEJ["abej"] -= WAMEI["amej"]*T[1]["bm"];
            //WABEJ["abej"] -= 0.5*WMNEF["mnef"]*T[3]["afbmnj"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WMBEJ3)
            WAMEI["amei"] += 0.5*WMNIE["nmie"]*T[1]["an"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WABCD)
            WABEF["abef"] -= WAMEF["amef"]*T[1]["bm"];
            WABEF["abef"] += 0.5*WMNEF["mnef"]*Tau["abmn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WAIBC)
            WAMEF["amef"] -= WMNEF["nmef"]*T[1]["an"];
            PROFILE_STOP
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * CCSDT Iteration
             */
            //PROFILE_SECTION(calc_T[2]_IN_T[3]_ABCI)
            //Z[3]["abcijk"] = WABEJ["bcek"]*T[2]["aeij"];
            //PROFILE_STOP

            PROFILE_SECTION(calc_T[2]_IN_T[3]_IJKA)
            //Z[3]["abcijk"] -= WAMIJ["cmkj"]*T[2]["abim"];
            Z[3]["abcijk"] = moints.getIJKA()["jkmc"]*T[2]["abim"];
            PROFILE_STOP

            //PROFILE_SECTION(calc_T[3]_IN_T[2]_ABCI)
            //Z[2]["abij"] += 0.5*WAMEF["bmef"]*T[3]["aefijm"];
            //PROFILE_STOP

            PROFILE_SECTION(calc_T[3]_IN_T[2]_IJKA)
            //Z[2]["abij"] -= 0.5*WMNIE["mnje"]*T[3]["abeimn"];
            Z[2]["abij"] -= 0.5*moints.getIJKA()["mnje"]*T[3]["abeimn"];
            PROFILE_STOP

            /*
            PROFILE_SECTION(calc_T[3]_IN_T[2]_FME)
            Z[2]["abij"] += FME["me"]*T[3]["abeijm"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[3]_IN_T[1])
            Z[1]["ai"] += 0.25*WMNEF["mnef"]*T[3]["aefimn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[3]_IN_T[3]_FAE)
            Z[3]["abcijk"] += FAE["ce"]*T[3]["abeijk"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[3]_IN_T[3]_FMI)
            Z[3]["abcijk"] -= FMI["mk"]*T[3]["abcijm"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[3]_IN_T[3]_ABCD)
            Z[3]["abcijk"] += 0.5*WABEF["abef"]*T[3]["efcijk"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[3]_IN_T[3]_IJKL)
            Z[3]["abcijk"] += 0.5*WMNIJ["mnij"]*T[3]["abcmnk"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T[3]_IN_T[3]_RING)
            Z[3]["abcijk"] -= WAMEI["amei"]*T[3]["ebcmjk"];
            PROFILE_STOP
            */
            /*
             **************************************************************************/

            PROFILE_SECTION(calc_EN)

            Z *= D;
            Z -= T;
            T += Z;

            Tau["abij"]   = T[2]["abij"];
            Tau["abij"]  += 0.5*T[1]["ai"]*T[1]["bj"];
            energy = scalar(moints.getAI()*T[1]) + 0.25*scalar(moints.getABIJ()*Tau);

            conv =                Z[1].getSpinCase(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv, Z[1].getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv, Z[2].getSpinCase(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv, Z[2].getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv, Z[2].getSpinCase(2).reduce(CTF_OP_MAXABS));
            conv = std::max(conv, Z[3].getSpinCase(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv, Z[3].getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv, Z[3].getSpinCase(2).reduce(CTF_OP_MAXABS));
            conv = std::max(conv, Z[3].getSpinCase(3).reduce(CTF_OP_MAXABS));

            diis->extrapolate(T, Z);

            PROFILE_STOP
        }
};

}
}

#endif
