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

#ifndef _AQUARIUS_CC_CCSD_HPP_
#define _AQUARIUS_CC_CCSD_HPP_

#include <cfloat>

#include "time/time.hpp"
#include "tensor/spinorbital.hpp"
#include "tensor/dist_tensor.hpp"
#include "scf/moints.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/exponentialoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "convergence/diis.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class CCSD : public Iterative, public op::ExponentialOperator<U,2>
{
    protected:
        op::ExponentialOperator<U,2>& T;
        op::ExcitationOperator<U,2> D, Z;
        op::TwoElectronOperator<U>& moints;
        convergence::DIIS< op::ExcitationOperator<U,2> > diis;

    public:
        CCSD(const input::Config& config, op::TwoElectronOperator<U>& moints)
        : Iterative(config), op::ExponentialOperator<U,2>(moints.getSCF()),
          T(*this), D(moints.getSCF()), Z(moints.getSCF()),
          moints(moints), diis(config.get("diis"))
        {
            D[0] = 1;
            D[1]["ai"]  = moints.getIJ()["ii"];
            D[1]["ai"] -= moints.getAB()["aa"];
            D[2]["abij"]  = moints.getIJ()["ii"];
            D[2]["abij"] += moints.getIJ()["jj"];
            D[2]["abij"] -= moints.getAB()["aa"];
            D[2]["abij"] -= moints.getAB()["bb"];

            int64_t size;
            U * data;
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

            Z[0] = 0;
            T[0] = 0;
            T[1] = moints.getAI()*D[1];
            T[2] = moints.getABIJ()*D[2];

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
            op::TwoElectronOperator<U> H(moints, op::TwoElectronOperator<U>::AB|
                                                 op::TwoElectronOperator<U>::IJ|
                                                 op::TwoElectronOperator<U>::IA|
                                                 op::TwoElectronOperator<U>::IJKL|
                                                 op::TwoElectronOperator<U>::IJKA|
                                                 op::TwoElectronOperator<U>::IAJK|
                                                 op::TwoElectronOperator<U>::AIBJ);

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
            FME["me"] += WMNEF["mnef"]*T[1]["fn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_FMI)
            FMI["mi"] += 0.5*WMNEF["mnef"]*T[2]["efin"];
            FMI["mi"] += FME["me"]*T[1]["ei"];
            FMI["mi"] += WMNIE["mnif"]*T[1]["fn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WIJKL)
            WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*Tau["efij"];
            WMNIJ["mnij"] += WMNIE["mnie"]*T[1]["ej"];
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

            PROFILE_SECTION(calc_T1_IN_T1_RING)
            Z[1]["ai"] -= T[1]["em"]*WAMEI["amei"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T1_ABCI)
            Z[1]["ai"] += 0.5*WAMEF["amef"]*Tau["efim"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T1_IJKA)
            Z[1]["ai"] -= 0.5*WMNIE["mnie"]*T[2]["aemn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T1_FME)
            Z[1]["ai"] += T[2]["aeim"]*FME["me"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T1_IN_T1_FAE)
            Z[1]["ai"] += T[1]["ei"]*FAE["ae"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T1_IN_T1_FMI)
            Z[1]["ai"] -= T[1]["am"]*FMI["mi"];
            PROFILE_STOP
            /*
             *************************************************************************/

            /**************************************************************************
             *
             * Intermediates for T[1]->T[2] and T[2]->T[2]
             */
            PROFILE_SECTION(calc_FAE)
            FAE["ae"] -= 0.5*WMNEF["mnef"]*T[2]["afmn"];
            FAE["ae"] -= FME["me"]*T[1]["am"];
            FAE["ae"] += WAMEF["amef"]*T[1]["fm"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WAIJK)
            WMBIJ["mbij"] += 0.5*WAMEF["bmfe"]*Tau["efij"];
            WMBIJ["mbij"] -= WAMEI["bmej"]*T[1]["ei"];
            PROFILE_STOP

            PROFILE_SECTION(calc_WMBEJ)
            WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T[2]["afin"];
            WAMEI["amei"] -= WAMEF["amfe"]*T[1]["fi"];
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

            PROFILE_SECTION(calc_T2_IN_T2_FAE)
            Z[2]["abij"] += FAE["af"]*T[2]["fbij"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T2_FMI)
            Z[2]["abij"] -= FMI["ni"]*T[2]["abnj"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T1_IN_T2_ABCI)
            Z[2]["abij"] += WABEJ["abej"]*T[1]["ei"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T1_IN_T2_IJKA)
            Z[2]["abij"] -= WMBIJ["mbij"]*T[1]["am"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T2_ABCD)
            Z[2]["abij"] += 0.5*WABEF["abef"]*Tau["efij"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T2_IJKL)
            Z[2]["abij"] += 0.5*WMNIJ["mnij"]*Tau["abmn"];
            PROFILE_STOP

            PROFILE_SECTION(calc_T2_IN_T2_RING)
            Z[2]["abij"] -= WAMEI["amei"]*T[2]["ebmj"];
            PROFILE_STOP
            /*
             *************************************************************************/

            PROFILE_SECTION(calc_EN)

            Z *= D;
            Z -= T;
            T += Z;

            Tau["abij"]  = T[2]["abij"];
            Tau["abij"] += 0.5*T[1]["ai"]*T[1]["bj"];
            energy = scalar(moints.getAI()*T[1]) + 0.25*scalar(moints.getABIJ()*Tau);

            conv =               Z[1].getSpinCase(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,Z[1].getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z[2].getSpinCase(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z[2].getSpinCase(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z[2].getSpinCase(2).reduce(CTF_OP_MAXABS));

            diis.extrapolate(T, Z);

            PROFILE_STOP
        }

        double getProjectedS2() const
        {
            int N = this->uhf.getMolecule().getNumOrbitals();
            int nI = this->uhf.getMolecule().getNumAlphaElectrons();
            int ni = this->uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;

            int shapeNN[] = {NS,NS};
            int shapeNNNN[] = {NS,NS,NS,NS};
            int sizeAI[] = {nA,nI};
            int sizeAi[] = {nA,ni};
            int sizeaI[] = {na,nI};
            int sizeai[] = {na,ni};
            int sizeIi[] = {nI,ni};
            int sizeIn[] = {nI,N};
            int sizein[] = {ni,N};
            int sizeAaIi[] = {nA,na,nI,ni};

            const tensor::DistTensor<U>& CA = this->uhf.getCA();
            const tensor::DistTensor<U>& Ca = this->uhf.getCa();
            const tensor::DistTensor<U>& CI = this->uhf.getCI();
            const tensor::DistTensor<U>& Ci = this->uhf.getCi();
            const tensor::DistTensor<U>& S = this->uhf.getOverlap();
            tensor::DistTensor<U> DAI(this->ctf, 2, sizeAI, shapeNN, false);
            tensor::DistTensor<U> DAi(this->ctf, 2, sizeAi, shapeNN, false);
            tensor::DistTensor<U> DaI(this->ctf, 2, sizeaI, shapeNN, false);
            tensor::DistTensor<U> Dai(this->ctf, 2, sizeai, shapeNN, false);
            tensor::DistTensor<U> DIj(this->ctf, 2, sizeIi, shapeNN, false);
            tensor::DistTensor<U> DAbIj(this->ctf, 4, sizeAaIi, shapeNNNN, false);
            tensor::DistTensor<U> tmp1(this->ctf, 2, sizeIn, shapeNN, false);
            tensor::DistTensor<U> tmp2(this->ctf, 2, sizein, shapeNN, false);

            tmp1["Iq"] = CI["pI"]*S["pq"];
            DIj["Ij"] = tmp1["Iq"]*Ci["qj"];
            DaI["aI"] = tmp1["Iq"]*Ca["qa"];

            tmp2["iq"] = Ci["pi"]*S["pq"];
            DAi["Ai"] = tmp2["iq"]*CA["qA"];

            DAI["AI"] = DAi["Aj"]*DIj["Ij"];
            Dai["ai"] = DaI["aJ"]*DIj["Ji"];
            DAbIj["AbIj"] = DAi["Aj"]*DaI["bI"];

            const tensor::DistTensor<U>& T1A = T[1].getSpinCase(0);
            const tensor::DistTensor<U>& T1B = T[1].getSpinCase(1);
            tensor::DistTensor<U> TauAB(T[2].getSpinCase(1));

            TauAB["AbIj"] += T1A["AI"]*T1B["bj"];

            double S2 = this->uhf.getS2();

            S2 -= tensor::scalar(DAI*T1A);
            S2 -= tensor::scalar(Dai*T1B);
            S2 -= tensor::scalar(DAbIj*TauAB);

            return S2;
        }

        double getProjectedMultiplicity() const
        {
            return sqrt(1+4*getProjectedS2());
        }
};

}
}

#endif
