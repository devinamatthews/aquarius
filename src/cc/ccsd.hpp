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

#include "time/time.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/exponentialoperator.hpp"
#include "operator/stexcitationoperator.hpp"
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
        op::TwoElectronOperator<U>& H;
        convergence::DIIS< op::ExcitationOperator<U,2> > diis;

    public:
        CCSD(const input::Config& config, op::TwoElectronOperator<U>& H)
        : Iterative(config), op::ExponentialOperator<U,2>(H.getSCF()),
          T(*this), D(H.getSCF()), Z(H.getSCF()),
          H(H), diis(config.get("diis"))
        {
            D(0) = 1;
            D(1)["ai"]  = H.getIJ()["ii"];
            D(1)["ai"] -= H.getAB()["aa"];
            D(2)["abij"]  = H.getIJ()["ii"];
            D(2)["abij"] += H.getIJ()["jj"];
            D(2)["abij"] -= H.getAB()["aa"];
            D(2)["abij"] -= H.getAB()["bb"];

            D = 1/D;

            Z(0) = 0;
            T(0) = 0;
            T(1) = H.getAI()*D(1);
            T(2) = H.getABIJ()*D(2);

            tensor::SpinorbitalTensor< tensor::DistTensor<U> > Tau(T(2));
            Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

            energy = scalar(H.getAI()*T(1)) + 0.25*scalar(H.getABIJ()*Tau);

            conv =               T(1)(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,T(1)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,T(2)(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,T(2)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,T(2)(2).reduce(CTF_OP_MAXABS));
        }

        void _iterate()
        {
            op::STExcitationOperator<U,2>::transform(H, T, Z);
            //Z(0) = 0;

            Z *= D;
            T += Z;

            tensor::SpinorbitalTensor<tensor::DistTensor<U> > Tau(T(2));
            Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

            energy = scalar(H.getAI()*T(1)) + 0.25*scalar(H.getABIJ()*Tau);

            conv =               Z(1)(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,Z(1)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(2).reduce(CTF_OP_MAXABS));

            diis.extrapolate(T, Z);
        }

        static double getProjectedS2(const scf::UHF<U>& uhf,
                                     const tensor::SpinorbitalTensor<tensor::DistTensor<U> >& T1,
                                     const tensor::SpinorbitalTensor<tensor::DistTensor<U> >& T2)
        {
            tCTF_World<U>& ctf = T1(0).ctf;

            int N = uhf.getMolecule().getNumOrbitals();
            int nI = uhf.getMolecule().getNumAlphaElectrons();
            int ni = uhf.getMolecule().getNumBetaElectrons();
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

            const tensor::DistTensor<U>& CA = uhf.getCA();
            const tensor::DistTensor<U>& Ca = uhf.getCa();
            const tensor::DistTensor<U>& CI = uhf.getCI();
            const tensor::DistTensor<U>& Ci = uhf.getCi();
            const tensor::DistTensor<U>& S = uhf.getOverlap();
            tensor::DistTensor<U> DAI(ctf, 2, sizeAI, shapeNN, false);
            tensor::DistTensor<U> DAi(ctf, 2, sizeAi, shapeNN, false);
            tensor::DistTensor<U> DaI(ctf, 2, sizeaI, shapeNN, false);
            tensor::DistTensor<U> Dai(ctf, 2, sizeai, shapeNN, false);
            tensor::DistTensor<U> DIj(ctf, 2, sizeIi, shapeNN, false);
            tensor::DistTensor<U> DAbIj(ctf, 4, sizeAaIi, shapeNNNN, false);
            tensor::DistTensor<U> tmp1(ctf, 2, sizeIn, shapeNN, false);
            tensor::DistTensor<U> tmp2(ctf, 2, sizein, shapeNN, false);

            tmp1["Iq"] = CI["pI"]*S["pq"];
            DIj["Ij"] = tmp1["Iq"]*Ci["qj"];
            DaI["aI"] = tmp1["Iq"]*Ca["qa"];

            tmp2["iq"] = Ci["pi"]*S["pq"];
            DAi["Ai"] = tmp2["iq"]*CA["qA"];

            DAI["AI"] = DAi["Aj"]*DIj["Ij"];
            Dai["ai"] = DaI["aJ"]*DIj["Ji"];
            DAbIj["AbIj"] = DAi["Aj"]*DaI["bI"];

            const tensor::DistTensor<U>& T1A = T1(1,0,0,1);
            const tensor::DistTensor<U>& T1B = T1(0,0,0,0);
            tensor::DistTensor<U> TauAB(T2(1,0,0,1));

            TauAB["AbIj"] += T1A["AI"]*T1B["bj"];

            U S2 = uhf.getS2();

            U S2T11 = -scalar(DAI*T1A);
            U S2T12 = -scalar(Dai*T1B);
            U S2T2 = -scalar(DAbIj*TauAB);

            return std::abs(S2+S2T11+S2T12+S2T2);
        }

        double getProjectedS2() const
        {
            return getProjectedS2(this->uhf, T(1), T(2));
        }

        static double getProjectedMultiplicity(const scf::UHF<U>& uhf,
                                               const tensor::SpinorbitalTensor<tensor::DistTensor<U> >& T1,
                                               const tensor::SpinorbitalTensor<tensor::DistTensor<U> >& T2)
        {
            return sqrt(1+4*getProjectedS2(uhf, T1, T2));
        }

        double getProjectedMultiplicity() const
        {
            return getProjectedMultiplicity(this->uhf, T(1), T(2));
        }
};

}
}

#endif
