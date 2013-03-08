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

#ifndef _AQUARIUS_SCF_MOINTS_HPP_
#define _AQUARIUS_SCF_MOINTS_HPP_

#include "util/distributed.hpp"
#include "tensor/dist_tensor.hpp"
#include "tensor/spinorbital.hpp"
#include "operator/2eoperator.hpp"

#include "scf.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class MOIntegrals : public op::TwoElectronOperator<T>
{
    protected:
        tensor::DistTensor<T> *fAB_;
        tensor::DistTensor<T> *fab_;
        tensor::DistTensor<T> *fAI_;
        tensor::DistTensor<T> *fai_;
        tensor::DistTensor<T> *fIJ_;
        tensor::DistTensor<T> *fij_;
        tensor::DistTensor<T> *IJKL_;
        tensor::DistTensor<T> *IjKl_;
        tensor::DistTensor<T> *ijkl_;
        tensor::DistTensor<T> *IJKA_;
        tensor::DistTensor<T> *IjKa_;
        tensor::DistTensor<T> *iJkA_;
        tensor::DistTensor<T> *ijka_;
        tensor::DistTensor<T> *ABIJ_;
        tensor::DistTensor<T> *AbIj_;
        tensor::DistTensor<T> *abij_;
        tensor::DistTensor<T> *AIBJ_;
        tensor::DistTensor<T> *AiBj_;
        tensor::DistTensor<T> *aIbJ_;
        tensor::DistTensor<T> *aibj_;
        tensor::DistTensor<T> *ABCI_;
        tensor::DistTensor<T> *AbCi_;
        tensor::DistTensor<T> *aBcI_;
        tensor::DistTensor<T> *abci_;
        tensor::DistTensor<T> *ABCD_;
        tensor::DistTensor<T> *AbCd_;
        tensor::DistTensor<T> *abcd_;

        MOIntegrals(const scf::UHF<T>& uhf)
        : op::TwoElectronOperator<T>(uhf)
        {
            int N = uhf.getMolecule().getNumOrbitals();
            int nI = uhf.getMolecule().getNumAlphaElectrons();
            int ni = uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;

            int sizeAA[] = {nA, nA};
            int sizeaa[] = {na, na};
            int sizeAI[] = {nA, nI};
            int sizeai[] = {na, ni};
            int sizeII[] = {nI, nI};
            int sizeii[] = {ni, ni};

            int sizeIIII[] = {nI, nI, nI, nI};
            int sizeIiIi[] = {nI, ni, nI, ni};
            int sizeiiii[] = {ni, ni, ni, ni};
            int sizeIIIA[] = {nI, nI, nI, nA};
            int sizeIiIa[] = {nI, ni, nI, na};
            int sizeiIiA[] = {ni, nI, ni, nA};
            int sizeiiia[] = {ni, ni, ni, na};
            int sizeAAII[] = {nA, nA, nI, nI};
            int sizeAaIi[] = {nA, na, nI, ni};
            int sizeaaii[] = {na, na, ni, ni};
            int sizeAIAI[] = {nA, nI, nA, nI};
            int sizeAiAi[] = {nA, ni, nA, ni};
            int sizeaIaI[] = {na, nI, na, nI};
            int sizeaiai[] = {na, ni, na, ni};
            int sizeAAAI[] = {nA, nA, nA, nI};
            int sizeAaAi[] = {nA, na, nA, ni};
            int sizeaAaI[] = {na, nA, na, nI};
            int sizeaaai[] = {na, na, na, ni};
            int sizeAAAA[] = {nA, nA, nA, nA};
            int sizeAaAa[] = {nA, na, nA, na};
            int sizeaaaa[] = {na, na, na, na};

            int shapeNN[] = {NS, NS};
            int shapeNNNN[] = {NS, NS, NS, NS};
            int shapeANNN[] = {AS, NS, NS, NS};
            int shapeANAN[] = {AS, NS, AS, NS};

            fAB_ = new tensor::DistTensor<T>(this->ctf, 2, sizeAA, shapeNN, true);
            fab_ = new tensor::DistTensor<T>(this->ctf, 2, sizeaa, shapeNN, true);
            fAI_ = new tensor::DistTensor<T>(this->ctf, 2, sizeAI, shapeNN, true);
            fai_ = new tensor::DistTensor<T>(this->ctf, 2, sizeai, shapeNN, true);
            fIJ_ = new tensor::DistTensor<T>(this->ctf, 2, sizeII, shapeNN, true);
            fij_ = new tensor::DistTensor<T>(this->ctf, 2, sizeii, shapeNN, true);
            IJKL_ = new tensor::DistTensor<T>(this->ctf, 4, sizeIIII, shapeANAN, false);
            IjKl_ = new tensor::DistTensor<T>(this->ctf, 4, sizeIiIi, shapeNNNN, false);
            ijkl_ = new tensor::DistTensor<T>(this->ctf, 4, sizeiiii, shapeANAN, false);
            IJKA_ = new tensor::DistTensor<T>(this->ctf, 4, sizeIIIA, shapeANNN, false);
            IjKa_ = new tensor::DistTensor<T>(this->ctf, 4, sizeIiIa, shapeNNNN, false);
            iJkA_ = new tensor::DistTensor<T>(this->ctf, 4, sizeiIiA, shapeNNNN, false);
            ijka_ = new tensor::DistTensor<T>(this->ctf, 4, sizeiiia, shapeANNN, false);
            ABIJ_ = new tensor::DistTensor<T>(this->ctf, 4, sizeAAII, shapeANAN, false);
            AbIj_ = new tensor::DistTensor<T>(this->ctf, 4, sizeAaIi, shapeNNNN, false);
            abij_ = new tensor::DistTensor<T>(this->ctf, 4, sizeaaii, shapeANAN, false);
            AIBJ_ = new tensor::DistTensor<T>(this->ctf, 4, sizeAIAI, shapeNNNN, false);
            AiBj_ = new tensor::DistTensor<T>(this->ctf, 4, sizeAiAi, shapeNNNN, false);
            aIbJ_ = new tensor::DistTensor<T>(this->ctf, 4, sizeaIaI, shapeNNNN, false);
            aibj_ = new tensor::DistTensor<T>(this->ctf, 4, sizeaiai, shapeNNNN, false);
            ABCI_ = new tensor::DistTensor<T>(this->ctf, 4, sizeAAAI, shapeANNN, false);
            AbCi_ = new tensor::DistTensor<T>(this->ctf, 4, sizeAaAi, shapeNNNN, false);
            aBcI_ = new tensor::DistTensor<T>(this->ctf, 4, sizeaAaI, shapeNNNN, false);
            abci_ = new tensor::DistTensor<T>(this->ctf, 4, sizeaaai, shapeANNN, false);
            ABCD_ = new tensor::DistTensor<T>(this->ctf, 4, sizeAAAA, shapeANAN, false);
            AbCd_ = new tensor::DistTensor<T>(this->ctf, 4, sizeAaAa, shapeNNNN, false);
            abcd_ = new tensor::DistTensor<T>(this->ctf, 4, sizeaaaa, shapeANAN, false);

            this->ab.addSpinCase(fAB_, "A,B", "AB");
            this->ab.addSpinCase(fab_, "a,b", "ab");

            this->ai.addSpinCase(fAI_, "A,I", "AI");
            this->ai.addSpinCase(fai_, "a,i", "ai");

            this->ij.addSpinCase(fIJ_, "I,J", "IJ");
            this->ij.addSpinCase(fij_, "i,j", "ij");

            this->ijkl.addSpinCase(IJKL_, "IJ,KL", "IJKL");
            this->ijkl.addSpinCase(IjKl_, "Ij,Kl", "IjKl");
            this->ijkl.addSpinCase(ijkl_, "ij,kl", "ijkl");

            this->ijka.addSpinCase(IJKA_, "IJ,KA", "IJKA");
            this->ijka.addSpinCase(IjKa_, "Ij,Ka", "IjKa");
            this->ijka.addSpinCase(iJkA_, "iJ,kA", "iJkA");
            this->ijka.addSpinCase(ijka_, "ij,ka", "ijka");

            this->abij.addSpinCase(ABIJ_, "AB,IJ", "ABIJ");
            this->abij.addSpinCase(AbIj_, "Ab,Ij", "AbIj");
            this->abij.addSpinCase(abij_, "ab,ij", "abij");

            this->aibj.addSpinCase(AIBJ_, "AI,BJ", "AIBJ");
            this->aibj.addSpinCase(AiBj_, "Ai,Bj", "AiBj");
            this->aibj.addSpinCase(aIbJ_, "aI,bJ", "aIbJ");
            this->aibj.addSpinCase(aibj_, "ai,bj", "aibj");
            this->aibj.addSpinCase(*AbIj_, "Ai,bJ", "AbJi", -1.0);
            this->aibj.addSpinCase(*AbIj_, "aI,Bj", "BaIj", -1.0);

            this->abci.addSpinCase(ABCI_, "AB,CI", "ABCI");
            this->abci.addSpinCase(AbCi_, "Ab,Ci", "AbCi");
            this->abci.addSpinCase(aBcI_, "aB,cI", "aBcI");
            this->abci.addSpinCase(abci_, "ab,ci", "abci");

            this->abcd.addSpinCase(ABCD_, "AB,CD", "ABCD");
            this->abcd.addSpinCase(AbCd_, "Ab,Cd", "AbCd");
            this->abcd.addSpinCase(abcd_, "ab,cd", "abcd");

            const T* ea = uhf.getAlphaEigenvalues();
            const T* eb = uhf.getBetaEigenvalues();
            int rank = this->comm.Get_rank();
            int np = this->comm.Get_size();

            {
                std::vector< tkv_pair<T> > pairs;

                for (int i = 0;i < nI;i++)
                {
                    if (i%np == rank) pairs.push_back(tkv_pair<T>(i+i*nI, ea[i]));
                }

                fIJ_->writeRemoteData(pairs.size(), pairs.data());
            }

            {
                std::vector< tkv_pair<T> > pairs;

                for (int i = 0;i < ni;i++)
                {
                    if (i%np == rank) pairs.push_back(tkv_pair<T>(i+i*ni, eb[i]));
                }

                fij_->writeRemoteData(pairs.size(), pairs.data());
            }

            {
                std::vector< tkv_pair<T> > pairs;

                for (int i = 0;i < nA;i++)
                {
                    if (i%np == rank) pairs.push_back(tkv_pair<T>(i+i*nA, ea[i+nI]));
                }

                fAB_->writeRemoteData(pairs.size(), pairs.data());
            }

            {
                std::vector< tkv_pair<T> > pairs;

                for (int i = 0;i < na;i++)
                {
                    if (i%np == rank) pairs.push_back(tkv_pair<T>(i+i*na, eb[i+ni]));
                }

                fab_->writeRemoteData(pairs.size(), pairs.data());
            }
        }
};

}
}

#endif
