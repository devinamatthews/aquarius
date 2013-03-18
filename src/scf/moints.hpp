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
        tensor::DistTensor<T> *AibJ_;
        tensor::DistTensor<T> *ABCI_;
        tensor::DistTensor<T> *AbCi_;
        tensor::DistTensor<T> *aBcI_;
        tensor::DistTensor<T> *abci_;
        tensor::DistTensor<T> *ABCD_;
        tensor::DistTensor<T> *AbCd_;
        tensor::DistTensor<T> *abcd_;

        MOIntegrals(const scf::UHF<T>& uhf)
        : op::TwoElectronOperator<T>(uhf, true)
        {
            int N = uhf.getMolecule().getNumOrbitals();
            int nI = uhf.getMolecule().getNumAlphaElectrons();
            int ni = uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;

            fAB_ = &this->ab(0);
            fab_ = &this->ab(1);

            fAI_ = &this->ai(0);
            fai_ = &this->ai(1);

            fIJ_ = &this->ij(0);
            fij_ = &this->ij(1);

            IJKL_ = &this->ijkl(0);
            IjKl_ = &this->ijkl(1);
            ijkl_ = &this->ijkl(2);

            IJKA_ = &this->ijka(0);
            IjKa_ = &this->ijka(1);
            iJkA_ = &this->ijka(2);
            ijka_ = &this->ijka(3);

            AIBJ_ = &this->aibj(0);
            AiBj_ = &this->aibj(1);
            aIbJ_ = &this->aibj(2);
            aibj_ = &this->aibj(3);
            AibJ_ = &this->aibj(4);

            ABIJ_ = &this->abij(0);
            AbIj_ = &this->abij(1);
            abij_ = &this->abij(2);

            ABCI_ = &this->abci(0);
            AbCi_ = &this->abci(1);
            aBcI_ = &this->abci(2);
            abci_ = &this->abci(3);

            ABCD_ = &this->abcd(0);
            AbCd_ = &this->abcd(1);
            abcd_ = &this->abcd(2);

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
