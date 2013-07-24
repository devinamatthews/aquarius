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

#include "moints.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::scf;
using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::slide;

template <typename T>
MOIntegrals<T>::MOIntegrals(const UHF<T>& uhf)
: TwoElectronOperator<T>(uhf, true)
{
    int N = uhf.getMolecule().getNumOrbitals();
    int nI = uhf.getMolecule().getNumAlphaElectrons();
    int ni = uhf.getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;

    fAB_ = &this->ab(1,0,1,0);
    fab_ = &this->ab(0,0,0,0);

    fAI_ = &this->ai(1,0,0,1);
    fai_ = &this->ai(0,0,0,0);

    fIJ_ = &this->ij(0,1,0,1);
    fij_ = &this->ij(0,0,0,0);

    IJKL_ = &this->ijkl(0,2,0,2);
    IjKl_ = &this->ijkl(0,1,0,1);
    ijkl_ = &this->ijkl(0,0,0,0);

    IJKA_ = &this->ijka(0,2,1,1);
    IjKa_ = &this->ijka(0,1,0,1);
    iJkA_ = &this->ijka(0,1,1,0);
    ijka_ = &this->ijka(0,0,0,0);

    AIBJ_ = &this->aibj(1,1,1,1);
    AiBj_ = &this->aibj(1,0,1,0);
    aIbJ_ = &this->aibj(0,1,0,1);
    aibj_ = &this->aibj(0,0,0,0);
    AibJ_ = &this->aibj(1,0,0,1);
    aIBj_ = &this->aibj(0,1,1,0);

    ABIJ_ = &this->abij(2,0,0,2);
    AbIj_ = &this->abij(1,0,0,1);
    abij_ = &this->abij(0,0,0,0);

    ABCI_ = &this->abci(2,0,1,1);
    AbCi_ = &this->abci(1,0,1,0);
    aBcI_ = &this->abci(1,0,0,1);
    abci_ = &this->abci(0,0,0,0);

    ABCD_ = &this->abcd(2,0,2,0);
    AbCd_ = &this->abcd(1,0,1,0);
    abcd_ = &this->abcd(0,0,0,0);

    const T* ea = uhf.getAlphaEigenvalues();
    const T* eb = uhf.getBetaEigenvalues();
    int rank = this->comm.Get_rank();
    int np = this->comm.Get_size();

    {
        vector< tkv_pair<T> > pairs;

        for (int i = 0;i < nI;i++)
        {
            if (i%np == rank) pairs.push_back(tkv_pair<T>(i+i*nI, ea[i]));
        }

        fIJ_->writeRemoteData(pairs.size(), pairs.data());
    }

    {
        vector< tkv_pair<T> > pairs;

        for (int i = 0;i < ni;i++)
        {
            if (i%np == rank) pairs.push_back(tkv_pair<T>(i+i*ni, eb[i]));
        }

        fij_->writeRemoteData(pairs.size(), pairs.data());
    }

    {
        vector< tkv_pair<T> > pairs;

        for (int i = 0;i < nA;i++)
        {
            if (i%np == rank) pairs.push_back(tkv_pair<T>(i+i*nA, ea[i+nI]));
        }

        fAB_->writeRemoteData(pairs.size(), pairs.data());
    }

    {
        vector< tkv_pair<T> > pairs;

        for (int i = 0;i < na;i++)
        {
            if (i%np == rank) pairs.push_back(tkv_pair<T>(i+i*na, eb[i+ni]));
        }

        fab_->writeRemoteData(pairs.size(), pairs.data());
    }
}

INSTANTIATE_SPECIALIZATIONS(MOIntegrals);
