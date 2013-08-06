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

#ifndef _AQUARIUS_OPERATOR_MOINTS_HPP_
#define _AQUARIUS_OPERATOR_MOINTS_HPP_

#include "2eoperator.hpp"

#include "scf/scf.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class MOIntegrals : public TwoElectronOperator<T>
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
        tensor::DistTensor<T> *aIBj_;
        tensor::DistTensor<T> *ABCI_;
        tensor::DistTensor<T> *AbCi_;
        tensor::DistTensor<T> *aBcI_;
        tensor::DistTensor<T> *abci_;
        tensor::DistTensor<T> *ABCD_;
        tensor::DistTensor<T> *AbCd_;
        tensor::DistTensor<T> *abcd_;

        MOIntegrals(const scf::UHF<T>& uhf);
};

}
}

#endif
