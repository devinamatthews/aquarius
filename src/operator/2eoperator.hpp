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

#ifndef _AQUARIUS_OPERATOR_2EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_2EOPERATOR_HPP_

#include "1eoperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class TwoElectronOperator : public OneElectronOperatorBase<T,TwoElectronOperator<T> >
{
    protected:
        tensor::SpinorbitalTensor<T>& ijkl;
        tensor::SpinorbitalTensor<T>& iajk;
        tensor::SpinorbitalTensor<T>& ijka;
        tensor::SpinorbitalTensor<T>& abij;
        tensor::SpinorbitalTensor<T>& ijab;
        tensor::SpinorbitalTensor<T>& aibj;
        tensor::SpinorbitalTensor<T>& aibc;
        tensor::SpinorbitalTensor<T>& abci;
        tensor::SpinorbitalTensor<T>& abcd;

        void initialize();

    public:
        enum
        {
            IJKL = 0x0010,
            IAJK = 0x0020,
            IJKA = 0x0040,
            ABIJ = 0x0080,
            IJAB = 0x0100,
            AIBJ = 0x0200,
            AIBC = 0x0400,
            ABCI = 0x0800,
            ABCD = 0x1000
        };

        TwoElectronOperator(const Arena& arena, const Space& occ, const Space& vrt, bool hermitian=true);

        TwoElectronOperator(OneElectronOperator<T>& other, int copy, bool breakhermicity=false);

        TwoElectronOperator(const OneElectronOperator<T>& other, bool breakhermicity=false);

        TwoElectronOperator(TwoElectronOperator<T>& other, int copy, bool breakhermicity=false);

        TwoElectronOperator(const TwoElectronOperator<T>& other, bool breakhermicity=false);

        T dot(bool conja, const TwoElectronOperator<T>& A, bool conjb) const;

        tensor::SpinorbitalTensor<T>& getIJKL() { return ijkl; }
        tensor::SpinorbitalTensor<T>& getIAJK() { return iajk; }
        tensor::SpinorbitalTensor<T>& getIJKA() { return ijka; }
        tensor::SpinorbitalTensor<T>& getABIJ() { return abij; }
        tensor::SpinorbitalTensor<T>& getIJAB() { return ijab; }
        tensor::SpinorbitalTensor<T>& getAIBJ() { return aibj; }
        tensor::SpinorbitalTensor<T>& getAIBC() { return aibc; }
        tensor::SpinorbitalTensor<T>& getABCI() { return abci; }
        tensor::SpinorbitalTensor<T>& getABCD() { return abcd; }

        const tensor::SpinorbitalTensor<T>& getIJKL() const { return ijkl; }
        const tensor::SpinorbitalTensor<T>& getIAJK() const { return iajk; }
        const tensor::SpinorbitalTensor<T>& getIJKA() const { return ijka; }
        const tensor::SpinorbitalTensor<T>& getABIJ() const { return abij; }
        const tensor::SpinorbitalTensor<T>& getIJAB() const { return ijab; }
        const tensor::SpinorbitalTensor<T>& getAIBJ() const { return aibj; }
        const tensor::SpinorbitalTensor<T>& getAIBC() const { return aibc; }
        const tensor::SpinorbitalTensor<T>& getABCI() const { return abci; }
        const tensor::SpinorbitalTensor<T>& getABCD() const { return abcd; }
};

}
}

#endif
