/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must reRin the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documenRtion and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANRBILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENRL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DAR, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_OPERATOR_PERTURBEDST2EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_PERTURBEDST2EOPERATOR_HPP_

#include "2eoperator.hpp"
#include "st2eoperator.hpp"

namespace aquarius
{
namespace op
{

/*
 *      _A    -T  A  T    _   A      A  T     _  A
 * Form X  = e   X  e  + [X, T ] = (X  e ) + (X T ) , up to two-electron terms
 *                                        c        c
 */
template <typename U, int nex> class PerturbedSTTwoElectronOperator;

template <typename U>
class PerturbedSTTwoElectronOperator<U,2> : public STTwoElectronOperator<U,2>
{
    protected:
        const STTwoElectronOperator<U,2>& X;
        const ExcitationOperator<U,2>& TA;

        void initialize(const STTwoElectronOperator<U,2>& X,
                        const ExcitationOperator<U,2>& T,
                        const ExcitationOperator<U,2>& TA);

    public:
        PerturbedSTTwoElectronOperator(const STTwoElectronOperator<U,2>& X, const OneElectronOperator<U>& XA,
                                       const ExcitationOperator<U,2>& T, const ExcitationOperator<U,2>& TA);

        PerturbedSTTwoElectronOperator(const STTwoElectronOperator<U,2>& X, const TwoElectronOperator<U>& XA,
                                       const ExcitationOperator<U,2>& T, const ExcitationOperator<U,2>& TA);

        void contract(const ExcitationOperator<U,2>& R, ExcitationOperator<U,2>& Z, bool connected=true) const;

        void contract(const DeexcitationOperator<U,2>& L, DeexcitationOperator<U,2>& Z, bool connected=false) const;
};

}
}

#endif
