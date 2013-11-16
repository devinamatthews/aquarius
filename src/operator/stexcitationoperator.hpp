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

#ifndef _AQUARIUS_OPERATOR_STEXCITATIONOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_STEXCITATIONOPERATOR_HPP_

#include "1eoperator.hpp"
#include "2eoperator.hpp"
#include "excitationoperator.hpp"

namespace aquarius
{
namespace op
{

/*
 *                                  _    -T   T       T
 * Form the pure excitation part of X = e  X e  = (X e )
 *                                                      c
 */
template <typename U, int nex> class STExcitationOperator;

template <typename U>
class STExcitationOperator<U,2> : public ExcitationOperator<U,2>
{
    public:
        STExcitationOperator(const std::string& name, const OneElectronOperator<U>& X, const ExcitationOperator<U,2>& T);

        STExcitationOperator(const std::string& name, const TwoElectronOperator<U>& X, const ExcitationOperator<U,2>& T);

        static void transform(const OneElectronOperator<U>& X, const ExcitationOperator<U,2>& T,
                              ExcitationOperator<U,2>& Z);

        static void transform(const OneElectronOperator<U>& X, const ExcitationOperator<U,2>& T,
                              ExcitationOperator<U,2>& Z, OneElectronOperator<U>& W);

        static void transform(const TwoElectronOperator<U>& X, const ExcitationOperator<U,2>& T,
                              ExcitationOperator<U,2>& Z);

        static void transform(const TwoElectronOperator<U>& X, const ExcitationOperator<U,2>& T,
                              tensor::SpinorbitalTensor<U>& Tau, ExcitationOperator<U,2>& Z, TwoElectronOperator<U>& W);
};

}
}

#endif
