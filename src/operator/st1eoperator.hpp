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

#ifndef _AQUARIUS_OPERATOR_ST1EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_ST1EOPERATOR_HPP_

#include "1eoperator.hpp"
#include "excitationoperator.hpp"

namespace aquarius
{
namespace op
{

/*
 *      _    -T   T       T
 * Form X = e  X e  = (X e ) , up to one-electron terms
 *                          c
 */
template <typename U>
class STOneElectronOperator<U> : public OneElectronOperator<U>
{
    public:
        template <int N>
        STOneElectronOperator(const std::string& name, const OneElectronOperator<U>& X, const ExcitationOperator<U,N>& T);

        template <int N>
        STOneElectronOperator(const std::string& name, const TwoElectronOperator<U>& X, const ExcitationOperator<U,N>& T);
};

template <typename U> template <int N>
STOneElectronOperator<U>::STOneElectronOperator(const std::string& name, const OneElectronOperator<U>& X,
                                                const ExcitationOperator<U,N>& T)
: OneElectronOperator<U>(name, X)
{
    assert(N >= 2 && N <= 4);

    this->ij["mi"] += this->ia["me"]*T(1)["ei"];

    this->ai["ai"] += T(2)["aeim"]*this->ia["me"];
    this->ai["ai"] += T(1)["ei"]*this->ab["ae"];
    this->ai["ai"] -= T(1)["am"]*this->ij["mi"];

    this->ab["ae"] -= this->ia["me"]*T(1)["am"];
}

template <typename U> template <int N>
STOneElectronOperator<U>::STOneElectronOperator(const std::string& name, const TwoElectronOperator<U>& X,
                                                const ExcitationOperator<U,N>& T)
: OneElectronOperator<U>(name, X)
{
    assert(N >= 2 && N <= 4);

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    this->ia["me"] = this->ijab["mnef"]*T(1)["fn"];

    this->ij["mi"] += 0.5*this->ijab["nmef"]*T(2)["efni"];
    this->ij["mi"] += this->ia["me"]*T(1)["ei"];
    this->ij["mi"] += this->ijak["nmfi"]*T(1)["fn"];

    this->ai["ai"] -= T(1)["em"]*this->aibj["amei"];
    this->ai["ai"] += 0.5*this->aibc["amef"]*Tau["efim"];
    this->ai["ai"] -= 0.5*this->ijak["mnei"]*T(2)["eamn"];
    this->ai["ai"] += T(2)["aeim"]*this->ia["me"];
    this->ai["ai"] += T(1)["ei"]*this->ab["ae"];
    this->ai["ai"] -= T(1)["am"]*this->ij["mi"];

    this->ab["ae"] -= 0.5*this->ijab["mnfe"]*T(2)["famn"];
    this->ab["ae"] -= this->ia["me"]*T(1)["am"];
    this->ab["ae"] += this->aibc["anef"]*T(1)["fn"];

    if (N > 2)
    {
        this->ai["ai"] += 0.25*this->ijab["mnef"]*T(3)["aefimn"];
    }
}

}
}

#endif
