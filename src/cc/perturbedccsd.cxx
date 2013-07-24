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

#include "perturbedccsd.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::scf;
using namespace aquarius::input;
using namespace aquarius::tensor;

template <typename U>
PerturbedCCSD<U>::PerturbedCCSD(const Config& config, const STTwoElectronOperator<U,2>& H,
                                const ExponentialOperator<U,2>& T, const OneElectronOperator<U>& A, const U omega)
: Iterative(config), ExcitationOperator<U,2>(T.getSCF()),
  H(H), omega(omega), TA(*this), D(this->uhf), Z(this->uhf), X(this->uhf),
  diis(config.get("diis"))
{
    D(0) = (U)1.0;
    D(1)["ai"]  = H.getIJ()["ii"];
    D(1)["ai"] -= H.getAB()["aa"];
    D(2)["abij"]  = H.getIJ()["ii"];
    D(2)["abij"] += H.getIJ()["jj"];
    D(2)["abij"] -= H.getAB()["aa"];
    D(2)["abij"] -= H.getAB()["bb"];

    D += omega;
    D = 1/D;

    STExcitationOperator<U,2>::transform(A, T, X);
    X(0) = (U)0.0;

    TA = X*D;
}

template <typename U>
PerturbedCCSD<U>::PerturbedCCSD(const Config& config, const STTwoElectronOperator<U,2>& H,
                                const ExponentialOperator<U,2>& T, const TwoElectronOperator<U>& A, const U omega)
: Iterative(config), ExcitationOperator<U,2>(T.getSCF()),
  H(H), omega(omega), TA(*this), D(this->uhf), Z(this->uhf), X(this->uhf),
  diis(config.get("diis"))
{
    D(0) = (U)1.0;
    D(1)["ai"]  = H.getIJ()["ii"];
    D(1)["ai"] -= H.getAB()["aa"];
    D(2)["abij"]  = H.getIJ()["ii"];
    D(2)["abij"] += H.getIJ()["jj"];
    D(2)["abij"] -= H.getAB()["aa"];
    D(2)["abij"] -= H.getAB()["bb"];

    D += omega;
    D = 1/D;

    STExcitationOperator<U,2>::transform(A, T, X);
    X(0) = (U)0.0;

    TA = X*D;
}

template <typename U>
void PerturbedCCSD<U>::_iterate()
{
    Z = X;
    H.contract(TA, Z);

     Z *=  D;
    // Z -= TA;
    TA +=  Z;

    conv =          Z(1)(0).reduce(CTF_OP_MAXABS);
    conv = max(conv,Z(1)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(0).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(2).reduce(CTF_OP_MAXABS));

    diis.extrapolate(TA, Z);
}

INSTANTIATE_SPECIALIZATIONS(PerturbedCCSD);
