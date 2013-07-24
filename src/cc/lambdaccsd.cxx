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

#include "lambdaccsd.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::scf;
using namespace aquarius::input;
using namespace aquarius::tensor;

template <typename U>
LambdaCCSD<U>::LambdaCCSD(const Config& config, const STTwoElectronOperator<U,2>& H,
                          const ExponentialOperator<U,2>& T, const double Ecc)
: Iterative(config), DeexcitationOperator<U,2>(H.getSCF()), L(*this),
  Z(this->uhf), D(this->uhf),
  H(H), T(T), Ecc(Ecc), diis(config.get("diis"))
{
    D(0) = (U)1.0;
    D(1)["ia"]  = H.getIJ()["ii"];
    D(1)["ia"] -= H.getAB()["aa"];
    D(2)["ijab"]  = H.getIJ()["ii"];
    D(2)["ijab"] += H.getIJ()["jj"];
    D(2)["ijab"] -= H.getAB()["aa"];
    D(2)["ijab"] -= H.getAB()["bb"];

    D = 1/D;

    L(0) = (U)1.0;
    L(1) = H.getIA()*D(1);
    L(2) = H.getIJAB()*D(2);
}

template <typename U>
void LambdaCCSD<U>::_iterate()
{
    Z = (U)0.0;
    H.contract(L, Z);

    energy = Ecc + real(scalar(Z*conj(L))/scalar(L*conj(L)));

    Z *= D;
    L += Z;

    conv =          Z(1)(0).reduce(CTF_OP_MAXABS);
    conv = max(conv,Z(1)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(0).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(2).reduce(CTF_OP_MAXABS));

    diis.extrapolate(L, Z);
}

INSTANTIATE_SPECIALIZATIONS(LambdaCCSD);
