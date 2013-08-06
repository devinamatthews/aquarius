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

#include "perturbedlambdaccsd.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::scf;
using namespace aquarius::input;
using namespace aquarius::tensor;

template <typename U>
PerturbedLambdaCCSD<U>::PerturbedLambdaCCSD(const Config& config, const STTwoElectronOperator<U,2>& H,
                                            const DeexcitationOperator<U,2>& L,
                                            const PerturbedSTTwoElectronOperator<U,2>& A, const U omega)
: Iterative(config), DeexcitationOperator<U,2>(H.getSCF()), LA(*this),
  Z(this->uhf), D(this->uhf), N(this->uhf), H(H),
  omega(omega), diis(config.get("diis"))
{
    D(0) = (U)1.0;
    D(1)["ia"]  = H.getIJ()["ii"];
    D(1)["ia"] -= H.getAB()["aa"];
    D(2)["ijab"]  = H.getIJ()["ii"];
    D(2)["ijab"] += H.getIJ()["jj"];
    D(2)["ijab"] -= H.getAB()["aa"];
    D(2)["ijab"] -= H.getAB()["bb"];

    D -= omega;
    D = 1/D;

    N(1) = (U)0.0;
    N(2) = (U)0.0;
    A.contract(L, N);
    N(0) = (U)0.0;

    LA = N*D;
}

template <typename U>
void PerturbedLambdaCCSD<U>::_iterate()
{
    Z = N;
    H.contract(LA, Z);

     Z *= D;
    LA += Z;

    conv =          Z(1)(0).norm(00);
    conv = max(conv,Z(1)(1).norm(00));
    conv = max(conv,Z(2)(0).norm(00));
    conv = max(conv,Z(2)(1).norm(00));
    conv = max(conv,Z(2)(2).norm(00));

    diis.extrapolate(LA, Z);
}

INSTANTIATE_SPECIALIZATIONS(PerturbedLambdaCCSD);
