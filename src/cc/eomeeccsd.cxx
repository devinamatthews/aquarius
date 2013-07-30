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

#include "eomeeccsd.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::scf;
using namespace aquarius::input;
using namespace aquarius::tensor;

template <typename U>
EOMEECCSD<U>::EOMEECCSD(const Config& config, const STTwoElectronOperator<U,2>& H,
                        const ExponentialOperator<U,2>& T)
: Iterative(config), ExcitationOperator<U,2>(H.getSCF()), R(*this),
  Z(this->uhf), D(this->uhf), H(H), T(T), davidson(config.get("davidson"))
{
    D(0) = (U)1.0;
    D(1)["ai"]  = H.getIJ()["ii"];
    D(1)["ai"] -= H.getAB()["aa"];
    //D(1)["ai"] -= H.getAIBJ()["aiai"];
    D(2)["abij"]  = H.getIJ()["ii"];
    D(2)["abij"] += H.getIJ()["jj"];
    D(2)["abij"] -= H.getAB()["aa"];
    D(2)["abij"] -= H.getAB()["bb"];
    //D(2)["abij"] -= H.getAIBJ()["aiai"];
    //D(2)["abij"] -= H.getAIBJ()["bjbj"];
    //D(2)["abij"] += H.getABCD()["abab"];
    //D(2)["abij"] += H.getIJKL()["ijij"];
    //D(2)["abij"] -= H.getIJAB()["imab"]*T(2)["abim"];
    //D(2)["abij"] += H.getIJAB()["ijae"]*T(2)["aeij"];

    //TODO: guess
}

template <typename U>
void EOMEECCSD<U>::_iterate()
{
    H.contract(R, Z);

    energy = davidson.extrapolate(R, Z, D);

    conv =          Z(1)(0).reduce(CTF_OP_MAXABS);
    conv = max(conv,Z(1)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(0).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(1).reduce(CTF_OP_MAXABS));
    conv = max(conv,Z(2)(2).reduce(CTF_OP_MAXABS));
}

INSTANTIATE_SPECIALIZATIONS(EOMEECCSD);
