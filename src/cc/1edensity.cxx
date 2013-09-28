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

#include "1edensity.hpp"

using namespace std;
using namespace aquarius::cc;
using namespace aquarius::op;
using namespace aquarius::tensor;

template <typename U>
OneElectronDensity<U>::OneElectronDensity(const MOSpace<U>& occ, const MOSpace<U>& vrt,
                                          const SymmetryBlockedTensor<U>& Da,
                                          const SymmetryBlockedTensor<U>& Db)
: OneElectronOperator<U>(occ, vrt, Da, Db) {}

/*
 * Form the unrelaxed CCSD density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const ExcitationOperator<U,2>& T)
: OneElectronOperator<U>(T.arena, T.occ, T.vrt)
{
    this->ai["ai"] = T(1)["ai"];
}

/*
 * Form the partial perturbed CCSD Density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const DeexcitationOperator<U,2>& L,
                                          const ExcitationOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: OneElectronOperator<U>(T.arena, T.occ, T.vrt)
{
    OneElectronOperator<U> I(this->arena, this->occ, this->vrt);

    SpinorbitalTensor<U>& IIJ = I.getIJ();
    SpinorbitalTensor<U>& IAB = I.getAB();

    IAB["ab"] = 0.5*T(2)["aemn"]*L(2)["mnbe"];
    IIJ["ij"] = 0.5*T(2)["efim"]*L(2)["jmef"];

    this->ab["ab"] += TA(1)["am"]*L(1)["mb"];
    this->ab["ab"] += 0.5*TA(2)["aemn"]*L(2)["mnbe"];

    this->ij["ij"] -= TA(1)["ei"]*L(1)["je"];
    this->ij["ij"] -= 0.5*TA(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] += TA(1)["ai"];
    this->ai["ai"] += TA(2)["aeim"]*L(1)["me"];
    this->ai["ai"] += this->ij["mi"]*T(1)["am"];
    this->ai["ai"] -= this->ab["ae"]*T(1)["ei"];
    this->ai["ai"] -= IIJ["mi"]*TA(1)["am"];
    this->ai["ai"] -= IAB["ae"]*TA(1)["ei"];
}

/*
 * Form the relaxed CCSD density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const DeexcitationOperator<U,2>& L,
                                          const ExcitationOperator<U,2>& T)
: OneElectronOperator<U>(T.arena, T.occ, T.vrt)
{
    this->ia["ia"] += L(1)["ia"];

    this->ab["ab"] += 0.5*T(2)["aemn"]*L(2)["mnbe"];

    this->ij["ij"] -= T(1)["ei"]*L(1)["je"];
    this->ij["ij"] -= 0.5*T(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] += T(1)["ai"];
    this->ai["ai"] += T(2)["aeim"]*L(1)["me"];
    this->ai["ai"] += this->ij["mi"]*T(1)["am"];
    this->ai["ai"] -= this->ab["ae"]*T(1)["ei"];

    this->ab["ab"] += T(1)["am"]*L(1)["mb"];
}

/*
 * Form the relaxed perturbed CCSD Density
 */
template <typename U>
OneElectronDensity<U>::OneElectronDensity(const DeexcitationOperator<U,2>& L,
                                          const DeexcitationOperator<U,2>& LA,
                                          const ExcitationOperator<U,2>& T,
                                          const ExcitationOperator<U,2>& TA)
: OneElectronOperator<U>(T.arena, T.occ, T.vrt)
{
    OneElectronOperator<U> I(this->arena, this->occ, this->vrt);

    SpinorbitalTensor<U>& IIJ = I.getIJ();
    SpinorbitalTensor<U>& IAB = I.getAB();

    this->ia["ia"] += LA(1)["ia"];

    this->ab["ab"] += 0.5*T(2)["aemn"]*LA(2)["mnbe"];

    this->ij["ij"] -= T(1)["ei"]*LA(1)["je"];
    this->ij["ij"] -= 0.5*T(2)["efim"]*LA(2)["jmef"];

    this->ai["ai"] += T(1)["ai"];
    this->ai["ai"] += T(2)["aeim"]*LA(1)["me"];
    this->ai["ai"] += this->ij["mi"]*T(1)["am"];
    this->ai["ai"] -= this->ab["ae"]*T(1)["ei"];

    this->ab["ab"] += T(1)["am"]*LA(1)["mb"];

    IAB["ab"]  = TA(1)["am"]*L(1)["mb"];
    IAB["ab"] += 0.5*TA(2)["aemn"]*L(2)["mnbe"];

    this->ab["ab"] += IAB["ab"];
    this->ij["ij"] -= IIJ["ij"];

    IIJ["ij"]  = TA(1)["ei"]*L(1)["je"];
    IIJ["ij"] += 0.5*TA(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] += TA(1)["ai"];
    this->ai["ai"] += TA(2)["aeim"]*L(1)["me"];
    this->ai["ai"] -= IIJ["mi"]*T(1)["am"];
    this->ai["ai"] -= IAB["ae"]*T(1)["ei"];

    IAB["ab"] = 0.5*T(2)["aemn"]*L(2)["mnbe"];
    IIJ["ij"] = 0.5*T(2)["efim"]*L(2)["jmef"];

    this->ai["ai"] -= IIJ["mi"]*TA(1)["am"];
    this->ai["ai"] -= IAB["ae"]*TA(1)["ei"];
}

INSTANTIATE_SPECIALIZATIONS(OneElectronDensity);
