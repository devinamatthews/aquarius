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
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;

template <typename U>
LambdaCCSD<U>::LambdaCCSD(const string& name, const Config& config)
: Iterative("lambdaccsd", name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("ccsd.T", "T"));
    addProduct(Product("double", "energy", reqs));
    addProduct(Product("double", "convergence", reqs));
    addProduct(Product("ccsd.L", "L", reqs));
}

template <typename U>
void LambdaCCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    put("L", new DeexcitationOperator<U,2>(arena, occ, vrt));
    puttmp("D", new Denominator<U>(H));
    puttmp("Z", new DeexcitationOperator<U,2>(arena, occ, vrt));

    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");
    DeexcitationOperator<U,2>& L = get<DeexcitationOperator<U,2> >("L");
    Denominator<U>& D = gettmp<Denominator<U> >("D");

    L(0) = (U)1.0;
    L(1) = H.getIA();
    L(2) = H.getIJAB();

    L.weight(D);

    {
        SpinorbitalTensor<U> mTau(T(2));
        mTau["abij"] -= 0.5*T(1)["ai"]*T(1)["bj"];
        Ecc = scalar(T(1)["ai"]*H.getIA()["ia"]) + 0.25*scalar(mTau["abij"]*H.getIJAB()["ijab"]);
    }

    Iterative::run(dag, arena);

    put("energy", new Scalar(arena, energy));
    put("convergence", new Scalar(arena, conv));
}

template <typename U>
void LambdaCCSD<U>::iterate()
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    DeexcitationOperator<U,2>& L = get<DeexcitationOperator<U,2> >("L");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    DeexcitationOperator<U,2>& Z = gettmp<DeexcitationOperator<U,2> >("Z");

    Z = (U)0.0;
    H.contract(L, Z);

    energy = Ecc + real(scalar(Z*conj(L))/scalar(L*conj(L)));

    Z.weight(D);
    L += Z;

    conv =          Z(1)(0).norm(00);
    conv = max(conv,Z(1)(1).norm(00));
    conv = max(conv,Z(2)(0).norm(00));
    conv = max(conv,Z(2)(1).norm(00));
    conv = max(conv,Z(2)(2).norm(00));

    diis.extrapolate(L, Z);
}

INSTANTIATE_SPECIALIZATIONS(LambdaCCSD);
REGISTER_TASK(LambdaCCSD<double>,"lambdaccsd");
