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
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;

template <typename U>
PerturbedCCSD<U>::PerturbedCCSD(const string& name, const Config& config)
: Iterative("perturbedccsd", name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.T", "T"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("1epert", "A"));
    reqs.push_back(Requirement("double", "omega"));
    addProduct(Product("double", "convergence", reqs));
    addProduct(Product("ccsd.TA", "TA", reqs));
}

template <typename U>
void PerturbedCCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const OneElectronOperator<U>& A = get<OneElectronOperator<U> >("A");
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    put("TA", new ExcitationOperator<U,2>("T^A", arena, occ, vrt));
    puttmp("D", new ExcitationOperator<U,2>("D", arena, occ, vrt));
    puttmp("X", new ExcitationOperator<U,2>("X", arena, occ, vrt));
    puttmp("Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));

    double omega = get<Scalar>("omega");

    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");
    ExcitationOperator<U,2>& TA = get<ExcitationOperator<U,2> >("TA");
    ExcitationOperator<U,2>& D = gettmp<ExcitationOperator<U,2> >("D");
    ExcitationOperator<U,2>& X = gettmp<ExcitationOperator<U,2> >("X");

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
void PerturbedCCSD<U>::iterate()
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    ExcitationOperator<U,2>& TA = get<ExcitationOperator<U,2> >("TA");
    ExcitationOperator<U,2>& D = gettmp<ExcitationOperator<U,2> >("D");
    ExcitationOperator<U,2>& X = gettmp<ExcitationOperator<U,2> >("X");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");

    Z = X;
    H.contract(TA, Z);

     Z *=  D;
    // Z -= TA;
    TA +=  Z;

    conv = Z.norm(00);

    diis.extrapolate(TA, Z);
}

INSTANTIATE_SPECIALIZATIONS(PerturbedCCSD);
REGISTER_TASK(PerturbedCCSD<double>, "perturbedccsd");
