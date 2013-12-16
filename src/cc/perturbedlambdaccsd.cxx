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
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;

template <typename U>
PerturbedLambdaCCSD<U>::PerturbedLambdaCCSD(const std::string& name, const Config& config)
: Iterative("perturbedlambdaccsd", name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.L", "L"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("1epert", "A"));
    reqs.push_back(Requirement("double", "omega"));
    addProduct(Product("double", "convergence", reqs));
    addProduct(Product("ccsd.LA", "LA", reqs));
}

template <typename U>
void PerturbedLambdaCCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const PerturbedSTTwoElectronOperator<U,2>& A = get<PerturbedSTTwoElectronOperator<U,2> >("A");
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    put("LA", new DeexcitationOperator<U,2>("L^A", arena, occ, vrt));
    puttmp("D", new DeexcitationOperator<U,2>("D", arena, occ, vrt));
    puttmp("N", new DeexcitationOperator<U,2>("N", arena, occ, vrt));
    puttmp("Z", new DeexcitationOperator<U,2>("Z", arena, occ, vrt));

    double omega = get<Scalar>("omega");

    DeexcitationOperator<U,2>& L = get<DeexcitationOperator<U,2> >("L");
    DeexcitationOperator<U,2>& LA = get<DeexcitationOperator<U,2> >("LA");
    DeexcitationOperator<U,2>& D = gettmp<DeexcitationOperator<U,2> >("D");
    DeexcitationOperator<U,2>& N = gettmp<DeexcitationOperator<U,2> >("N");

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

    Iterative::run(dag, arena);

    put("convergence", new Scalar(arena, conv));
}

template <typename U>
void PerturbedLambdaCCSD<U>::iterate(const Arena& arena)
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    DeexcitationOperator<U,2>& LA = get<DeexcitationOperator<U,2> >("LA");
    DeexcitationOperator<U,2>& D = gettmp<DeexcitationOperator<U,2> >("D");
    DeexcitationOperator<U,2>& N = gettmp<DeexcitationOperator<U,2> >("N");
    DeexcitationOperator<U,2>& Z = gettmp<DeexcitationOperator<U,2> >("Z");

    Z = N;
    H.contract(LA, Z);

     Z *= D;
    LA += Z;

    conv = Z.norm(00);

    diis.extrapolate(LA, Z);
}

INSTANTIATE_SPECIALIZATIONS(PerturbedLambdaCCSD);
REGISTER_TASK(PerturbedLambdaCCSD<double>, "perturbedlambdaccsd");
