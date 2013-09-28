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
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;

template <typename U>
EOMEECCSD<U>::EOMEECCSD(const std::string& name, const Config& config)
: Iterative("eomeeccsd", name, config), davidson(config.get("davidson"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.T", "T"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    addProduct(Product("double", "energy", reqs));
    addProduct(Product("double", "convergence", reqs));
    addProduct(Product("eomeeccsd.R", "R", reqs));
}

template <typename U>
void EOMEECCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    put("R", new ExcitationOperator<U,2>(arena, occ, vrt));
    puttmp("D", new ExcitationOperator<U,2>(arena, occ, vrt));
    puttmp("Z", new ExcitationOperator<U,2>(arena, occ, vrt));

    ExcitationOperator<U,2>& R = get<ExcitationOperator<U,2> >("R");
    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");
    ExcitationOperator<U,2>& D = gettmp<ExcitationOperator<U,2> >("D");

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
void EOMEECCSD<U>::iterate()
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    ExcitationOperator<U,2>& R = get<ExcitationOperator<U,2> >("R");
    ExcitationOperator<U,2>& D = gettmp<ExcitationOperator<U,2> >("D");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");

    H.contract(R, Z);

    energy = davidson.extrapolate(R, Z, D);

    conv = Z.norm(00);
}

INSTANTIATE_SPECIALIZATIONS(EOMEECCSD);
REGISTER_TASK(EOMEECCSD<double>, "eomeeccsd");
