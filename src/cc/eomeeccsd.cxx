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
#include "util/lapack.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;

template <typename U>
EOMEECCSD<U>::EOMEECCSD(const std::string& name, const Config& config)
: Iterative("eomeeccsd", name, config), nroot(config.get<int>("nroot")), davidson(config.get("davidson"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.T", "T"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("tda.TDAevals", "TDAevals"));
    reqs.push_back(Requirement("tda.TDAevecs", "TDAevecs"));
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

    int mysize = occ.nalpha[0] * vrt.nalpha[0];

    put("R", new ExcitationOperator<U,2>("R", arena, occ, vrt));
    puttmp("D", new Denominator<U>(H));
    puttmp("Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));

    ExcitationOperator<U,2>& R = get<ExcitationOperator<U,2> >("R");
    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");

    CTFTensor<U>& TDAevecs = get<CTFTensor<U> >("TDAevecs");
    vector<U> TDAevecsVec;
    TDAevecs.getAllData(TDAevecsVec);

    vector<U> zerovec(mysize*mysize);
    // cout << zerovec << endl;
    while (TDAevecsVec == zerovec)
    {
        TDAevecs.getAllData(TDAevecsVec);
        cout << TDAevecsVec << endl;
    }
    // cout << TDAevecsVec << endl;

    

    
    R(0) = 0; R(1) = 0; R(2) = 0;
    double thisval;

    for (int ii=0; ii<mysize; ii++)
    {
        thisval = TDAevecsVec[nroot*mysize+ii];
        if (R(1).arena.rank == 0)
            R(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval/sqrt(2))));
        else
            R(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData();

        if (R(1).arena.rank == 0)
            R(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval/sqrt(2))));
        else
            R(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData();
    }

    Iterative::run(dag, arena);
    put("energy", new Scalar(arena, energy));
    put("convergence", new Scalar(arena, conv));
}

template <typename U>
void EOMEECCSD<U>::iterate(const Arena& arena)
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    int mysize = occ.nalpha[0] * vrt.nalpha[0];

    ExcitationOperator<U,2>& R = get<ExcitationOperator<U,2> >("R");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");

    CTFTensor<U>& TDAevals = get<CTFTensor<U> >("TDAevals");
    vector<U> TDAevalsVec;
    TDAevals.getAllData(TDAevalsVec);

    vector<U> zerovec(mysize);
    // cout << zerovec << endl;
    while (TDAevalsVec == zerovec)
    {
        TDAevals.getAllData(TDAevalsVec);
        // cout << TDAevalsVec << endl;
    }
    // cout << TDAevalsVec << endl;

    Z = 0;
    H.contract(R, Z);

    energy = davidson.extrapolate(R, Z, D, TDAevalsVec[nroot]);

    conv = Z.norm(00);
}

INSTANTIATE_SPECIALIZATIONS(EOMEECCSD);
REGISTER_TASK(EOMEECCSD<double>, "eomeeccsd");
