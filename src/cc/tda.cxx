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

#include "tda.hpp"
#include "util/lapack.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;

template <typename U>
TDA<U>::TDA(const std::string& name, const Config& config)
: Task("tda", name)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    addProduct(Product("tda.TDAevals", "TDAevals", reqs));
    addProduct(Product("tda.TDAevecs", "TDAevecs", reqs));
}

template <typename U>
void TDA<U>::run(TaskDAG& dag, const Arena& arena)
{
    TwoElectronOperator<U>& W = get<TwoElectronOperator<U> >("H");
    const Space& occ = W.occ;
    const Space& vrt = W.vrt;
    int mysize = occ.nalpha[0] * vrt.nalpha[0];

    put("TDAevals", new CTFTensor<U>("TDAevals", arena, 1, vec(mysize), vec(NS), true));
    put("TDAevecs", new CTFTensor<U>("TDAevecs", arena, 2, vec(mysize,mysize), vec(NS,NS), true));

    CTFTensor<U>& TDAevals = get<CTFTensor<U> >("TDAevals");
    CTFTensor<U>& TDAevecs = get<CTFTensor<U> >("TDAevecs");

    SpinorbitalTensor<U> Hguess("Hguess", W.getAIBJ());
    Hguess = 0;

    SpinorbitalTensor<U>& FAB = W.getAB();
    SpinorbitalTensor<U>& FIJ = W.getIJ();
    SpinorbitalTensor<U>& WAIBJ = W.getAIBJ();

    Hguess["aiaj"] -= FIJ["ij"];
    Hguess["aibi"] += FAB["ab"];
    Hguess["aibj"] -= WAIBJ["aibj"];

    vector<U> data;

    CTFTensor<U>& AAAA = Hguess(vec(0,0),vec(0,0))(vec(0,0,0,0));
    CTFTensor<U>& AABB = Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0));
    CTFTensor<U>& BBAA = Hguess(vec(0,1),vec(1,0))(vec(0,0,0,0));
    CTFTensor<U>& BBBB = Hguess(vec(1,1),vec(1,1))(vec(0,0,0,0));
    CTFTensor<U> H(AAAA);
    H["ajbi"]  = 0.5*AAAA["aibj"];
    H["ajbi"] += 0.5*AABB["aibj"];
    H["ajbi"] += 0.5*BBAA["aibj"];
    H["ajbi"] += 0.5*BBBB["aibj"];

    H.getAllData(data);
    printf("H: %18.15f\n", H.norm(2));

    vector<U> w(mysize);
    heev('V','U',mysize,data.data(),mysize,w.data());

    vector<int> evalorder = range(mysize);

    cosort(w.begin(), w.end(),
           evalorder.begin(), evalorder.end());
    printf("%18.15f\n", w[0]);
    printf("%18.15f\n", w[1]);

    vector<tkv_pair<U> > pairs;

    for (int i = 0;i < mysize;i++)
    {
        pairs.push_back(tkv_pair<U>(i, w[evalorder[i]]));
    }

    if (arena.rank == 0)
        TDAevals.writeRemoteData(pairs);
    else
        TDAevals.writeRemoteData();

    pairs.clear();

    for (int i = 0;i < mysize;i++)
    {
        for (int j = 0;j < mysize;j++)
        {
            pairs.push_back(kv_pair(j+i*mysize,data[evalorder[i]*mysize + j]));
            // cout << j+i*mysize << " " << data[evalorder[i]*mysize + j] << endl;
        }
    }

    if (arena.rank == 0)
        TDAevecs.writeRemoteData(pairs);
    else
        TDAevecs.writeRemoteData();

    // vector<U> data1;
    // TDAevecs.getAllData(data1);
    // cout << "data1.size() = " << data1.size() << endl;

}

INSTANTIATE_SPECIALIZATIONS(TDA);
REGISTER_TASK(TDA<double>, "tda");
