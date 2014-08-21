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
: NonIterative("tda", name, config)
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

    put("TDAevals", new CTFTensor<U>("TDAevals", arena, 1, vec(mysize), vec(1), true));
    put("TDAevecs", new CTFTensor<U>("TDAevecs", arena, 1, vec(mysize*mysize), vec(1), true));

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

    vector<U> data1;
    vector<U> data2;

    CTFTensor<U> AAAA = Hguess(vec(0,0),vec(0,0))(vec(0,0,0,0));
    CTFTensor<U> AABB = Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0));
    CTFTensor<U> AAAAtrans = Hguess(vec(0,0),vec(0,0))(vec(0,0,0,0));
    CTFTensor<U> AABBtrans = Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0));
    AAAAtrans = 0;
    AABBtrans = 0;
    AAAAtrans["ajbi"] += AAAA["aibj"];
    AABBtrans["ajbi"] += AABB["aibj"];

    AABBtrans.getAllData(data1);
    AAAAtrans.getAllData(data2);

    for (int ii = 0; ii < data1.size(); ii ++)
    {
        data1[ii] += data2[ii];
    }

    vector<complex<double> > w(mysize);
    vector<U> vl(mysize*mysize);
    vector<U> vr(mysize*mysize);
    geev('N','V',mysize,data1.data(),mysize,w.data(),vl.data(),mysize,vr.data(),mysize);

    vector<int> evalorder(mysize);
    complex<double> lowesteval = w[0]*0.0 + 100000000.0;
    int lowestindex = mysize+1;
    for (int ii = 0; ii<mysize; ii++)
    {
        if (abs(w[ii]) < abs(lowesteval))
        {
            lowestindex = ii;
            lowesteval = w[ii];
        }
    }
    evalorder[0] = lowestindex;

    for (int jj = 0; jj<mysize-1; jj++)
    {
        lowestindex = mysize+1;
        lowesteval = w[0]*0.0 + 100000000.0;

        for (int ii = 0; ii<mysize; ii++)
        {
            if (abs(w[ii]) < abs(lowesteval) and abs(w[ii]) > abs(w[evalorder[jj]]))
            {
                lowestindex = ii;
                lowesteval = w[ii];
            }
        }

        evalorder[jj+1] = lowestindex;
    }

    for (int ii = 0; ii<mysize; ii++)
    {
        if (std::abs(std::imag(w[evalorder[ii]])) > 1e-5)
           throw std::runtime_error("TDA: complex eigenvalue");

        // TDAevals.writeRemoteData(vec(kv_pair(ii,real(w[evalorder[ii]]))));
        if (arena.rank == 0)
            TDAevals.writeRemoteData(vec(kv_pair(ii,real(w[evalorder[ii]]))));
        else
            TDAevals.writeRemoteData();

        for (int jj = 0; jj<mysize; jj++)
        {
            // TDAevecs.writeRemoteData(vec(kv_pair(jj+ii*10,vr[evalorder[ii]*mysize + jj])));
            if (arena.rank == 0)
                TDAevecs.writeRemoteData(vec(kv_pair(jj+ii*10,vr[evalorder[ii]*mysize + jj])));
            else
                TDAevecs.writeRemoteData();
        }
    }
    // vector<U> test1;
    // vector<U> test2;
    // TDAevals.getAllData(test1);
    // TDAevecs.getAllData(test2);
    // cout << test1 << endl;
    // cout << test2 << endl;
    
}

INSTANTIATE_SPECIALIZATIONS(TDA);
REGISTER_TASK(TDA<double>, "tda");
