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

    put("R", new ExcitationOperator<U,2>("R", arena, occ, vrt));
    puttmp("D", new ExcitationOperator<U,2>("D", arena, occ, vrt));
    puttmp("Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));

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

    //SpinorbitalTensor<U>& R1 = R(1);
    //SymmetryBlockedTensor<U>& R11 = R1(vec(1,0),vec(0,1));
    //CTFTensor<U>& R11CTF = R11(vec(0,0));
    //R11CTF.writeRemoteData(vec(tkv_pair<U>(76,1)));
    //Iterative::run(dag, arena);
    //put("energy", new Scalar(arena, energy));
    //put("convergence", new Scalar(arena, conv));

    vector<vector<U> > cols;

    for (int nex = 1;nex <= 2;nex++)
    {  
        SpinorbitalTensor<U>& Rso = R(nex);
        SymmetryBlockedTensor<U>& Rsb = Rso(vec(1,0),vec(0,1));
        CTFTensor<U>& Rctf = Rsb(vector<int>(nex*2,0));
        int max_key;
        if (nex == 1)
        {
            max_key = 10;
        }
        else
        {
            max_key = 100;
        }
        for (int key = 0;key < max_key;key++)
        {
            cols.push_back(vector<U>());
            if (arena.rank == 0)
            {
                //cout << "key = " << key << endl;
                Rctf.writeRemoteData(vec(tkv_pair<U>(key,1.0)));
            }
            else
            {
                Rctf.writeRemoteData();
            }
            // do Hbar*R product into Z
            ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");
            H.contractsam(R, Z);
            for (int nex2 = 1;nex2 <= 2;nex2++)
            {
                SpinorbitalTensor<U>& Zso = Z(nex2);
                SymmetryBlockedTensor<U>& Zsb = Zso(vec(1,0),vec(0,1));
                CTFTensor<U>& Zctf = Zsb(vector<int>(nex2*2,0));
                vector<U> data;
                Zctf.getAllData(data);
                cols.back().insert(cols.back().end(), data.begin(), data.end());
            }
        }
    }
    cout << "cols.size() = " << cols.size() << endl;
    cout << "cols[0].size() = " << cols[0].size() << endl;

    vector<U> matrix;
    for (int i=0;i<cols.size();i++)
    {
        matrix.insert(matrix.end(),cols[i].begin(),cols[i].end());
    }

    //cout << "matrix.size() = " << matrix.size() << endl;

    int mysize = cols.size();
    vector<complex<double> > w(mysize);
    vector<U> vl(mysize);
    vector<U> vr(mysize);
    geev('N','N',mysize,matrix.data(),mysize,w.data(),vl.data(),mysize,vr.data(),mysize);
    
    cout << "w.size() = " << w.size() << endl;
    cout << w << endl;
}

template <typename U>
void EOMEECCSD<U>::iterate(const Arena& arena)
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
