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
    reqs.push_back(Requirement("moints", "H"));
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
    puttmp("D", new Denominator<U>(H));
    puttmp("Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));

    ExcitationOperator<U,2>& R = get<ExcitationOperator<U,2> >("R");
    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");

    //TODO: guess

    // R(0) = 0;R(1) = 0;R(2) = 0;

    // cout << occ.nalpha[0] << ' ' << vrt.nalpha[0] << ' ' << (occ.nalpha[0]-1)*vrt.nalpha[0] << endl;

    // if (R(1).arena.rank == 0)
    //     R(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData(vec(kv_pair((occ.nalpha[0]-1)*vrt.nalpha[0],1/sqrt(2))));
    // else
    //     R(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData();

    // if (R(1).arena.rank == 0)
    //     R(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData(vec(kv_pair((occ.nalpha[0]-1)*vrt.nalpha[0],1/sqrt(2))));
    // else
    //     R(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData();

    // Iterative::run(dag, arena);
    // put("energy", new Scalar(arena, energy));
    // put("convergence", new Scalar(arena, conv));

    
    TwoElectronOperator<U>& W = get<TwoElectronOperator<U> >("H");
    // TwoElectronOperator<U>& W = get<TwoElectronOperator<U> >("Hbar");
    // STTwoElectronOperator<U,2>& W = get<STTwoElectronOperator<U,2> >("Hbar");
    SpinorbitalTensor<U> Hguess("Hguess", W.getAIBJ()); // TDA
    // SpinorbitalTensor<U> Hguess("Hguess", W.getABIJ()); // Diags of Hbar?
    Hguess = 0;

    // SpinorbitalTensor<U> Hnew("Hnew", W.getAIBJ());
    // Hnew = 0;

    // SymmetryBlockedTensor<U> Hguess_sb1 = Hguess(vec(1,0),vec(0,1));
    // CTFTensor<U> Hguess_ctf1 = Hguess_sb1(vec(0,0,0,0));
    // vector<U> data1;
    // Hguess_ctf1.getAllData(data1);
    // cout << data1 << endl;

    // cout << "occ: " << occ.nalpha[0] << endl;
    // cout << "vrt: " << vrt.nalpha[0] << endl;


    SpinorbitalTensor<U>& FAB = W.getAB();
    SpinorbitalTensor<U>& FIJ = W.getIJ();
    SpinorbitalTensor<U>& WAIBJ = W.getAIBJ(); // TDA
    // Hguess["abij"] += FAB["aa"]; // Diags of Hbar?
    // Hguess["abij"] += FAB["bb"]; // Diags of Hbar?
    // Hguess["abij"] -= FIJ["ii"]; // Diags of Hbar?
    // Hguess["abij"] -= FIJ["jj"]; // Diags of Hbar?

    Hguess["aiaj"] -= FIJ["ij"];
    Hguess["aibi"] += FAB["ab"];
    Hguess["aibj"] -= WAIBJ["aibj"];  

    

    // Hguess["aiaj"] -= 0.5*FIJ["ij"];
    // Hguess["aiaj"] -= 0.5*FIJ["ji"]; 
    // Hguess["aibi"] += 0.5*FAB["ab"]; 
    // Hguess["aibi"] += 0.5*FAB["ba"];
    // Hguess["aibj"] -= 0.25*WAIBJ["aibj"]; 
    // Hguess["aibj"] -= 0.25*WAIBJ["biaj"]; 
    // Hguess["aibj"] -= 0.25*WAIBJ["ajbi"]; 
    // Hguess["aibj"] -= 0.25*WAIBJ["bjai"];

    // Hnew["ajbi"] += Hguess["aibj"]; 


    // Hguess(vec(0,0),vec(0,0))(vec(0,0)) = Hguess(vec(1,0),vec(0,1))(vec(0,0));
    // Hguess(vec(0,0),vec(0,0))(vec(0,0,0,0)) = 0.5*Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0));
    // Hguess(vec(2,0),vec(0,2))(vec(0,0,0,0)) = 0.5*Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0));

    // SpinorbitalTensor<U> SIJ = W.getIJ(); 
    // SIJ=0; 
    // cout << "2norm of SIJ: " << SIJ.norm(2) << endl;
    // Hguess["aiaj"] += FIJ["ij"]; 
    // SIJ["ij"] += Hguess["aibj"];
    // SIJ["ij"] -= 2*FIJ["ij"]; // where n is the dimension of a, b
    // cout << "2norm of FIJ: " << FIJ.norm(2) << endl;
    // cout << "2norm of Hguess: " << Hguess.norm(2) << endl;
    // cout << "2norm of SIJ: " << SIJ.norm(2) << endl;

    // SymmetryBlockedTensor<U> Hguess_sb = Hguess(vec(1,0),vec(1,0));
    // CTFTensor<U> Hguess_ctf = Hguess_sb(vec(0,0,0,0));
    // vector<U> data;
    // Hguess_ctf.getAllData(data);
    vector<U> data1;
    vector<U> data2;
    vector<U> data3;
    vector<U> data4;
    vector<U> data5;
    vector<U> data6;

    CTFTensor<U> AAAA = Hguess(vec(0,0),vec(0,0))(vec(0,0,0,0));
    CTFTensor<U> AABB = Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0));
    CTFTensor<U> AAAAtrans = Hguess(vec(0,0),vec(0,0))(vec(0,0,0,0));
    CTFTensor<U> AABBtrans = Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0));
    AAAAtrans = 0;
    AABBtrans = 0;
    AAAAtrans["ajbi"] += AAAA["aibj"];
    AABBtrans["ajbi"] += AABB["aibj"];

    Hguess(vec(1,0),vec(1,0))(vec(0,0,0,0)).getAllData(data1); // TDA
    // Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0)).getAllData(data2); // TDA
    AABBtrans.getAllData(data2);
    Hguess(vec(0,1),vec(1,0))(vec(0,0,0,0)).getAllData(data3); // TDA
    Hguess(vec(0,1),vec(0,1))(vec(0,0,0,0)).getAllData(data4); // TDA
    // Hguess(vec(0,0),vec(0,0))(vec(0,0,0,0)).getAllData(data5); // TDA
    AAAAtrans.getAllData(data5);
    Hguess(vec(1,1),vec(1,1))(vec(0,0,0,0)).getAllData(data6); // TDA
    // Hguess(vec(0,1),vec(1,0))(vec(0,0,0,0)).getAllData(data3); // TDA
    // Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0)).getAllData(data1); // Diags of Hbar
    // Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0)).getAllData(data2);
    // Hguess(vec(0,0),vec(0,0))(vec(0,0,0,0)).getAllData(data3);
    // Hguess(vec(0,0),vec(0,0))(vec(0,0)).getAllData(data4);

    // cout << data1.size() << endl;
    // cout << data1 << endl;
    // cout << data2.size() << endl;
    // cout << data2 << endl;
    // cout << data3.size() << endl;
    // cout << data3 << endl;
    // cout << data4.size() << endl;
    // cout << data4 << endl;
    // cout << data5.size() << endl;
    // cout << data5 << endl;
    // cout << data6.size() << endl;
    // cout << data6 << endl;

    for (int ii = 0; ii < data1.size(); ii ++)
    {
        data5[ii] += data2[ii];
    }

    cout << "data5" << endl;
    cout << data5.size() << endl;
    cout << data5 << endl;

    int mysize = occ.nalpha[0] * vrt.nalpha[0];
    cout << "mysize = " << mysize << endl;
    vector<complex<double> > w(mysize);
    vector<U> vl(mysize*mysize);
    vector<U> vr(mysize*mysize);
    geev('N','V',mysize,data5.data(),mysize,w.data(),vl.data(),mysize,vr.data(),mysize);

    // for (int ii = 0; ii<10; ii++)
    // {
    //     cout << w[ii] << endl;
    //     cout << vr[ii*10] << " " << vr[ii*10 + 1] << " " << vr[ii*10 + 2] << " " << vr[ii*10 + 3] << " " << vr[ii*10 + 4] << " " << vr[ii*10 + 5] << " " << vr[ii*10 + 6] << " " << vr[ii*10 + 7] << " " << vr[ii*10 + 8] << " " << vr[ii*10 + 9] << " " << endl;
    // }
    // cout << "vr = " << vr << endl;
    // cout << setprecision(9) << w << endl;

    vector<int> evalorder(mysize);
    complex<double> lowesteval = w[0]*0.0 + 100000000.0;
    int lowestindex = mysize+1;
    for (int ii = 0; ii<mysize; ii++)
    {
        // cout << abs(w[ii]) << endl;
        // cout << abs(lowesteval) << endl;
        if (abs(w[ii]) < abs(lowesteval))
        {
            lowestindex = ii;
            lowesteval = w[ii];
        }
    }
    evalorder[0] = lowestindex;

    // cout << "li = " << lowestindex << endl;
    // cout << "le = " << lowesteval << endl;

    for (int jj = 0; jj<mysize-1; jj++)
    {
        lowestindex = mysize+1;
        lowesteval = w[0]*0.0 + 100000000.0;

        for (int ii = 0; ii<mysize; ii++)
        {
            // cout << abs(w[ii]) << endl;
            // cout << abs(lowesteval) << endl;
            if (abs(w[ii]) < abs(lowesteval) and abs(w[ii]) > abs(w[evalorder[jj]]))
            {
                lowestindex = ii;
                lowesteval = w[ii];
            }
        }

        evalorder[jj+1] = lowestindex;
        // cout << "li = " << lowestindex << endl;
        // cout << "le = " << lowesteval << endl;
    }

    // int myroot = nroot;
    // cout << myroot << endl;
    vector<double> myval(mysize);
    double totval;
    double thisval;
    for (int ii=0; ii<mysize; ii++)
    {
        myval[ii] = vr[evalorder[nroot]*mysize+ii];
        totval += myval[ii];
    }
    // cout << "totval = " << totval << endl;
    // cout << "eval = " << w[evalorder[nroot]] << endl;
    // cout << "eval = " << setprecision(9) << w[evalorder[0]] << endl;
    // cout << "eval = " << setprecision(9) << w[evalorder[1]] << endl;
    // cout << "eval = " << setprecision(9) << w[evalorder[2]] << endl;
    // cout << "eval = " << setprecision(9) << w[evalorder[3]] << endl;
    // cout << "eval = " << setprecision(9) << w[evalorder[4]] << endl;
    // cout << "eval = " << setprecision(9) << w[evalorder[5]] << endl;
    for (int ii = 0; ii<10; ii++)
    {
        cout << "eval" << ii+1 << " = " << setprecision(9) << w[evalorder[ii]] << endl;
        cout << "evec" << ii+1 << " = " << setprecision(3) << vr[evalorder[ii]*10] << " " << vr[evalorder[ii]*10 + 1] << " " << vr[evalorder[ii]*10 + 2] << " " << vr[evalorder[ii]*10 + 3] << " " << vr[evalorder[ii]*10 + 4] << " " << vr[evalorder[ii]*10 + 5] << " " << vr[evalorder[ii]*10 + 6] << " " << vr[evalorder[ii]*10 + 7] << " " << vr[evalorder[ii]*10 + 8] << " " << vr[evalorder[ii]*10 + 9] << " " << endl;
        cout << "" << endl;
    }

    R(0) = 0; R(1) = 0; R(2) = 0;

    for (int ii=0; ii<mysize; ii++)
    {
        // thisval = myval[ii]/totval;
        thisval = myval[ii];
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

    

    /*

    int nex = 2;
    SpinorbitalTensor<U>& Rso = R(nex);
    SymmetryBlockedTensor<U>& Rsb = Rso(vec(1,0),vec(0,1));
    CTFTensor<U>& Rctf = Rsb(vector<int>(nex*2,0));
    int max_key = mysize*mysize;
    for (int key = 0;key < max_key;key++)
    {
        Rctf.writeRemoteData(vec(tkv_pair<U>(key,1.0)));
    }
    R(2)(vec(0,0),vec(0,0))(vec(0,0,0,0)) = 0.5*R(2)(vec(1,0),vec(0,1))(vec(0,0,0,0));
    R(2)(vec(2,0),vec(0,2))(vec(0,0,0,0)) = 0.5*R(2)(vec(1,0),vec(0,1))(vec(0,0,0,0));
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");
    Z(0) = 0; Z(1) = 0; Z(2) = 0;
    H.contract(R, Z);
    SpinorbitalTensor<U>& Zso = Z(nex);
    SymmetryBlockedTensor<U>& Zsb = Zso(vec(1,0),vec(0,1));
    CTFTensor<U>& Zctf = Zsb(vector<int>(nex*2,0));
    vector<U> data;
    Zctf.getAllData(data);

    cout << "data = " << data << endl;

    // int mysize = occ.nalpha[0] * vrt.nalpha[0];
    // cout << "mysize = " << mysize << endl;
    vector<complex<double> > w1(mysize);
    vector<U> vl1(mysize*mysize);
    vector<U> vr1(mysize*mysize);
    geev('N','V',mysize,data.data(),mysize,w1.data(),vl1.data(),mysize,vr1.data(),mysize);

    for (int ii = 0; ii<10; ii++)
    {
        cout << w1[ii] << endl;
        cout << vr1[ii*10] << " " << vr1[ii*10 + 1] << " " << vr1[ii*10 + 2] << " " << vr1[ii*10 + 3] << " " << vr1[ii*10 + 4] << " " << vr1[ii*10 + 5] << " " << vr1[ii*10 + 6] << " " << vr1[ii*10 + 7] << " " << vr1[ii*10 + 8] << " " << vr1[ii*10 + 9] << " " << endl;
    }

    */
    // cout << "vr1 = " << vr1 << endl;
    // cout << setprecision(9) << w << endl;

    /*
    vector<vector<U> > cols;

    for (int nex = 1;nex <= 1;nex++)
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
            R(0) = 0; R(1) = 0; R(2) = 0;
            if (arena.rank == 0)
            {
                //cout << "key = " << key << endl;
                Rctf.writeRemoteData(vec(tkv_pair<U>(key,1.0)));
            }
            else
            {
                Rctf.writeRemoteData();
            }

            R(1)(vec(0,0),vec(0,0))(vec(0,0)) = R(1)(vec(1,0),vec(0,1))(vec(0,0));
            // R(2)(vec(0,0),vec(0,0))(vec(0,0,0,0)) = 0.5*R(2)(vec(1,0),vec(0,1))(vec(0,0,0,0));
            // R(2)(vec(2,0),vec(0,2))(vec(0,0,0,0)) = 0.5*R(2)(vec(1,0),vec(0,1))(vec(0,0,0,0));

            // do Hbar*R product into Z
            ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");
            Z(0) = 0; Z(1) = 0; Z(2) = 0;
            H.contract(R, Z);
            for (int nex2 = 1;nex2 <= 1;nex2++)
            {
                SpinorbitalTensor<U>& Zso = Z(nex2);
                SymmetryBlockedTensor<U>& Zsb = Zso(vec(1,0),vec(0,1));
                CTFTensor<U>& Zctf = Zsb(vector<int>(nex2*2,0));
                vector<U> data;
                // cout << "first " << data << endl;
                Zctf.getAllData(data);
                cout << "second " << data << endl;
                cols.back().insert(cols.back().end(), data.begin(), data.end());
                data.assign(data.size(),0.0);
            }
        }
    }
    cout << "cols.size() = " << cols.size() << endl;
    cout << "cols[0].size() = " << cols[0].size() << endl;

    vector<U> matrix;
    for (int i=0;i<cols.size();i++)
    {
        //cout << setprecision(9) << cols[0][i] << endl;
        //cout << i << endl;
        //cout << cols[i] << endl;
        // if (i > 0)
        // {
        //     for (int jj=0;jj<cols.size();jj++)
        //     {
        //         cout << cols[i][jj] - cols[i-1][jj] << endl;
        //     }
        // }
        matrix.insert(matrix.end(),cols[i].begin(),cols[i].end());
    }

    //cout << "matrix.size() = " << matrix.size() << endl;

    int mysize = cols.size();
    vector<complex<double> > w(mysize);
    vector<U> vl(mysize);
    vector<U> vr(mysize);
    geev('N','N',mysize,matrix.data(),mysize,w.data(),vl.data(),mysize,vr.data(),mysize);

    // cout << "vl = " << vl << endl;
    // cout << "vr = " << vr << endl;
    cout << "w.size() = " << w.size() << endl;
    cout << setprecision(9) << w << endl;
    */

}

template <typename U>
void EOMEECCSD<U>::iterate(const Arena& arena)
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    ExcitationOperator<U,2>& R = get<ExcitationOperator<U,2> >("R");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");

    Z = 0;
    H.contract(R, Z);

    energy = davidson.extrapolate(R, Z, D);

    conv = Z.norm(00);
}

INSTANTIATE_SPECIALIZATIONS(EOMEECCSD);
REGISTER_TASK(EOMEECCSD<double>, "eomeeccsd");
