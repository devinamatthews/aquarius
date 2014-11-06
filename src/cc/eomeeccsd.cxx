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
using namespace aquarius::convergence;

template <typename U>
EOMEECCSD<U>::EOMEECCSD(const std::string& name, const Config& config)
: MultiIterative<U>("eomeeccsd", name, config), nroot(config.get<int>("nroot"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.T", "T"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("tda.TDAevals", "TDAevals"));
    reqs.push_back(Requirement("tda.TDAevecs", "TDAevecs"));
    addProduct(Product("eomeeccsd.energy", "energy", reqs));
    addProduct(Product("double", "convergence", reqs));

    assert(nroot > 0);
    addProduct(Product("eomeeccsd.R", "R1", reqs));
    if (nroot > 1)
    {
        addProduct(Product("eomeeccsd.R", "R2", reqs));
    }
    if (nroot > 2)
    {
        addProduct(Product("eomeeccsd.R", "R3", reqs));
    }
    if (nroot > 3)
    {
        addProduct(Product("eomeeccsd.R", "R4", reqs));
    }
    if (nroot > 4)
    {
        addProduct(Product("eomeeccsd.R", "R5", reqs));
    }
    if (nroot > 5)
    {
        cout << "Too many roots! Not yet available." << endl;
    }

    puttmp("Davidson", new Davidson<ExcitationOperator<U,2> >(config.get("davidson"), nroot));
}

template <typename U>
void EOMEECCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const STTwoElectronOperator<U,2>& H = this->template get<STTwoElectronOperator<U,2> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    int mysize = occ.nalpha[0] * vrt.nalpha[0];
    puttmp("D", new Denominator<U>(H));
    ExcitationOperator<U,2>& T = this->template get<ExcitationOperator<U,2> >("T");
    CTFTensor<U>& TDAevecs = this->template get<CTFTensor<U> >("TDAevecs");
    vector<U> TDAevecsVec;
    TDAevecs.getAllData(TDAevecsVec);

    put("R1", new ExcitationOperator<U,2>("R1", arena, occ, vrt));
    puttmp("Z1", new ExcitationOperator<U,2>("Z1", arena, occ, vrt));
    ExcitationOperator<U,2>& R1 = this->template get<ExcitationOperator<U,2> >("R1");
    R1(0) = 0; R1(1) = 0; R1(2) = 0;
    for (int ii=0; ii<mysize; ii++)
    {
        double thisval1 = TDAevecsVec[ii];

        if (R1(1).arena.rank == 0)
            R1(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval1/sqrt(2))));
        else
            R1(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData();

        if (R1(1).arena.rank == 0)
            R1(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval1/sqrt(2))));
        else
            R1(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData();
    }

    if (nroot > 1)
    {
        put("R2", new ExcitationOperator<U,2>("R2", arena, occ, vrt));
        puttmp("Z2", new ExcitationOperator<U,2>("Z2", arena, occ, vrt));
        ExcitationOperator<U,2>& R2 = this->template get<ExcitationOperator<U,2> >("R2");
        R2(0) = 0; R2(1) = 0; R2(2) = 0;
        double thisval2;
        for (int ii=0; ii<mysize; ii++)
        {
            thisval2 = TDAevecsVec[mysize+ii];
            if (R2(1).arena.rank == 0)
                R2(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval2/sqrt(2))));
            else
                R2(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData();

            if (R2(1).arena.rank == 0)
                R2(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval2/sqrt(2))));
            else
                R2(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData();
        }
    }
    if (nroot > 2)
    {
        put("R3", new ExcitationOperator<U,2>("R3", arena, occ, vrt));
        puttmp("Z3", new ExcitationOperator<U,2>("Z3", arena, occ, vrt));
        ExcitationOperator<U,2>& R3 = this->template get<ExcitationOperator<U,2> >("R3");
        R3(0) = 0; R3(1) = 0; R3(2) = 0;
        double thisval3;
        for (int ii=0; ii<mysize; ii++)
        {
            thisval3 = TDAevecsVec[2*mysize+ii];
            if (R3(1).arena.rank == 0)
                R3(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval3/sqrt(2))));
            else
                R3(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData();

            if (R3(1).arena.rank == 0)
                R3(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval3/sqrt(2))));
            else
                R3(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData();
        }
    }
    if (nroot > 3)
    {
        put("R4", new ExcitationOperator<U,2>("R4", arena, occ, vrt));
        puttmp("Z4", new ExcitationOperator<U,2>("Z4", arena, occ, vrt));
        ExcitationOperator<U,2>& R4 = this->template get<ExcitationOperator<U,2> >("R4");
        R4(0) = 0; R4(1) = 0; R4(2) = 0;
        double thisval4;
        for (int ii=0; ii<mysize; ii++)
        {
            thisval4 = TDAevecsVec[3*mysize+ii];
            if (R4(1).arena.rank == 0)
                R4(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval4/sqrt(2))));
            else
                R4(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData();

            if (R4(1).arena.rank == 0)
                R4(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval4/sqrt(2))));
            else
                R4(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData();
        }
    }
    if (nroot > 4)
    {
        put("R5", new ExcitationOperator<U,2>("R5", arena, occ, vrt));
        puttmp("Z5", new ExcitationOperator<U,2>("Z5", arena, occ, vrt));
        ExcitationOperator<U,2>& R5 = this->template get<ExcitationOperator<U,2> >("R5");
        R5(0) = 0; R5(1) = 0; R5(2) = 0;
        double thisval5;
        for (int ii=0; ii<mysize; ii++)
        {
            thisval5 = TDAevecsVec[4*mysize+ii];
            if (R5(1).arena.rank == 0)
                R5(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval5/sqrt(2))));
            else
                R5(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData();

            if (R5(1).arena.rank == 0)
                R5(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval5/sqrt(2))));
            else
                R5(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData();
        }
    }

    // Iterative::run(dag, arena);
    // put("TDAevals", new CTFTensor<U>("TDAevals", arena, 1, vec(mysize), vec(NS), true));
    // put("energy", new Scalar(arena, energy));
    // put("convergence", new Scalar(arena, conv));

    put("energy", new CTFTensor<U>("energy", arena, 1, vec(nroot), vec(NS), true));
    MultiIterative<U>::run(dag, arena, nroot);
    put("convergence", new U(conv));
}

template <typename U>
void EOMEECCSD<U>::iterate(const Arena& arena)
{
    const STTwoElectronOperator<U,2>& H = this->template get<STTwoElectronOperator<U,2> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    Davidson<ExcitationOperator<U,2> >& davidson =
        this->template gettmp<Davidson<ExcitationOperator<U,2> > >("Davidson");

    Denominator<U>& D = this->template gettmp<Denominator<U> >("D");
    std::vector<ExcitationOperator<U,2>*> R;
    std::vector<ExcitationOperator<U,2>*> Z;

    ExcitationOperator<U,2>& R1 = this->template get<ExcitationOperator<U,2> >("R1");
    ExcitationOperator<U,2>& Z1 = this->template gettmp<ExcitationOperator<U,2> >("Z1");
    Z1 = 0;
    H.contract(R1, Z1);
    R.push_back(&R1);
    Z.push_back(&Z1);

    if (nroot > 1)
    {
        ExcitationOperator<U,2>& R2 = this->template get<ExcitationOperator<U,2> >("R2");
        ExcitationOperator<U,2>& Z2 = this->template gettmp<ExcitationOperator<U,2> >("Z2");
        Z2 = 0;
        H.contract(R2, Z2);
        R.push_back(&R2);
        Z.push_back(&Z2);
    }

    if (nroot > 2)
    {
        ExcitationOperator<U,2>& R3 = this->template get<ExcitationOperator<U,2> >("R3");
        ExcitationOperator<U,2>& Z3 = this->template gettmp<ExcitationOperator<U,2> >("Z3");
        Z3 = 0;
        H.contract(R3, Z3);
        R.push_back(&R3);
        Z.push_back(&Z3);
    }

    if (nroot > 3)
    {
        ExcitationOperator<U,2>& R4 = this->template get<ExcitationOperator<U,2> >("R4");
        ExcitationOperator<U,2>& Z4 = this->template gettmp<ExcitationOperator<U,2> >("Z4");
        Z4 = 0;
        H.contract(R4, Z4);
        R.push_back(&R4);
        Z.push_back(&Z4);
    }

    if (nroot > 4)
    {
        ExcitationOperator<U,2>& R5 = this->template get<ExcitationOperator<U,2> >("R5");
        ExcitationOperator<U,2>& Z5 = this->template gettmp<ExcitationOperator<U,2> >("Z5");
        Z5 = 0;
        H.contract(R5, Z5);
        R.push_back(&R5);
        Z.push_back(&Z5);
    }

    vector<U> energyvec = davidson.extrapolate(R, Z, D);
    cout << "energyvec[0] = " << energyvec[0] << endl;
    conv = Z[0]->norm(00);
    cout << "conv = " << conv << endl;
    // for (int i = 1;i < nroot;i++)
    //     conv = max(conv, Z[i]->norm(00));

    // cout << "energyvec = " << energyvec << endl;

    vector<tkv_pair<U> > pairs;
    for (int ii = 0; ii < nroot; ii++)
    {
        pairs.push_back(kv_pair(ii, energyvec[ii]));
        // cout << "pair = " << ii << ", " << energyvec[ii] << endl;
    }

    CTFTensor<U>& energy = this->template gettmp<CTFTensor<U> >("energy");
    if (arena.rank == 0)
        energy.writeRemoteData(pairs);
    else
        energy.writeRemoteData();

    // vector<U> data1;
    // energy.getAllData(data1);
    // cout << "data1.size() = " << data1.size() << endl;

}

INSTANTIATE_SPECIALIZATIONS(EOMEECCSD);
REGISTER_TASK(EOMEECCSD<double>, "eomeeccsd");
