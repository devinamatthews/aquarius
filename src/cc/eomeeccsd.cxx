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
: MultiIterative<U>("eomeeccsd", name, config), nroot(config.get<int>("nroot")), multiroot(config.get<int>("multiroot")), davidson(config.get("davidson"),multiroot*nroot+1)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("ccsd.T", "T"));
    reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
    reqs.push_back(Requirement("tda.TDAevals", "TDAevals"));
    reqs.push_back(Requirement("tda.TDAevecs", "TDAevecs"));
    addProduct(Product("eomeeccsd.energy", "energy", reqs));
    addProduct(Product("double", "convergence", reqs));
    if (multiroot == 0)
    {
        addProduct(Product("eomeeccsd.R", "R", reqs));
    }
    else
    {
        assert(multiroot == 1);
        assert(nroot != 0);
        addProduct(Product("eomeeccsd.R", "R0", reqs));
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

    }
}

template <typename U>
void EOMEECCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    int mysize = occ.nalpha[0] * vrt.nalpha[0];
    puttmp("D", new Denominator<U>(H));
    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");
    CTFTensor<U>& TDAevecs = get<CTFTensor<U> >("TDAevecs");
    vector<U> TDAevecsVec;
    TDAevecs.getAllData(TDAevecsVec);


    if (multiroot == 0)
    {
        put("R", new ExcitationOperator<U,2>("R", arena, occ, vrt));
        puttmp("Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));
        ExcitationOperator<U,2>& R = get<ExcitationOperator<U,2> >("R");
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

        put("energy", new CTFTensor<U>("energy", arena, 1, vec(multiroot*nroot+1), vec(NS), true));
        MultiIterative::run(dag, arena);
        put("convergence", new Scalar(arena, conv));
    }
    else
    {
        put("R0", new ExcitationOperator<U,2>("R0", arena, occ, vrt));
        puttmp("Z0", new ExcitationOperator<U,2>("Z0", arena, occ, vrt));
        ExcitationOperator<U,2>& R0 = get<ExcitationOperator<U,2> >("R0");
        put("R1", new ExcitationOperator<U,2>("R1", arena, occ, vrt));
        puttmp("Z1", new ExcitationOperator<U,2>("Z1", arena, occ, vrt));
        ExcitationOperator<U,2>& R1 = get<ExcitationOperator<U,2> >("R1");
        R0(0) = 0; R0(1) = 0; R0(2) = 0;
        R1(0) = 0; R1(1) = 0; R1(2) = 0;
        double thisval0;
        double thisval1;

        for (int ii=0; ii<mysize; ii++)
        {
            thisval0 = TDAevecsVec[ii];
            thisval1 = TDAevecsVec[mysize+ii];
            if (R0(1).arena.rank == 0)
                R0(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval0/sqrt(2))));
            else
                R0(1)(vec(1,0),vec(0,1))(vec(0,0)).writeRemoteData();

            if (R0(1).arena.rank == 0)
                R0(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData(vec(kv_pair(ii,thisval0/sqrt(2))));
            else
                R0(1)(vec(0,0),vec(0,0))(vec(0,0)).writeRemoteData();

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
            ExcitationOperator<U,2>& R2 = get<ExcitationOperator<U,2> >("R2");
            R2(0) = 0; R2(1) = 0; R2(2) = 0;
            double thisval2;
            for (int ii=0; ii<mysize; ii++)
            {
                thisval2 = TDAevecsVec[2*mysize+ii];
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
            ExcitationOperator<U,2>& R3 = get<ExcitationOperator<U,2> >("R3");
            R3(0) = 0; R3(1) = 0; R3(2) = 0;
            double thisval3;
            for (int ii=0; ii<mysize; ii++)
            {
                thisval3 = TDAevecsVec[3*mysize+ii];
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
            ExcitationOperator<U,2>& R4 = get<ExcitationOperator<U,2> >("R4");
            R4(0) = 0; R4(1) = 0; R4(2) = 0;
            double thisval4;
            for (int ii=0; ii<mysize; ii++)
            {
                thisval4 = TDAevecsVec[4*mysize+ii];
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
            ExcitationOperator<U,2>& R5 = get<ExcitationOperator<U,2> >("R5");
            R5(0) = 0; R5(1) = 0; R5(2) = 0;
            double thisval5;
            for (int ii=0; ii<mysize; ii++)
            {
                thisval5 = TDAevecsVec[5*mysize+ii];
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
    }

}

template <typename U>
void EOMEECCSD<U>::iterate(const Arena& arena)
{
    const STTwoElectronOperator<U,2>& H = get<STTwoElectronOperator<U,2> >("Hbar");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    int mysize = occ.nalpha[0] * vrt.nalpha[0];

    CTFTensor<U>& TDAevals = get<CTFTensor<U> >("TDAevals");
    vector<U> TDAevalsVec;
    TDAevals.getAllData(TDAevalsVec);
    Denominator<U>& D = gettmp<Denominator<U> >("D");

    vector<double> energyvec(multiroot*nroot+1);
    vector<tkv_pair<U> > pairs;

    CTFTensor<U>& energy = gettmp<CTFTensor<U> >("energy");

    if (multiroot == 0)
    {
        ExcitationOperator<U,2>& R = get<ExcitationOperator<U,2> >("R");
        ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");
        Z = 0;
        H.contract(R, Z);
        energyvec = davidson.extrapolate(R, Z, D, TDAevalsVec[nroot]);
        cout << "evec = " << energyvec << endl;
        conv = Z.norm(00);
        pairs.push_back(tkv_pair<U>(0, energyvec[0]));
    }
    else
    {
        ExcitationOperator<U,2>& R0 = get<ExcitationOperator<U,2> >("R0");
        ExcitationOperator<U,2>& Z0 = gettmp<ExcitationOperator<U,2> >("Z0");
        Z0 = 0;
        H.contract(R0, Z0);
        ExcitationOperator<U,2>& R1 = get<ExcitationOperator<U,2> >("R1");
        ExcitationOperator<U,2>& Z1 = gettmp<ExcitationOperator<U,2> >("Z1");
        Z1 = 0;
        H.contract(R1, Z1);


        energyvec = davidson.extrapolate(R0, R1, Z0, Z1, D, TDAevalsVec[0], TDAevalsVec[1]);
        conv = Z0.norm(00);

        for (int ii = 0; ii < nroot; ii++)
        {
            pairs.push_back(tkv_pair<U>(ii, energyvec[ii]));
        }
    }
    if (arena.rank == 0)
        energy.writeRemoteData(pairs);
    else
        energy.writeRemoteData();

}

INSTANTIATE_SPECIALIZATIONS(EOMEECCSD);
REGISTER_TASK(EOMEECCSD<double>, "eomeeccsd");
