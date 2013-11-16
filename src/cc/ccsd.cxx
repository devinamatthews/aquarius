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

#include "ccsd.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;

template <typename U>
CCSD<U>::CCSD(const std::string& name, const Config& config)
: Iterative("ccsd", name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    addProduct(Product("double", "mp2", reqs));
    addProduct(Product("double", "energy", reqs));
    addProduct(Product("double", "convergence", reqs));
    addProduct(Product("double", "S2", reqs));
    addProduct(Product("double", "multiplicity", reqs));
    addProduct(Product("ccsd.T", "T", reqs));
    addProduct(Product("ccsd.Hbar", "Hbar", reqs));
}

template <typename U>
void CCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = get<TwoElectronOperator<U> >("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    put("T", new ExcitationOperator<U,2>("T", arena, occ, vrt));
    puttmp("Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));
    puttmp("Tau", new SpinorbitalTensor<U>("Tau", H.getABIJ()));
    puttmp("D", new Denominator<U>(H));
    puttmp("W", new TwoElectronOperator<U>("W", const_cast<TwoElectronOperator<U>&>(H),
                                           TwoElectronOperator<U>::AB|
                                           TwoElectronOperator<U>::IJ|
                                           TwoElectronOperator<U>::IA|
                                           TwoElectronOperator<U>::IJKL|
                                           TwoElectronOperator<U>::IJAK|
                                           TwoElectronOperator<U>::AIJK|
                                           TwoElectronOperator<U>::AIBJ));

    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");
    SpinorbitalTensor<U>& Tau = gettmp<SpinorbitalTensor<U> >("Tau");

    Z(0) = (U)0.0;
    T(0) = (U)0.0;
    T(1) = H.getAI();
    T(2) = H.getABIJ();

    T.weight(D);

    //cout << setprecision(15) << real(scalar(T(1)(0)*T(1)(0))) << " " << pow(T(1)(0).norm(2),2) << endl;
    //cout << setprecision(15) << real(scalar(T(1)(1)*T(1)(1))) << " " << pow(T(1)(1).norm(2),2) << endl;
    //cout << setprecision(15) << real(scalar(T(2)(0)*T(2)(0))) << " " << pow(T(2)(0).norm(2),2) << endl;
    //cout << setprecision(15) << real(scalar(T(2)(1)*T(2)(1))) << " " << pow(T(2)(1).norm(2),2) << endl;
    //cout << setprecision(15) << real(scalar(T(2)(2)*T(2)(2))) << " " << pow(T(2)(2).norm(2),2) << endl;
    //cout << setprecision(15) << real(scalar(T(1)*T(1))) << " " << pow(T(1).norm(2),2) << endl;
    //cout << setprecision(15) << real(scalar(T(2)*T(2))) << " " << pow(T(2).norm(2),2) << endl;

    //exit(1);

    Tau = T(2);
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    energy = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));

    conv = T.norm(00);

    Logger::log(arena) << "MP2 energy = " << setprecision(15) << energy << endl;
    put("mp2", new Scalar(arena, energy));

    Iterative::run(dag, arena);

    put("energy", new Scalar(arena, energy));
    put("convergence", new Scalar(arena, conv));

    /*
    if (isUsed("S2") || isUsed("multiplicity"))
    {
        double s2 = getProjectedS2(occ, vrt, T(1), T(2));
        double mult = sqrt(4*s2+1);

        put("S2", new Scalar(arena, s2));
        put("multiplicity", new Scalar(arena, mult));
    }
    */

    if (isUsed("Hbar"))
    {
        put("Hbar", new STTwoElectronOperator<U,2>("Hbar", H, T, true));
    }
}

template <typename U>
void CCSD<U>::iterate()
{
    const TwoElectronOperator<U>& H = get<TwoElectronOperator<U> >("H");

    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");
    TwoElectronOperator<U>& W = gettmp<TwoElectronOperator<U> >("W");
    SpinorbitalTensor<U>& Tau = gettmp<SpinorbitalTensor<U> >("Tau");

    STExcitationOperator<U,2>::transform(H, T, Tau, Z, W);

    Z.weight(D);
    T += Z;

    Tau = T(2);
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    energy = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));

    conv = Z.norm(00);

    diis.extrapolate(T, Z);
}

/*
template <typename U>
double CCSD<U>::getProjectedS2(const MOSpace<U>& occ, const MOSpace<U>& vrt,
                               const SpinorbitalTensor<U>& T1,
                               const SpinorbitalTensor<U>& T2)
{
    const Arena& arena = occ.arena;

    int N = occ.nao;
    int nI = occ.nalpha;
    int ni = occ.nbeta;
    int nA = vrt.nalpha;
    int na = vrt.nbeta;

    vector<int> shapeNN = vec(NS,NS);
    vector<int> shapeNNNN = vec(NS,NS,NS,NS);
    vector<int> sizeAI = vec(nA,nI);
    vector<int> sizeAi = vec(nA,ni);
    vector<int> sizeaI = vec(na,nI);
    vector<int> sizeai = vec(na,ni);
    vector<int> sizeIi = vec(nI,ni);
    vector<int> sizeIn = vec(nI,N);
    vector<int> sizein = vec(ni,N);
    vector<int> sizeAaIi = vec(nA,na,nI,ni);

    const CTFTensor<U>& CA = vrt.Calpha;
    const CTFTensor<U>& Ca = vrt.Cbeta;
    const CTFTensor<U>& CI = occ.Calpha;
    const CTFTensor<U>& Ci = occ.Cbeta;

    //TODO
    CTFTensor<U> S(arena, 2, vec(N,N), shapeNN, true);

    CTFTensor<U> DAI(arena, 2, sizeAI, shapeNN, false);
    CTFTensor<U> DAi(arena, 2, sizeAi, shapeNN, false);
    CTFTensor<U> DaI(arena, 2, sizeaI, shapeNN, false);
    CTFTensor<U> Dai(arena, 2, sizeai, shapeNN, false);
    CTFTensor<U> DIj(arena, 2, sizeIi, shapeNN, false);
    CTFTensor<U> DAbIj(arena, 4, sizeAaIi, shapeNNNN, false);
    CTFTensor<U> tmp1(arena, 2, sizeIn, shapeNN, false);
    CTFTensor<U> tmp2(arena, 2, sizein, shapeNN, false);

    tmp1["Iq"] = CI["pI"]*S["pq"];
    DIj["Ij"] = tmp1["Iq"]*Ci["qj"];
    DaI["aI"] = tmp1["Iq"]*Ca["qa"];

    tmp2["iq"] = Ci["pi"]*S["pq"];
    DAi["Ai"] = tmp2["iq"]*CA["qA"];

    DAI["AI"] = DAi["Aj"]*DIj["Ij"];
    Dai["ai"] = DaI["aJ"]*DIj["Ji"];
    DAbIj["AbIj"] = DAi["Aj"]*DaI["bI"];

    const CTFTensor<U>& T1A = T1(1,0,0,1);
    const CTFTensor<U>& T1B = T1(0,0,0,0);
    CTFTensor<U> TauAB(T2(1,0,0,1));

    TauAB["AbIj"] += T1A["AI"]*T1B["bj"];

    U S2 = (U)0;

    U S2T11 = -scalar(DAI*T1A);
    U S2T12 = -scalar(Dai*T1B);
    U S2T2 = -scalar(DAbIj*TauAB);

    return abs(S2+S2T11+S2T12+S2T2);
}
*/

INSTANTIATE_SPECIALIZATIONS(CCSD);
REGISTER_TASK(CCSD<double>,"ccsd");
