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

#include "lccd.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::time;
using namespace aquarius::task;

template <typename U>
LCCD<U>::LCCD(const string& name, const Config& config)
: Iterative("lccd", name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    addProduct(Product("double", "mp2", reqs));
    addProduct(Product("double", "energy", reqs));
    addProduct(Product("double", "convergence", reqs));
    addProduct(Product("double", "S2", reqs));
    addProduct(Product("double", "multiplicity", reqs));
    addProduct(Product("lccd.T", "T", reqs));
    addProduct(Product("lccd.Hbar", "Hbar", reqs));
}

template <typename U>
void LCCD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = get<TwoElectronOperator<U> >("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    put("T", new ExcitationOperator<U,2>("T", arena, occ, vrt));
    puttmp("D", new Denominator<U>(H));
    puttmp("Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));

    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");

    Z(0) = (U)0.0;
    T(0) = (U)0.0;
    T(1) = (U)0.0;
    T(2) = H.getABIJ();

    T.weight(D);

    energy = 0.25*real(scalar(H.getABIJ()*T(2)));

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
void LCCD<U>::iterate(const Arena& arena)
{
    TwoElectronOperator<U>& H_ = get<TwoElectronOperator<U> >("H");

    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");

    TwoElectronOperator<U> H("W", H_, TwoElectronOperator<U>::AB|
                                 TwoElectronOperator<U>::IJ|
                                 TwoElectronOperator<U>::IJKL|
                                 TwoElectronOperator<U>::AIBJ);

    SpinorbitalTensor<U>& FBC = H.getAB();
    SpinorbitalTensor<U>& FKJ = H.getIJ();
    SpinorbitalTensor<U>& WABIJ = H.getABIJ();
    SpinorbitalTensor<U>& WABCD = H.getABCD();
    SpinorbitalTensor<U>& WKLIJ = H.getIJKL();
    SpinorbitalTensor<U>& WBKCJ = H.getAIBJ();

//    sched.set_max_partitions(1);
    /**************************************************************************
     *
     * Intermediates
     */
    // FMI["mi"] += 0.5*WMNEF["mnef"]*T(2)["efin"];


    // WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*T(2)["efij"];
    // FAE["ae"] -= 0.5*WMNEF["mnef"]*T(2)["afmn"];
    // WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * T(1)->T(2) and T(2)->T(2), now aligned with Shavitt and Bartlett 9.126
     */
    Z(2)["abij"] = WABIJ["abij"];
    Z(2)["abij"] += FBC["bc"]*T(2)["acij"];
    Z(2)["abij"] -= FKJ["kj"]*T(2)["abik"];
    Z(2)["abij"] += 0.5*WABCD["abcd"]*T(2)["cdij"];
    Z(2)["abij"] += 0.5*WKLIJ["klij"]*T(2)["abkl"];
    Z(2)["abij"] -= WBKCJ["bkcj"]*T(2)["acik"]; /* minus sign because getAIBJ instead of getIABJ */
    /*
     *************************************************************************/

    Z.weight(D);
    T += Z;

    energy = 0.25*real(scalar(H.getABIJ()*T(2)));

    conv = Z.norm(00);

    diis.extrapolate(T, Z);
}

INSTANTIATE_SPECIALIZATIONS(LCCD);
REGISTER_TASK(LCCD<double>,"lccd");
