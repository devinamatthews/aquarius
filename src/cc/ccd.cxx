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

#include "ccd.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::time;
using namespace aquarius::task;

template <typename U>
CCD<U>::CCD(const string& name, const Config& config)
: Iterative("ccd", name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    addProduct(Product("double", "mp2", reqs));
    addProduct(Product("double", "energy", reqs));
    addProduct(Product("double", "convergence", reqs));
    addProduct(Product("double", "S2", reqs));
    addProduct(Product("double", "multiplicity", reqs));
    addProduct(Product("ccd.T", "T", reqs));
    addProduct(Product("ccd.Hbar", "Hbar", reqs));
}

template <typename U>
void CCD<U>::run(TaskDAG& dag, const Arena& arena)
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
    put("mp2", new U(energy));

    Iterative::run(dag, arena);

    put("energy", new U(energy));
    put("convergence", new U(conv));

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
void CCD<U>::iterate(const Arena& arena)
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
    SpinorbitalTensor<U>& WIJAB = H.getIJAB();
    SpinorbitalTensor<U>& WABCD = H.getABCD();
    SpinorbitalTensor<U>& WKLIJ = H.getIJKL();
    SpinorbitalTensor<U>& WBKCJ = H.getAIBJ();

//    sched.set_max_partitions(1);
    /**************************************************************************
     *
     * Intermediates, now aligned with Shavitt and Bartlett, 9.126
     */
    FKJ["kj"] += 0.5*WIJAB["klcd"]*T(2)["dclj"]; /* 9, through 3 */
    WKLIJ["klij"] += 0.5*WIJAB["klcd"]*T(2)["cdij"]; /* 7, through 5 */
    FBC["bc"] -= 0.5*WIJAB["klcd"]*T(2)["dblk"]; /* 10, through 2 */
    WBKCJ["bkcj"] -= 0.5*WIJAB["klcd"]*T(2)["bdjl"]; /* 8, through 6, Minus sign cancels that in 6*/
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * T(1)->T(2) and T(2)->T(2), now aligned with Shavitt and Bartlett 9.126
     */
    Z(2)["abij"] = WIJAB["ijab"]; /* 1, ijab instead of abij because getIJAB instead of getABIJ. Makes intermeds 7, 8 nicer. */
    Z(2)["abij"] += FBC["bc"]*T(2)["acij"]; /* 2 */
    Z(2)["abij"] -= FKJ["kj"]*T(2)["abik"]; /* 3 */
    Z(2)["abij"] += 0.5*WABCD["abcd"]*T(2)["cdij"]; /* 4 */
    Z(2)["abij"] += 0.5*WKLIJ["klij"]*T(2)["abkl"]; /* 5 */
    Z(2)["abij"] -= WBKCJ["bkcj"]*T(2)["acik"]; /* 6, minus sign because getAIBJ instead of getIABJ, which necessitates -bkcj instead of +kbcj */
    /*
     *************************************************************************/

    Z.weight(D);
    T += Z;

    energy = 0.25*real(scalar(H.getABIJ()*T(2)));

    conv = Z.norm(00);

    diis.extrapolate(T, Z);
}

INSTANTIATE_SPECIALIZATIONS(CCD);
REGISTER_TASK(CCD<double>,"ccd");
