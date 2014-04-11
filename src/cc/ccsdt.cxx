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

#include "ccsdt.hpp"

using namespace std;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;

template <typename U>
CCSDT<U>::CCSDT(const string& name, const Config& config)
: Iterative("ccsdt", name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    addProduct(Product("double", "mp2", reqs));
    addProduct(Product("double", "energy", reqs));
    addProduct(Product("double", "convergence", reqs));
    addProduct(Product("double", "S2", reqs));
    addProduct(Product("double", "multiplicity", reqs));
    addProduct(Product("ccsdt.T", "T", reqs));
    addProduct(Product("ccsdt.Hbar", "Hbar", reqs));
}

template <typename U>
void CCSDT<U>::run(task::TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = get<TwoElectronOperator<U> >("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    put("T", new ExcitationOperator<U,3>("T", arena, occ, vrt));
    puttmp("D", new Denominator<U>(H));
    puttmp("Z", new ExcitationOperator<U,3>("Z", arena, occ, vrt));

    ExcitationOperator<U,3>& T = get<ExcitationOperator<U,3> >("T");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,3>& Z = gettmp<ExcitationOperator<U,3> >("Z");

    Z(0) = (U)0.0;
    T(0) = (U)0.0;
    T(1) = H.getAI();
    T(2) = H.getABIJ();
    T(3) = (U)0.0;

    T.weight(D);

    SpinorbitalTensor<U> Tau("Tau", T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    energy = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));

    conv = T.norm(00);

    Logger::log(arena) << "MP2 energy = " << setprecision(15) << energy << endl;
    put("mp2", new Scalar(arena, energy));

    CTF_Timer_epoch ep(name.c_str());
    ep.begin();
    Iterative::run(dag, arena);
    ep.end();

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

    if (isUsed("Hbar"))
    {
        put("Hbar", new STTwoElectronOperator<U,3>(H, T, true));
    }
    */
}

template <typename U>
void CCSDT<U>::iterate(const Arena& arena)
{
    TwoElectronOperator<U>& H = get<TwoElectronOperator<U> >("H");

    ExcitationOperator<U,3>& T = get<ExcitationOperator<U,3> >("T");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,3>& Z = gettmp<ExcitationOperator<U,3> >("Z");

    TwoElectronOperator<U> W("W", H, TwoElectronOperator<U>::AB|
                                TwoElectronOperator<U>::IJ|
                                TwoElectronOperator<U>::IA|
                                TwoElectronOperator<U>::AIBC|
                                TwoElectronOperator<U>::ABCI|
                                TwoElectronOperator<U>::ABCD|
                                TwoElectronOperator<U>::IJKL|
                                TwoElectronOperator<U>::IJAK|
                                TwoElectronOperator<U>::AIJK|
                                TwoElectronOperator<U>::AIBJ);

    SpinorbitalTensor<U>& FAI = W.getAI();
    SpinorbitalTensor<U>& FME = W.getIA();
    SpinorbitalTensor<U>& FAE = W.getAB();
    SpinorbitalTensor<U>& FMI = W.getIJ();
    SpinorbitalTensor<U>& WABIJ = W.getABIJ();
    SpinorbitalTensor<U>& WMNEF = W.getIJAB();
    SpinorbitalTensor<U>& WAMEF = W.getAIBC();
    SpinorbitalTensor<U>& WABEJ = W.getABCI();
    SpinorbitalTensor<U>& WABEF = W.getABCD();
    SpinorbitalTensor<U>& WMNIJ = W.getIJKL();
    SpinorbitalTensor<U>& WMNEJ = W.getIJAK();
    SpinorbitalTensor<U>& WAMIJ = W.getAIJK();
    SpinorbitalTensor<U>& WAMEI = W.getAIBJ();

    SpinorbitalTensor<U> Tau("Tau", T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    /**************************************************************************
     *
     * Intermediates for T(1)->T(1) and T(2)->T(1)
     */
    FME["me"] += WMNEF["mnef"]*T(1)["fn"];

    FMI["mi"] += 0.5*WMNEF["mnef"]*T(2)["efin"];
    FMI["mi"] += FME["me"]*T(1)["ei"];
    FMI["mi"] += WMNEJ["nmfi"]*T(1)["fn"];

    WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*Tau["efij"];
    WMNIJ["mnij"] += WMNEJ["mnej"]*T(1)["ei"];

    WMNEJ["mnej"] += WMNEF["mnef"]*T(1)["fj"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * T(1)->T(1) and T(2)->T(1)
     */
    Z(1)["ai"]  = FAI["ai"];
    Z(1)["ai"] -= T(1)["em"]*WAMEI["amei"];
    Z(1)["ai"] += 0.5*WAMEF["amef"]*Tau["efim"];
    Z(1)["ai"] -= 0.5*WMNEJ["mnei"]*T(2)["eamn"];
    Z(1)["ai"] += T(2)["aeim"]*FME["me"];
    Z(1)["ai"] += T(1)["ei"]*FAE["ae"];
    Z(1)["ai"] -= T(1)["am"]*FMI["mi"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * Intermediates for T(1)->T(2) and T(2)->T(2)
     */
    FAE["ae"] -= 0.5*WMNEF["mnef"]*T(2)["afmn"];
    FAE["ae"] -= FME["me"]*T(1)["am"];
    FAE["ae"] += WAMEF["amef"]*T(1)["fm"];

    WAMIJ["amij"] += 0.5*WAMEF["amef"]*Tau["efij"];
    WAMIJ["amij"] += WAMEI["amej"]*T(1)["ei"];

    WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
    WAMEI["amei"] -= WAMEF["amfe"]*T(1)["fi"];
    WAMEI["amei"] -= WMNEJ["nmei"]*T(1)["an"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * T(1)->T(2) and T(2)->T(2)
     */
    Z(2)["abij"]  = WABIJ["abij"];
    Z(2)["abij"] += FAE["af"]*T(2)["fbij"];
    Z(2)["abij"] -= FMI["ni"]*T(2)["abnj"];
    Z(2)["abij"] += WABEJ["abej"]*T(1)["ei"];
    Z(2)["abij"] -= WAMIJ["amij"]*T(1)["bm"];
    Z(2)["abij"] += 0.5*WABEF["abef"]*Tau["efij"];
    Z(2)["abij"] += 0.5*WMNIJ["mnij"]*Tau["abmn"];
    Z(2)["abij"] -= WAMEI["amei"]*T(2)["ebmj"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * Intermediates for CCSDT
     */
    WAMIJ["amij"] += WMNEJ["nmej"]*T(2)["aein"];
    WAMIJ["amij"] -= WMNIJ["nmij"]*T(1)["an"];
    WAMIJ["amij"] += FME["me"]*T(2)["aeij"];
    WAMIJ["amij"] += 0.5*WMNEF["mnef"]*T(3)["aefijn"];

    WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
    WAMEI["amei"] += 0.5*WMNEJ["nmei"]*T(1)["an"];

    WABEJ["abej"] += WAMEF["amef"]*T(2)["fbmj"];
    WABEJ["abej"] += 0.5*WMNEJ["mnej"]*T(2)["abmn"];
    WABEJ["abej"] += WABEF["abef"]*T(1)["fj"];
    WABEJ["abej"] -= WAMEI["amej"]*T(1)["bm"];
    WABEJ["abej"] -= 0.5*WMNEF["mnef"]*T(3)["afbmnj"];

    WAMEI["amei"] -= 0.5*WMNEJ["nmei"]*T(1)["an"];

    WABEF["abef"] -= WAMEF["amef"]*T(1)["bm"];
    WABEF["abef"] += 0.5*WMNEF["mnef"]*Tau["abmn"];

    WAMEF["amef"] -= WMNEF["nmef"]*T(1)["an"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * CCSDT Iteration
     */
    Z(1)["ai"] += 0.25*WMNEF["mnef"]*T(3)["aefimn"];

    Z(2)["abij"] += 0.5*WAMEF["bmef"]*T(3)["aefijm"];
    Z(2)["abij"] -= 0.5*WMNEJ["mnej"]*T(3)["abeinm"];
    Z(2)["abij"] += FME["me"]*T(3)["abeijm"];

    Z(3)["abcijk"]  = WABEJ["bcek"]*T(2)["aeij"];
    Z(3)["abcijk"] -= WAMIJ["bmjk"]*T(2)["acim"];
    Z(3)["abcijk"] += FAE["ce"]*T(3)["abeijk"];
    Z(3)["abcijk"] -= FMI["mk"]*T(3)["abcijm"];
    Z(3)["abcijk"] += 0.5*WABEF["abef"]*T(3)["efcijk"];
    Z(3)["abcijk"] += 0.5*WMNIJ["mnij"]*T(3)["abcmnk"];
    Z(3)["abcijk"] -= WAMEI["amei"]*T(3)["ebcmjk"];
    /*
     **************************************************************************/

    //cout << Z(3)(vec(0,0),vec(0,0)).norm(00) << endl;
    //cout << Z(3)(vec(1,0),vec(0,1)).norm(00) << endl;
    //cout << Z(3)(vec(2,0),vec(0,2)).norm(00) << endl;
    //cout << Z(3)(vec(3,0),vec(0,3)).norm(00) << endl;

    Z.weight(D);
    T += Z;

    Tau["abij"]  = T(2)["abij"];
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];
    energy = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));

    conv = Z.norm(00);

    diis.extrapolate(T, Z);
}

/*
template <typename U>
double CCSDT<U>::getProjectedS2() const
{
    return CCSD<U>::getProjectedS2(this->occ, this->vrt, T(1), T(2));
}

template <typename U>
double CCSDT<U>::getProjectedMultiplicity() const
{
    return CCSD<U>::getProjectedS2(this->occ, this->vrt, T(1), T(2));
}
*/

INSTANTIATE_SPECIALIZATIONS(CCSDT);
REGISTER_TASK(CCSDT<double>,"ccsdt");
