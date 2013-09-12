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

    put("T", new ExcitationOperator<U,3>(arena, occ, vrt));
    puttmp("D", new ExcitationOperator<U,3>(arena, occ, vrt, 0, SY));
    puttmp("Z", new ExcitationOperator<U,3>(arena, occ, vrt));

    ExcitationOperator<U,3>& T = get<ExcitationOperator<U,3> >("T");
    ExcitationOperator<U,3>& D = gettmp<ExcitationOperator<U,3> >("D");
    ExcitationOperator<U,3>& Z = gettmp<ExcitationOperator<U,3> >("Z");

    D(0) = (U)1.0;
    D(1)["ai"]  = H.getIJ()["ii"];
    D(1)["ai"] -= H.getAB()["aa"];
    D(2)["abij"]  = H.getIJ()["ii"];
    D(2)["abij"] -= H.getAB()["aa"];
    D(3)["abcijk"]  = H.getIJ()["ii"];
    D(3)["abcijk"] -= H.getAB()["aa"];

    D(2)(1,0,0,1) *= 2;
    D(3)(1,0,0,1) *= 3;
    D(3)(2,0,0,2) *= 3;

    D = 1/D;

    Z(0) = (U)0.0;
    T(0) = (U)0.0;
    T(1) = H.getAI()*D(1);
    T(2) = H.getABIJ()*D(2);
    T(3) = (U)0.0;

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    energy = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));

    conv =          T(1)(0).norm(00);
    conv = max(conv,T(1)(1).norm(00));
    conv = max(conv,T(2)(0).norm(00));
    conv = max(conv,T(2)(1).norm(00));
    conv = max(conv,T(2)(2).norm(00));

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

    if (isUsed("Hbar"))
    {
        put("Hbar", new STTwoElectronOperator<U,3>(H, T, true));
    }
    */
}

template <typename U>
void CCSDT<U>::iterate()
{
    const TwoElectronOperator<U>& H = get<TwoElectronOperator<U> >("H");

    ExcitationOperator<U,3>& T = get<ExcitationOperator<U,3> >("T");
    ExcitationOperator<U,3>& D = gettmp<ExcitationOperator<U,3> >("D");
    ExcitationOperator<U,3>& Z = gettmp<ExcitationOperator<U,3> >("Z");

    TwoElectronOperator<U> W(H, TwoElectronOperator<U>::AB|
                                TwoElectronOperator<U>::IJ|
                                TwoElectronOperator<U>::IA|
                                TwoElectronOperator<U>::AIBC|
                                TwoElectronOperator<U>::ABCI|
                                TwoElectronOperator<U>::ABCD|
                                TwoElectronOperator<U>::IJKL|
                                TwoElectronOperator<U>::IJKA|
                                TwoElectronOperator<U>::IAJK|
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
    SpinorbitalTensor<U>& WMNIE = W.getIJKA();
    SpinorbitalTensor<U>& WMBIJ = W.getIAJK();
    SpinorbitalTensor<U>& WAMEI = W.getAIBJ();

    //FAE["aa"] = (U)0.0;
    //FMI["ii"] = (U)0.0;

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    /**************************************************************************
     *
     * Intermediates for T(1)->T(1) and T(2)->T(1)
     */
    FME["me"] += WMNEF["mnef"]*T(1)["fn"];

    FMI["mi"] += 0.5*WMNEF["mnef"]*T(2)["efin"];
    FMI["mi"] += FME["me"]*T(1)["ei"];
    FMI["mi"] += WMNIE["mnif"]*T(1)["fn"];

    WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*Tau["efij"];
    WMNIJ["mnij"] += WMNIE["mnie"]*T(1)["ej"];

    WMNIE["mnie"] += WMNEF["mnfe"]*T(1)["fi"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * T(1)->T(1) and T(2)->T(1)
     */
    Z(1)["ai"]  = FAI["ai"];
    Z(1)["ai"] -= T(1)["em"]*WAMEI["amei"];
    Z(1)["ai"] += 0.5*WAMEF["amef"]*Tau["efim"];
    Z(1)["ai"] -= 0.5*WMNIE["mnie"]*T(2)["aemn"];
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

    WMBIJ["mbij"] += 0.5*WAMEF["bmfe"]*Tau["efij"];
    WMBIJ["mbij"] -= WAMEI["bmej"]*T(1)["ei"];

    WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
    WAMEI["amei"] -= WAMEF["amfe"]*T(1)["fi"];
    WAMEI["amei"] += WMNIE["nmie"]*T(1)["an"];
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
    Z(2)["abij"] -= WMBIJ["mbij"]*T(1)["am"];
    Z(2)["abij"] += 0.5*WABEF["abef"]*Tau["efij"];
    Z(2)["abij"] += 0.5*WMNIJ["mnij"]*Tau["abmn"];
    Z(2)["abij"] -= WAMEI["amei"]*T(2)["ebmj"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * Intermediates for CCSDT
     */
    WMBIJ["mbij"] += WMNIE["mnie"]*T(2)["bejn"];
    WMBIJ["mbij"] -= WMNIJ["mnij"]*T(1)["bn"];
    WMBIJ["mbij"] += FME["me"]*T(2)["ebij"];
    WMBIJ["mbij"] += 0.5*WMNEF["mnef"]*T(3)["efbinj"];

    WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
    WAMEI["amei"] -= 0.5*WMNIE["nmie"]*T(1)["an"];

    WABEJ["abej"] += WAMEF["amef"]*T(2)["fbmj"];
    WABEJ["abej"] += 0.5*WMNIE["nmje"]*T(2)["abmn"];
    WABEJ["abej"] += WABEF["abef"]*T(1)["fj"];
    WABEJ["abej"] -= WAMEI["amej"]*T(1)["bm"];
    WABEJ["abej"] -= 0.5*WMNEF["mnef"]*T(3)["afbmnj"];

    WAMEI["amei"] += 0.5*WMNIE["nmie"]*T(1)["an"];

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
    Z(2)["abij"] -= 0.5*WMNIE["mnje"]*T(3)["abeimn"];
    Z(2)["abij"] += FME["me"]*T(3)["abeijm"];

    Z(3)["abcijk"]  = WABEJ["bcek"]*T(2)["aeij"];
    Z(3)["abcijk"] -= WMBIJ["mcjk"]*T(2)["abim"];
    Z(3)["abcijk"] += FAE["ce"]*T(3)["abeijk"];
    Z(3)["abcijk"] -= FMI["mk"]*T(3)["abcijm"];
    Z(3)["abcijk"] += 0.5*WABEF["abef"]*T(3)["efcijk"];
    Z(3)["abcijk"] += 0.5*WMNIJ["mnij"]*T(3)["abcmnk"];
    Z(3)["abcijk"] -= WAMEI["amei"]*T(3)["ebcmjk"];
    /*
     **************************************************************************/

    Z *= D;
    //Z -= T;
    T += Z;

    Tau["abij"]  = T(2)["abij"];
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];
    energy = real(scalar(H.getAI()*T(1))) + 0.25*real(scalar(H.getABIJ()*Tau));

    conv =           Z(1)(0).norm(00);
    conv = max(conv, Z(1)(1).norm(00));
    conv = max(conv, Z(2)(0).norm(00));
    conv = max(conv, Z(2)(1).norm(00));
    conv = max(conv, Z(2)(2).norm(00));
    conv = max(conv, Z(3)(0).norm(00));
    conv = max(conv, Z(3)(1).norm(00));
    conv = max(conv, Z(3)(2).norm(00));
    conv = max(conv, Z(3)(3).norm(00));

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
