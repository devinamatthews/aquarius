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
using namespace aquarius::task;
using namespace aquarius::time;

template <typename U>
LCCD<U>::LCCD(const std::string& name, const Config& config)
: Iterative<U>("lccd", name, config), diis(config.get("diis"))
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    this->addProduct(Product("double", "energy", reqs));
    this->addProduct(Product("double", "convergence", reqs));
    this->addProduct(Product("double", "S2", reqs));
    this->addProduct(Product("double", "multiplicity", reqs));
    this->addProduct(Product("lccd.T", "T", reqs));
    this->addProduct(Product("lccd.Hbar", "Hbar", reqs));
}

template <typename U>
void LCCD<U>::run(TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U> >("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    this->put   (  "T", new ExcitationOperator<U,2>("T", arena, occ, vrt));
    this->puttmp(  "Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));
    this->puttmp(  "D", new Denominator<U>(H));

    ExcitationOperator<U,2>& T = this->template get   <ExcitationOperator<U,2> >(  "T");
    Denominator<U>&          D = this->template gettmp<Denominator<U> >         (  "D");
    ExcitationOperator<U,2>& Z = this->template gettmp<ExcitationOperator<U,2> >(  "Z");

    Z(0) = 0;
    T(0) = 0;
    T(1) = 0;
    T(2) = H.getABIJ();

    T.weight(D);

    CTF_Timer_epoch ep(this->name.c_str());
    ep.begin();
    Iterative<U>::run(dag, arena);
    ep.end();

    this->put("energy", new U(this->energy()));
    this->put("convergence", new U(this->conv()));

    /*
    if (isUsed("S2") || isUsed("multiplicity"))
    {
        double s2 = this->template getProjectedS2(occ, vrt, T(1), T(2));
        double mult = sqrt(4*s2+1);

        this->put("S2", new Scalar(arena, s2));
        this->put("multiplicity", new Scalar(arena, mult));
    }
    */

    if (this->isUsed("Hbar"))
    {
        this->put("Hbar", new STTwoElectronOperator<U>("Hbar", H, T, true));
    }
}

template <typename U>
void LCCD<U>::iterate(const Arena& arena)
{
    const TwoElectronOperator<U>& H = this->template get<TwoElectronOperator<U> >("H");

    const SpinorbitalTensor<U>&   fAE =   H.getAB();
    const SpinorbitalTensor<U>&   fMI =   H.getIJ();
    const SpinorbitalTensor<U>& VABIJ = H.getABIJ();
    const SpinorbitalTensor<U>& VABEF = H.getABCD();
    const SpinorbitalTensor<U>& VMNIJ = H.getIJKL();
    const SpinorbitalTensor<U>& VAMEI = H.getAIBJ();

    ExcitationOperator<U,2>& T = this->template get   <ExcitationOperator<U,2> >(  "T");
    Denominator<U>&          D = this->template gettmp<Denominator<U>          >(  "D");
    ExcitationOperator<U,2>& Z = this->template gettmp<ExcitationOperator<U,2> >(  "Z");

    /**************************************************************************
     *
     * LCCD Iteration
     */
    Z(2)["abij"]  =     VABIJ["abij"];
    Z(2)["abij"] +=       fAE[  "af"]*T(2)["fbij"];
    Z(2)["abij"] -=       fMI[  "ni"]*T(2)["abnj"];
    Z(2)["abij"] += 0.5*VABEF["abef"]*T(2)["efij"];
    Z(2)["abij"] += 0.5*VMNIJ["mnij"]*T(2)["abmn"];
    Z(2)["abij"] +=     VAMEI["amei"]*T(2)["ebjm"];
    /*
     *************************************************************************/

    Z.weight(D);
    T += Z;

    this->energy() = 0.25*real(scalar(H.getABIJ()*T(2)));
    this->conv() = Z.norm(00);

    diis.extrapolate(T, Z);
}

INSTANTIATE_SPECIALIZATIONS(LCCD);
REGISTER_TASK(LCCD<double>,"lccd");
