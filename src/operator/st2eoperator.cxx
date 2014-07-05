/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must reRin the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documenRtion and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANRBILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENRL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DAR, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "st2eoperator.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::tensor;

template <typename U>
STTwoElectronOperator<U,2>::STTwoElectronOperator(const std::string& name, const OneElectronOperator<U>& X,
                                                  const ExcitationOperator<U,2>& T)
: TwoElectronOperator<U>(name, X), T(T)
{
    this->ij["mi"] += this->ia["me"]*T(1)["ei"];

    this->ai["ai"] += T(2)["aeim"]*this->ia["me"];
    this->ai["ai"] += T(1)["ei"]*this->ab["ae"];
    this->ai["ai"] -= T(1)["am"]*this->ij["mi"];

    this->ab["ae"] -= this->ia["me"]*T(1)["am"];

    this->abij["abij"] += this->ab["af"]*T(2)["fbij"];
    this->abij["abij"] -= this->ij["ni"]*T(2)["abnj"];

    this->aijk["amij"] += this->ia["me"]*T(2)["aeij"];

    this->abci["abej"] -= this->ia["me"]*T(2)["abmj"];
}

template <typename U>
STTwoElectronOperator<U,2>::STTwoElectronOperator(const std::string& name, const TwoElectronOperator<U>& X,
                                                  const ExcitationOperator<U,2>& T,
                                                  bool isHbar)
: TwoElectronOperator<U>(name, X), T(T)
{
    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    this->ia["me"] = this->ijab["mnef"]*T(1)["fn"];

    this->ij["mi"] += 0.5*this->ijab["nmef"]*T(2)["efni"];
    this->ij["mi"] += this->ia["me"]*T(1)["ei"];
    this->ij["mi"] += this->ijak["nmfi"]*T(1)["fn"];

    this->ijkl["mnij"] += 0.5*this->ijab["mnef"]*Tau["efij"];
    this->ijkl["mnij"] += this->ijak["nmei"]*T(1)["ej"];

    this->ijak["mnej"] += this->ijab["mnef"]*T(1)["fj"];

    if (!isHbar)
    {
        this->ai["ai"] -= T(1)["em"]*this->aibj["amei"];
        this->ai["ai"] += 0.5*this->aibc["amef"]*Tau["efim"];
        this->ai["ai"] -= 0.5*this->ijak["mnei"]*T(2)["eamn"];
        this->ai["ai"] += T(2)["aeim"]*this->ia["me"];
        this->ai["ai"] += T(1)["ei"]*this->ab["ae"];
        this->ai["ai"] -= T(1)["am"]*this->ij["mi"];
    }
    else
    {
        this->ai = 0;
    }

    this->ab["ae"] -= 0.5*this->ijab["mnfe"]*T(2)["famn"];
    this->ab["ae"] -= this->ia["me"]*T(1)["am"];
    this->ab["ae"] += this->aibc["anef"]*T(1)["fn"];

    this->aijk["amij"] += 0.5*this->aibc["amef"]*Tau["efij"];
    this->aijk["amij"] += this->aibj["amej"]*T(1)["ei"];

    this->aibj["amei"] -= 0.5*this->ijab["mnef"]*T(2)["afin"];
    this->aibj["amei"] -= this->aibc["amfe"]*T(1)["fi"];
    this->aibj["amei"] -= this->ijak["nmei"]*T(1)["an"];

    if (!isHbar)
    {
        this->abij["abij"] += this->ab["af"]*T(2)["fbij"];
        this->abij["abij"] -= this->ij["ni"]*T(2)["abnj"];
        this->abij["abij"] += this->abci["abej"]*T(1)["ei"];
        this->abij["abij"] -= this->aijk["amij"]*T(1)["bm"];
        this->abij["abij"] += 0.5*this->abcd["abef"]*Tau["efij"];
        this->abij["abij"] += 0.5*this->ijkl["mnij"]*Tau["abmn"];
        this->abij["abij"] -= this->aibj["amei"]*T(2)["ebmj"];
    }
    else
    {
        this->abij = 0;
    }

    this->aijk["amij"] += this->ijak["nmej"]*T(2)["aein"];
    this->aijk["amij"] -= this->ijkl["nmij"]*T(1)["an"];
    this->aijk["amij"] += this->ia["me"]*T(2)["aeij"];

    this->aibj["amei"] -= 0.5*this->ijab["mnef"]*T(2)["afin"];
    this->aibj["amei"] += 0.5*this->ijak["nmei"]*T(1)["an"];

    this->abci["abej"] += 0.5*this->ijak["mnej"]*T(2)["abmn"];
    this->abci["abej"] -= this->aibj["amej"]*T(1)["bm"];
    this->abci["abej"] += this->aibc["amef"]*T(2)["fbmj"];
    this->abci["abej"] += this->abcd["abef"]*T(1)["fj"];
    this->abci["abej"] -= this->ia["me"]*T(2)["abmj"];

    this->aibj["amei"] -= 0.5*this->ijak["nmei"]*T(1)["an"];

    this->abcd["abef"] += 0.5*this->ijab["mnef"]*Tau["abmn"];
    this->abcd["abef"] -= this->aibc["amef"]*T(1)["bm"];

    this->aibc["amef"] -= this->ijab["nmef"]*T(1)["an"];
}

template <typename U>
void STTwoElectronOperator<U,2>::contract(const ExcitationOperator<U,2>& R,
                                               ExcitationOperator<U,2>& Z,
                                         bool connected) const
{
   OneElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

   SpinorbitalTensor<U>& IMI = I.getIJ();
   SpinorbitalTensor<U>& IAE = I.getAB();

   IMI["mi"]  = this->ijak["nmei"]*R(1)["en"];
   IMI["mi"] += 0.5*this->ijab["mnef"]*R(2)["efin"];
   IAE["ae"]  = this->aibc["amef"]*R(1)["fm"];
   IAE["ae"] -= 0.5*this->ijab["mnef"]*R(2)["afmn"];
   

   Z(1)["ai"] += this->ab["ae"]*R(1)["ei"];
   Z(1)["ai"] -= this->ij["mi"]*R(1)["am"];
   Z(1)["ai"] -= this->aibj["amei"]*R(1)["em"];
   Z(1)["ai"] += this->ia["me"]*R(2)["aeim"];
   Z(1)["ai"] += 0.5*this->aibc["amef"]*R(2)["efim"];
   Z(1)["ai"] -= 0.5*this->ijak["mnei"]*R(2)["eamn"];

   

   Z(2)["abij"] += this->ab["ae"]*R(2)["ebij"];
   Z(2)["abij"] -= this->ij["mi"]*R(2)["abmj"];
   Z(2)["abij"] += IAE["ae"]*T(2)["ebij"];
   Z(2)["abij"] -= IMI["mi"]*T(2)["abmj"];
   Z(2)["abij"] += this->abci["abej"]*R(1)["ei"];
   Z(2)["abij"] -= this->aijk["amij"]*R(1)["bm"];
   Z(2)["abij"] += 0.5*this->abcd["abef"]*R(2)["efij"];
   Z(2)["abij"] += 0.5*this->ijkl["mnij"]*R(2)["abmn"];
   Z(2)["abij"] -= this->aibj["amei"]*R(2)["ebmj"];

   if (!connected)
   {
       Z(1) += this->ai*R(0);
       Z(2) += this->abij*R(0);

       Z(2)["abij"] += this->ai["ai"]*R(1)["bj"];
   }
}

template <typename U>
void STTwoElectronOperator<U,2>::contractsam(const ExcitationOperator<U,2>& R,
                                                ExcitationOperator<U,2>& Z,
                                          bool connected) const
{
    OneElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

    SpinorbitalTensor<U>& IMI = I.getIJ();
    SpinorbitalTensor<U>& IAE = I.getAB();


    // Z(1)["ai"] += this->ab["ae"]*R(1)["ei"];
    // Z(1)["ai"] -= this->ij["mi"]*R(1)["am"];
    // Z(1)["ai"] -= this->aibj["amei"]*R(1)["em"];
    // Z(1)["ai"] += this->ia["me"]*R(2)["aeim"];
    // Z(1)["ai"] += 0.5*this->aibc["amef"]*R(2)["efim"];
    // Z(1)["ai"] -= 0.5*this->ijak["mnei"]*R(2)["eamn"]; // Orig
    Z(1)["ai"] += this->ab["ad"]*R(1)["di"];
    Z(1)["ai"] -= this->ij["li"]*R(1)["al"];
    Z(1)["ai"] -= this->aibj["aldi"]*R(1)["dl"];
    Z(1)["ai"] += this->ia["ld"]*R(2)["adil"];
    Z(1)["ai"] += 0.5*this->aibc["alde"]*R(2)["deil"];
    Z(1)["ai"] += 0.5*this->ijak["lmdi"]*R(2)["adlm"]; // No Change. 



    // Z(2)["abij"] += this->ab["ae"]*R(2)["ebij"];
    // Z(2)["abij"] -= this->ij["mi"]*R(2)["abmj"]; // Orig
    Z(2)["abij"] += this->ab["bd"]*R(2)["adij"];
    Z(2)["abij"] -= this->ij["lj"]*R(2)["abil"]; // No Change. 


    //Z(2)["abij"] += this->abci["abej"]*R(1)["ei"];
    //Z(2)["abij"] -= this->aijk["amij"]*R(1)["bm"]; // Orig
    Z(2)["abij"] += this->abci["abdj"]*R(1)["di"];
    Z(2)["abij"] += this->aijk["blij"]*R(1)["al"]; // No Change

    // Z(2)["abij"] += 0.5*this->abcd["abef"]*R(2)["efij"];
    // Z(2)["abij"] += 0.5*this->ijkl["mnij"]*R(2)["abmn"]; // Orig
    Z(2)["abij"] += 0.5*this->abcd["abde"]*R(2)["deij"];
    Z(2)["abij"] += 0.5*this->ijkl["lmij"]*R(2)["ablm"]; // No Change

    //Z(2)["abij"] -= this->aibj["amei"]*R(2)["ebmj"]; // Orig
    Z(2)["abij"] -= 2.0*this->aibj["bldj"]*R(2)["adil"]; // No Change
    // With both of these ^^^^ turned on accidentally, I got a lot more real evals. 2.0* same thing.

    TwoElectronOperator<U> V("H", this->arena, this->occ, this->vrt);
    SpinorbitalTensor<U>& VIJAB = V.getIJAB();
    SpinorbitalTensor<U>& VIJAK = V.getIJAK();

    IMI["mi"]  = this->ijak["nmei"]*R(1)["en"];
    IMI["mi"] += 0.5*this->ijab["mnef"]*R(2)["efin"];
    IAE["ae"]  = this->aibc["amef"]*R(1)["fm"];
    IAE["ae"] -= 0.5*this->ijab["mnef"]*R(2)["afmn"]; // Orig

    //IMI["mi"]  = -VIJAK["mldi"]*R(1)["dl"];
    //IMI["mi"] += 0.5*VIJAB["mlde"]*R(2)["deil"];
    ////IAE["ae"]  = this->aibc["amef"]*R(1)["fm"];
    //IAE["ae"]  = -0.5*VIJAB["mlde"]*R(2)["adlm"];


    Z(2)["abij"] -= IMI["mi"]*T(2)["abmj"]; // Orig
    Z(2)["abij"] += IAE["ae"]*T(2)["ebij"]; // Orig
    


    if (!connected)
    {
        Z(1) += this->ai*R(0);
        Z(2) += this->abij*R(0);

        Z(2)["abij"] += this->ai["ai"]*R(1)["bj"];
    }
}

template <typename U>
void STTwoElectronOperator<U,2>::contract(const DeexcitationOperator<U,2>& L,
                                                DeexcitationOperator<U,2>& Z,
                                          bool connected) const
{
    OneElectronOperator<U> I("I", this->arena, this->occ, this->vrt);

    SpinorbitalTensor<U>& IMN = I.getIJ();
    SpinorbitalTensor<U>& IEF = I.getAB();

    SpinorbitalTensor<U> Tau(T(2));
    Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

    IMN["mn"] = 0.5*T(2)["efno"]*L(2)["moef"];
    IEF["ef"] = 0.5*L(2)["mnfg"]*T(2)["egmn"];

    Z(1)["ia"] -= L(1)["me"]*this->aibj["eiam"];
    Z(1)["ia"] += 0.5*this->abci["efam"]*L(2)["imef"];
    Z(1)["ia"] -= 0.5*this->aijk["eimn"]*L(2)["mnea"];
    Z(1)["ia"] += L(1)["ie"]*this->ab["ea"];
    Z(1)["ia"] -= L(1)["ma"]*this->ij["im"];
    Z(1)["ia"] -= IMN["mn"]*this->ijak["inam"];
    Z(1)["ia"] += IEF["ef"]*this->aibc["fiea"];

    Z(2)["ijab"] += this->ab["ea"]*L(2)["ijeb"];
    Z(2)["ijab"] -= this->ij["im"]*L(2)["mjab"];
    Z(2)["ijab"] += this->aibc["ejab"]*L(1)["ie"];
    Z(2)["ijab"] -= this->ijak["jibm"]*L(1)["ma"];
    Z(2)["ijab"] += 0.5*this->abcd["efab"]*L(2)["ijef"];
    Z(2)["ijab"] += 0.5*this->ijkl["ijmn"]*L(2)["mnab"];
    Z(2)["ijab"] -= this->aibj["eiam"]*L(2)["mjeb"];
    Z(2)["ijab"] -= IMN["im"]*this->ijab["mjab"];
    Z(2)["ijab"] -= IEF["ae"]*this->ijab["ijeb"];

    if (!connected)
    {
        Z(1)["ia"] += this->ia["ia"]*L(0)[""];
        Z(2)["ijab"] += this->ijab["ijab"]*L(0)[""];

        Z(1)["ia"] += this->ia["ia"];
        Z(2)["ijab"] += this->ijab["ijab"];

        Z(2)["ijab"] += this->ia["ia"]*L(1)["jb"];
    }
}

INSTANTIATE_SPECIALIZATIONS_2(STTwoElectronOperator,2);
