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

#include "choleskyscf.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::scf;
using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::slide;
using namespace aquarius::op;

template <typename T>
CholeskyUHF<T>::CholeskyUHF(const Config& config, const CholeskyIntegrals<T>& chol)
: UHF<T>(chol.arena, config, chol.molecule), chol(chol)
{
    vector<int> shapeN = vec(NS);
    vector<int> shapeNNN = vec(NS,NS,NS);
    vector<int> sizer = vec(chol.getRank());
    vector<int> sizenOr = vec(this->norb,this->nalpha,chol.getRank());
    vector<int> sizenor = vec(this->norb,this->nbeta,chol.getRank());

    J = new DistTensor<T>(this->arena, 1, sizer, shapeN, false);
    JD = new DistTensor<T>(this->arena, 1, sizer, shapeN, false);
    La_occ = new DistTensor<T>(this->arena, 3, sizenOr, shapeNNN, false);
    Lb_occ = new DistTensor<T>(this->arena, 3, sizenor, shapeNNN, false);
    LDa_occ = new DistTensor<T>(this->arena, 3, sizenOr, shapeNNN, false);
    LDb_occ = new DistTensor<T>(this->arena, 3, sizenor, shapeNNN, false);
}

template <typename T>
CholeskyUHF<T>::~CholeskyUHF()
{
    delete J;
    delete JD;
    delete La_occ;
    delete Lb_occ;
    delete LDa_occ;
    delete LDb_occ;
}

template <typename T>
void CholeskyUHF<T>::buildFock()
{
    /*
     * Coulomb contribution:
     *
     * F[ab] = (Da[cd]+Db[cd])*(ab|cd)
     *
     *       = (Da[cd]+Db[cd])*L[abJ]*D[J]*L[cdJ]
     *
     *       = L(abJ)*{D[J]*{(Da[cd]+Db[cd])*L[cdJ]}}
     *
     *       = L(abJ)*{D[J]*J[J]}
     */
    (*this->Da) += (*this->Db);
    (*J)["J"] = chol.getL()["cdJ"]*(*this->Da)["cd"];
    (*JD)["J"] = chol.getD()["J"]*(*J)["J"];
    (*this->Da) -= (*this->Db);
    (*this->Fa)["ab"] = (*JD)["J"]*chol.getL()["abJ"];

    /*
     * Core contribution:
     *
     * F += H
     *
     * Up though this point, Fa = Fb
     */
    (*this->Fa) += (*this->H);
    (*this->Fb) = (*this->Fa);

    /*
     * Exchange contribution:
     *
     * Fa[ab] -= Da[cd]*(ac|bd)
     *
     *         = C[ci]*C[di]*L[acJ]*D[J]*L[bdJ]
     *
     *         = {C[ci]*L[acJ]}*{D[J]*{C[di]*L[bdJ]}}
     *
     *         = L[aiJ]*{D[J]*L[biJ]}
     */
    (*La_occ)["aiJ"] = chol.getL()["acJ"]*(*this->Ca_occ)["ci"];
    (*LDa_occ)["aiJ"] = chol.getD()["J"]*(*La_occ)["aiJ"];
    (*this->Fa)["ab"] -= (*LDa_occ)["aiJ"]*(*La_occ)["biJ"];

    (*Lb_occ)["aiJ"] = chol.getL()["acJ"]*(*this->Cb_occ)["ci"];
    (*LDb_occ)["aiJ"] = chol.getD()["J"]*(*Lb_occ)["aiJ"];
    (*this->Fb)["ab"] -= (*LDb_occ)["aiJ"]*(*Lb_occ)["biJ"];
}

INSTANTIATE_SPECIALIZATIONS(CholeskyUHF);
