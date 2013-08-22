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
using namespace aquarius::integrals;
using namespace aquarius::task;

template <typename T>
CholeskyUHF<T>::CholeskyUHF(const string& name, const Config& config)
: UHF<T>("choleskyscf", name, config)
{
    for (vector<Product>::iterator i = this->products.begin();i != this->products.end();++i)
    {
        i->addRequirement(Requirement("cholesky", "cholesky"));
    }
}

template <typename T>
void CholeskyUHF<T>::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = this->template get<Molecule>("molecule");
    const CholeskyIntegrals<T>& chol = this->template get<CholeskyIntegrals<T> >("cholesky");

    int norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumAlphaElectrons();

    vector<int> shapeN = vec(NS);
    vector<int> shapeNNN = vec(NS,NS,NS);
    vector<int> sizer = vec(chol.getRank());
    vector<int> sizenOr = vec(norb,nalpha,chol.getRank());
    vector<int> sizenor = vec(norb,nbeta,chol.getRank());

    this->puttmp("J", new DistTensor<T>(arena, 1, sizer, shapeN, false));
    this->puttmp("JD", new DistTensor<T>(arena, 1, sizer, shapeN, false));
    this->puttmp("La_occ", new DistTensor<T>(arena, 3, sizenOr, shapeNNN, false));
    this->puttmp("Lb_occ", new DistTensor<T>(arena, 3, sizenor, shapeNNN, false));
    this->puttmp("LDa_occ", new DistTensor<T>(arena, 3, sizenOr, shapeNNN, false));
    this->puttmp("LDb_occ", new DistTensor<T>(arena, 3, sizenor, shapeNNN, false));

    UHF<T>::run(dag, arena);
}

template <typename T>
void CholeskyUHF<T>::buildFock()
{
    const CholeskyIntegrals<T>& chol = this->template get<CholeskyIntegrals<T> >("cholesky");

    DistTensor<T>& H = this->template get<DistTensor<T> >("H");
    DistTensor<T>& Da = this->template get<DistTensor<T> >("Da");
    DistTensor<T>& Db = this->template get<DistTensor<T> >("Db");
    DistTensor<T>& Fa = this->template get<DistTensor<T> >("Fa");
    DistTensor<T>& Fb = this->template get<DistTensor<T> >("Fb");

    DistTensor<T>& Ca_occ = this->template gettmp<DistTensor<T> >("Ca_occ");
    DistTensor<T>& Cb_occ = this->template gettmp<DistTensor<T> >("Cb_occ");

    DistTensor<T>& J = this->template gettmp<DistTensor<T> >("J");
    DistTensor<T>& JD = this->template gettmp<DistTensor<T> >("JD");
    DistTensor<T>& La_occ = this->template gettmp<DistTensor<T> >("La_occ");
    DistTensor<T>& Lb_occ = this->template gettmp<DistTensor<T> >("Lb_occ");
    DistTensor<T>& LDa_occ = this->template gettmp<DistTensor<T> >("LDa_occ");
    DistTensor<T>& LDb_occ = this->template gettmp<DistTensor<T> >("LDb_occ");

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
    Da += Db;
    J["J"] = chol.getL()["cdJ"]*Da["cd"];
    JD["J"] = chol.getD()["J"]*J["J"];
    Da -= Db;
    Fa["ab"] = JD["J"]*chol.getL()["abJ"];

    /*
     * Core contribution:
     *
     * F += H
     *
     * Up though this point, Fa = Fb
     */
    Fa += H;
    Fb  = Fa;

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
    La_occ["aiJ"] = chol.getL()["acJ"]*Ca_occ["ci"];
    LDa_occ["aiJ"] = chol.getD()["J"]*La_occ["aiJ"];
    Fa["ab"] -= LDa_occ["aiJ"]*La_occ["biJ"];

    Lb_occ["aiJ"] = chol.getL()["acJ"]*Cb_occ["ci"];
    LDb_occ["aiJ"] = chol.getD()["J"]*Lb_occ["aiJ"];
    Fb["ab"] -= LDb_occ["aiJ"]*Lb_occ["biJ"];
}

INSTANTIATE_SPECIALIZATIONS(CholeskyUHF);
REGISTER_TASK(CholeskyUHF<double>, "choleskyscf");
