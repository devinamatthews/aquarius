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
using namespace aquarius::symmetry;

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

    const PointGroup& group = molecule.getGroup();

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumAlphaElectrons();

    vector<int> shapeN{NS};
    vector<int> shapeNNN{NS,NS,NS};
    vector<vector<int> > sizer{{chol.getRank()}};
    vector<vector<int> > sizenOr{norb,{nalpha},{chol.getRank()}};
    vector<vector<int> > sizenor{norb,{nbeta},{chol.getRank()}};

    this->puttmp("J", new SymmetryBlockedTensor<T>("J", arena, group, 1, sizer, shapeN, false));
    this->puttmp("JD", new SymmetryBlockedTensor<T>("JD", arena, group, 1, sizer, shapeN, false));
    this->puttmp("La_occ", new SymmetryBlockedTensor<T>("LpI", arena, group, 3, sizenOr, shapeNNN, false));
    this->puttmp("Lb_occ", new SymmetryBlockedTensor<T>("Lpi", arena, group, 3, sizenor, shapeNNN, false));
    this->puttmp("LDa_occ", new SymmetryBlockedTensor<T>("LDpI", arena, group, 3, sizenOr, shapeNNN, false));
    this->puttmp("LDb_occ", new SymmetryBlockedTensor<T>("LDpi", arena, group, 3, sizenor, shapeNNN, false));

    UHF<T>::run(dag, arena);
}

template <typename T>
void CholeskyUHF<T>::buildFock()
{
    const Molecule& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumAlphaElectrons();

    const CholeskyIntegrals<T>& chol = this->template get<CholeskyIntegrals<T> >("cholesky");

    SymmetryBlockedTensor<T>& H = this->template get<SymmetryBlockedTensor<T> >("H");
    SymmetryBlockedTensor<T>& Da = this->template get<SymmetryBlockedTensor<T> >("Da");
    SymmetryBlockedTensor<T>& Db = this->template get<SymmetryBlockedTensor<T> >("Db");
    SymmetryBlockedTensor<T>& Fa = this->template get<SymmetryBlockedTensor<T> >("Fa");
    SymmetryBlockedTensor<T>& Fb = this->template get<SymmetryBlockedTensor<T> >("Fb");

    SymmetryBlockedTensor<T> Ca_occ("CI", this->template gettmp<SymmetryBlockedTensor<T> >("Ca"),
                                    {{0},{0}}, {norb,{nalpha}});
    SymmetryBlockedTensor<T> Cb_occ("Ci", this->template gettmp<SymmetryBlockedTensor<T> >("Cb"),
                                    {{0},{0}}, {norb,{nbeta}});

    SymmetryBlockedTensor<T>& J = this->template gettmp<SymmetryBlockedTensor<T> >("J");
    SymmetryBlockedTensor<T>& JD = this->template gettmp<SymmetryBlockedTensor<T> >("JD");
    SymmetryBlockedTensor<T>& La_occ = this->template gettmp<SymmetryBlockedTensor<T> >("La_occ");
    SymmetryBlockedTensor<T>& Lb_occ = this->template gettmp<SymmetryBlockedTensor<T> >("Lb_occ");
    SymmetryBlockedTensor<T>& LDa_occ = this->template gettmp<SymmetryBlockedTensor<T> >("LDa_occ");
    SymmetryBlockedTensor<T>& LDb_occ = this->template gettmp<SymmetryBlockedTensor<T> >("LDb_occ");

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
