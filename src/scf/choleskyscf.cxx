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
 * ARE DISCLAIMED. IN NO EVENT SHALL EDGAR SOLOMONIK BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "choleskyscf.hpp"

#include <algorithm>
#include <limits>

using namespace std;
using namespace libtensor;
using namespace aquarius::input;

namespace aquarius
{
namespace scf
{

CholeskyUHF::CholeskyUHF(DistWorld* dw, const CholeskyIntegrals& chol, const Config& config)
: UHF(dw, chol.getMolecule(), config), chol(chol)
{
    int shapeN[] = {NS};
    int shapeNNN[] = {NS,NS,NS};

    int sizer[] = {chol.getRank()};
    J = new DistTensor(1, sizer, shapeN, dw);
    JD = new DistTensor(1, sizer, shapeN, dw);
    int sizenOr[] = {norb,nalpha,chol.getRank()};
    int sizenor[] = {norb,nbeta,chol.getRank()};
    La_occ = new DistTensor(3, sizenOr, shapeNNN, dw);
    Lb_occ = new DistTensor(3, sizenor, shapeNNN, dw);
    LDa_occ = new DistTensor(3, sizenOr, shapeNNN, dw);
    LDb_occ = new DistTensor(3, sizenor, shapeNNN, dw);
}

CholeskyUHF::~CholeskyUHF()
{
    delete J;
    delete JD;
    delete La_occ;
    delete Lb_occ;
    delete LDa_occ;
    delete LDb_occ;
}

void CholeskyUHF::buildFock()
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
    (*Da)["cd"] += (*Db)["cd"];
    (*J)["J"] = chol.getL()["cdJ"]*(*Da)["cd"];
    (*JD)["J"] = chol.getD()["J"]*(*J)["J"];
    (*Da)["cd"] -= (*Db)["cd"];
    (*Fa)["ab"] = (*JD)["J"]*chol.getL()["abJ"];

    /*
     * Core contribution:
     *
     * F += H
     *
     * Up though this point, Fa = Fb
     */
    (*Fa)["ab"] += (*H)["ab"];
    (*Fb)["ab"] = (*Fa)["ab"];

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
    (*La_occ)["aiJ"] = chol.getL()["acJ"]*(*Ca_occ)["ci"];
    (*LDa_occ)["aiJ"] = chol.getD()["J"]*(*La_occ)["aiJ"];
    (*Fa)["ab"] -= (*LDa_occ)["aiJ"]*(*La_occ)["biJ"];

    (*Lb_occ)["aiJ"] = chol.getL()["acJ"]*(*Cb_occ)["ci"];
    (*LDb_occ)["aiJ"] = chol.getD()["J"]*(*Lb_occ)["aiJ"];
    (*Fb)["ab"] -= (*LDb_occ)["aiJ"]*(*Lb_occ)["biJ"];
}

}
}
