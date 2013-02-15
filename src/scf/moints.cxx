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

#include "moints.hpp"

using namespace std;

namespace aquarius
{
namespace scf
{

MOIntegrals::MOIntegrals(DistWorld *dw, const UHF& uhf)
: uhf(uhf),
  fab("a,b"),
  fai("a,i"),
  fij("i,j"),
  ijkl("ij,kl"),
  ijka("ij,ka"),
  abij("ab,ij"),
  aibj("ai,bj"),
  abci("ab,ci"),
  abcd("ab,cd"),
  comm(dw->comm),
  dw(dw)
{
    int N = uhf.getMolecule().getNumOrbitals();
    int nI = uhf.getMolecule().getNumAlphaElectrons();
    int ni = uhf.getMolecule().getNumBetaElectrons();
    int nA = N-nI;
    int na = N-ni;

    int sizeAA[] = {nA, nA};
    int sizeaa[] = {na, na};
    int sizeAI[] = {nA, nI};
    int sizeai[] = {na, ni};
    int sizeII[] = {nI, nI};
    int sizeii[] = {ni, ni};

    int sizeIIII[] = {nI, nI, nI, nI};
    int sizeIiIi[] = {nI, ni, nI, ni};
    int sizeiiii[] = {ni, ni, ni, ni};
    int sizeIIIA[] = {nI, nI, nI, nA};
    int sizeIiIa[] = {nI, ni, nI, na};
    int sizeiIiA[] = {ni, nI, ni, nA};
    int sizeiiia[] = {ni, ni, ni, na};
    int sizeAAII[] = {nA, nA, nI, nI};
    int sizeAaIi[] = {nA, na, nI, ni};
    int sizeaaii[] = {na, na, ni, ni};
    int sizeAIAI[] = {nA, nI, nA, nI};
    int sizeAiAi[] = {nA, ni, nA, ni};
    int sizeaIaI[] = {na, nI, na, nI};
    int sizeaiai[] = {na, ni, na, ni};
    int sizeAAAI[] = {nA, nA, nA, nI};
    int sizeAaAi[] = {nA, na, nA, ni};
    int sizeaAaI[] = {na, nA, na, nI};
    int sizeaaai[] = {na, na, na, ni};
    int sizeAAAA[] = {nA, nA, nA, nA};
    int sizeAaAa[] = {nA, na, nA, na};
    int sizeaaaa[] = {na, na, na, na};

    int shapeNN[] = {NS, NS};
    int shapeNNNN[] = {NS, NS, NS, NS};
    int shapeANNN[] = {AS, NS, NS, NS};
    int shapeANAN[] = {AS, NS, AS, NS};

    fAB_ = new DistTensor(2, sizeAA, shapeNN, dw, true);
    fab_ = new DistTensor(2, sizeaa, shapeNN, dw, true);
    fAI_ = new DistTensor(2, sizeAI, shapeNN, dw, true);
    fai_ = new DistTensor(2, sizeai, shapeNN, dw, true);
    fIJ_ = new DistTensor(2, sizeII, shapeNN, dw, true);
    fij_ = new DistTensor(2, sizeii, shapeNN, dw, true);
    IJKL_ = new DistTensor(4, sizeIIII, shapeANAN, dw, false);
    IjKl_ = new DistTensor(4, sizeIiIi, shapeNNNN, dw, false);
    ijkl_ = new DistTensor(4, sizeiiii, shapeANAN, dw, false);
    IJKA_ = new DistTensor(4, sizeIIIA, shapeANNN, dw, false);
    IjKa_ = new DistTensor(4, sizeIiIa, shapeNNNN, dw, false);
    iJkA_ = new DistTensor(4, sizeiIiA, shapeNNNN, dw, false);
    ijka_ = new DistTensor(4, sizeiiia, shapeANNN, dw, false);
    ABIJ_ = new DistTensor(4, sizeAAII, shapeANAN, dw, false);
    AbIj_ = new DistTensor(4, sizeAaIi, shapeNNNN, dw, false);
    abij_ = new DistTensor(4, sizeaaii, shapeANAN, dw, false);
    AIBJ_ = new DistTensor(4, sizeAIAI, shapeNNNN, dw, false);
    AiBj_ = new DistTensor(4, sizeAiAi, shapeNNNN, dw, false);
    aIbJ_ = new DistTensor(4, sizeaIaI, shapeNNNN, dw, false);
    aibj_ = new DistTensor(4, sizeaiai, shapeNNNN, dw, false);
    ABCI_ = new DistTensor(4, sizeAAAI, shapeANNN, dw, false);
    AbCi_ = new DistTensor(4, sizeAaAi, shapeNNNN, dw, false);
    aBcI_ = new DistTensor(4, sizeaAaI, shapeNNNN, dw, false);
    abci_ = new DistTensor(4, sizeaaai, shapeANNN, dw, false);
    ABCD_ = new DistTensor(4, sizeAAAA, shapeANAN, dw, false);
    AbCd_ = new DistTensor(4, sizeAaAa, shapeNNNN, dw, false);
    abcd_ = new DistTensor(4, sizeaaaa, shapeANAN, dw, false);

    fab.addSpinCase(fAB_, "A,B", "AB");
    fab.addSpinCase(fab_, "a,b", "ab");

    fai.addSpinCase(fAI_, "A,I", "AI");
    fai.addSpinCase(fai_, "a,i", "ai");

    fij.addSpinCase(fIJ_, "I,J", "IJ");
    fij.addSpinCase(fij_, "i,j", "ij");

    ijkl.addSpinCase(IJKL_, "IJ,KL", "IJKL");
    ijkl.addSpinCase(IjKl_, "Ij,Kl", "IjKl");
    ijkl.addSpinCase(ijkl_, "ij,kl", "ijkl");

    ijka.addSpinCase(IJKA_, "IJ,KA", "IJKA");
    ijka.addSpinCase(IjKa_, "Ij,Ka", "IjKa");
    ijka.addSpinCase(iJkA_, "iJ,kA", "iJkA");
    ijka.addSpinCase(ijka_, "ij,ka", "ijka");

    abij.addSpinCase(ABIJ_, "AB,IJ", "ABIJ");
    abij.addSpinCase(AbIj_, "Ab,Ij", "AbIj");
    abij.addSpinCase(abij_, "ab,ij", "abij");

    aibj.addSpinCase(AIBJ_, "AI,BJ", "AIBJ");
    aibj.addSpinCase(AiBj_, "Ai,Bj", "AiBj");
    aibj.addSpinCase(aIbJ_, "aI,bJ", "aIbJ");
    aibj.addSpinCase(aibj_, "ai,bj", "aibj");
    aibj.addSpinCase(*AbIj_, "Ai,bJ", "AbJi", -1.0);
    aibj.addSpinCase(*AbIj_, "aI,Bj", "BaIj", -1.0);

    abci.addSpinCase(ABCI_, "AB,CI", "ABCI");
    abci.addSpinCase(AbCi_, "Ab,Ci", "AbCi");
    abci.addSpinCase(aBcI_, "aB,cI", "aBcI");
    abci.addSpinCase(abci_, "ab,ci", "abci");

    abcd.addSpinCase(ABCD_, "AB,CD", "ABCD");
    abcd.addSpinCase(AbCd_, "Ab,Cd", "AbCd");
    abcd.addSpinCase(abcd_, "ab,cd", "abcd");

    const double* ea = uhf.getAlphaEigenvalues();
    const double* eb = uhf.getBetaEigenvalues();
    int rank, np;
    MPI_Comm_rank(dw->comm, &rank);
    MPI_Comm_size(dw->comm, &np);

    {
        vector<kv_pair> pairs;

        for (int i = 0;i < nI;i++)
        {
            if (i%np == rank) pairs.push_back(kv_pair(i+i*nI, ea[i]));
        }

        fIJ_->writeRemoteData(pairs.size(), pairs.data());
    }

    {
        vector<kv_pair> pairs;

        for (int i = 0;i < ni;i++)
        {
            if (i%np == rank) pairs.push_back(kv_pair(i+i*ni, eb[i]));
        }

        fij_->writeRemoteData(pairs.size(), pairs.data());
    }

    {
        vector<kv_pair> pairs;

        for (int i = 0;i < nA;i++)
        {
            if (i%np == rank) pairs.push_back(kv_pair(i+i*nA, ea[i+nI]));
        }

        fAB_->writeRemoteData(pairs.size(), pairs.data());
    }

    {
        vector<kv_pair> pairs;

        for (int i = 0;i < na;i++)
        {
            if (i%np == rank) pairs.push_back(kv_pair(i+i*na, eb[i+ni]));
        }

        fab_->writeRemoteData(pairs.size(), pairs.data());
    }
}

}
}
