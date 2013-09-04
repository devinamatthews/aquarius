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

#include "2eoperator.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::tensor;

template <typename T>
void TwoElectronOperator<T>::initialize()
{
    int nI = this->occ.nalpha;
    int ni = this->occ.nbeta;
    int nA = this->vrt.nalpha;
    int na = this->vrt.nbeta;

    vector<int> sizeIIII = vec(nI, nI, nI, nI);
    vector<int> sizeIiIi = vec(nI, ni, nI, ni);
    vector<int> sizeiiii = vec(ni, ni, ni, ni);
    vector<int> sizeIIIA = vec(nI, nI, nI, nA);
    vector<int> sizeIiIa = vec(nI, ni, nI, na);
    vector<int> sizeiIiA = vec(ni, nI, ni, nA);
    vector<int> sizeiiia = vec(ni, ni, ni, na);
    vector<int> sizeAAII = vec(nA, nA, nI, nI);
    vector<int> sizeAaIi = vec(nA, na, nI, ni);
    vector<int> sizeaaii = vec(na, na, ni, ni);
    vector<int> sizeAIAI = vec(nA, nI, nA, nI);
    vector<int> sizeAiAi = vec(nA, ni, nA, ni);
    vector<int> sizeaIaI = vec(na, nI, na, nI);
    vector<int> sizeaiai = vec(na, ni, na, ni);
    vector<int> sizeAAAI = vec(nA, nA, nA, nI);
    vector<int> sizeAaAi = vec(nA, na, nA, ni);
    vector<int> sizeaAaI = vec(na, nA, na, nI);
    vector<int> sizeaaai = vec(na, na, na, ni);
    vector<int> sizeAAAA = vec(nA, nA, nA, nA);
    vector<int> sizeAaAa = vec(nA, na, nA, na);
    vector<int> sizeaaaa = vec(na, na, na, na);

    vector<int> shapeNNNN = vec(NS, NS, NS, NS);
    vector<int> shapeANNN = vec(AS, NS, NS, NS);
    vector<int> shapeANAN = vec(AS, NS, AS, NS);

    ijkl.addSpinCase(new DistTensor<T>(this->arena, 4, sizeIIII, shapeANAN, true), "IJ,KL", "IJKL");
    ijkl.addSpinCase(new DistTensor<T>(this->arena, 4, sizeIiIi, shapeNNNN, true), "Ij,Kl", "IjKl");
    ijkl.addSpinCase(new DistTensor<T>(this->arena, 4, sizeiiii, shapeANAN, true), "ij,kl", "ijkl");

    iajk.addSpinCase(new DistTensor<T>(this->arena, 4, sizeIIIA, shapeANNN, true), "IA,JK", "JKIA");
    iajk.addSpinCase(new DistTensor<T>(this->arena, 4, sizeIiIa, shapeNNNN, true), "Ia,Jk", "JkIa");
    iajk.addSpinCase(new DistTensor<T>(this->arena, 4, sizeiIiA, shapeNNNN, true), "iA,jK", "jKiA");
    iajk.addSpinCase(new DistTensor<T>(this->arena, 4, sizeiiia, shapeANNN, true), "ia,jk", "jkia");

    if (this->hermitian)
    {
        ijka.addSpinCase(iajk(0), "IJ,KA", "IJKA");
        ijka.addSpinCase(iajk(1), "Ij,Ka", "IjKa");
        ijka.addSpinCase(iajk(2), "iJ,kA", "iJkA");
        ijka.addSpinCase(iajk(3), "ij,ka", "ijka");
    }
    else
    {
        ijka.addSpinCase(new DistTensor<T>(this->arena, 4, sizeIIIA, shapeANNN, true), "IJ,KA", "IJKA");
        ijka.addSpinCase(new DistTensor<T>(this->arena, 4, sizeIiIa, shapeNNNN, true), "Ij,Ka", "IjKa");
        ijka.addSpinCase(new DistTensor<T>(this->arena, 4, sizeiIiA, shapeNNNN, true), "iJ,kA", "iJkA");
        ijka.addSpinCase(new DistTensor<T>(this->arena, 4, sizeiiia, shapeANNN, true), "ij,ka", "ijka");
    }

    abij.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAAII, shapeANAN, true), "AB,IJ", "ABIJ");
    abij.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAaIi, shapeNNNN, true), "Ab,Ij", "AbIj");
    abij.addSpinCase(new DistTensor<T>(this->arena, 4, sizeaaii, shapeANAN, true), "ab,ij", "abij");

    if (this->hermitian)
    {
        ijab.addSpinCase(abij(0), "IJ,AB", "ABIJ");
        ijab.addSpinCase(abij(1), "Ij,Ab", "AbIj");
        ijab.addSpinCase(abij(2), "ij,ab", "abij");
    }
    else
    {
        ijab.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAAII, shapeANAN, true), "IJ,AB", "ABIJ");
        ijab.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAaIi, shapeNNNN, true), "Ij,Ab", "AbIj");
        ijab.addSpinCase(new DistTensor<T>(this->arena, 4, sizeaaii, shapeANAN, true), "ij,ab", "abij");
    }

    aibj.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAIAI, shapeNNNN, true), "AI,BJ", "AIBJ");
    aibj.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAiAi, shapeNNNN, true), "Ai,Bj", "AiBj");
    aibj.addSpinCase(new DistTensor<T>(this->arena, 4, sizeaIaI, shapeNNNN, true), "aI,bJ", "aIbJ");
    aibj.addSpinCase(new DistTensor<T>(this->arena, 4, sizeaiai, shapeNNNN, true), "ai,bj", "aibj");
    aibj.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAaIi, shapeNNNN, true), "Ai,bJ", "AbJi");
    aibj.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAaIi, shapeNNNN, true), "aI,Bj", "BaIj");

    aibc.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAAAI, shapeANNN, true), "AI,BC", "BCAI");
    aibc.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAaAi, shapeNNNN, true), "Ai,Bc", "BcAi");
    aibc.addSpinCase(new DistTensor<T>(this->arena, 4, sizeaAaI, shapeNNNN, true), "aI,bC", "bCaI");
    aibc.addSpinCase(new DistTensor<T>(this->arena, 4, sizeaaai, shapeANNN, true), "ai,bc", "bcai");

    if (this->hermitian)
    {
        abci.addSpinCase(aibc(0), "AB,CI", "ABCI");
        abci.addSpinCase(aibc(1), "Ab,Ci", "AbCi");
        abci.addSpinCase(aibc(2), "aB,cI", "aBcI");
        abci.addSpinCase(aibc(3), "ab,ci", "abci");
    }
    else
    {
        abci.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAAAI, shapeANNN, true), "AB,CI", "ABCI");
        abci.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAaAi, shapeNNNN, true), "Ab,Ci", "AbCi");
        abci.addSpinCase(new DistTensor<T>(this->arena, 4, sizeaAaI, shapeNNNN, true), "aB,cI", "aBcI");
        abci.addSpinCase(new DistTensor<T>(this->arena, 4, sizeaaai, shapeANNN, true), "ab,ci", "abci");
    }

    abcd.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAAAA, shapeANAN, true), "AB,CD", "ABCD");
    abcd.addSpinCase(new DistTensor<T>(this->arena, 4, sizeAaAa, shapeNNNN, true), "Ab,Cd", "AbCd");
    abcd.addSpinCase(new DistTensor<T>(this->arena, 4, sizeaaaa, shapeANAN, true), "ab,cd", "abcd");
}

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(const Arena& arena, const Space& occ, const Space& vrt, bool hermitian)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(arena, occ, vrt, hermitian),
  ijkl(this->addTensor(new SpinorbitalTensor<T>("ij,kl"))),
  iajk(this->addTensor(new SpinorbitalTensor<T>("ia,jk"))),
  ijka(this->addTensor(new SpinorbitalTensor<T>("ij,ka"))),
  abij(this->addTensor(new SpinorbitalTensor<T>("ab,ij"))),
  ijab(this->addTensor(new SpinorbitalTensor<T>("ij,ab"))),
  aibj(this->addTensor(new SpinorbitalTensor<T>("ai,bj"))),
  aibc(this->addTensor(new SpinorbitalTensor<T>("ai,bc"))),
  abci(this->addTensor(new SpinorbitalTensor<T>("ab,ci"))),
  abcd(this->addTensor(new SpinorbitalTensor<T>("ab,cd")))
{
    initialize();
}

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(OneElectronOperator<T>& other, int copy, bool breakhermicity)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(other, copy, breakhermicity),
  ijkl(this->addTensor(new SpinorbitalTensor<T>("ij,kl"))),
  iajk(this->addTensor(new SpinorbitalTensor<T>("ia,jk"))),
  ijka(this->addTensor(new SpinorbitalTensor<T>("ij,ka"))),
  abij(this->addTensor(new SpinorbitalTensor<T>("ab,ij"))),
  ijab(this->addTensor(new SpinorbitalTensor<T>("ij,ab"))),
  aibj(this->addTensor(new SpinorbitalTensor<T>("ai,bj"))),
  aibc(this->addTensor(new SpinorbitalTensor<T>("ai,bc"))),
  abci(this->addTensor(new SpinorbitalTensor<T>("ab,ci"))),
  abcd(this->addTensor(new SpinorbitalTensor<T>("ab,cd")))
{
  initialize();
}

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(const OneElectronOperator<T>& other, bool breakhermicity)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(other),
  ijkl(this->addTensor(new SpinorbitalTensor<T>("ij,kl"))),
  iajk(this->addTensor(new SpinorbitalTensor<T>("ia,jk"))),
  ijka(this->addTensor(new SpinorbitalTensor<T>("ij,ka"))),
  abij(this->addTensor(new SpinorbitalTensor<T>("ab,ij"))),
  ijab(this->addTensor(new SpinorbitalTensor<T>("ij,ab"))),
  aibj(this->addTensor(new SpinorbitalTensor<T>("ai,bj"))),
  aibc(this->addTensor(new SpinorbitalTensor<T>("ai,bc"))),
  abci(this->addTensor(new SpinorbitalTensor<T>("ab,ci"))),
  abcd(this->addTensor(new SpinorbitalTensor<T>("ab,cd")))
{
  initialize();
}

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(TwoElectronOperator<T>& other, int copy, bool breakhermicity)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(other, copy, breakhermicity),
  ijkl(this->addTensor(new SpinorbitalTensor<T>("ij,kl"))),
  iajk(this->addTensor(new SpinorbitalTensor<T>("ia,jk"))),
  ijka(this->addTensor(new SpinorbitalTensor<T>("ij,ka"))),
  abij(this->addTensor(new SpinorbitalTensor<T>("ab,ij"))),
  ijab(this->addTensor(new SpinorbitalTensor<T>("ij,ab"))),
  aibj(this->addTensor(new SpinorbitalTensor<T>("ai,bj"))),
  aibc(this->addTensor(new SpinorbitalTensor<T>("ai,bc"))),
  abci(this->addTensor(new SpinorbitalTensor<T>("ab,ci"))),
  abcd(this->addTensor(new SpinorbitalTensor<T>("ab,cd")))
{
    if (copy&IJKL)
    {
        ijkl.addSpinCase(new DistTensor<T>(other.getIJKL()(0)), "IJ,KL", "IJKL");
        ijkl.addSpinCase(new DistTensor<T>(other.getIJKL()(1)), "Ij,Kl", "IjKl");
        ijkl.addSpinCase(new DistTensor<T>(other.getIJKL()(2)), "ij,kl", "ijkl");
    }
    else
    {
        ijkl.addSpinCase(other.getIJKL()(0), "IJ,KL", "IJKL");
        ijkl.addSpinCase(other.getIJKL()(1), "Ij,Kl", "IjKl");
        ijkl.addSpinCase(other.getIJKL()(2), "ij,kl", "ijkl");
    }

    if (copy&IAJK)
    {
        iajk.addSpinCase(new DistTensor<T>(other.getIAJK()(0)), "IA,JK", "JKIA");
        iajk.addSpinCase(new DistTensor<T>(other.getIAJK()(1)), "Ia,Jk", "JkIa");
        iajk.addSpinCase(new DistTensor<T>(other.getIAJK()(2)), "iA,jK", "jKiA");
        iajk.addSpinCase(new DistTensor<T>(other.getIAJK()(3)), "ia,jk", "jkia");
    }
    else
    {
        iajk.addSpinCase(other.getIAJK()(0), "IA,JK", "JKIA");
        iajk.addSpinCase(other.getIAJK()(1), "Ia,Jk", "JkIa");
        iajk.addSpinCase(other.getIAJK()(2), "iA,jK", "jKiA");
        iajk.addSpinCase(other.getIAJK()(3), "ia,jk", "jkia");
    }

    if (copy&IJKA)
    {
        if (!this->hermitian || !(copy&IAJK))
        {
            ijka.addSpinCase(new DistTensor<T>(other.getIJKA()(0)), "IJ,KA", "IJKA");
            ijka.addSpinCase(new DistTensor<T>(other.getIJKA()(1)), "Ij,Ka", "IjKa");
            ijka.addSpinCase(new DistTensor<T>(other.getIJKA()(2)), "iJ,kA", "iJkA");
            ijka.addSpinCase(new DistTensor<T>(other.getIJKA()(3)), "ij,ka", "ijka");
        }
        else
        {
            ijka.addSpinCase(this->iajk(0), "IJ,KA", "IJKA");
            ijka.addSpinCase(this->iajk(1), "Ij,Ka", "IjKa");
            ijka.addSpinCase(this->iajk(2), "iJ,kA", "iJkA");
            ijka.addSpinCase(this->iajk(3), "ij,ka", "ijka");
        }
    }
    else
    {
        ijka.addSpinCase(other.getIJKA()(0), "IJ,KA", "IJKA");
        ijka.addSpinCase(other.getIJKA()(1), "Ij,Ka", "IjKa");
        ijka.addSpinCase(other.getIJKA()(2), "iJ,kA", "iJkA");
        ijka.addSpinCase(other.getIJKA()(3), "ij,ka", "ijka");
    }

    if (copy&ABIJ)
    {
        abij.addSpinCase(new DistTensor<T>(other.getABIJ()(0)), "AB,IJ", "ABIJ");
        abij.addSpinCase(new DistTensor<T>(other.getABIJ()(1)), "Ab,Ij", "AbIj");
        abij.addSpinCase(new DistTensor<T>(other.getABIJ()(2)), "ab,ij", "abij");
    }
    else
    {
        abij.addSpinCase(other.getABIJ()(0), "AB,IJ", "ABIJ");
        abij.addSpinCase(other.getABIJ()(1), "Ab,Ij", "AbIj");
        abij.addSpinCase(other.getABIJ()(2), "ab,ij", "abij");
    }

    if (copy&IJAB)
    {
        if (!this->hermitian || !(copy&ABIJ))
        {
            ijab.addSpinCase(new DistTensor<T>(other.getIJAB()(0)), "IJ,AB", "ABIJ");
            ijab.addSpinCase(new DistTensor<T>(other.getIJAB()(1)), "Ij,Ab", "AbIj");
            ijab.addSpinCase(new DistTensor<T>(other.getIJAB()(2)), "ij,ab", "abij");
        }
        else
        {
            ijab.addSpinCase(this->abij(0), "IJ,AB", "ABIJ");
            ijab.addSpinCase(this->abij(1), "Ij,Ab", "AbIj");
            ijab.addSpinCase(this->abij(2), "ij,ab", "abij");
        }
    }
    else
    {
        ijab.addSpinCase(other.getIJAB()(0), "IJ,AB", "ABIJ");
        ijab.addSpinCase(other.getIJAB()(1), "Ij,Ab", "AbIj");
        ijab.addSpinCase(other.getIJAB()(2), "ij,ab", "abij");
    }

    if (copy&AIBJ)
    {
        aibj.addSpinCase(new DistTensor<T>(other.getAIBJ()(0)), "AI,BJ", "AIBJ");
        aibj.addSpinCase(new DistTensor<T>(other.getAIBJ()(1)), "Ai,Bj", "AiBj");
        aibj.addSpinCase(new DistTensor<T>(other.getAIBJ()(2)), "aI,bJ", "aIbJ");
        aibj.addSpinCase(new DistTensor<T>(other.getAIBJ()(3)), "ai,bj", "aibj");
        aibj.addSpinCase(new DistTensor<T>(other.getABIJ()(1)), "aI,Bj", "BaIj", -1.0);
        aibj.addSpinCase(new DistTensor<T>(other.getABIJ()(1)), "Ai,bJ", "AbJi", -1.0);
    }
    else
    {
        aibj.addSpinCase(other.getAIBJ()(0), "AI,BJ", "AIBJ");
        aibj.addSpinCase(other.getAIBJ()(1), "Ai,Bj", "AiBj");
        aibj.addSpinCase(other.getAIBJ()(2), "aI,bJ", "aIbJ");
        aibj.addSpinCase(other.getAIBJ()(3), "ai,bj", "aibj");
        aibj.addSpinCase(other.getABIJ()(1), "aI,Bj", "BaIj", -1.0);
        aibj.addSpinCase(other.getABIJ()(1), "Ai,bJ", "AbJi", -1.0);
    }

    if (copy&AIBC)
    {
        aibc.addSpinCase(new DistTensor<T>(other.getAIBC()(0)), "AI,BC", "BCAI");
        aibc.addSpinCase(new DistTensor<T>(other.getAIBC()(1)), "Ai,Bc", "BcAi");
        aibc.addSpinCase(new DistTensor<T>(other.getAIBC()(2)), "aI,bC", "bCaI");
        aibc.addSpinCase(new DistTensor<T>(other.getAIBC()(3)), "ai,bc", "bcai");
    }
    else
    {
        aibc.addSpinCase(other.getAIBC()(0), "AI,BC", "BCAI");
        aibc.addSpinCase(other.getAIBC()(1), "Ai,Bc", "BcAi");
        aibc.addSpinCase(other.getAIBC()(2), "aI,bC", "bCaI");
        aibc.addSpinCase(other.getAIBC()(3), "ai,bc", "bcai");
    }

    if (copy&ABCI)
    {
        if (!this->hermitian || !(copy&AIBC))
        {
            abci.addSpinCase(new DistTensor<T>(other.getABCI()(0)), "AB,CI", "ABCI");
            abci.addSpinCase(new DistTensor<T>(other.getABCI()(1)), "Ab,Ci", "AbCi");
            abci.addSpinCase(new DistTensor<T>(other.getABCI()(2)), "aB,cI", "aBcI");
            abci.addSpinCase(new DistTensor<T>(other.getABCI()(3)), "ab,ci", "abci");
        }
        else
        {
            abci.addSpinCase(this->aibc(0), "AB,CI", "ABCI");
            abci.addSpinCase(this->aibc(1), "Ab,Ci", "AbCi");
            abci.addSpinCase(this->aibc(2), "aB,cI", "aBcI");
            abci.addSpinCase(this->aibc(3), "ab,ci", "abci");
        }
    }
    else
    {
        abci.addSpinCase(other.getABCI()(0), "AB,CI", "ABCI");
        abci.addSpinCase(other.getABCI()(1), "Ab,Ci", "AbCi");
        abci.addSpinCase(other.getABCI()(2), "aB,cI", "aBcI");
        abci.addSpinCase(other.getABCI()(3), "ab,ci", "abci");
    }

    if (copy&ABCD)
    {
        abcd.addSpinCase(new DistTensor<T>(other.getABCD()(0)), "AB,CD", "ABCD");
        abcd.addSpinCase(new DistTensor<T>(other.getABCD()(1)), "Ab,Cd", "AbCd");
        abcd.addSpinCase(new DistTensor<T>(other.getABCD()(2)), "ab,cd", "abcd");
    }
    else
    {
        abcd.addSpinCase(other.getABCD()(0), "AB,CD", "ABCD");
        abcd.addSpinCase(other.getABCD()(1), "Ab,Cd", "AbCd");
        abcd.addSpinCase(other.getABCD()(2), "ab,cd", "abcd");
    }
}

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(const TwoElectronOperator<T>& other, bool breakhermicity)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(other, breakhermicity),
  ijkl(this->addTensor(new SpinorbitalTensor<T>("ij,kl"))),
  iajk(this->addTensor(new SpinorbitalTensor<T>("ia,jk"))),
  ijka(this->addTensor(new SpinorbitalTensor<T>("ij,ka"))),
  abij(this->addTensor(new SpinorbitalTensor<T>("ab,ij"))),
  ijab(this->addTensor(new SpinorbitalTensor<T>("ij,ab"))),
  aibj(this->addTensor(new SpinorbitalTensor<T>("ai,bj"))),
  aibc(this->addTensor(new SpinorbitalTensor<T>("ai,bc"))),
  abci(this->addTensor(new SpinorbitalTensor<T>("ab,ci"))),
  abcd(this->addTensor(new SpinorbitalTensor<T>("ab,cd")))
{
    ijkl.addSpinCase(new DistTensor<T>(other.getIJKL()(0)), "IJ,KL", "IJKL");
    ijkl.addSpinCase(new DistTensor<T>(other.getIJKL()(1)), "Ij,Kl", "IjKl");
    ijkl.addSpinCase(new DistTensor<T>(other.getIJKL()(2)), "ij,kl", "ijkl");

    iajk.addSpinCase(new DistTensor<T>(other.getIAJK()(0)), "IA,JK", "JKIA");
    iajk.addSpinCase(new DistTensor<T>(other.getIAJK()(1)), "Ia,Jk", "JkIa");
    iajk.addSpinCase(new DistTensor<T>(other.getIAJK()(2)), "iA,jK", "jKiA");
    iajk.addSpinCase(new DistTensor<T>(other.getIAJK()(3)), "ia,jk", "jkia");

    if (this->hermitian)
    {
        ijka.addSpinCase(this->iajk(0), "IJ,KA", "IJKA");
        ijka.addSpinCase(this->iajk(1), "Ij,Ka", "IjKa");
        ijka.addSpinCase(this->iajk(2), "iJ,kA", "iJkA");
        ijka.addSpinCase(this->iajk(3), "ij,ka", "ijka");
    }
    else
    {
        ijka.addSpinCase(new DistTensor<T>(other.getIJKA()(0)), "IJ,KA", "IJKA");
        ijka.addSpinCase(new DistTensor<T>(other.getIJKA()(1)), "Ij,Ka", "IjKa");
        ijka.addSpinCase(new DistTensor<T>(other.getIJKA()(2)), "iJ,kA", "iJkA");
        ijka.addSpinCase(new DistTensor<T>(other.getIJKA()(3)), "ij,ka", "ijka");
    }

    abij.addSpinCase(new DistTensor<T>(other.getABIJ()(0)), "AB,IJ", "ABIJ");
    abij.addSpinCase(new DistTensor<T>(other.getABIJ()(1)), "Ab,Ij", "AbIj");
    abij.addSpinCase(new DistTensor<T>(other.getABIJ()(2)), "ab,ij", "abij");

    if (this->hermitian)
    {
        ijab.addSpinCase(this->abij(0), "IJ,AB", "ABIJ");
        ijab.addSpinCase(this->abij(1), "Ij,Ab", "AbIj");
        ijab.addSpinCase(this->abij(2), "ij,ab", "abij");
    }
    else
    {
        ijab.addSpinCase(new DistTensor<T>(other.getIJAB()(0)), "IJ,AB", "ABIJ");
        ijab.addSpinCase(new DistTensor<T>(other.getIJAB()(1)), "Ij,Ab", "AbIj");
        ijab.addSpinCase(new DistTensor<T>(other.getIJAB()(2)), "ij,ab", "abij");
    }

    aibj.addSpinCase(new DistTensor<T>(other.getAIBJ()(0)), "AI,BJ", "AIBJ");
    aibj.addSpinCase(new DistTensor<T>(other.getAIBJ()(1)), "Ai,Bj", "AiBj");
    aibj.addSpinCase(new DistTensor<T>(other.getAIBJ()(2)), "aI,bJ", "aIbJ");
    aibj.addSpinCase(new DistTensor<T>(other.getAIBJ()(3)), "ai,bj", "aibj");
    aibj.addSpinCase(new DistTensor<T>(other.getABIJ()(1)), "aI,Bj", "BaIj", -1.0);
    aibj.addSpinCase(new DistTensor<T>(other.getABIJ()(1)), "Ai,bJ", "AbJi", -1.0);

    aibc.addSpinCase(new DistTensor<T>(other.getAIBC()(0)), "AI,BC", "BCAI");
    aibc.addSpinCase(new DistTensor<T>(other.getAIBC()(1)), "Ai,Bc", "BcAi");
    aibc.addSpinCase(new DistTensor<T>(other.getAIBC()(2)), "aI,bC", "bCaI");
    aibc.addSpinCase(new DistTensor<T>(other.getAIBC()(3)), "ai,bc", "bcai");

    if (this->hermitian)
    {
        abci.addSpinCase(this->aibc(0), "AB,CI", "ABCI");
        abci.addSpinCase(this->aibc(1), "Ab,Ci", "AbCi");
        abci.addSpinCase(this->aibc(2), "aB,cI", "aBcI");
        abci.addSpinCase(this->aibc(3), "ab,ci", "abci");
    }
    else
    {
        abci.addSpinCase(new DistTensor<T>(other.getABCI()(0)), "AB,CI", "ABCI");
        abci.addSpinCase(new DistTensor<T>(other.getABCI()(1)), "Ab,Ci", "AbCi");
        abci.addSpinCase(new DistTensor<T>(other.getABCI()(2)), "aB,cI", "aBcI");
        abci.addSpinCase(new DistTensor<T>(other.getABCI()(3)), "ab,ci", "abci");
    }

    abcd.addSpinCase(new DistTensor<T>(other.getABCD()(0)), "AB,CD", "ABCD");
    abcd.addSpinCase(new DistTensor<T>(other.getABCD()(1)), "Ab,Cd", "AbCd");
    abcd.addSpinCase(new DistTensor<T>(other.getABCD()(2)), "ab,cd", "abcd");
}

template <typename T>
T TwoElectronOperator<T>::dot(bool conja, const TwoElectronOperator<T>& A, bool conjb) const
{
    T sum = 0;

    sum += this->ab.dot(conja, A.ab, conjb);
    sum += this->ai.dot(conja, A.ai, conjb);
    sum += this->ia.dot(conja, A.ia, conjb);
    sum += this->ij.dot(conja, A.ij, conjb);

    sum += (T)0.25*ijkl.dot(conja, A.ijkl, conjb);
    sum += (T)0.25*abcd.dot(conja, A.abcd, conjb);
    sum += (T)0.25*abij.dot(conja, A.abij, conjb);
    sum += (T)0.25*ijab.dot(conja, A.ijab, conjb);
    sum +=  (T)0.5*abci.dot(conja, A.abci, conjb);
    sum +=  (T)0.5*aibc.dot(conja, A.aibc, conjb);
    sum +=  (T)0.5*ijka.dot(conja, A.ijka, conjb);
    sum +=  (T)0.5*iajk.dot(conja, A.iajk, conjb);
    sum +=         aibj.dot(conja, A.aibj, conjb);

    return sum;
}

INSTANTIATE_SPECIALIZATIONS(TwoElectronOperator);
