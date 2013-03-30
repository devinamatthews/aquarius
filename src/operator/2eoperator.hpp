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

#ifndef _AQUARIUS_OPERATOR_2EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_2EOPERATOR_HPP_

#include "1eoperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class TwoElectronOperator : public OneElectronOperator<T>
{
    protected:
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& ijkl;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& iajk;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& ijka;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& abij;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& ijab;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& aibj;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& aibc;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& abci;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& abcd;

    public:
        enum
        {
            IJKL = 0x0010,
            IAJK = 0x0020,
            IJKA = 0x0040,
            ABIJ = 0x0080,
            IJAB = 0x0100,
            AIBJ = 0x0200,
            AIBC = 0x0400,
            ABCI = 0x0800,
            ABCD = 0x1000
        };

        TwoElectronOperator(const scf::UHF<T>& uhf, const bool hermitian=true)
        : OneElectronOperator<T>(uhf, hermitian),
          ijkl(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ij,kl"))),
          iajk(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ia,jk"))),
          ijka(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ij,ka"))),
          abij(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ab,ij"))),
          ijab(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ij,ab"))),
          aibj(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ai,bj"))),
          aibc(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ai,bc"))),
          abci(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ab,ci"))),
          abcd(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ab,cd")))
        {
            int N = uhf.getMolecule().getNumOrbitals();
            int nI = uhf.getMolecule().getNumAlphaElectrons();
            int ni = uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;

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

            int shapeNNNN[] = {NS, NS, NS, NS};
            int shapeANNN[] = {AS, NS, NS, NS};
            int shapeANAN[] = {AS, NS, AS, NS};

            ijkl.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeIIII, shapeANAN, true), "IJ,KL", "IJKL");
            ijkl.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeIiIi, shapeNNNN, true), "Ij,Kl", "IjKl");
            ijkl.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeiiii, shapeANAN, true), "ij,kl", "ijkl");

            iajk.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeIIIA, shapeANNN, true), "IA,JK", "JKIA");
            iajk.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeIiIa, shapeNNNN, true), "Ia,Jk", "JkIa");
            iajk.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeiIiA, shapeNNNN, true), "iA,jK", "jKiA");
            iajk.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeiiia, shapeANNN, true), "ia,jk", "jkia");

            if (hermitian)
            {
                ijka.addSpinCase(iajk(0), "IJ,KA", "IJKA");
                ijka.addSpinCase(iajk(1), "Ij,Ka", "IjKa");
                ijka.addSpinCase(iajk(2), "iJ,kA", "iJkA");
                ijka.addSpinCase(iajk(3), "ij,ka", "ijka");
            }
            else
            {
                ijka.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeIIIA, shapeANNN, true), "IJ,KA", "IJKA");
                ijka.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeIiIa, shapeNNNN, true), "Ij,Ka", "IjKa");
                ijka.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeiIiA, shapeNNNN, true), "iJ,kA", "iJkA");
                ijka.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeiiia, shapeANNN, true), "ij,ka", "ijka");
            }

            abij.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAAII, shapeANAN, true), "AB,IJ", "ABIJ");
            abij.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaIi, shapeNNNN, true), "Ab,Ij", "AbIj");
            abij.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaaii, shapeANAN, true), "ab,ij", "abij");

            if (hermitian)
            {
                ijab.addSpinCase(abij(0), "IJ,AB", "ABIJ");
                ijab.addSpinCase(abij(1), "Ij,Ab", "AbIj");
                ijab.addSpinCase(abij(2), "ij,ab", "abij");
            }
            else
            {
                ijab.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAAII, shapeANAN, true), "IJ,AB", "ABIJ");
                ijab.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaIi, shapeNNNN, true), "Ij,Ab", "AbIj");
                ijab.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaaii, shapeANAN, true), "ij,ab", "abij");
            }

            aibj.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAIAI, shapeNNNN, true), "AI,BJ", "AIBJ");
            aibj.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAiAi, shapeNNNN, true), "Ai,Bj", "AiBj");
            aibj.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaIaI, shapeNNNN, true), "aI,bJ", "aIbJ");
            aibj.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaiai, shapeNNNN, true), "ai,bj", "aibj");
            aibj.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaIi, shapeNNNN, true), "Ai,bJ", "AbJi");

            if (hermitian)
            {
                aibj.addSpinCase(aibj(4), "aI,Bj", "BaIj");
            }
            else
            {
                aibj.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaIi, shapeNNNN, true), "aI,Bj", "BaIj");
            }

            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAAAI, shapeANNN, true), "AI,BC", "BCAI");
            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaAi, shapeNNNN, true), "Ai,Bc", "BcAi");
            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaAaI, shapeNNNN, true), "aI,bC", "bCaI");
            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaaai, shapeANNN, true), "ai,bc", "bcai");

            if (hermitian)
            {
                abci.addSpinCase(aibc(0), "AB,CI", "ABCI");
                abci.addSpinCase(aibc(1), "Ab,Ci", "AbCi");
                abci.addSpinCase(aibc(2), "aB,cI", "aBcI");
                abci.addSpinCase(aibc(3), "ab,ci", "abci");
            }
            else
            {
                abci.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAAAI, shapeANNN, true), "AB,CI", "ABCI");
                abci.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaAi, shapeNNNN, true), "Ab,Ci", "AbCi");
                abci.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaAaI, shapeNNNN, true), "aB,cI", "aBcI");
                abci.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaaai, shapeANNN, true), "ab,ci", "abci");
            }

            abcd.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAAAA, shapeANAN, true), "AB,CD", "ABCD");
            abcd.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaAa, shapeNNNN, true), "Ab,Cd", "AbCd");
            abcd.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaaaa, shapeANAN, true), "ab,cd", "abcd");
        }

        TwoElectronOperator(OneElectronOperator<T>& other, int copy)
        : OneElectronOperator<T>(other, copy) {}

        TwoElectronOperator(const OneElectronOperator<T>& other)
        : OneElectronOperator<T>(other) {}

        TwoElectronOperator(TwoElectronOperator<T>& other, int copy)
        : OneElectronOperator<T>(other, copy),
          ijkl(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ij,kl"))),
          iajk(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ia,jk"))),
          ijka(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ij,ka"))),
          abij(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ab,ij"))),
          ijab(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ij,ab"))),
          aibj(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ai,bj"))),
          aibc(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ai,bc"))),
          abci(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ab,ci"))),
          abcd(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ab,cd")))
        {
            if (copy&IJKL)
            {
                ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL()(0)), "IJ,KL", "IJKL");
                ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL()(1)), "Ij,Kl", "IjKl");
                ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL()(2)), "ij,kl", "ijkl");
            }
            else
            {
                ijkl.addSpinCase(other.getIJKL()(0), "IJ,KL", "IJKL");
                ijkl.addSpinCase(other.getIJKL()(1), "Ij,Kl", "IjKl");
                ijkl.addSpinCase(other.getIJKL()(2), "ij,kl", "ijkl");
            }

            if (copy&IAJK)
            {
                iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(0)), "IA,JK", "JKIA");
                iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(1)), "Ia,Jk", "JkIa");
                iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(2)), "iA,jK", "jKiA");
                iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(3)), "ia,jk", "jkia");
            }
            else
            {
                iajk.addSpinCase(other.getIAJK()(0), "IA,JK", "JKIA");
                iajk.addSpinCase(other.getIAJK()(1), "Ia,Jk", "JkIa");
                iajk.addSpinCase(other.getIAJK()(2), "iA,jK", "jKiA");
                iajk.addSpinCase(other.getIAJK()(3), "ia,jk", "jkia");
            }

            if (copy&IAJK)
            {
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(0)), "IJ,KA", "IJKA");
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(1)), "Ij,Ka", "IjKa");
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(2)), "iJ,kA", "iJkA");
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(3)), "ij,ka", "ijka");
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
                abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(0)), "AB,IJ", "ABIJ");
                abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(1)), "Ab,Ij", "AbIj");
                abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(2)), "ab,ij", "abij");
            }
            else
            {
                abij.addSpinCase(other.getABIJ()(0), "AB,IJ", "ABIJ");
                abij.addSpinCase(other.getABIJ()(1), "Ab,Ij", "AbIj");
                abij.addSpinCase(other.getABIJ()(2), "ab,ij", "abij");
            }

            if (copy&IJAB)
            {
                ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB()(0)), "IJ,AB", "ABIJ");
                ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB()(1)), "Ij,Ab", "AbIj");
                ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB()(2)), "ij,ab", "abij");
            }
            else
            {
                ijab.addSpinCase(other.getIJAB()(0), "IJ,AB", "ABIJ");
                ijab.addSpinCase(other.getIJAB()(1), "Ij,Ab", "AbIj");
                ijab.addSpinCase(other.getIJAB()(2), "ij,ab", "abij");
            }

            if (copy&AIBJ)
            {
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ()(0)), "AI,BJ", "AIBJ");
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ()(1)), "Ai,Bj", "AiBj");
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ()(2)), "aI,bJ", "aIbJ");
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ()(3)), "ai,bj", "aibj");
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(1)), "aI,Bj", "BaIj", -1.0);
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(1)), "Ai,bJ", "AbJi", -1.0);
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
                aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC()(0)), "AI,BC", "BCAI");
                aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC()(1)), "Ai,Bc", "BcAi");
                aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC()(2)), "aI,bC", "bCaI");
                aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC()(3)), "ai,bc", "bcai");
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
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(0)), "AB,CI", "ABCI");
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(1)), "Ab,Ci", "AbCi");
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(2)), "aB,cI", "aBcI");
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(3)), "ab,ci", "abci");
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
                abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD()(0)), "AB,CD", "ABCD");
                abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD()(1)), "Ab,Cd", "AbCd");
                abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD()(2)), "ab,cd", "abcd");
            }
            else
            {
                abcd.addSpinCase(other.getABCD()(0), "AB,CD", "ABCD");
                abcd.addSpinCase(other.getABCD()(1), "Ab,Cd", "AbCd");
                abcd.addSpinCase(other.getABCD()(2), "ab,cd", "abcd");
            }
        }

        TwoElectronOperator(const TwoElectronOperator<T>& other)
        : OneElectronOperator<T>(other),
          ijkl(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ij,kl"))),
          iajk(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ia,jk"))),
          ijka(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ij,ka"))),
          abij(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ab,ij"))),
          ijab(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ij,ab"))),
          aibj(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ai,bj"))),
          aibc(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ai,bc"))),
          abci(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ab,ci"))),
          abcd(this->addTensor(new tensor::SpinorbitalTensor<tensor::DistTensor<T> >("ab,cd")))
        {
            ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL()(0)), "IJ,KL", "IJKL");
            ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL()(1)), "Ij,Kl", "IjKl");
            ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL()(2)), "ij,kl", "ijkl");

            iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(0)), "IA,JK", "JKIA");
            iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(1)), "Ia,Jk", "JkIa");
            iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(2)), "iA,jK", "jKiA");
            iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(3)), "ia,jk", "jkia");

            ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(0)), "IJ,KA", "IJKA");
            ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(1)), "Ij,Ka", "IjKa");
            ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(2)), "iJ,kA", "iJkA");
            ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(3)), "ij,ka", "ijka");

            abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(0)), "AB,IJ", "ABIJ");
            abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(1)), "Ab,Ij", "AbIj");
            abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(2)), "ab,ij", "abij");

            ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB()(0)), "IJ,AB", "ABIJ");
            ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB()(1)), "Ij,Ab", "AbIj");
            ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB()(2)), "ij,ab", "abij");

            aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ()(0)), "AI,BJ", "AIBJ");
            aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ()(1)), "Ai,Bj", "AiBj");
            aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ()(2)), "aI,bJ", "aIbJ");
            aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ()(3)), "ai,bj", "aibj");
            aibj.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(1)), "aI,Bj", "BaIj", -1.0);
            aibj.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(1)), "Ai,bJ", "AbJi", -1.0);

            aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC()(0)), "AI,BC", "BCAI");
            aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC()(1)), "Ai,Bc", "BcAi");
            aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC()(2)), "aI,bC", "bCaI");
            aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC()(3)), "ai,bc", "bcai");

            abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(0)), "AB,CI", "ABCI");
            abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(1)), "Ab,Ci", "AbCi");
            abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(2)), "aB,cI", "aBcI");
            abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(3)), "ab,ci", "abci");

            abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD()(0)), "AB,CD", "ABCD");
            abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD()(1)), "Ab,Cd", "AbCd");
            abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD()(2)), "ab,cd", "abcd");
        }

        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getIJKL() { return ijkl; }
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getIAJK() { return iajk; }
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getIJKA() { return ijka; }
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getABIJ() { return abij; }
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getIJAB() { return ijab; }
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getAIBJ() { return aibj; }
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getAIBC() { return aibc; }
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getABCI() { return abci; }
        tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getABCD() { return abcd; }

        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getIJKL() const { return ijkl; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getIAJK() const { return iajk; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getIJKA() const { return ijka; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getABIJ() const { return abij; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getIJAB() const { return ijab; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getAIBJ() const { return aibj; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getAIBC() const { return aibc; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getABCI() const { return abci; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<T> >& getABCD() const { return abcd; }
};

}
}

#endif
