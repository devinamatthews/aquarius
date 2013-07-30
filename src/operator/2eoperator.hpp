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
class TwoElectronOperator : public OneElectronOperatorBase<T,TwoElectronOperator<T> >
{
    protected:
        tensor::SpinorbitalTensor<T>& ijkl;
        tensor::SpinorbitalTensor<T>& iajk;
        tensor::SpinorbitalTensor<T>& ijka;
        tensor::SpinorbitalTensor<T>& abij;
        tensor::SpinorbitalTensor<T>& ijab;
        tensor::SpinorbitalTensor<T>& aibj;
        tensor::SpinorbitalTensor<T>& aibc;
        tensor::SpinorbitalTensor<T>& abci;
        tensor::SpinorbitalTensor<T>& abcd;

        void initialize()
        {
            int N = this->uhf.getMolecule().getNumOrbitals();
            int nI = this->uhf.getMolecule().getNumAlphaElectrons();
            int ni = this->uhf.getMolecule().getNumBetaElectrons();
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

            if (this->hermitian)
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

            if (this->hermitian)
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
            aibj.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaIi, shapeNNNN, true), "aI,Bj", "BaIj");

            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAAAI, shapeANNN, true), "AI,BC", "BCAI");
            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaAi, shapeNNNN, true), "Ai,Bc", "BcAi");
            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaAaI, shapeNNNN, true), "aI,bC", "bCaI");
            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaaai, shapeANNN, true), "ai,bc", "bcai");

            if (this->hermitian)
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
        : OneElectronOperatorBase<T,TwoElectronOperator<T> >(uhf, hermitian),
          ijkl(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,kl"))),
          iajk(this->addTensor(new tensor::SpinorbitalTensor<T>("ia,jk"))),
          ijka(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ka"))),
          abij(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ij"))),
          ijab(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ab"))),
          aibj(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bj"))),
          aibc(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bc"))),
          abci(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ci"))),
          abcd(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,cd")))
        {
            initialize();
        }

        TwoElectronOperator(OneElectronOperator<T>& other, int copy)
        : OneElectronOperatorBase<T,TwoElectronOperator<T> >(other, copy),
          ijkl(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,kl"))),
          iajk(this->addTensor(new tensor::SpinorbitalTensor<T>("ia,jk"))),
          ijka(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ka"))),
          abij(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ij"))),
          ijab(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ab"))),
          aibj(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bj"))),
          aibc(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bc"))),
          abci(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ci"))),
          abcd(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,cd")))
        {
          initialize();
        }

        TwoElectronOperator(const OneElectronOperator<T>& other)
        : OneElectronOperatorBase<T,TwoElectronOperator<T> >(other),
          ijkl(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,kl"))),
          iajk(this->addTensor(new tensor::SpinorbitalTensor<T>("ia,jk"))),
          ijka(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ka"))),
          abij(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ij"))),
          ijab(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ab"))),
          aibj(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bj"))),
          aibc(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bc"))),
          abci(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ci"))),
          abcd(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,cd")))
        {
          initialize();
        }

        TwoElectronOperator(TwoElectronOperator<T>& other, int copy)
        : OneElectronOperatorBase<T,TwoElectronOperator<T> >(other, copy),
          ijkl(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,kl"))),
          iajk(this->addTensor(new tensor::SpinorbitalTensor<T>("ia,jk"))),
          ijka(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ka"))),
          abij(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ij"))),
          ijab(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ab"))),
          aibj(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bj"))),
          aibc(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bc"))),
          abci(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ci"))),
          abcd(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,cd")))
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

            if (copy&IJKA)
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
        : OneElectronOperatorBase<T,TwoElectronOperator<T> >(other),
          ijkl(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,kl"))),
          iajk(this->addTensor(new tensor::SpinorbitalTensor<T>("ia,jk"))),
          ijka(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ka"))),
          abij(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ij"))),
          ijab(this->addTensor(new tensor::SpinorbitalTensor<T>("ij,ab"))),
          aibj(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bj"))),
          aibc(this->addTensor(new tensor::SpinorbitalTensor<T>("ai,bc"))),
          abci(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,ci"))),
          abcd(this->addTensor(new tensor::SpinorbitalTensor<T>("ab,cd")))
        {
            ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL()(0)), "IJ,KL", "IJKL");
            ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL()(1)), "Ij,Kl", "IjKl");
            ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL()(2)), "ij,kl", "ijkl");

            iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(0)), "IA,JK", "JKIA");
            iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(1)), "Ia,Jk", "JkIa");
            iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(2)), "iA,jK", "jKiA");
            iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK()(3)), "ia,jk", "jkia");

            if (this->hermitian)
            {
                ijka.addSpinCase(this->iajk(0), "IJ,KA", "IJKA");
                ijka.addSpinCase(this->iajk(1), "Ij,Ka", "IjKa");
                ijka.addSpinCase(this->iajk(2), "iJ,kA", "iJkA");
                ijka.addSpinCase(this->iajk(3), "ij,ka", "ijka");
            }
            else
            {
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(0)), "IJ,KA", "IJKA");
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(1)), "Ij,Ka", "IjKa");
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(2)), "iJ,kA", "iJkA");
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA()(3)), "ij,ka", "ijka");
            }

            abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(0)), "AB,IJ", "ABIJ");
            abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(1)), "Ab,Ij", "AbIj");
            abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ()(2)), "ab,ij", "abij");

            if (this->hermitian)
            {
                ijab.addSpinCase(this->abij(0), "IJ,AB", "ABIJ");
                ijab.addSpinCase(this->abij(1), "Ij,Ab", "AbIj");
                ijab.addSpinCase(this->abij(2), "ij,ab", "abij");
            }
            else
            {
                ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB()(0)), "IJ,AB", "ABIJ");
                ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB()(1)), "Ij,Ab", "AbIj");
                ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB()(2)), "ij,ab", "abij");
            }

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

            if (this->hermitian)
            {
                abci.addSpinCase(this->aibc(0), "AB,CI", "ABCI");
                abci.addSpinCase(this->aibc(1), "Ab,Ci", "AbCi");
                abci.addSpinCase(this->aibc(2), "aB,cI", "aBcI");
                abci.addSpinCase(this->aibc(3), "ab,ci", "abci");
            }
            else
            {
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(0)), "AB,CI", "ABCI");
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(1)), "Ab,Ci", "AbCi");
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(2)), "aB,cI", "aBcI");
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI()(3)), "ab,ci", "abci");
            }

            abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD()(0)), "AB,CD", "ABCD");
            abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD()(1)), "Ab,Cd", "AbCd");
            abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD()(2)), "ab,cd", "abcd");
        }

        T dot(bool conja, const TwoElectronOperator<T>& A, bool conjb) const
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

        tensor::SpinorbitalTensor<T>& getIJKL() { return ijkl; }
        tensor::SpinorbitalTensor<T>& getIAJK() { return iajk; }
        tensor::SpinorbitalTensor<T>& getIJKA() { return ijka; }
        tensor::SpinorbitalTensor<T>& getABIJ() { return abij; }
        tensor::SpinorbitalTensor<T>& getIJAB() { return ijab; }
        tensor::SpinorbitalTensor<T>& getAIBJ() { return aibj; }
        tensor::SpinorbitalTensor<T>& getAIBC() { return aibc; }
        tensor::SpinorbitalTensor<T>& getABCI() { return abci; }
        tensor::SpinorbitalTensor<T>& getABCD() { return abcd; }

        const tensor::SpinorbitalTensor<T>& getIJKL() const { return ijkl; }
        const tensor::SpinorbitalTensor<T>& getIAJK() const { return iajk; }
        const tensor::SpinorbitalTensor<T>& getIJKA() const { return ijka; }
        const tensor::SpinorbitalTensor<T>& getABIJ() const { return abij; }
        const tensor::SpinorbitalTensor<T>& getIJAB() const { return ijab; }
        const tensor::SpinorbitalTensor<T>& getAIBJ() const { return aibj; }
        const tensor::SpinorbitalTensor<T>& getAIBC() const { return aibc; }
        const tensor::SpinorbitalTensor<T>& getABCI() const { return abci; }
        const tensor::SpinorbitalTensor<T>& getABCD() const { return abcd; }
};

}

/*
template <typename T>
T scalar(const tensor::TensorMult<op::TwoElectronOperator<T>,T>& tm)
{
    std::cout << "here2" << std::endl;
    return tm.factor_*tm.B_.tensor_.dot(tm.A_.conj_, tm.A_.tensor_, tm.B_.conj_);
}
*/

}

#endif
