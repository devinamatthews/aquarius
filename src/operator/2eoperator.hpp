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

#include "tensor/dist_tensor.hpp"
#include "tensor/spinorbital.hpp"
#include "scf/scf.hpp"

#include "1eoperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class TwoElectronOperator : public OneElectronOperator<T>
{
    protected:
        tensor::SpinorbitalTensor< tensor::DistTensor<T> > ijkl;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> > iajk;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> > ijka;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> > abij;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> > ijab;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> > aibj;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> > aibc;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> > abci;
        tensor::SpinorbitalTensor< tensor::DistTensor<T> > abcd;

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
          ijkl("ij,kl"),
          iajk("ai,jk"),
          ijka("ij,ka"),
          abij("ab,ij"),
          ijab("ij,ab"),
          aibj("ai,bj"),
          aibc("ai,bc"),
          abci("ab,ci"),
          abcd("ab,cd")
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
                ijka.addSpinCase(iajk.getSpinCase(0), "IJ,KA", "IJKA");
                ijka.addSpinCase(iajk.getSpinCase(1), "Ij,Ka", "IjKa");
                ijka.addSpinCase(iajk.getSpinCase(2), "iJ,kA", "iJkA");
                ijka.addSpinCase(iajk.getSpinCase(3), "ij,ka", "ijka");
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
                ijab.addSpinCase(abij.getSpinCase(0), "IJ,AB", "ABIJ");
                ijab.addSpinCase(abij.getSpinCase(1), "Ij,Ab", "AbIj");
                ijab.addSpinCase(abij.getSpinCase(2), "ij,ab", "abij");
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
            aibj.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaIi, shapeNNNN, true), "aI,Bj", "BaIj");
            aibj.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaIi, shapeNNNN, true), "Ai,bJ", "AbJi");

            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAAAI, shapeANNN, true), "AI,BC", "BCAI");
            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaAi, shapeNNNN, true), "Ai,Bc", "BcAi");
            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaAaI, shapeNNNN, true), "aI,bC", "bCaI");
            aibc.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaaai, shapeANNN, true), "ai,bc", "bcai");

            if (hermitian)
            {
                abci.addSpinCase(aibc.getSpinCase(0), "AB,CI", "ABCI");
                abci.addSpinCase(aibc.getSpinCase(1), "Ab,Ci", "AbCi");
                abci.addSpinCase(aibc.getSpinCase(2), "aB,cI", "aBcI");
                abci.addSpinCase(aibc.getSpinCase(3), "aB,ci", "abci");
            }
            else
            {
                abci.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAAAI, shapeANNN, true), "AB,CI", "ABCI");
                abci.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaAi, shapeNNNN, true), "Ab,Ci", "AbCi");
                abci.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaAaI, shapeNNNN, true), "aB,cI", "aBcI");
                abci.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaaai, shapeANNN, true), "aB,ci", "abci");
            }

            abcd.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAAAA, shapeANAN, true), "AB,CD", "ABCD");
            abcd.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeAaAa, shapeNNNN, true), "Ab,Cd", "AbCd");
            abcd.addSpinCase(new tensor::DistTensor<T>(this->ctf, 4, sizeaaaa, shapeANAN, true), "ab,cd", "abcd");
        }

        TwoElectronOperator(OneElectronOperator<T>& other, int copy)
        : OneElectronOperator<T>(other, copy) {}

        TwoElectronOperator(TwoElectronOperator<T>& other, int copy)
        : OneElectronOperator<T>(other, copy),
          ijkl("ij,kl"),
          iajk("ai,jk"),
          ijka("ij,ka"),
          abij("ab,ij"),
          ijab("ij,ab"),
          aibj("ai,bj"),
          aibc("ai,bc"),
          abci("ab,ci"),
          abcd("ab,cd")
        {
            if (copy&IJKL)
            {
                ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL().getSpinCase(0)), "IJ,KL", "IJKL");
                ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL().getSpinCase(1)), "Ij,Kl", "IjKl");
                ijkl.addSpinCase(new tensor::DistTensor<T>(other.getIJKL().getSpinCase(2)), "ij,kl", "ijkl");
            }
            else
            {
                ijkl.addSpinCase(other.getIJKL().getSpinCase(0), "IJ,KL", "IJKL");
                ijkl.addSpinCase(other.getIJKL().getSpinCase(1), "Ij,Kl", "IjKl");
                ijkl.addSpinCase(other.getIJKL().getSpinCase(2), "ij,kl", "ijkl");
            }

            if (copy&IAJK)
            {
                iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK().getSpinCase(0)), "IA,JK", "JKIA");
                iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK().getSpinCase(1)), "Ia,Jk", "JkIa");
                iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK().getSpinCase(2)), "iA,jK", "jKiA");
                iajk.addSpinCase(new tensor::DistTensor<T>(other.getIAJK().getSpinCase(3)), "ia,jk", "jkia");
            }
            else
            {
                iajk.addSpinCase(other.getIAJK().getSpinCase(0), "IA,JK", "JKIA");
                iajk.addSpinCase(other.getIAJK().getSpinCase(1), "Ia,Jk", "JkIa");
                iajk.addSpinCase(other.getIAJK().getSpinCase(2), "iA,jK", "jKiA");
                iajk.addSpinCase(other.getIAJK().getSpinCase(3), "ia,jk", "jkia");
            }

            if (copy&IAJK)
            {
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA().getSpinCase(0)), "IJ,KA", "IJKA");
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA().getSpinCase(1)), "Ij,Ka", "IjKa");
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA().getSpinCase(2)), "iJ,kA", "iJkA");
                ijka.addSpinCase(new tensor::DistTensor<T>(other.getIJKA().getSpinCase(3)), "ij,ka", "ijka");
            }
            else
            {
                ijka.addSpinCase(other.getIJKA().getSpinCase(0), "IJ,KA", "IJKA");
                ijka.addSpinCase(other.getIJKA().getSpinCase(1), "Ij,Ka", "IjKa");
                ijka.addSpinCase(other.getIJKA().getSpinCase(2), "iJ,kA", "iJkA");
                ijka.addSpinCase(other.getIJKA().getSpinCase(3), "ij,ka", "ijka");
            }

            if (copy&ABIJ)
            {
                abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ().getSpinCase(0)), "AB,IJ", "ABIJ");
                abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ().getSpinCase(1)), "Ab,Ij", "AbIj");
                abij.addSpinCase(new tensor::DistTensor<T>(other.getABIJ().getSpinCase(2)), "ab,ij", "abij");
            }
            else
            {
                abij.addSpinCase(other.getABIJ().getSpinCase(0), "AB,IJ", "ABIJ");
                abij.addSpinCase(other.getABIJ().getSpinCase(1), "Ab,Ij", "AbIj");
                abij.addSpinCase(other.getABIJ().getSpinCase(2), "ab,ij", "abij");
            }

            if (copy&IJAB)
            {
                ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB().getSpinCase(0)), "IJ,AB", "ABIJ");
                ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB().getSpinCase(1)), "Ij,Ab", "AbIj");
                ijab.addSpinCase(new tensor::DistTensor<T>(other.getIJAB().getSpinCase(2)), "ij,ab", "abij");
            }
            else
            {
                ijab.addSpinCase(other.getIJAB().getSpinCase(0), "IJ,AB", "ABIJ");
                ijab.addSpinCase(other.getIJAB().getSpinCase(1), "Ij,Ab", "AbIj");
                ijab.addSpinCase(other.getIJAB().getSpinCase(2), "ij,ab", "abij");
            }

            if (copy&AIBJ)
            {
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ().getSpinCase(0)), "AI,BJ", "AIBJ");
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ().getSpinCase(1)), "Ai,Bj", "AiBj");
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ().getSpinCase(2)), "aI,bJ", "aIbJ");
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getAIBJ().getSpinCase(3)), "ai,bj", "aibj");
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getABIJ().getSpinCase(1)), "aI,Bj", "BaIj", -1.0);
                aibj.addSpinCase(new tensor::DistTensor<T>(other.getABIJ().getSpinCase(1)), "Ai,bJ", "AbJi", -1.0);
            }
            else
            {
                aibj.addSpinCase(other.getAIBJ().getSpinCase(0), "AI,BJ", "AIBJ");
                aibj.addSpinCase(other.getAIBJ().getSpinCase(1), "Ai,Bj", "AiBj");
                aibj.addSpinCase(other.getAIBJ().getSpinCase(2), "aI,bJ", "aIbJ");
                aibj.addSpinCase(other.getAIBJ().getSpinCase(3), "ai,bj", "aibj");
                aibj.addSpinCase(other.getABIJ().getSpinCase(1), "aI,Bj", "BaIj", -1.0);
                aibj.addSpinCase(other.getABIJ().getSpinCase(1), "Ai,bJ", "AbJi", -1.0);
            }

            if (copy&AIBC)
            {
                aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC().getSpinCase(0)), "AI,BC", "BCAI");
                aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC().getSpinCase(1)), "Ai,Bc", "BcAi");
                aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC().getSpinCase(2)), "aI,bC", "bCaI");
                aibc.addSpinCase(new tensor::DistTensor<T>(other.getAIBC().getSpinCase(3)), "ai,bc", "bcai");
            }
            else
            {
                aibc.addSpinCase(other.getAIBC().getSpinCase(0), "AI,BC", "BCAI");
                aibc.addSpinCase(other.getAIBC().getSpinCase(1), "Ai,Bc", "BcAi");
                aibc.addSpinCase(other.getAIBC().getSpinCase(2), "aI,bC", "bCaI");
                aibc.addSpinCase(other.getAIBC().getSpinCase(3), "ai,bc", "bcai");
            }

            if (copy&ABCI)
            {
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI().getSpinCase(0)), "AB,CI", "ABCI");
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI().getSpinCase(1)), "Ab,Ci", "AbCi");
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI().getSpinCase(2)), "aB,cI", "aBcI");
                abci.addSpinCase(new tensor::DistTensor<T>(other.getABCI().getSpinCase(3)), "aB,ci", "abci");
            }
            else
            {
                abci.addSpinCase(other.getABCI().getSpinCase(0), "AB,CI", "ABCI");
                abci.addSpinCase(other.getABCI().getSpinCase(1), "Ab,Ci", "AbCi");
                abci.addSpinCase(other.getABCI().getSpinCase(2), "aB,cI", "aBcI");
                abci.addSpinCase(other.getABCI().getSpinCase(3), "aB,ci", "abci");
            }

            if (copy&ABCD)
            {
                abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD().getSpinCase(0)), "AB,CD", "ABCD");
                abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD().getSpinCase(1)), "Ab,Cd", "AbCd");
                abcd.addSpinCase(new tensor::DistTensor<T>(other.getABCD().getSpinCase(2)), "ab,cd", "abcd");
            }
            else
            {
                abcd.addSpinCase(other.getABCD().getSpinCase(0), "AB,CD", "ABCD");
                abcd.addSpinCase(other.getABCD().getSpinCase(1), "Ab,Cd", "AbCd");
                abcd.addSpinCase(other.getABCD().getSpinCase(2), "ab,cd", "abcd");
            }
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
