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

#ifndef _AQUARIUS_OPERATOR_1EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_1EOPERATOR_HPP_

#include "tensor/spinorbital_tensor.hpp"

#include "mooperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T, typename Derived>
class OneElectronOperatorBase : public MOOperator<T>,
    public tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>
{
    INHERIT_FROM_COMPOSITE_TENSOR(Derived,tensor::SpinorbitalTensor<T>,T)

    protected:
        const bool hermitian;
        tensor::SpinorbitalTensor<T>& ab;
        tensor::SpinorbitalTensor<T>& ij;
        tensor::SpinorbitalTensor<T>& ai;
        tensor::SpinorbitalTensor<T>& ia;

    public:
        enum
        {
            NONE = 0x0000,
            AB   = 0x0001,
            IJ   = 0x0002,
            AI   = 0x0004,
            IA   = 0x0008,
            ALL  = 0xffffffff
        };

        OneElectronOperatorBase(const scf::UHF<T>& uhf, const bool hermitian=true)
        : MOOperator<T>(uhf),
          hermitian(hermitian),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>("a,b"))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>("i,j"))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>("a,i"))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>("i,a")))
        {
            int N = uhf.getMolecule().getNumOrbitals();
            int nI = uhf.getMolecule().getNumAlphaElectrons();
            int ni = uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;

            std::vector<int> sizeAA = vec(nA, nA);
            std::vector<int> sizeaa = vec(na, na);
            std::vector<int> sizeAI = vec(nA, nI);
            std::vector<int> sizeai = vec(na, ni);
            std::vector<int> sizeII = vec(nI, nI);
            std::vector<int> sizeii = vec(ni, ni);

            std::vector<int> shapeNN = vec(NS, NS);

            ab.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeAA, shapeNN, true), "A,B", "AB");
            ab.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeaa, shapeNN, true), "a,b", "ab");

            ij.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeII, shapeNN, true), "I,J", "IJ");
            ij.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeii, shapeNN, true), "i,j", "ij");

            ai.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeAI, shapeNN, true), "A,I", "AI");
            ai.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeai, shapeNN, true), "a,i", "ai");

            if (hermitian)
            {
                ia.addSpinCase(ai(0), "I,A", "AI");
                ia.addSpinCase(ai(1), "i,a", "ai");
            }
            else
            {
                ia.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeAI, shapeNN, true), "I,A", "AI");
                ia.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeai, shapeNN, true), "i,a", "ai");
            }
        }

        OneElectronOperatorBase(const scf::UHF<T>& uhf, const tensor::DistTensor<T>& ao, const bool hermitian=true)
        : MOOperator<T>(uhf),
          hermitian(hermitian),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>("a,b"))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>("i,j"))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>("a,i"))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>("i,a")))
        {
            const tensor::DistTensor<T>& cA = uhf.getCA();
            const tensor::DistTensor<T>& ca = uhf.getCa();
            const tensor::DistTensor<T>& cI = uhf.getCI();
            const tensor::DistTensor<T>& ci = uhf.getCi();

            int N = uhf.getMolecule().getNumOrbitals();
            int nI = uhf.getMolecule().getNumAlphaElectrons();
            int ni = uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;


            std::vector<int> sizeAA = vec(nA, nA);
            std::vector<int> sizeaa = vec(na, na);
            std::vector<int> sizeAI = vec(nA, nI);
            std::vector<int> sizeai = vec(na, ni);
            std::vector<int> sizeII = vec(nI, nI);
            std::vector<int> sizeii = vec(ni, ni);
            std::vector<int> sizeAN = vec(nA, N);
            std::vector<int> sizeaN = vec(na, N);
            std::vector<int> sizeIN = vec(nI, N);
            std::vector<int> sizeiN = vec(ni, N);

            std::vector<int>tor<int> shapeNN = vec(NS, NS);

            tensor::DistTensor<T> Aq(this->arena, 2, sizeAN, shapeNN, false);
            tensor::DistTensor<T> aq(this->arena, 2, sizeaN, shapeNN, false);
            tensor::DistTensor<T> Iq(this->arena, 2, sizeIN, shapeNN, false);
            tensor::DistTensor<T> iq(this->arena, 2, sizeiN, shapeNN, false);

            Aq["Aq"] = cA["pA"]*ao["pq"];
            aq["aq"] = ca["pa"]*ao["pq"];
            Iq["Iq"] = cI["pI"]*ao["pq"];
            iq["iq"] = ci["pi"]*ao["pq"];

            ab.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeAA, shapeNN, true), "A,B", "AB");
            ab.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeaa, shapeNN, true), "a,b", "ab");

            ab(1,0,1,0)["AB"] = Aq["Aq"]*cA["qB"];
            ab(0,0,0,0)["ab"] = aq["aq"]*ca["qb"];

            ij.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeII, shapeNN, true), "I,J", "IJ");
            ij.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeii, shapeNN, true), "i,j", "ij");

            ij(0,1,0,1)["IJ"] = Iq["Iq"]*cI["qJ"];
            ij(0,0,0,0)["ij"] = iq["iq"]*ci["qj"];

            ai.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeAI, shapeNN, true), "A,I", "AI");
            ai.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeai, shapeNN, true), "a,i", "ai");

            ai(1,0,0,1)["AI"] = Aq["Aq"]*cI["qI"];
            ai(0,0,0,0)["ai"] = aq["aq"]*ci["qi"];

            if (hermitian)
            {
                ia.addSpinCase(ai(0), "I,A", "AI");
                ia.addSpinCase(ai(1), "i,a", "ai");
            }
            else
            {
                ia.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeAI, shapeNN, true), "I,A", "AI");
                ia.addSpinCase(new tensor::DistTensor<T>(this->arena, 2, sizeai, shapeNN, true), "i,a", "ai");

                ia(0,1,1,0)["IA"] = Iq["Iq"]*cA["qA"];
                ia(0,0,0,0)["ia"] = iq["iq"]*ca["qa"];
            }
        }

        template <typename otherDerived>
        OneElectronOperatorBase(OneElectronOperatorBase<T,otherDerived>& other, int copy)
        : MOOperator<T>(other),
          hermitian(copy == NONE && other.isHermitian()),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>("a,b"))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>("i,j"))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>("a,i"))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>("i,a")))
        {
            if (copy&AB)
            {
                ab.addSpinCase(new tensor::DistTensor<T>(other.getAB()(0)), "A,B", "AB");
                ab.addSpinCase(new tensor::DistTensor<T>(other.getAB()(1)), "a,b", "ab");
            }
            else
            {
                ab.addSpinCase(other.getAB()(0), "A,B", "AB");
                ab.addSpinCase(other.getAB()(1), "a,b", "ab");
            }

            if (copy&IJ)
            {
                ij.addSpinCase(new tensor::DistTensor<T>(other.getIJ()(0)), "I,J", "IJ");
                ij.addSpinCase(new tensor::DistTensor<T>(other.getIJ()(1)), "i,j", "ij");
            }
            else
            {
                ij.addSpinCase(other.getIJ()(0), "I,J", "IJ");
                ij.addSpinCase(other.getIJ()(1), "i,j", "ij");
            }

            if (copy&AI)
            {
                ai.addSpinCase(new tensor::DistTensor<T>(other.getAI()(0)), "A,I", "AI");
                ai.addSpinCase(new tensor::DistTensor<T>(other.getAI()(1)), "a,i", "ai");
            }
            else
            {
                ai.addSpinCase(other.getAI()(0), "A,I", "AI");
                ai.addSpinCase(other.getAI()(1), "a,i", "ai");
            }

            if (copy&IA)
            {
                ia.addSpinCase(new tensor::DistTensor<T>(other.getIA()(0)), "I,A", "AI");
                ia.addSpinCase(new tensor::DistTensor<T>(other.getIA()(1)), "i,a", "ai");
            }
            else
            {
                ia.addSpinCase(other.getIA()(0), "I,A", "AI");
                ia.addSpinCase(other.getIA()(1), "i,a", "ai");
            }
        }

        template <typename otherDerived>
        OneElectronOperatorBase(const OneElectronOperatorBase<T,otherDerived>& other)
        : MOOperator<T>(other),
          hermitian(other.isHermitian()),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>("a,b"))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>("i,j"))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>("a,i"))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>("i,a")))
        {
            ab.addSpinCase(new tensor::DistTensor<T>(other.getAB()(0)), "A,B", "AB");
            ab.addSpinCase(new tensor::DistTensor<T>(other.getAB()(1)), "a,b", "ab");

            ij.addSpinCase(new tensor::DistTensor<T>(other.getIJ()(0)), "I,J", "IJ");
            ij.addSpinCase(new tensor::DistTensor<T>(other.getIJ()(1)), "i,j", "ij");

            ai.addSpinCase(new tensor::DistTensor<T>(other.getAI()(0)), "A,I", "AI");
            ai.addSpinCase(new tensor::DistTensor<T>(other.getAI()(1)), "a,i", "ai");

            ia.addSpinCase(new tensor::DistTensor<T>(other.getIA()(0)), "I,A", "AI");
            ia.addSpinCase(new tensor::DistTensor<T>(other.getIA()(1)), "i,a", "ai");
        }

        bool isHermitian() const { return hermitian; }

        tensor::SpinorbitalTensor<T>& getAB()   { return ab; }
        tensor::SpinorbitalTensor<T>& getIJ()   { return ij; }
        tensor::SpinorbitalTensor<T>& getAI()   { return ai; }
        tensor::SpinorbitalTensor<T>& getIA()   { return ia; }

        const tensor::SpinorbitalTensor<T>& getAB()   const { return ab; }
        const tensor::SpinorbitalTensor<T>& getIJ()   const { return ij; }
        const tensor::SpinorbitalTensor<T>& getAI()   const { return ai; }
        const tensor::SpinorbitalTensor<T>& getIA()   const { return ia; }
};

template <typename T>
class OneElectronOperator : public OneElectronOperatorBase<T,OneElectronOperator<T> >
{
    public:
        OneElectronOperator(const scf::UHF<T>& uhf, const bool hermitian=true)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(uhf, hermitian) {}

        OneElectronOperator(const scf::UHF<T>& uhf, const tensor::DistTensor<T>& ao, const bool hermitian=true)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(uhf, ao, hermitian) {}

        template <typename Derived>
        OneElectronOperator(OneElectronOperatorBase<T,Derived>& other, int copy)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(other, copy) {}

        template <typename Derived>
        OneElectronOperator(const OneElectronOperatorBase<T,Derived>& other)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(other) {}
};

}
}

#endif
