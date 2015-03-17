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
#include "util/stl_ext.hpp"

#include "mooperator.hpp"
#include "space.hpp"

namespace aquarius
{
namespace op
{

template <typename T, typename Derived>
class OneElectronOperatorBase : public MOOperator,
    public tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>
{
    INHERIT_FROM_COMPOSITE_TENSOR(Derived,tensor::SpinorbitalTensor<T>,T)

    protected:
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
            ALL  = 0xFFFF
        };

        OneElectronOperatorBase(const std::string& name, const Arena& arena, const Space& occ, const Space& vrt)
        : MOOperator(arena, occ, vrt), tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>(name),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>(name, arena, occ.group, std::vec(vrt, occ), std::vec(1,0), std::vec(1,0)))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>(name, arena, occ.group, std::vec(vrt, occ), std::vec(0,1), std::vec(0,1)))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>(name, arena, occ.group, std::vec(vrt, occ), std::vec(1,0), std::vec(0,1)))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>(name, arena, occ.group, std::vec(vrt, occ), std::vec(0,1), std::vec(1,0)))) {}

        OneElectronOperatorBase(const std::string& name, const MOSpace<T>& occ, const MOSpace<T>& vrt,
                                const tensor::SymmetryBlockedTensor<T>& aoa,
                                const tensor::SymmetryBlockedTensor<T>& aob)
        : MOOperator(occ.arena, occ, vrt), tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>(name),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>(name, occ.arena, occ.group, std::vec<Space>(vrt, occ), std::vec(1,0), std::vec(1,0)))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>(name, occ.arena, occ.group, std::vec<Space>(vrt, occ), std::vec(0,1), std::vec(0,1)))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>(name, occ.arena, occ.group, std::vec<Space>(vrt, occ), std::vec(1,0), std::vec(0,1)))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>(name, occ.arena, occ.group, std::vec<Space>(vrt, occ), std::vec(0,1), std::vec(1,0))))
        {
            const tensor::SymmetryBlockedTensor<T>& cA = vrt.Calpha;
            const tensor::SymmetryBlockedTensor<T>& ca = vrt.Cbeta;
            const tensor::SymmetryBlockedTensor<T>& cI = occ.Calpha;
            const tensor::SymmetryBlockedTensor<T>& ci = occ.Cbeta;

            const std::vector<int>& N = occ.nao;
            const std::vector<int>& nI = occ.nalpha;
            const std::vector<int>& ni = occ.nbeta;
            const std::vector<int>& nA = vrt.nalpha;
            const std::vector<int>& na = vrt.nbeta;

            std::vector<std::vector<int> > sizeAA = std::vec(nA, nA);
            std::vector<std::vector<int> > sizeaa = std::vec(na, na);
            std::vector<std::vector<int> > sizeAI = std::vec(nA, nI);
            std::vector<std::vector<int> > sizeai = std::vec(na, ni);
            std::vector<std::vector<int> > sizeII = std::vec(nI, nI);
            std::vector<std::vector<int> > sizeii = std::vec(ni, ni);
            std::vector<std::vector<int> > sizeAN = std::vec(nA, N);
            std::vector<std::vector<int> > sizeaN = std::vec(na, N);
            std::vector<std::vector<int> > sizeIN = std::vec(nI, N);
            std::vector<std::vector<int> > sizeiN = std::vec(ni, N);

            std::vector<int> shapeNN = std::vec((int)NS, (int)NS);

            tensor::SymmetryBlockedTensor<T> Aq("Aq", this->arena, occ.group, 2, sizeAN, shapeNN, false);
            tensor::SymmetryBlockedTensor<T> aq("aq", this->arena, occ.group, 2, sizeaN, shapeNN, false);
            tensor::SymmetryBlockedTensor<T> Iq("Iq", this->arena, occ.group, 2, sizeIN, shapeNN, false);
            tensor::SymmetryBlockedTensor<T> iq("iq", this->arena, occ.group, 2, sizeiN, shapeNN, false);

            Aq["Aq"] = cA["pA"]*aoa["pq"];
            aq["aq"] = ca["pa"]*aob["pq"];
            Iq["Iq"] = cI["pI"]*aoa["pq"];
            iq["iq"] = ci["pi"]*aob["pq"];

            ab(std::vec(1,0), std::vec(1,0))["AB"] = Aq["Aq"]*cA["qB"];
            ab(std::vec(0,0), std::vec(0,0))["ab"] = aq["aq"]*ca["qb"];

            ij(std::vec(0,1), std::vec(0,1))["IJ"] = Iq["Iq"]*cI["qJ"];
            ij(std::vec(0,0), std::vec(0,0))["ij"] = iq["iq"]*ci["qj"];

            ai(std::vec(1,0), std::vec(0,1))["AI"] = Aq["Aq"]*cI["qI"];
            ai(std::vec(0,0), std::vec(0,0))["ai"] = aq["aq"]*ci["qi"];

            ia(std::vec(0,1), std::vec(1,0))["IA"] = Iq["Iq"]*cA["qA"];
            ia(std::vec(0,0), std::vec(0,0))["ia"] = iq["iq"]*ca["qa"];
        }

        template <typename otherDerived>
        OneElectronOperatorBase(const std::string& name, OneElectronOperatorBase<T,otherDerived>& other, int copy)
        : MOOperator(other), tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>(name),
          ab(copy&AB ? this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getAB())) : this->addTensor(other.getAB())),
          ij(copy&IJ ? this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getIJ())) : this->addTensor(other.getIJ())),
          ai(copy&AI ? this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getAI())) : this->addTensor(other.getAI())),
          ia(copy&IA ? this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getIA())) : this->addTensor(other.getIA())) {}

        template <typename otherDerived>
        OneElectronOperatorBase(const OneElectronOperatorBase<T,otherDerived>& other)
        : MOOperator(other), tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>(other.name),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>(other.getAB()))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>(other.getIJ()))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>(other.getAI()))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>(other.getIA()))) {}

        template <typename otherDerived>
        OneElectronOperatorBase(const std::string& name, const OneElectronOperatorBase<T,otherDerived>& other)
        : MOOperator(other), tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>(name),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getAB()))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getIJ()))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getAI()))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getIA()))) {}

        OneElectronOperatorBase(const OneElectronOperatorBase<T,Derived>& other)
        : MOOperator(other), tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>(other.name),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>(other.getAB()))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>(other.getIJ()))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>(other.getAI()))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>(other.getIA()))) {}

        OneElectronOperatorBase(const std::string& name, const OneElectronOperatorBase<T,Derived>& other)
        : MOOperator(other), tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>(name),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getAB()))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getIJ()))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getAI()))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>(name, other.getIA()))) {}

        tensor::SpinorbitalTensor<T>& getAB() { return ab; }
        tensor::SpinorbitalTensor<T>& getIJ() { return ij; }
        tensor::SpinorbitalTensor<T>& getAI() { return ai; }
        tensor::SpinorbitalTensor<T>& getIA() { return ia; }

        const tensor::SpinorbitalTensor<T>& getAB() const { return ab; }
        const tensor::SpinorbitalTensor<T>& getIJ() const { return ij; }
        const tensor::SpinorbitalTensor<T>& getAI() const { return ai; }
        const tensor::SpinorbitalTensor<T>& getIA() const { return ia; }
};

template <typename T>
class OneElectronOperator : public OneElectronOperatorBase<T,OneElectronOperator<T> >
{
    public:
        OneElectronOperator(const std::string& name, const Arena& arena, const Space& occ, const Space& vrt)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, arena, occ, vrt) {}

        OneElectronOperator(const std::string& name, const MOSpace<T>& occ, const MOSpace<T>& vrt,
                            const tensor::SymmetryBlockedTensor<T>& aoa,
                            const tensor::SymmetryBlockedTensor<T>& aob)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, occ, vrt, aoa, aob) {}

        template <typename Derived>
        OneElectronOperator(const std::string& name, OneElectronOperatorBase<T,Derived>& other, int copy)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, other, copy) {}

        template <typename Derived>
        OneElectronOperator(const OneElectronOperatorBase<T,Derived>& other)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(other) {}

        template <typename Derived>
        OneElectronOperator(const std::string& name, const OneElectronOperatorBase<T,Derived>& other)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, other) {}

        OneElectronOperator(const OneElectronOperator<T>& other)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(other) {}

        OneElectronOperator(const std::string& name, const OneElectronOperator<T>& other)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, other) {}
};

}
}

#endif
