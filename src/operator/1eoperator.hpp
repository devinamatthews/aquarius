#ifndef _AQUARIUS_OPERATOR_1EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_1EOPERATOR_HPP_

#include "util/global.hpp"

#include "tensor/spinorbital_tensor.hpp"

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

        OneElectronOperatorBase(const string& name, const Arena& arena, const Space& occ, const Space& vrt)
        : MOOperator(arena, occ, vrt), tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>(name),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>(name, arena, occ.group, {vrt, occ}, {1,0}, {1,0}))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>(name, arena, occ.group, {vrt, occ}, {0,1}, {0,1}))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>(name, arena, occ.group, {vrt, occ}, {1,0}, {0,1}))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>(name, arena, occ.group, {vrt, occ}, {0,1}, {1,0}))) {}

        OneElectronOperatorBase(const string& name, const MOSpace<T>& occ, const MOSpace<T>& vrt,
                                const tensor::SymmetryBlockedTensor<T>& aoa,
                                const tensor::SymmetryBlockedTensor<T>& aob)
        : MOOperator(occ.arena, occ, vrt), tensor::CompositeTensor<Derived,tensor::SpinorbitalTensor<T>,T>(name),
          ab(this->addTensor(new tensor::SpinorbitalTensor<T>(name, occ.arena, occ.group, {vrt, occ}, {1,0}, {1,0}))),
          ij(this->addTensor(new tensor::SpinorbitalTensor<T>(name, occ.arena, occ.group, {vrt, occ}, {0,1}, {0,1}))),
          ai(this->addTensor(new tensor::SpinorbitalTensor<T>(name, occ.arena, occ.group, {vrt, occ}, {1,0}, {0,1}))),
          ia(this->addTensor(new tensor::SpinorbitalTensor<T>(name, occ.arena, occ.group, {vrt, occ}, {0,1}, {1,0})))
        {
            const tensor::SymmetryBlockedTensor<T>& cA = vrt.Calpha;
            const tensor::SymmetryBlockedTensor<T>& ca = vrt.Cbeta;
            const tensor::SymmetryBlockedTensor<T>& cI = occ.Calpha;
            const tensor::SymmetryBlockedTensor<T>& ci = occ.Cbeta;

            const vector<int>& N = occ.nao;
            const vector<int>& nI = occ.nalpha;
            const vector<int>& ni = occ.nbeta;
            const vector<int>& nA = vrt.nalpha;
            const vector<int>& na = vrt.nbeta;

            vector<vector<int> > sizeAA = {nA, nA};
            vector<vector<int> > sizeaa = {na, na};
            vector<vector<int> > sizeAI = {nA, nI};
            vector<vector<int> > sizeai = {na, ni};
            vector<vector<int> > sizeII = {nI, nI};
            vector<vector<int> > sizeii = {ni, ni};
            vector<vector<int> > sizeAN = {nA, N};
            vector<vector<int> > sizeaN = {na, N};
            vector<vector<int> > sizeIN = {nI, N};
            vector<vector<int> > sizeiN = {ni, N};

            vector<int> shapeNN = {NS, NS};

            tensor::SymmetryBlockedTensor<T> Aq("Aq", this->arena, occ.group, 2, sizeAN, shapeNN, false);
            tensor::SymmetryBlockedTensor<T> aq("aq", this->arena, occ.group, 2, sizeaN, shapeNN, false);
            tensor::SymmetryBlockedTensor<T> Iq("Iq", this->arena, occ.group, 2, sizeIN, shapeNN, false);
            tensor::SymmetryBlockedTensor<T> iq("iq", this->arena, occ.group, 2, sizeiN, shapeNN, false);

            Aq["Aq"] = cA["pA"]*aoa["pq"];
            aq["aq"] = ca["pa"]*aob["pq"];
            Iq["Iq"] = cI["pI"]*aoa["pq"];
            iq["iq"] = ci["pi"]*aob["pq"];

            ab({1,0}, {1,0})["AB"] = Aq["Aq"]*cA["qB"];
            ab({0,0}, {0,0})["ab"] = aq["aq"]*ca["qb"];

            ij({0,1}, {0,1})["IJ"] = Iq["Iq"]*cI["qJ"];
            ij({0,0}, {0,0})["ij"] = iq["iq"]*ci["qj"];

            ai({1,0}, {0,1})["AI"] = Aq["Aq"]*cI["qI"];
            ai({0,0}, {0,0})["ai"] = aq["aq"]*ci["qi"];

            ia({0,1}, {1,0})["IA"] = Iq["Iq"]*cA["qA"];
            ia({0,0}, {0,0})["ia"] = iq["iq"]*ca["qa"];
        }

        template <typename otherDerived>
        OneElectronOperatorBase(const string& name, OneElectronOperatorBase<T,otherDerived>& other, int copy)
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
        OneElectronOperatorBase(const string& name, const OneElectronOperatorBase<T,otherDerived>& other)
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

        OneElectronOperatorBase(const string& name, const OneElectronOperatorBase<T,Derived>& other)
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
        OneElectronOperator(const string& name, const Arena& arena, const Space& occ, const Space& vrt)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, arena, occ, vrt) {}

        OneElectronOperator(const string& name, const MOSpace<T>& occ, const MOSpace<T>& vrt,
                            const tensor::SymmetryBlockedTensor<T>& aoa,
                            const tensor::SymmetryBlockedTensor<T>& aob)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, occ, vrt, aoa, aob) {}

        template <typename Derived>
        OneElectronOperator(const string& name, OneElectronOperatorBase<T,Derived>& other, int copy)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, other, copy) {}

        template <typename Derived>
        OneElectronOperator(const OneElectronOperatorBase<T,Derived>& other)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(other) {}

        template <typename Derived>
        OneElectronOperator(const string& name, const OneElectronOperatorBase<T,Derived>& other)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, other) {}

        OneElectronOperator(const OneElectronOperator<T>& other)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(other) {}

        OneElectronOperator(const string& name, const OneElectronOperator<T>& other)
        : OneElectronOperatorBase<T,OneElectronOperator<T> >(name, other) {}
};

}
}

#endif
