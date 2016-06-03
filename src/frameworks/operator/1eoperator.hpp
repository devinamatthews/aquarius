#ifndef _AQUARIUS_FRAMEWORKS_OPERATOR_1EOPERATOR_HPP_
#define _AQUARIUS_FRAMEWORKS_OPERATOR_1EOPERATOR_HPP_

#include "frameworks/util.hpp"
#include "frameworks/tensor.hpp"

namespace aquarius
{
namespace op
{

class OneElectronOperator
{
    protected:
        tensor::Tensor<SPINORBITAL> ab;
        tensor::Tensor<SPINORBITAL> ij;
        tensor::Tensor<SPINORBITAL> ai;
        tensor::Tensor<SPINORBITAL> ia;

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

        template <capability_type C, typename=enable_if_t<IS_SUPERSET_OF(C,SPINORBITAL)>>
        OneElectronOperator(const tensor::TensorInitializerList<C>& ilist)
        : ab(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as<>().name+"(AB)") <<
                                          tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {1,0}, {1,0}))),
          ij(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as<>().name+"(AB)") <<
                                          tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {0,1}, {0,1}))),
          ai(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as<>().name+"(AB)") <<
                                          tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {1,0}, {0,1}))),
          ia(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as<>().name+"(AB)") <<
                                          tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {0,1}, {1,0}))) {}

        OneElectronOperator(OneElectronOperator& other, int copy)
        : ab(copy&AB ? other.getAB().construct() : other.getAB()),
          ij(copy&IJ ? other.getIJ().construct() : other.getIJ()),
          ai(copy&AI ? other.getAI().construct() : other.getAI()),
          ia(copy&IA ? other.getIA().construct() : other.getIA())
        {
            if (copy&AB) getAB() = other.getAB();
            if (copy&IJ) getIJ() = other.getIJ();
            if (copy&AI) getAI() = other.getAI();
            if (copy&IA) getIA() = other.getIA();
        }

        OneElectronOperator(const string& name, OneElectronOperator& other, int copy)
        : ab(copy&AB ? other.getAB().construct(name+"(AB)") : other.getAB()),
          ij(copy&IJ ? other.getIJ().construct(name+"(IJ)") : other.getIJ()),
          ai(copy&AI ? other.getAI().construct(name+"(AI)") : other.getAI()),
          ia(copy&IA ? other.getIA().construct(name+"(IA)") : other.getIA())
        {
            if (copy&AB) getAB() = other.getAB();
            if (copy&IJ) getIJ() = other.getIJ();
            if (copy&AI) getAI() = other.getAI();
            if (copy&IA) getIA() = other.getIA();
        }

        OneElectronOperator(const OneElectronOperator& other)
        : ab(other.getAB().construct()),
          ij(other.getIJ().construct()),
          ai(other.getAI().construct()),
          ia(other.getIA().construct())
        {
            getAB() = other.getAB();
            getIJ() = other.getIJ();
            getAI() = other.getAI();
            getIA() = other.getIA();
        }

        OneElectronOperator(const string& name, const OneElectronOperator& other)
        : ab(other.getAB().construct(name+"(AB)")),
          ij(other.getIJ().construct(name+"(IJ)")),
          ai(other.getAI().construct(name+"(AI)")),
          ia(other.getIA().construct(name+"(IA)"))
        {
            getAB() = other.getAB();
            getIJ() = other.getIJ();
            getAI() = other.getAI();
            getIA() = other.getIA();
        }

        tensor::Tensor<SPINORBITAL|PGSYMMETRIC>& getAB() { return ab; }
        tensor::Tensor<SPINORBITAL|PGSYMMETRIC>& getIJ() { return ij; }
        tensor::Tensor<SPINORBITAL|PGSYMMETRIC>& getAI() { return ai; }
        tensor::Tensor<SPINORBITAL|PGSYMMETRIC>& getIA() { return ia; }

        const tensor::Tensor<SPINORBITAL|PGSYMMETRIC>& getAB() const { return ab; }
        const tensor::Tensor<SPINORBITAL|PGSYMMETRIC>& getIJ() const { return ij; }
        const tensor::Tensor<SPINORBITAL|PGSYMMETRIC>& getAI() const { return ai; }
        const tensor::Tensor<SPINORBITAL|PGSYMMETRIC>& getIA() const { return ia; }
};

}
}

#endif
