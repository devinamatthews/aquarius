#ifndef _AQUARIUS_FRAMEWORKS_OPERATOR_2EOPERATOR_HPP_
#define _AQUARIUS_FRAMEWORKS_OPERATOR_2EOPERATOR_HPP_

#include "frameworks/util.hpp"
#include "frameworks/tensor.hpp"

#include "1eoperator.hpp"

namespace aquarius
{
namespace op
{

class TwoElectronOperator : public OneElectronOperator
{
    protected:
        tensor::Tensor<SPINORBITAL> abcd;
        tensor::Tensor<SPINORBITAL> abci;
        tensor::Tensor<SPINORBITAL> aibc;
        tensor::Tensor<SPINORBITAL> abij;
        tensor::Tensor<SPINORBITAL> ijab;
        tensor::Tensor<SPINORBITAL> aibj;
        tensor::Tensor<SPINORBITAL> aijk;
        tensor::Tensor<SPINORBITAL> ijak;
        tensor::Tensor<SPINORBITAL> ijkl;

    public:
        enum
        {
            ABCD = 0x0010,
            ABCI = 0x0020,
            AIBC = 0x0040,
            ABIJ = 0x0080,
            IJAB = 0x0100,
            AIBJ = 0x0200,
            AIJK = 0x0400,
            IJAK = 0x0800,
            IJKL = 0x1000
        };

        template <capability_type C, typename=enable_if_t<IS_SUPERSET_OF(C,SPINORBITAL)>>
        TwoElectronOperator(const string& name, const tensor::TensorInitializerList<C>& ilist)
        : abcd(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(name+"(ABCD)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {2,0}, {2,0}))),
          abci(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(name+"(ABCI)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {2,0}, {1,1}))),
          aibc(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(name+"(AIBC)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {1,1}, {2,0}))),
          abij(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(name+"(ABIJ)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {2,0}, {0,2}))),
          ijab(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(name+"(IJAB)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {0,2}, {2,0}))),
          aibj(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(name+"(AIBJ)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {1,1}, {1,1}))),
          aijk(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(name+"(AIJK)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {1,1}, {0,2}))),
          ijak(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(name+"(IJAK)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {0,2}, {1,1}))),
          ijkl(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(name+"(IJKL)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {0,2}, {0,2}))) {}

        template <capability_type C, typename=enable_if_t<IS_SUPERSET_OF(C,SPINORBITAL)>>
        TwoElectronOperator(const tensor::TensorInitializerList<C>& ilist)
        : abcd(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as().name+"(ABCD)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {2,0}, {2,0}))),
          abci(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as().name+"(ABCI)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {2,0}, {1,1}))),
          aibc(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as().name+"(AIBC)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {1,1}, {2,0}))),
          abij(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as().name+"(ABIJ)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {2,0}, {0,2}))),
          ijab(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as().name+"(IJAB)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {0,2}, {2,0}))),
          aibj(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as().name+"(AIBJ)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {1,1}, {1,1}))),
          aijk(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as().name+"(AIJK)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {1,1}, {0,2}))),
          ijak(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as().name+"(IJAK)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {0,2}, {1,1}))),
          ijkl(tensor::Tensor<C>::construct(ilist << tensor::TensorInitializer<>(ilist.as().name+"(IJKL)") <<
                                            tensor::TensorInitializer<C&(SPINORBITAL|PGSYMMETRIC)>(ilist, {0,2}, {0,2}))) {}

        TwoElectronOperator(TwoElectronOperator& other, int copy)
        : OneElectronOperator(other, copy),
          abcd(copy&ABCD ? other.getABCD().construct() : other.getABCD()),
          abci(copy&ABCI ? other.getABCI().construct() : other.getABCI()),
          aibc(copy&AIBC ? other.getAIBC().construct() : other.getAIBC()),
          abij(copy&ABIJ ? other.getABIJ().construct() : other.getABIJ()),
          ijab(copy&IJAB ? other.getIJAB().construct() : other.getIJAB()),
          aibj(copy&AIBJ ? other.getAIBJ().construct() : other.getAIBJ()),
          aijk(copy&AIJK ? other.getAIJK().construct() : other.getAIJK()),
          ijak(copy&IJAK ? other.getIJAK().construct() : other.getIJAK()),
          ijkl(copy&IJKL ? other.getIJKL().construct() : other.getIJKL())
        {
            if (copy&ABCD) getABCD() = other.getABCD();
            if (copy&ABCI) getABCI() = other.getABCI();
            if (copy&AIBC) getAIBC() = other.getAIBC();
            if (copy&ABIJ) getABIJ() = other.getABIJ();
            if (copy&IJAB) getIJAB() = other.getIJAB();
            if (copy&AIBJ) getAIBJ() = other.getAIBJ();
            if (copy&AIJK) getAIJK() = other.getAIJK();
            if (copy&IJAK) getIJAK() = other.getIJAK();
            if (copy&IJKL) getIJKL() = other.getIJKL();
        }

        TwoElectronOperator(const string& name, TwoElectronOperator& other, int copy)
        : OneElectronOperator(name, other, copy),
          abcd(copy&ABCD ? other.getABCD().construct(name+"(ABCD)") : other.getABCD()),
          abci(copy&ABCI ? other.getABCI().construct(name+"(ABCI)") : other.getABCI()),
          aibc(copy&AIBC ? other.getAIBC().construct(name+"(AIBC)") : other.getAIBC()),
          abij(copy&ABIJ ? other.getABIJ().construct(name+"(ABIJ)") : other.getABIJ()),
          ijab(copy&IJAB ? other.getIJAB().construct(name+"(IJAB)") : other.getIJAB()),
          aibj(copy&AIBJ ? other.getAIBJ().construct(name+"(AIBJ)") : other.getAIBJ()),
          aijk(copy&AIJK ? other.getAIJK().construct(name+"(AIJK)") : other.getAIJK()),
          ijak(copy&IJAK ? other.getIJAK().construct(name+"(IJAK)") : other.getIJAK()),
          ijkl(copy&IJKL ? other.getIJKL().construct(name+"(IJKL)") : other.getIJKL())
        {
            if (copy&ABCD) getABCD() = other.getABCD();
            if (copy&ABCI) getABCI() = other.getABCI();
            if (copy&AIBC) getAIBC() = other.getAIBC();
            if (copy&ABIJ) getABIJ() = other.getABIJ();
            if (copy&IJAB) getIJAB() = other.getIJAB();
            if (copy&AIBJ) getAIBJ() = other.getAIBJ();
            if (copy&AIJK) getAIJK() = other.getAIJK();
            if (copy&IJAK) getIJAK() = other.getIJAK();
            if (copy&IJKL) getIJKL() = other.getIJKL();
        }

        TwoElectronOperator(const TwoElectronOperator& other)
        : OneElectronOperator(other),
          abcd(other.getABCD().construct()),
          abci(other.getABCI().construct()),
          aibc(other.getAIBC().construct()),
          abij(other.getABIJ().construct()),
          ijab(other.getIJAB().construct()),
          aibj(other.getAIBJ().construct()),
          aijk(other.getAIJK().construct()),
          ijak(other.getIJAK().construct()),
          ijkl(other.getIJKL().construct())
        {
            getABCD() = other.getABCD();
            getABCI() = other.getABCI();
            getAIBC() = other.getAIBC();
            getABIJ() = other.getABIJ();
            getIJAB() = other.getIJAB();
            getAIBJ() = other.getAIBJ();
            getAIJK() = other.getAIJK();
            getIJAK() = other.getIJAK();
            getIJKL() = other.getIJKL();
        }

        TwoElectronOperator(const string& name, const TwoElectronOperator& other)
        : OneElectronOperator(name, other),
          abcd(other.getABCD().construct(name+"(ABCD)")),
          abci(other.getABCI().construct(name+"(ABCI)")),
          aibc(other.getAIBC().construct(name+"(AIBC)")),
          abij(other.getABIJ().construct(name+"(ABIJ)")),
          ijab(other.getIJAB().construct(name+"(IJAB)")),
          aibj(other.getAIBJ().construct(name+"(AIBJ)")),
          aijk(other.getAIJK().construct(name+"(AIJK)")),
          ijak(other.getIJAK().construct(name+"(IJAK)")),
          ijkl(other.getIJKL().construct(name+"(IJKL)"))
        {
            getABCD() = other.getABCD();
            getABCI() = other.getABCI();
            getAIBC() = other.getAIBC();
            getABIJ() = other.getABIJ();
            getIJAB() = other.getIJAB();
            getAIBJ() = other.getAIBJ();
            getAIJK() = other.getAIJK();
            getIJAK() = other.getIJAK();
            getIJKL() = other.getIJKL();
        }

        tensor::Tensor<SPINORBITAL>& getABCD() { return abcd; }
        tensor::Tensor<SPINORBITAL>& getABCI() { return abci; }
        tensor::Tensor<SPINORBITAL>& getAIBC() { return aibc; }
        tensor::Tensor<SPINORBITAL>& getABIJ() { return abij; }
        tensor::Tensor<SPINORBITAL>& getIJAB() { return ijab; }
        tensor::Tensor<SPINORBITAL>& getAIBJ() { return aibj; }
        tensor::Tensor<SPINORBITAL>& getAIJK() { return aijk; }
        tensor::Tensor<SPINORBITAL>& getIJAK() { return ijak; }
        tensor::Tensor<SPINORBITAL>& getIJKL() { return ijkl; }

        const tensor::Tensor<SPINORBITAL>& getABCD() const { return abcd; }
        const tensor::Tensor<SPINORBITAL>& getABCI() const { return abci; }
        const tensor::Tensor<SPINORBITAL>& getAIBC() const { return aibc; }
        const tensor::Tensor<SPINORBITAL>& getABIJ() const { return abij; }
        const tensor::Tensor<SPINORBITAL>& getIJAB() const { return ijab; }
        const tensor::Tensor<SPINORBITAL>& getAIBJ() const { return aibj; }
        const tensor::Tensor<SPINORBITAL>& getAIJK() const { return aijk; }
        const tensor::Tensor<SPINORBITAL>& getIJAK() const { return ijak; }
        const tensor::Tensor<SPINORBITAL>& getIJKL() const { return ijkl; }
};

}
}

#endif
