#ifndef _AQUARIUS_OPERATOR_2EOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_2EOPERATOR_HPP_

#include "util/global.hpp"

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
        tensor::SpinorbitalTensor<T>& aijk;
        tensor::SpinorbitalTensor<T>& ijak;
        tensor::SpinorbitalTensor<T>& abij;
        tensor::SpinorbitalTensor<T>& ijab;
        tensor::SpinorbitalTensor<T>& aibj;
        tensor::SpinorbitalTensor<T>& aibc;
        tensor::SpinorbitalTensor<T>& abci;
        tensor::SpinorbitalTensor<T>& abcd;

    public:
        enum
        {
            IJKL = 0x0010,
            AIJK = 0x0020,
            IJAK = 0x0040,
            ABIJ = 0x0080,
            IJAB = 0x0100,
            AIBJ = 0x0200,
            AIBC = 0x0400,
            ABCI = 0x0800,
            ABCD = 0x1000
        };

        TwoElectronOperator(const string& name, const Arena& arena, const Space& occ, const Space& vrt);

        TwoElectronOperator(const string& name, OneElectronOperator<T>& other, int copy);

        TwoElectronOperator(const OneElectronOperator<T>& other);

        TwoElectronOperator(const string& name, const OneElectronOperator<T>& other);

        TwoElectronOperator(const string& name, TwoElectronOperator<T>& other, int copy);

        TwoElectronOperator(const TwoElectronOperator<T>& other);

        TwoElectronOperator(const string& name, const TwoElectronOperator<T>& other);

        T dot(bool conja, const TwoElectronOperator<T>& A, bool conjb) const;

        tensor::SpinorbitalTensor<T>& getIJKL() { return ijkl; }
        tensor::SpinorbitalTensor<T>& getAIJK() { return aijk; }
        tensor::SpinorbitalTensor<T>& getIJAK() { return ijak; }
        tensor::SpinorbitalTensor<T>& getABIJ() { return abij; }
        tensor::SpinorbitalTensor<T>& getIJAB() { return ijab; }
        tensor::SpinorbitalTensor<T>& getAIBJ() { return aibj; }
        tensor::SpinorbitalTensor<T>& getAIBC() { return aibc; }
        tensor::SpinorbitalTensor<T>& getABCI() { return abci; }
        tensor::SpinorbitalTensor<T>& getABCD() { return abcd; }

        const tensor::SpinorbitalTensor<T>& getIJKL() const { return ijkl; }
        const tensor::SpinorbitalTensor<T>& getAIJK() const { return aijk; }
        const tensor::SpinorbitalTensor<T>& getIJAK() const { return ijak; }
        const tensor::SpinorbitalTensor<T>& getABIJ() const { return abij; }
        const tensor::SpinorbitalTensor<T>& getIJAB() const { return ijab; }
        const tensor::SpinorbitalTensor<T>& getAIBJ() const { return aibj; }
        const tensor::SpinorbitalTensor<T>& getAIBC() const { return aibc; }
        const tensor::SpinorbitalTensor<T>& getABCI() const { return abci; }
        const tensor::SpinorbitalTensor<T>& getABCD() const { return abcd; }
};

}
}

#endif
