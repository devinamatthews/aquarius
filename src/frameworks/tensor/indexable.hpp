#ifndef _AQUARIUS_FRAMEWORKS_TENSOR_INDEXABLE_HPP_
#define _AQUARIUS_FRAMEWORKS_TENSOR_INDEXABLE_HPP_

#include "frameworks/util.hpp"
#include "frameworks/tensor.hpp"

namespace aquarius
{
namespace tensor
{

class ConstIndexedTensor;
class IndexedTensor;
class IndexedTensorMult;

namespace detail
{
    inline string implicit(int ndim)
    {
        string inds(ndim, ' ');
        for (int i = 0;i < ndim;i++) inds[i] = (char)('A'+i);
        return inds;
    }
}

template <> class TensorInitializer<INDEXABLE>
{
    public:
        int ndim;

        TensorInitializer(int ndim) : ndim(ndim) {}
};

TENSOR_INTERFACE(INDEXABLE)
{
    public:
        virtual int getDimension() const = 0;

        /**********************************************************************
         *
         * Binary tensor operations (multiplication)
         *
         *********************************************************************/
        virtual void mult(const Scalar& alpha, bool conja, const TensorImplementation<>& A, const string& idxA,
                                               bool conjb, const TensorImplementation<>& B, const string& idxB,
                          const Scalar& beta,                                               const string& idxC) = 0;

        /**********************************************************************
         *
         * Unary tensor operations (summation)
         *
         *********************************************************************/
        virtual void sum(const Scalar& alpha, bool conja, const TensorImplementation<>& A, const string& idxA,
                         const Scalar& beta,                                               const string& idxB) = 0;

        /**********************************************************************
         *
         * Scalar operations
         *
         *********************************************************************/
        virtual void scale(const Scalar& alpha, const string& idxA) = 0;

        virtual Scalar dot(bool conja, const TensorImplementation<>& A, const string& idxA,
                           bool conjb,                                  const string& idxB) const = 0;
};

TENSOR_DEFINITION(INDEXABLE)
{
    public:
        TensorDefinition(const INITIALIZER_TYPE(C)& ilist) : Base(ilist) {}

        int getDimension() const { return this->ndim; }

        void mult(const Scalar& alpha, bool conja, const TensorImplementation<>& A,
                                       bool conjb, const TensorImplementation<>& B, const Scalar& beta)
        {
            this->template as<INDEXABLE_>().mult(alpha, conja, A, detail::implicit(A.as<INDEXABLE_>().getDimension()),
                                                        conjb, B, detail::implicit(B.as<INDEXABLE_>().getDimension()),
                                                  beta,           detail::implicit(                   getDimension()));
        }

        void sum(const Scalar& alpha, bool conja, const TensorImplementation<>& A, const Scalar& beta)
        {
            this->template as<INDEXABLE_>().sum(alpha, conja, A, detail::implicit(A.as<INDEXABLE_>().getDimension()),
                                                 beta,           detail::implicit(                   getDimension()));
        }

        void scale(const Scalar& alpha)
        {
            this->template as<INDEXABLE_>().scale(alpha, detail::implicit(getDimension()));
        }

        Scalar dot(bool conja, const TensorImplementation<>& A, bool conjb) const
        {
            return this->template as<INDEXABLE_>().dot(conja, A, detail::implicit(A.as<INDEXABLE_>().getDimension()),
                                                       conjb,    detail::implicit(                   getDimension()));
        }
};

class IndexedTensorMult
{
    private:
        const IndexedTensorMult& operator=(const IndexedTensorMult& other);

    public:
        const TensorImplementation<>& A;
        const TensorImplementation<>& B;
        string idxa;
        string idxb;
        bool conja;
        bool conjb;
        Scalar factor;

        IndexedTensorMult(const TensorImplementation<>& A,
                          const TensorImplementation<>& B,
                          string idxa, string idxb, bool conja, bool conjb,
                          const Scalar& factor)
        : A(A), B(B), idxa(idxa), idxb(idxb), conja(conja), conjb(conjb), factor(factor) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/

        IndexedTensorMult operator-() const
        {
            IndexedTensorMult ret(*this);
            ret.factor = -ret.factor;
            return ret;
        }

        IndexedTensorMult conj() const
        {
            IndexedTensorMult ret(*this);
            ret.conja = !ret.conja;
            ret.conjb = !ret.conjb;
            return ret;
        }

        friend IndexedTensorMult conj(const IndexedTensorMult& other)
        {
            return other.conj();
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/

        IndexedTensorMult operator*(const Scalar& factor) const
        {
            IndexedTensorMult ret(*this);
            ret.factor *= factor;
            return ret;
        }

        IndexedTensorMult operator/(const Scalar& factor) const
        {
            IndexedTensorMult ret(*this);
            ret.factor /= factor;
            return ret;
        }

        friend IndexedTensorMult operator*(const Scalar& factor, const IndexedTensorMult& other)
        {
            return other*factor;
        }
};

class ConstIndexedTensor
{
    public:
        TensorImplementation<>& tensor;
        string idx;
        Scalar factor;
        bool conj_;

        ConstIndexedTensor(const TensorImplementation<>& tensor,
                           const string& idx,
                           Scalar factor = 1,
                           bool conj=false)
        : tensor(const_cast<TensorImplementation<>&>(tensor)),
          idx(idx), factor(factor), conj_(conj)
        {
            assert(idx.size() == tensor.as<INDEXABLE>().getDimension());
        }

        ConstIndexedTensor& operator=(const ConstIndexedTensor& other) = delete;

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/

        ConstIndexedTensor operator-() const
        {
            ConstIndexedTensor ret(*this);
            ret.factor = -ret.factor;
            return ret;
        }

        ConstIndexedTensor conj() const
        {
            ConstIndexedTensor ret(*this);
            ret.conj_ = !ret.conj_;
            return ret;
        }

        friend ConstIndexedTensor conj(const ConstIndexedTensor& other)
        {
            return other.conj();
        }

        /**********************************************************************
         *
         * Binary tensor operations (multiplication)
         *
         *********************************************************************/

        IndexedTensorMult operator*(const ConstIndexedTensor& other) const
        {
            return IndexedTensorMult(tensor,
                                     other.tensor,
                                     idx,
                                     other.idx,
                                     conj_,
                                     other.conj_,
                                     factor*other.factor);
        }

        template <capability_type C>
        IndexedTensorMult operator*(const ConstTensor<C>& other) const
        {
            static_assert(C&INDEXABLE, "The operand must be INDEXABLE.");
            return IndexedTensorMult(tensor,
                                     other.impl(),
                                     idx,
                                     detail::implicit(other.impl().template as<INDEXABLE>().getDimension()),
                                     conj_,
                                     false,
                                     factor);
        }

        IndexedTensorMult operator*(const ConstScaledTensor& other) const
        {
            assert(other.tensor.C&INDEXABLE);
            return IndexedTensorMult(tensor,
                                     other.tensor,
                                     idx,
                                     detail::implicit(other.tensor.template as<INDEXABLE>().getDimension()),
                                     conj_,
                                     other.conj_,
                                     factor*other.factor);
        }

        template <capability_type C>
        friend IndexedTensorMult operator*(const ConstTensor<C>& t1,
                                           const ConstIndexedTensor& t2)
        {
            static_assert(C&INDEXABLE, "The operand must be INDEXABLE.");
            return IndexedTensorMult(
                ConstIndexedTensor(t1.impl(),
                                   t2.tensor,
                                   detail::implicit(t1.impl().template as<INDEXABLE>().getDimension())),
                                   t2.idx,
                                   false,
                                   t2.conj_,
                                   t2.factor);
        }

        friend IndexedTensorMult operator*(const ConstScaledTensor& t1,
                                           const ConstIndexedTensor& t2)
        {
            assert(t1.tensor.C&INDEXABLE);
            return IndexedTensorMult(t1.tensor,
                                     t2.tensor,
                                     detail::implicit(t1.tensor.template as<INDEXABLE>().getDimension()),
                                     t2.idx,
                                     t1.conj_,
                                     t2.conj_,
                                     t1.factor*t2.factor);
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/

        ConstIndexedTensor operator*(const Scalar& factor) const
        {
            ConstIndexedTensor it(*this);
            it.factor *= factor;
            return it;
        }

        friend ConstIndexedTensor operator*(const Scalar& factor, const ConstIndexedTensor& other)
        {
            return other*factor;
        }
};

class IndexedTensor : public ConstIndexedTensor
{
    public:
        IndexedTensor(TensorImplementation<>& tensor,
                      const string& idx,
                      Scalar factor = 1,
                      bool conj=false)
        : ConstIndexedTensor(tensor, idx, factor, conj) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/

        using ConstIndexedTensor::operator-;
        using ConstIndexedTensor::conj;

        IndexedTensor operator-()
        {
            IndexedTensor ret(*this);
            ret.factor = -ret.factor;
            return ret;
        }

        IndexedTensor conj()
        {
            IndexedTensor ret(*this);
            ret.conj_ = !ret.conj_;
            return ret;
        }

        friend IndexedTensor conj(IndexedTensor& other)
        {
            return other.conj();
        }

        friend IndexedTensor conj(IndexedTensor&& other)
        {
            return other.conj();
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/

        using ConstIndexedTensor::operator*;

        IndexedTensor operator*(const Scalar& factor) const
        {
            IndexedTensor it(*this);
            it.factor *= factor;
            return it;
        }

        friend IndexedTensor operator*(const Scalar& factor, const IndexedTensor& other)
        {
            return other*factor;
        }

        /**********************************************************************
         *
         * Unary tensor operations (summation)
         *
         *********************************************************************/

        IndexedTensor& operator=(const IndexedTensor& other)
        {
            tensor.as<INDEXABLE>().sum(other.factor, conj_ != other.conj_, other.tensor, other.idx,
                                                  0,                                           idx);
            return *this;
        }

        IndexedTensor& operator=(const ConstIndexedTensor& other)
        {
            tensor.as<INDEXABLE>().sum(other.factor, conj_ != other.conj_, other.tensor, other.idx,
                                                  0,                                           idx);
            return *this;
        }

        IndexedTensor& operator+=(const ConstIndexedTensor& other)
        {
            tensor.as<INDEXABLE>().sum(other.factor, conj_ != other.conj_, other.tensor, other.idx,
                                             factor,                                           idx);
            return *this;
        }

        IndexedTensor& operator-=(const ConstIndexedTensor& other)
        {
            tensor.as<INDEXABLE>().sum(-other.factor, conj_ != other.conj_, other.tensor, other.idx,
                                              factor,                                           idx);
            return *this;
        }

        /**********************************************************************
         *
         * Binary tensor operations (multiplication)
         *
         *********************************************************************/

        IndexedTensor& operator=(const IndexedTensorMult& other)
        {
            tensor.as<INDEXABLE>().mult(other.factor, conj_ != other.conja, other.A, other.idxa,
                                                      conj_ != other.conjb, other.B, other.idxb,
                                                   0,                                      idx);
            return *this;
        }

        IndexedTensor& operator+=(const IndexedTensorMult& other)
        {
            tensor.as<INDEXABLE>().mult(other.factor, conj_ != other.conja, other.A, other.idxa,
                                                      conj_ != other.conjb, other.B, other.idxb,
                                              factor,                                      idx);
            return *this;
        }

        IndexedTensor& operator-=(const IndexedTensorMult& other)
        {
            tensor.as<INDEXABLE>().mult(-other.factor, conj_ != other.conja, other.A, other.idxa,
                                                       conj_ != other.conjb, other.B, other.idxb,
                                               factor,                                      idx);
            return *this;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/

        IndexedTensor& operator*=(const Scalar& factor)
        {
            tensor.as<INDEXABLE>().scale(factor, idx);
            return *this;
        }

        IndexedTensor& operator/=(const Scalar& factor)
        {
            tensor.as<INDEXABLE>().scale(1/factor, idx);
            return *this;
        }
};

TENSOR_WRAPPER(INDEXABLE)
{
    public:
        int getDimension() const
        {
            return this->template impl<INDEXABLE>().getDimension();
        }

        /**********************************************************************
         *
         * Explicit indexing operations
         *
         *********************************************************************/

        ConstIndexedTensor operator[](const string& idx) const
        {
            return ConstIndexedTensor(this->impl(), idx);
        }

        /**********************************************************************
         *
         * Explicit indexing operations
         *
         *********************************************************************/

        template <capability_type C_=C>
        enable_if_t<!(C_&CONST_),IndexedTensor> operator[](const string& idx)
        {
            return IndexedTensor(this->impl(), idx);
        }

        /**********************************************************************
         *
         * Implicitly indexed binary operations (inner product, trace, and weighting)
         *
         *********************************************************************/

        Tensor<C>& operator=(const IndexedTensorMult& other)
        {
            static_assert(!(C&CONST_), "The LHS must not be const.");
            (*this)[detail::implicit(this->getDimension())] = other;
            return *this;
        }

        Tensor<C>& operator+=(const IndexedTensorMult& other)
        {
            static_assert(!(C&CONST_), "The LHS must not be const.");
            (*this)[detail::implicit(this->getDimension())] += other;
            return *this;
        }

        Tensor<C>& operator-=(const IndexedTensorMult& other)
        {
            static_assert(!(C&CONST_), "The LHS must not be const.");
            (*this)[detail::implicit(this->getDimension())] -= other;
            return *this;
        }

        /**********************************************************************
         *
         * Implicitly indexed unary operations (assignment and summation)
         *
         *********************************************************************/

        Tensor<C>& operator=(const ConstIndexedTensor& other)
        {
            static_assert(!(C&CONST_), "The LHS must not be const.");
            (*this)[detail::implicit(this->getDimension())] = other;
            return *this;
        }

        Tensor<C>& operator+=(const ConstIndexedTensor& other)
        {
            static_assert(!(C&CONST_), "The LHS must not be const.");
            (*this)[detail::implicit(this->getDimension())] += other;
            return *this;
        }

        Tensor<C>& operator-=(const ConstIndexedTensor& other)
        {
            static_assert(!(C&CONST_), "The LHS must not be const.");
            (*this)[detail::implicit(this->getDimension())] -= other;
            return *this;
        }
};

/**************************************************************************
 *
 * Tensor to scalar operations
 *
 *************************************************************************/

inline Scalar scalar(const IndexedTensorMult& itm)
{
    return itm.factor*itm.B.as<INDEXABLE_>().dot(itm.conja, itm.A, itm.idxa,
                                                 itm.conjb,        itm.idxb);
}

}
}

#endif
