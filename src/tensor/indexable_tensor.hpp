#ifndef _AQUARIUS_INDEXABLE_TENSOR_HPP_
#define _AQUARIUS_INDEXABLE_TENSOR_HPP_

#include "util/global.hpp"

#include "tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <class Derived, class T> class IndexableTensor;
template <class Derived, class T> class IndexedTensor;
template <class Derived, class T> class IndexedTensorMult;

#define INHERIT_FROM_INDEXABLE_TENSOR(Derived,T) \
    protected: \
        using aquarius::tensor::IndexableTensor< Derived, T >::ndim; \
    public: \
        using aquarius::tensor::IndexableTensor< Derived, T >::mult; \
        using aquarius::tensor::IndexableTensor< Derived, T >::sum; \
        using aquarius::tensor::IndexableTensor< Derived, T >::scale; \
        using aquarius::tensor::IndexableTensor< Derived, T >::dot; \
        using aquarius::tensor::IndexableTensor< Derived, T >::operator=; \
        using aquarius::tensor::IndexableTensor< Derived, T >::operator+=; \
        using aquarius::tensor::IndexableTensor< Derived, T >::operator-=; \
        using aquarius::tensor::IndexableTensor< Derived, T >::operator[]; \
        using aquarius::tensor::Tensor< Derived,T >::getDerived; \
        using aquarius::tensor::Tensor< Derived,T >::operator*=; \
        using aquarius::tensor::Tensor< Derived,T >::operator/=; \
        using aquarius::tensor::Tensor< Derived,T >::operator*; \
        using aquarius::tensor::Tensor< Derived,T >::operator/; \
        Derived & operator=(const Derived & other) \
        { \
            sum((T)1, false, other, (T)0); \
            return *this; \
        } \
    private:

template <class Derived, typename T>
class IndexableTensorBase
{
    protected:
        int ndim;

    public:
        IndexableTensorBase(const int ndim = 0) : ndim(ndim) {}

        virtual ~IndexableTensorBase() {}

        Derived& getDerived() { return static_cast<Derived&>(*this); }

        const Derived& getDerived() const { return static_cast<const Derived&>(*this); }

        int getDimension() const { return ndim; }

        string implicit() const
        {
            string inds(ndim, ' ');
            for (int i = 0;i < ndim;i++) inds[i] = (char)('A'+i);
            return inds;
        }

        /**********************************************************************
         *
         * Explicit indexing operations
         *
         *********************************************************************/
        IndexedTensor<Derived,T> operator[](const string& idx)
        {
            return IndexedTensor<Derived,T>(getDerived(), idx);
        }

        IndexedTensor<const Derived,T> operator[](const string& idx) const
        {
            return IndexedTensor<const Derived,T>(getDerived(), idx);
        }

        /**********************************************************************
         *
         * Implicitly indexed binary operations (inner product, trace, and weighting)
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator=(const IndexedTensorMult<cvDerived,T>& other)
        {
            (*this)[implicit()] = other;
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator+=(const IndexedTensorMult<cvDerived,T>& other)
        {
            (*this)[implicit()] += other;
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator-=(const IndexedTensorMult<cvDerived,T>& other)
        {
            (*this)[implicit()] -= other;
            return getDerived();
        }

        /**********************************************************************
         *
         * Implicitly indexed unary operations (assignment and summation)
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator=(const IndexedTensor<cvDerived,T>& other)
        {
            (*this)[implicit()] = other;
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator+=(const IndexedTensor<cvDerived,T>& other)
        {
            (*this)[implicit()] += other;
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator-=(const IndexedTensor<cvDerived,T>& other)
        {
            (*this)[implicit()] -= other;
            return getDerived();
        }

        /**********************************************************************
         *
         * Binary tensor operations (multiplication)
         *
         *********************************************************************/
        virtual void mult(const T alpha,  bool conja, const Derived& A, const string& idx_A,
                                          bool conjb, const Derived& B, const string& idx_B,
                          const T beta,                                 const string& idx_C) = 0;


        /**********************************************************************
         *
         * Unary tensor operations (summation)
         *
         *********************************************************************/
        virtual void sum(const T alpha, bool conja, const Derived& A, const string& idx_A,
                         const T beta,                                const string& idx_B) = 0;


        /**********************************************************************
         *
         * Scalar operations
         *
         *********************************************************************/
        virtual void scale(const T alpha, const string& idx_A) = 0;

        virtual T dot(bool conja, const Derived& A, const string& idx_A,
                      bool conjb,                   const string& idx_B) const = 0;
};

template <class Derived, typename T>
class IndexableTensor : public IndexableTensorBase<Derived,T>, public Tensor<Derived,T>
{
    INHERIT_FROM_TENSOR(Derived,T)

    protected:
        using IndexableTensorBase<Derived,T>::ndim;

    public:
        using IndexableTensorBase<Derived,T>::scale;
        using IndexableTensorBase<Derived,T>::dot;
        using IndexableTensorBase<Derived,T>::mult;
        using IndexableTensorBase<Derived,T>::sum;
        using IndexableTensorBase<Derived,T>::implicit;

        IndexableTensor(const string& name, const int ndim = 0)
        : IndexableTensorBase<Derived,T>(ndim), Tensor<Derived,T>(name) {}

        virtual ~IndexableTensor() {}

        /**********************************************************************
         *
         * Binary tensor operations (multiplication)
         *
         *********************************************************************/
        void mult(const T alpha)
        {
            scale(alpha);
        }

        void mult(const T alpha, bool conja, const Derived& A,
                                 bool conjb, const Derived& B,
                  const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim != A.getDimension() || ndim != B_.getDimension()) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            mult(alpha, conja, A, A.implicit(),
                        conjb, B, B.implicit(),
                  beta,             implicit());
        }

        /**********************************************************************
         *
         * Unary tensor operations (summation)
         *
         *********************************************************************/
        void sum(const T alpha, const T beta)
        {
            Derived tensor("alpha", getDerived(), alpha);
            beta*(*this)[implicit()] = tensor[""];
        }

        void sum(const T alpha, bool conja, const Derived& A, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim != A.getDimension()) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            sum(alpha, conja, A, A.implicit(),
                 beta,             implicit());
        }

        /**********************************************************************
         *
         * Scalar operations
         *
         *********************************************************************/
        void scale(const T alpha)
        {
            scale(alpha, implicit());
        }

        T dot(bool conja, const Derived& A, bool conjb) const
        {
            #ifdef VALIDATE_INPUTS
            if (ndim != A.getDimension()) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            return dot(conja, A, A.implicit(),
                       conjb,      implicit());
        }
};

template <class Derived, typename T>
class IndexedTensor
{
    public:
        Derived& tensor_;
        string idx_;
        T factor_;
        bool conj_;

        template <typename cvDerived>
        IndexedTensor(const IndexedTensor<cvDerived,T>& other)
        : tensor_(other.tensor_), idx_(other.idx_), factor_(other.factor_), conj_(other.conj_) {}

        IndexedTensor(Derived& tensor, const string& idx, const T factor=(T)1, const bool conj=false)
        : tensor_(tensor), idx_(idx), factor_(factor), conj_(conj)
        {
            if (idx.size() != tensor.getDimension()) throw InvalidNdimError();
        }

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        IndexedTensor<Derived,T> operator-() const
        {
            IndexedTensor<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return ret;
        }

        friend IndexedTensor<const Derived,T> conj(const IndexedTensor<Derived,T>& other)
        {
            IndexedTensor<Derived,T> ret(other);
            ret.conj_ = !ret.conj_;
            return ret;
        }

        /**********************************************************************
         *
         * Unary tensor operations (summation)
         *
         *********************************************************************/
        IndexedTensor<Derived,T>& operator=(const IndexedTensor<Derived,T>& other)
        {
            tensor_.sum(other.factor_, other.conj_, other.tensor_, other.idx_, (T)0, idx_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexedTensor<Derived,T>&))
        operator=(const IndexedTensor<cvDerived,T>& other)
        {
            tensor_.sum(other.factor_, other.conj_, other.tensor_, other.idx_, (T)0, idx_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexedTensor<Derived,T>&))
        operator+=(const IndexedTensor<cvDerived,T>& other)
        {
            tensor_.sum(other.factor_, other.conj_, other.tensor_, other.idx_, factor_, idx_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexedTensor<Derived,T>&))
        operator-=(const IndexedTensor<cvDerived,T>& other)
        {
            tensor_.sum(-other.factor_, other.conj_, other.tensor_, other.idx_, factor_, idx_);
            return *this;
        }

        /**********************************************************************
         *
         * Binary tensor operations (multiplication)
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexedTensor<Derived,T>&))
        operator=(const IndexedTensorMult<cvDerived,T>& other)
        {
            tensor_.mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.A_.idx_,
                                        other.B_.conj_, other.B_.tensor_, other.B_.idx_,
                                  (T)0,                                            idx_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexedTensor<Derived,T>&))
        operator+=(const IndexedTensorMult<cvDerived,T>& other)
        {
            tensor_.mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.A_.idx_,
                                        other.B_.conj_, other.B_.tensor_, other.B_.idx_,
                               factor_,                                            idx_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexedTensor<Derived,T>&))
        operator-=(const IndexedTensorMult<cvDerived,T>& other)
        {
            tensor_.mult(-other.factor_, other.A_.conj_, other.A_.tensor_, other.A_.idx_,
                                         other.B_.conj_, other.B_.tensor_, other.B_.idx_,
                                factor_,                                            idx_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexedTensorMult<Derived,T>))
        operator*(const IndexedTensor<cvDerived,T>& other) const
        {
            return IndexedTensorMult<Derived,T>(*this, other);
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexedTensorMult<Derived,T>))
        operator*(const ScaledTensor<cvDerived,T>& other) const
        {
            cvDerived& B = other.tensor_.getDerived();

            if (other.conj_)
            {
                return IndexedTensorMult<Derived,T>(*this, B[B.implicit()]*other.factor_);
            }
            else
            {
                return IndexedTensorMult<Derived,T>(*this, conj(B[B.implicit()])*other.factor_);
            }
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexedTensorMult<Derived,T>))
        operator*(const IndexableTensor<cvDerived,T>& other) const
        {
            return IndexedTensorMult<Derived,T>(*this, other[other.implicit()]);
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        IndexedTensor<Derived,T> operator*(const T factor) const
        {
            IndexedTensor<Derived,T> it(*this);
            it.factor_ *= factor;
            return it;
        }

        friend IndexedTensor<Derived,T> operator*(const T factor, const IndexedTensor<Derived,T>& other)
        {
            return other*factor;
        }

        IndexedTensor<Derived,T>& operator*=(const T factor)
        {
            tensor_.scale(factor, idx_);
            return *this;
        }

        IndexedTensor<Derived,T>& operator=(const T val)
        {
            Derived tensor("alpha", tensor_, val);
            *this = tensor[""];
            return *this;
        }

        IndexedTensor<Derived,T>& operator+=(const T val)
        {
            Derived tensor("alpha", tensor_, val);
            *this += tensor[""];
            return *this;
        }

        IndexedTensor<Derived,T>& operator-=(const T val)
        {
            Derived tensor("alpha", tensor_, val);
            *this -= tensor[""];
            return *this;
        }
};

template <class Derived1, class Derived2, class T>
//enable_if_t<is_same<const Derived1, const Derived2>::value,IndexedTensorMult<Derived1,T>>
IndexedTensorMult<Derived1,T>
operator*(const IndexableTensorBase<Derived1,T>& t1, const IndexedTensor<Derived2,T>& t2)
{
    return IndexedTensorMult<Derived1,T>(t1[t1.implicit()], t2);
}

template <class Derived1, class Derived2, class T>
//enable_if_t<is_same<const Derived1, const Derived2>::value,IndexedTensorMult<Derived1,T>>
IndexedTensorMult<Derived1,T>
operator*(const ScaledTensor<Derived1,T>& t1, const IndexedTensor<Derived2,T>& t2)
{
    Derived1& A = t1.tensor_.getDerived();

    if (t1.conj_)
    {
        return IndexedTensorMult<Derived1,T>(conj(A[A.implicit()])*t1.factor_, t2);
    }
    else
    {
        return IndexedTensorMult<Derived1,T>(A[A.implicit()]*t1.factor_, t2);
    }
}

template <class Derived, typename T>
class IndexedTensorMult
{
    private:
        const IndexedTensorMult& operator=(const IndexedTensorMult<Derived,T>& other);

    public:
        IndexedTensor<const Derived,T> A_;
        IndexedTensor<const Derived,T> B_;
        T factor_;

        template <class Derived1, class Derived2>
        IndexedTensorMult(const IndexedTensor<Derived1,T>& A, const IndexedTensor<Derived2,T>& B)
        : A_(A), B_(B), factor_(A.factor_*B.factor_) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        IndexedTensorMult<Derived,T> operator-() const
        {
            IndexedTensorMult<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return ret;
        }

        friend IndexedTensorMult<Derived,T> conj(const IndexedTensorMult<Derived,T>& other)
        {
            IndexedTensorMult<Derived,T> ret(other);
            ret.A_.conj_ = !ret.A_.conj_;
            ret.B_.conj_ = !ret.B_.conj_;
            return ret;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        IndexedTensorMult<Derived,T> operator*(const T factor) const
        {
            IndexedTensorMult<Derived,T> ret(*this);
            ret.factor_ *= factor;
            return ret;
        }

        IndexedTensorMult<Derived,T> operator/(const T factor) const
        {
            IndexedTensorMult<Derived,T> ret(*this);
            ret.factor_ /= factor;
            return ret;
        }

        friend IndexedTensorMult<Derived,T> operator*(const T factor, const IndexedTensorMult<Derived,T>& other)
        {
            return other*factor;
        }
};

}

/**************************************************************************
 *
 * Tensor to scalar operations
 *
 *************************************************************************/
template <class Derived, typename T>
T scalar(const tensor::IndexedTensorMult<Derived,T>& itm)
{
    return itm.factor_*itm.B_.tensor_.dot(itm.A_.conj_, itm.A_.tensor_, itm.A_.idx_,
                                          itm.B_.conj_,                 itm.B_.idx_);
}

}

#endif
