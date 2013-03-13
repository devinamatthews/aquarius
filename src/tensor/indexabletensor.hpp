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

#ifndef _AQUARIUS_INDEXABLE_TENSOR_HPP_
#define _AQUARIUS_INDEXABLE_TENSOR_HPP_

#include <vector>
#include <string>
#include <algorithm>

#include "stl_ext/stl_ext.hpp"

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
        using aquarius::tensor::IndexableTensor< Derived, T >::ndim_; \
    public: \
        using aquarius::tensor::IndexableTensor< Derived, T >::mult; \
        using aquarius::tensor::IndexableTensor< Derived, T >::sum; \
        using aquarius::tensor::IndexableTensor< Derived, T >::scale; \
        using aquarius::tensor::IndexableTensor< Derived, T >::operator=; \
        using aquarius::tensor::IndexableTensor< Derived, T >::operator+=; \
        using aquarius::tensor::IndexableTensor< Derived, T >::operator-=; \
    INHERIT_FROM_TENSOR(CONCAT(Derived),T)

template <class Derived, typename T>
class IndexableTensor : virtual public Tensor<Derived,T>
{
    protected:
        using Tensor<Derived,T>::derived;
        int ndim_;

    public:
        IndexableTensor(const int ndim = 0)
        : Tensor<Derived,T>(*this), ndim_(ndim) {}

        virtual ~IndexableTensor() {}

        int getDimension() const { return ndim_; }

        std::string implicit() const
        {
            std::string inds(ndim_, ' ');
            for (int i = 0;i < ndim_;i++) inds[i] = (char)(i+1);
            return inds;
        }

        /**********************************************************************
         *
         * Explicit indexing operations
         *
         *********************************************************************/
        IndexedTensor<Derived,T> operator[](const std::string& idx)
        {
            return IndexedTensor<Derived,T>(derived, idx);
        }

        IndexedTensor<const Derived,T> operator[](const std::string& idx) const
        {
            return IndexedTensor<const Derived,T>(derived, idx);
        }

        /**********************************************************************
         *
         * Implicitly indexed binary operations (inner product, trace, and weighting)
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexableTensor<Derived,T>&))
        operator=(const IndexedTensorMult<cvDerived,T>& other)
        {
            (*this)[implicit()] = other;
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexableTensor<Derived,T>&))
        operator+=(const IndexedTensorMult<cvDerived,T>& other)
        {
            (*this)[implicit()] += other;
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexableTensor<Derived,T>&))
        operator-=(const IndexedTensorMult<cvDerived,T>& other)
        {
            (*this)[implicit()] -= other;
            return *this;
        }

        /**********************************************************************
         *
         * Implicitly indexed unary operations (assignment and summation)
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexableTensor<Derived,T>&))
        operator=(const IndexedTensor<cvDerived,T>& other)
        {
            (*this)[implicit()] = other;
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexableTensor<Derived,T>&))
        operator+=(const IndexedTensor<cvDerived,T>& other)
        {
            (*this)[implicit()] += other;
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(IndexableTensor<Derived,T>&))
        operator-=(const IndexedTensor<cvDerived,T>& other)
        {
            (*this)[implicit()] -= other;
            return *this;
        }

        /**********************************************************************
         *
         * Binary tensor operations (multiplication)
         *
         *********************************************************************/
        void mult(const T alpha)
        {
            scale(alpha);
        }

        void mult(const T alpha, bool conja, const Tensor<Derived,T>& A_,
                                 bool conjb, const Tensor<Derived,T>& B_,
                  const T beta)
        {
            const Derived& A = A_.getDerived();
            const Derived& B = B_.getDerived();

            #ifdef VALIDATE_INPUTS
            if (ndim_ != A.ndim_ || ndim_ != B_.ndim_) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            mult(alpha, conja, A, A.implicit(),
                        conjb, B, B.implicit(),
                  beta,             implicit());
        }

        void mult(const T alpha, bool conja, const IndexableTensor<Derived,T>& A, const std::string& idx_A,
                                 bool conjb, const IndexableTensor<Derived,T>& B, const std::string& idx_B,
                  const T beta,                                                   const std::string& idx_C)
        {
            std::vector<int> idx_A_(A.ndim_);
            std::vector<int> idx_B_(B.ndim_);
            std::vector<int> idx_C_(ndim_);

            for (int i = 0;i < A.ndim_;i++) idx_A_[i] = idx_A[i];
            for (int i = 0;i < B.ndim_;i++) idx_B_[i] = idx_B[i];
            for (int i = 0;i < ndim_;i++)   idx_C_[i] = idx_C[i];

            mult(alpha, conja, A, idx_A_.data(),
                        conjb, B, idx_B_.data(),
                  beta,           idx_C_.data());
        }

        virtual void mult(const T alpha,  bool conja, const IndexableTensor<Derived,T>& A, const int* idx_A,
                                          bool conjb, const IndexableTensor<Derived,T>& B, const int* idx_B,
                          const T beta,                                                    const int* idx_C) = 0;


        /**********************************************************************
         *
         * Unary tensor operations (summation)
         *
         *********************************************************************/
        void sum(const T alpha, const T beta)
        {
            Derived tensor(derived, alpha);
            beta*(*this)[implicit()] = tensor[""];
        }

        void sum(const T alpha, bool conja, const Tensor<Derived,T>& A_, const T beta)
        {
            const Derived& A = A_.getDerived();

            #ifdef VALIDATE_INPUTS
            if (ndim_ != A.ndim_) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            sum(alpha, conja, A, A.implicit(),
                 beta,             implicit());
        }

        void sum(const T alpha, bool conja, const IndexableTensor<Derived,T>& A, const std::string& idx_A,
                 const T beta,                                                   const std::string& idx_B)
        {
            std::vector<int> idx_A_(A.ndim_);
            std::vector<int> idx_B_(ndim_);

            for (int i = 0;i < A.ndim_;i++) idx_A_[i] = idx_A[i];
            for (int i = 0;i < ndim_;i++) idx_B_[i] = idx_B[i];

            sum(alpha, conja, A, idx_A_.data(),
                 beta,           idx_B_.data());
        }

        virtual void sum(const T alpha, bool conja, const IndexableTensor<Derived,T>& A, const int* idx_A,
                         const T beta,                                                   const int* idx_B) = 0;


        /**********************************************************************
         *
         * Scalar operations
         *
         *********************************************************************/
        void scale(const T alpha)
        {
            scale(alpha, implicit());
        }

        void scale(const T alpha, const std::string& idx_A)
        {
            std::vector<int> idx_A_(ndim_);
            for (int i = 0;i < ndim_;i++) idx_A_[i] = idx_A[i];
            scale(alpha, idx_A_.data());
        }

        virtual void scale(const T alpha, const int* idx_A) = 0;
};

/**************************************************************************
 *
 * Tensor to scalar operations
 *
 *************************************************************************/
template <class Derived, class foo=void> struct Scalar {};

inline float scalar(float val) { return val; }
inline double scalar(double val) { return val; }
inline scomplex scalar(scomplex val) { return val; }
inline dcomplex scalar(dcomplex val) { return val; }

template <class Derived, typename T>
struct Scalar<TensorMult<Derived,T>, typename std::enable_if<std::is_base_of<IndexableTensor<Derived,T>,Derived>::value>::type>
{
    static T value(const TensorMult<Derived,T>& other)
    {
        const Derived& A = other.A_.tensor_.getDerived();
        const Derived& B = other.B_.tensor_.getDerived();

        #ifdef VALIDATE_INPUTS
        if (A.getDimension() != B.getDimension()) throw InvalidNdimError();
        #endif //VALIDATE_INPUTS

        std::string implicit(A.getDimension(),' ');
        for (int i = 0;i < A.getDimension();i++) implicit[i] = (char)(i+1);

        if (other.A_.conj_)
        {
            if (other.B_.conj_)
            {
                return scalar(conj(A[implicit])*
                              conj(B[implicit]))*other.factor_;
            }
            else
            {
                return scalar(conj(A[implicit])*
                                   B[implicit])*other.factor_;
            }
        }
        else
        {
            if (other.B_.conj_)
            {
                return scalar(     A[implicit]*
                              conj(B[implicit]))*other.factor_;
            }
            else
            {
                return scalar(A[implicit]*
                              B[implicit])*other.factor_;
            }
        }
    }
};

template <class Derived, class T>
T scalar(const TensorMult<Derived,T>& d)
{
    return Scalar< TensorMult<Derived,T> >::value(d);
}

template <class Derived, class T>
T scalar(const IndexedTensorMult<Derived,T>& d)
{
    return Scalar< IndexedTensorMult<Derived,T> >::value(d);
}

template <class Derived, typename T>
class IndexedTensor
{
    public:
        Derived& tensor_;
        std::string idx_;
        T factor_;
        bool conj_;

        template <typename cvDerived>
        IndexedTensor(const IndexedTensor<cvDerived,T>& other)
        : tensor_(other.tensor_), idx_(other.idx_), factor_(other.factor_), conj_(other.conj_) {}

        IndexedTensor(Derived& tensor, const std::string& idx, const T factor=(T)1, const bool conj=false)
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
            Derived tensor(tensor_, val);
            *this = tensor[""];
            return *this;
        }

        IndexedTensor<Derived,T>& operator+=(const T val)
        {
            Derived tensor(tensor_, val);
            *this += tensor[""];
            return *this;
        }

        IndexedTensor<Derived,T>& operator-=(const T val)
        {
            Derived tensor(tensor_, val);
            *this -= tensor[""];
            return *this;
        }
};

template <class Derived1, class Derived2, class T>
typename std::enable_if<std::is_same<const Derived1, const Derived2>::value,IndexedTensorMult<Derived1,T> >::type
operator*(const IndexableTensor<Derived1,T>& t1, const IndexedTensor<Derived2,T>& t2)
{
    return IndexedTensorMult<Derived1,T>(t1[t1.implicit()], t2);
}

template <class Derived1, class Derived2, class T>
typename std::enable_if<std::is_same<const Derived1, const Derived2>::value,IndexedTensorMult<Derived1,T> >::type
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
        : A_(A), B_(B), factor_(A_.factor_*B_.factor_)  {}

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
}

#endif
