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

#ifndef _AQUAQRIUS_TENSOR_HPP_
#define _AQUAQRIUS_TENSOR_HPP_

#include <ostream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstdio>
#include <stdint.h>
#include <cstring>

namespace aquarius
{
namespace tensor
{

template <class Derived> class IndexableTensor;
template <class Derived> class IndexedTensor;
template <class Derived> class IndexedTensorMult;

class TensorError;
class OutOfBoundsError;
class LengthMismatchError;
class IndexMismatchError;
class InvalidNdimError;
class InvalidLengthError;
class InvalidLdError;
class LdTooSmallError;
class SymmetryMismatchError;
class InvalidSymmetryError;
class InvalidStartError;

#define INHERIT_FROM_INDEXABLE_TENSOR(Derived) \
    friend class aquarius::tensor::IndexableTensor< Derived >; \
    protected: \
        using aquarius::tensor::IndexableTensor< Derived >::ndim_; \
    public: \
        using aquarius::tensor::IndexableTensor< Derived >::mult; \
        using aquarius::tensor::IndexableTensor< Derived >::sum; \
        using aquarius::tensor::IndexableTensor< Derived >::scale; \
        using aquarius::tensor::IndexableTensor< Derived >::operator=; \
        using aquarius::tensor::IndexableTensor< Derived >::operator+=; \
        using aquarius::tensor::IndexableTensor< Derived >::operator-=; \
        Derived & operator=(const Derived & other) \
        { \
            static_cast< aquarius::tensor::IndexableTensor< Derived >* >(this)->operator=(other); \
            return *this; \
        } \
    private:

template <class Derived>
class IndexableTensor
{
    protected:
        int ndim_;

    public:
        IndexableTensor(const int ndim = 0)
        : ndim_(ndim) {}

        virtual ~IndexableTensor() {}

        int getDimension() const { return ndim_; }

        /**********************************************************************
         *
         * Explicit indexing operations
         *
         *********************************************************************/
        IndexedTensor<Derived> operator[](const char* idx)
        {
            return IndexedTensor<Derived>(*(Derived*)this, idx);
        }

        const IndexedTensor<Derived> operator[](const char* idx) const
        {
            return IndexedTensor<Derived>(*(Derived*)this, idx);
        }

        /**********************************************************************
         *
         * Operators with scalars
         *
         *********************************************************************/
        virtual Derived& operator=(const double val) = 0;

        Derived& operator*=(const double val)
        {
            scale(val);
            return static_cast<Derived&>(*this);
        }

        IndexedTensor<Derived> operator*(const double factor)
        {
            std::vector<char> idx(ndim_+1);
            for (int i = 0;i < ndim_;i++) idx[i] = (char)(i+1);
            idx[ndim_] = 0;

            return IndexedTensor<Derived>(*(Derived*)this, idx.data())*factor;
        }

        const IndexedTensor<Derived> operator*(const double factor) const
        {
            std::vector<char> idx(ndim_+1);
            for (int i = 0;i < ndim_;i++) idx[i] = (char)(i+1);
            idx[ndim_] = 0;

            return IndexedTensor<Derived>(*(Derived*)this, idx.data())*factor;
        }

        /**********************************************************************
         *
         * Implicitly indexed binary operations (inner product only)
         *
         *********************************************************************/
        Derived& operator=(const IndexedTensorMult<Derived>& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != 0) throw InvalidNdimError();
            if (strcmp(other.A_.idx_, other.B_.idx_) != 0) throw IndexMismatchError();
            #endif //VALIDATE_INPUTS

            (*this)[""] = other;

            return static_cast<Derived&>(*this);
        }

        Derived& operator+=(const IndexedTensorMult<Derived>& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != 0) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            (*this)[""] += other;

            return static_cast<Derived&>(*this);
        }

        Derived& operator-=(const IndexedTensorMult<Derived>& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != 0) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            (*this)[""] -= other;

            return static_cast<Derived&>(*this);
        }

        /**********************************************************************
         *
         * Implicitly indexed unary operations (assignment and summation)
         *
         *********************************************************************/
        Derived& operator=(const Derived& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            std::vector<char> idx(ndim_+1);
            for (int i = 0;i < ndim_;i++) idx[i] = (char)(i+1);
            idx[ndim_] = 0;

            (*this)[idx.data()] = other[idx.data()];

            return static_cast<Derived&>(*this);
        }

        Derived& operator+=(const Derived& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            std::vector<char> idx(ndim_+1);
            for (int i = 0;i < ndim_;i++) idx[i] = (char)(i+1);
            idx[ndim_] = 0;

            (*this)[idx.data()] += other[idx.data()];

            return static_cast<Derived&>(*this);
        }

        Derived& operator-=(const Derived& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            std::vector<char> idx(ndim_+1);
            for (int i = 0;i < ndim_;i++) idx[i] = (char)(i+1);
            idx[ndim_] = 0;

            (*this)[idx.data()] -= other[idx.data()];

            return static_cast<Derived&>(*this);
        }

        Derived& operator=(const IndexedTensor<Derived>& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.tensor_.ndim_) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            std::vector<char> idx(ndim_+1);
            for (int i = 0;i < ndim_;i++) idx[i] = (char)(i+1);
            idx[ndim_] = 0;

            #ifdef VALIDATE_INPUTS
            if (strcmp(idx.data(), other.idx_) != 0) throw IndexMismatchError();
            #endif //VALIDATE_INPUTS

            (*this)[idx.data()] = other;

            return static_cast<Derived&>(*this);
        }

        Derived& operator+=(const IndexedTensor<Derived>& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.tensor_.ndim_) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            std::vector<char> idx(ndim_+1);
            for (int i = 0;i < ndim_;i++) idx[i] = (char)(i+1);
            idx[ndim_] = 0;

            #ifdef VALIDATE_INPUTS
            if (strcmp(idx.data(), other.idx_) != 0) throw IndexMismatchError();
            #endif //VALIDATE_INPUTS

            (*this)[idx.data()] += other;

            return static_cast<Derived&>(*this);
        }

        Derived& operator-=(const IndexedTensor<Derived>& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.tensor_.ndim_) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            std::vector<char> idx(ndim_+1);
            for (int i = 0;i < ndim_;i++) idx[i] = (char)(i+1);
            idx[ndim_] = 0;

            #ifdef VALIDATE_INPUTS
            if (strcmp(idx.data(), other.idx_) != 0) throw IndexMismatchError();
            #endif //VALIDATE_INPUTS

            (*this)[idx.data()] -= other;

            return static_cast<Derived&>(*this);
        }

        /**********************************************************************
         *
         * Implicit dot product
         *
         *********************************************************************/
        IndexedTensorMult<Derived> operator*(const Derived& other) const
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            char* idx = new char[ndim_+1];
            for (int i = 0;i < ndim_;i++) idx[i] = (char)(i+1);
            idx[ndim_] = 0;

            IndexedTensorMult<Derived> itm(IndexedTensor<Derived>(*(Derived*)this, idx),
                                           IndexedTensor<Derived>(other, idx));

            delete[] idx;

            return itm;
        }

        /**********************************************************************
         *
         * Binary tensor operations (multiplication)
         *
         *********************************************************************/
        void mult(const double alpha, const Derived& A, const char* idx_A,
                                      const Derived& B, const char* idx_B,
                  const double beta,                    const char* idx_C)
        {
            int *idx_A_ = new int[A.ndim_];
            int *idx_B_ = new int[B.ndim_];
            int *idx_C_ = new int[ndim_];

            for (int i = 0;i < A.ndim_;i++) idx_A_[i] = idx_A[i];
            for (int i = 0;i < B.ndim_;i++) idx_B_[i] = idx_B[i];
            for (int i = 0;i < ndim_;i++)   idx_C_[i] = idx_C[i];

            mult(alpha, A, idx_A_,
                        B, idx_B_,
                  beta,    idx_C_);

            delete[] idx_A_;
            delete[] idx_B_;
            delete[] idx_C_;
        }

        virtual void mult(const double alpha, const Derived& A, const int* idx_A,
                                              const Derived& B, const int* idx_B,
                          const double beta,                    const int* idx_C) = 0;


        /**********************************************************************
         *
         * Unary tensor operations (summation)
         *
         *********************************************************************/
        void sum(const double alpha, const Derived& A, const char* idx_A,
                 const double beta,                    const char* idx_B)
        {
            int *idx_A_ = new int[A.ndim_];
            int *idx_B_ = new int[ndim_];

            for (int i = 0;i < A.ndim_;i++) idx_A_[i] = idx_A[i];
            for (int i = 0;i < ndim_;i++) idx_B_[i] = idx_B[i];

            sum(alpha, A, idx_A_,
                 beta,    idx_B_);

            delete[] idx_A_;
            delete[] idx_B_;
        }

        virtual void sum(const double alpha, const Derived& A, const int* idx_A,
                         const double beta,                    const int* idx_B) = 0;


        /**********************************************************************
         *
         * Scalar operations
         *
         *********************************************************************/
        void scale(const double alpha)
        {
            int* idx = new int[ndim_];
            for (int i = 0;i < ndim_;i++) idx[i] = i;

            scale(alpha, idx);

            delete[] idx;
        }

        void scale(const double alpha, const char* idx_A)
        {
            int *idx_A_ = new int[ndim_];

            for (int i = 0;i < ndim_;i++) idx_A_[i] = idx_A[i];

            scale(alpha, idx_A_);

            delete[] idx_A_;
        }

        virtual void scale(const double alpha, const int* idx_A) = 0;
};

/**************************************************************************
 *
 * Tensor to scalar operations
 *
 *************************************************************************/
inline double scalar(const double other) { return other; }

template<class Derived>
double scalar(const IndexedTensor<Derived>& other);

template<class Derived>
double scalar(const IndexedTensorMult<Derived>& other);

template <class Derived>
class IndexedTensor
{
    friend class IndexableTensor<Derived>;
    friend class IndexedTensorMult<Derived>;

    /**************************************************************************
     *
     * Reversed scalar operations
     *
     *************************************************************************/
    friend IndexedTensor<Derived> operator*(const double factor, const IndexedTensor<Derived>& other)
    {
        return other*factor;
    }

    friend IndexedTensor<Derived> operator*(const double factor, Derived& d)
    {
        std::vector<char> idx(d.getDimension()+1);
        for (int i = 0;i < d.getDimension();i++) idx[i] = (char)(i+1);
        idx[d.getDimension()] = 0;

        return IndexedTensor<Derived>(d, idx.data())*factor;
    }

    friend const IndexedTensor<Derived> operator*(const double factor, const Derived& d)
    {
        std::vector<char> idx(d.getDimension()+1);
        for (int i = 0;i < d.getDimension();i++) idx[i] = (char)(i+1);
        idx[d.getDimension()] = 0;

        return IndexedTensor<Derived>(d, idx.data())*factor;
    }

    public:
        double factor_;
        char* idx_;
        Derived& tensor_;

        IndexedTensor(const Derived& tensor, const char* idx)
        : factor_(1.0), tensor_(const_cast<Derived&>(tensor))
        {
            if (strlen(idx) != tensor.getDimension()) throw InvalidNdimError();
            idx_ = new char[tensor.getDimension()+1];
            strcpy(idx_, idx);
        }

        IndexedTensor(const IndexedTensor<Derived>& other)
        : factor_(other.factor_), tensor_(const_cast<Derived&>(other.tensor_))
        {
            idx_ = new char[tensor_.getDimension()+1];
            strcpy(idx_, other.idx_);
        }

        ~IndexedTensor()
        {
            delete[] idx_;
        }

        /**********************************************************************
         *
         * Unary negation
         *
         *********************************************************************/
        IndexedTensor<Derived>& operator-()
        {
            factor_ = -factor_;
            return *this;
        }

        /**********************************************************************
         *
         * Unary tensor operations (summation)
         *
         *********************************************************************/
        IndexedTensor<Derived>& operator=(const IndexedTensor<Derived>& other)
        {
            tensor_.sum(other.factor_, other.tensor_, other.idx_, 0.0, idx_);
            return *this;
        }

        IndexedTensor<Derived>& operator+=(const IndexedTensor<Derived>& other)
        {
            tensor_.sum(other.factor_, other.tensor_, other.idx_, factor_, idx_);
            return *this;
        }

        IndexedTensor<Derived>& operator-=(const IndexedTensor<Derived>& other)
        {
            tensor_.sum(-other.factor_, other.tensor_, other.idx_, factor_, idx_);
            return *this;
        }

        IndexedTensor<Derived>& operator=(const IndexedTensorMult<Derived>& other)
        {
            tensor_.mult(other.factor_*other.A_.factor_*other.B_.factor_, other.A_.tensor_, other.A_.idx_,
                                                                          other.B_.tensor_, other.B_.idx_,
                         0.0,                                                               idx_);
            return *this;
        }

        /**********************************************************************
         *
         * Binary tensor operations (multiplication)
         *
         *********************************************************************/
        IndexedTensor<Derived>& operator+=(const IndexedTensorMult<Derived>& other)
        {
            tensor_.mult(other.factor_*other.A_.factor_*other.B_.factor_, other.A_.tensor_, other.A_.idx_,
                                                                          other.B_.tensor_, other.B_.idx_,
                         factor_,                                                           idx_);
            return *this;
        }

        IndexedTensor<Derived>& operator-=(const IndexedTensorMult<Derived>& other)
        {
            tensor_.mult(-other.factor_*other.A_.factor_*other.B_.factor_, other.A_.tensor_, other.A_.idx_,
                                                                           other.B_.tensor_, other.B_.idx_,
                         factor_,                                                            idx_);
            return *this;
        }

        IndexedTensorMult<Derived> operator*(const IndexedTensor<Derived>& other) const
        {
            return IndexedTensorMult<Derived>(*this, other);
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        IndexedTensor<Derived> operator*(const double factor) const
        {
            IndexedTensor<Derived> it(*this);
            it.factor_ *= factor;
            return it;
        }

        IndexedTensor<Derived>& operator*=(const double factor)
        {
            tensor_.scale(factor, idx_);
            return *this;
        }
};

template <class Derived>
class IndexedTensorMult
{
    friend class IndexableTensor<Derived>;
    friend class IndexedTensor<Derived>;

    /**************************************************************************
     *
     * Reversed scalar operations
     *
     *************************************************************************/
    friend IndexedTensorMult<Derived> operator*(const double factor, const IndexedTensorMult<Derived>& other)
    {
        return other*factor;
    }

    private:
        const IndexedTensorMult& operator=(const IndexedTensorMult& other);

    public:
        double factor_;
        const IndexedTensor<Derived> A_;
        const IndexedTensor<Derived> B_;

        IndexedTensorMult(const IndexedTensor<Derived>& A, const IndexedTensor<Derived>& B)
        : factor_(1.0), A_(A), B_(B) {}

        IndexedTensorMult(const IndexedTensorMult<Derived>& other)
        : factor_(other.factor_), A_(other.A_), B_(other.B_) {}

        /**********************************************************************
         *
         * Unary negation
         *
         *********************************************************************/
        IndexedTensorMult<Derived>& operator-()
        {
            factor_ = -factor_;
            return *this;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        IndexedTensorMult<Derived> operator*(const double factor) const
        {
            IndexedTensorMult<Derived> itm(*this);
            itm.factor_ *= factor;
            return itm;
        }
};

class TensorError : public std::exception
{
    public:
        virtual const char* what() const throw() = 0;
};

class OutOfBoundsError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "out-of-bounds read or write"; }
};

class LengthMismatchError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "length mismatch error"; }
};

class IndexMismatchError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "index mismatch error"; }
};

class InvalidNdimError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid number of dimensions"; }
};

class InvalidLengthError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid length"; }
};

class InvalidLdError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid leading dimension"; }
};

class LdTooSmallError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "leading dimension is too small"; }
};

class SymmetryMismatchError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "symmetry mismatch error"; }
};

class InvalidSymmetryError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid symmetry value"; }
};

class InvalidStartError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid start value"; }
};

}
}

#endif
