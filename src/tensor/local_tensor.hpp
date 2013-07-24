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

#ifndef _AQUARIUS_TENSOR_LOCAL_TENSOR_HPP_
#define _AQUARIUS_TENSOR_LOCAL_TENSOR_HPP_

#include <ostream>
#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cstring>
#include <cassert>
#include <string>
#include <algorithm>
#include <cfloat>

#include "memory/memory.h"

#include "tensor.h"
#include "indexable_tensor.hpp"

namespace aquarius
{
namespace tensor
{

#define INHERIT_FROM_LOCAL_TENSOR(Derived, T) \
    protected: \
        using aquarius::tensor::LocalTensor< Derived, T >::len_; \
        using aquarius::tensor::LocalTensor< Derived, T >::ld_; \
        using aquarius::tensor::LocalTensor< Derived, T >::size_; \
        using aquarius::tensor::LocalTensor< Derived, T >::data_; \
    public: \
        using aquarius::tensor::LocalTensor< Derived, T >::CLONE; \
        using aquarius::tensor::LocalTensor< Derived, T >::REFERENCE; \
        using aquarius::tensor::LocalTensor< Derived, T >::REPLACE; \
        using aquarius::tensor::LocalTensor< Derived, T >::getSize; \
    INHERIT_FROM_INDEXABLE_TENSOR(Derived, T) \
    friend class aquarius::tensor::LocalTensor< Derived, T >;

#define CHECK_RETURN_VALUE(ret) \
switch (ret) \
{ \
    case TENSOR_OUT_OF_BOUNDS: \
        throw OutOfBoundsError(); \
        break; \
    case TENSOR_LENGTH_MISMATCH: \
        throw LengthMismatchError(); \
        break; \
    case TENSOR_INDEX_MISMATCH: \
        throw IndexMismatchError(); \
        break; \
    case TENSOR_INVALID_NDIM: \
        throw InvalidNdimError(); \
        break; \
    case TENSOR_INVALID_LENGTH: \
        throw InvalidLengthError(); \
        break; \
    case TENSOR_INVALID_LD: \
        throw InvalidLdError(); \
        break; \
    case TENSOR_LD_TOO_SMALL: \
        throw LdTooSmallError(); \
        break; \
    case TENSOR_SYMMETRY_MISMATCH: \
        throw SymmetryMismatchError(); \
        break; \
    case TENSOR_INVALID_SYMMETRY: \
        throw InvalidSymmetryError(); \
        break; \
    case TENSOR_INVALID_START: \
        throw InvalidStartError(); \
        break; \
    default: \
        break; \
}

#define VALIDATE_TENSOR_THROW(ndim,len,ld,sym) CHECK_RETURN_VALUE(validate_tensor(ndim,len,ld,sym));

template <class Derived, class T=double>
class LocalTensor : public IndexableTensor<Derived,T>
{
    INHERIT_FROM_INDEXABLE_TENSOR(Derived,T)

    protected:
        int* len_;
        int* ld_;
        uint64_t size_;
        T* data_;
        bool isAlloced;

    public:
        enum CopyType {CLONE, REFERENCE, REPLACE};

        LocalTensor(const LocalTensor& t, const T val)
        : IndexableTensor<Derived,T>(0), size_(1)
        {
            len_ = SAFE_MALLOC(int, ndim_);
            ld_ = NULL;
            data_ = SAFE_MALLOC(T, size_);
            isAlloced = true;

            data_[0] = val;
        }

        LocalTensor(const Derived& A, const CopyType type=CLONE)
        : IndexableTensor<Derived,T>(A.ndim_), size_(A.size_)
        {
            len_ = SAFE_MALLOC(int, ndim_);
            std::copy(A.len_, A.len_+ndim_, len_);

            if (A.ld_ == NULL)
            {
                ld_ = NULL;
            }
            else
            {
                ld_ = SAFE_MALLOC(int, ndim_);
                std::copy(A.ld_, A.ld_+ndim_, ld_);
            }

            switch (type)
            {
                case CLONE:
                    data_ = SAFE_MALLOC(T, size_);
                    std::copy(A.data_, A.data_+size_, data_);
                    isAlloced = true;
                    break;
                case REFERENCE:
                    data_ = A.data_;
                    isAlloced = false;
                    break;
                case REPLACE:
                    data_ = A.data_;
                    isAlloced = A.isAlloced;
                    const_cast<Derived&>(A).isAlloced = false;
                    break;
            }
        }

        LocalTensor(const int ndim, const int *len, const int *ld, uint64_t size, T* data, const bool zero = false)
        : IndexableTensor<Derived,T>(ndim), size_(size)
        {
            #ifdef VALIDATE_INPUTS
            if (validate_tensor(ndim,len,ld,NULL) != TENSOR_SUCCESS)
            {
                throw std::runtime_error("not a valid tensor");
            }
            #endif //VALIDATE_INPUTS

            len_ = SAFE_MALLOC(int, ndim_);
            std::copy(len, len+ndim_, len_);

            if (ld == NULL)
            {
                ld_ = NULL;
            }
            else
            {
                ld_ = SAFE_MALLOC(int, ndim_);
                std::copy(ld, ld+ndim_, ld_);
            }

            data_ = data;
            isAlloced = false;
            if (zero) std::fill(data_, data_+size_, (T)0);
        }

        LocalTensor(int ndim, const int *len, const int *ld, uint64_t size, bool zero = true)
        : IndexableTensor<Derived,T>(ndim), size_(size)
        {
            #ifdef VALIDATE_INPUTS
            if (validate_tensor(ndim,len,ld,NULL) != TENSOR_SUCCESS)
            {
                throw std::runtime_error("not a valid tensor");
            }
            #endif //VALIDATE_INPUTS

            len_ = SAFE_MALLOC(int, ndim_);
            std::copy(len, len+ndim_, len_);

            if (ld == NULL)
            {
                ld_ = NULL;
            }
            else
            {
                ld_ = SAFE_MALLOC(int, ndim_);
                std::copy(ld, ld+ndim_, ld_);
            }

            data_ = SAFE_MALLOC(T, size_);
            isAlloced = true;
            if (zero) std::fill(data_, data_+size_, (T)0);
        }

        ~LocalTensor()
        {
            FREE(len_);
            if (ld_ != NULL) FREE(ld_);
            if (isAlloced) FREE(data_);
        }

        const int* getLengths() const { return len_; }

        const int* getLeadingDims() const { return ld_; }

        uint64_t getSize() const { return size_; }

        void div(const T alpha, bool conja, const Derived& A,
                                bool conjb, const Derived& B, const T beta)
        {
            assert(size_ == A.size_ && size_ == B.size_);

            if (conja)
            {
                if (conjb)
                {
                    for (uint64_t i = 0;i < size_;i++)
                    {
                        if (std::abs(B.data_[i]) > DBL_MIN)
                        {
                            data_[i] = beta*data_[i] + alpha*conj(A.data_[i])/conj(B.data_[i]);
                        }
                    }
                }
                else
                {
                    for (uint64_t i = 0;i < size_;i++)
                    {
                        if (std::abs(B.data_[i]) > DBL_MIN)
                        {
                            data_[i] = beta*data_[i] + alpha*conj(A.data_[i])/B.data_[i];
                        }
                    }
                }
            }
            else
            {
                if (conjb)
                {
                    for (uint64_t i = 0;i < size_;i++)
                    {
                        if (std::abs(B.data_[i]) > DBL_MIN)
                        {
                            data_[i] = beta*data_[i] + alpha*A.data_[i]/conj(B.data_[i]);
                        }
                    }
                }
                else
                {
                    for (uint64_t i = 0;i < size_;i++)
                    {
                        if (std::abs(B.data_[i]) > DBL_MIN)
                        {
                            data_[i] = beta*data_[i] + alpha*A.data_[i]/B.data_[i];
                        }
                    }
                }
            }
        }

        void invert(const T alpha, bool conja, const Derived& A, const T beta)
        {
            assert(size_ == A.size_);

            if (conja)
            {
                for (uint64_t i = 0;i < size_;i++)
                {
                    if (std::abs(A.data_[i]) > DBL_MIN)
                    {
                        data_[i] = beta*data_[i] + alpha/conj(A.data_[i]);
                    }
                }
            }
            else
            {
                for (uint64_t i = 0;i < size_;i++)
                {
                    if (std::abs(A.data_[i]) > DBL_MIN)
                    {
                        data_[i] = beta*data_[i] + alpha/A.data_[i];
                    }
                }
            }
        }

        virtual void print(FILE* fp) const = 0;

        virtual void print(std::ostream& stream) const = 0;

        T* getData() { return data_; };

        const T* getData() const { return data_; };

        T dot(bool conja, const Derived& A, const int* idx_A,
              bool conjb,                   const int* idx_B) const
        {
            Derived dt(A, 0);
            dt.mult(1, conja,            A, idx_A,
                       conjb, getDerived(), idx_B,
                    0,                       NULL);
            return dt.getData()[0];
        }
};

}
}

#endif
