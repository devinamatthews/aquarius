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
        using aquarius::tensor::LocalTensor< Derived, T >::len; \
        using aquarius::tensor::LocalTensor< Derived, T >::ld; \
        using aquarius::tensor::LocalTensor< Derived, T >::size; \
        using aquarius::tensor::LocalTensor< Derived, T >::data; \
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

template <class Derived, class T>
class LocalTensor : public IndexableTensor<Derived,T>
{
    INHERIT_FROM_INDEXABLE_TENSOR(Derived,T)

    protected:
        std::vector<int> len;
        std::vector<int> ld;
        uint64_t size;
        T* data;
        bool isAlloced;

    public:
        enum CopyType {CLONE, REFERENCE, REPLACE};

        LocalTensor(const T val = (T)0)
        : IndexableTensor<Derived,T>(0), len(0), ld(0), size(1)
        {
            data = SAFE_MALLOC(T, 1);
            isAlloced = true;

            data[0] = val;
        }

        LocalTensor(const Derived& A)
        : IndexableTensor<Derived,T>(A.ndim), len(A.len), ld(A.ld), size(A.size)
        {
            data = SAFE_MALLOC(T, size);
            std::copy(A.data, A.data+size, data);
            isAlloced = true;
        }

        LocalTensor(Derived& A, const CopyType type=CLONE)
        : IndexableTensor<Derived,T>(A.ndim), len(A.len), ld(A.ld), size(A.size)
        {
            switch (type)
            {
                case CLONE:
                    data = SAFE_MALLOC(T, size);
                    std::copy(A.data, A.data+size, data);
                    isAlloced = true;
                    break;
                case REFERENCE:
                    data = A.data;
                    isAlloced = false;
                    break;
                case REPLACE:
                    data = A.data;
                    isAlloced = A.isAlloced;
                    A.isAlloced = false;
                    break;
            }
        }

        LocalTensor(const int ndim, const std::vector<int>& len, const std::vector<int>& ld_, uint64_t size, T* data, bool zero=false)
        : IndexableTensor<Derived,T>(ndim), len(len), ld(ld_), size(size)
        {
            assert(len.size() == ndim);

            if (ld.size() != ndim)
            {
                ld.resize(ndim);
                ld[0] = 1;
                for (int i = 1;i < ndim;i++) ld[i] = ld[i-1]*len[i-1];
            }

            #ifdef VALIDATE_INPUTS
            if (validate_tensor(ndim,len.data(),ld.data(),NULL) != TENSOR_SUCCESS)
            {
                throw std::runtime_error("not a valid tensor");
            }
            #endif //VALIDATE_INPUTS

            data = data;
            isAlloced = false;
            if (zero) std::fill(data, data+size, (T)0);
        }

        LocalTensor(int ndim, const std::vector<int>& len, const std::vector<int>& ld_, uint64_t size, bool zero=true)
        : IndexableTensor<Derived,T>(ndim), len(len), ld(ld_), size(size)
        {
            assert(len.size() == ndim);

            if (ld.size() != ndim)
            {
                ld.resize(ndim);
                ld[0] = 1;
                for (int i = 1;i < ndim;i++) ld[i] = ld[i-1]*len[i-1];
            }

            #ifdef VALIDATE_INPUTS
            if (validate_tensor(ndim,len.data(),ld.data(),NULL) != TENSOR_SUCCESS)
            {
                throw std::runtime_error("not a valid tensor");
            }
            #endif //VALIDATE_INPUTS

            data = SAFE_MALLOC(T, size);
            isAlloced = true;
            if (zero) std::fill(data, data+size, (T)0);
        }

        ~LocalTensor()
        {
            if (isAlloced) FREE(data);
        }

        const std::vector<int>& getLengths() const { return len; }

        const std::vector<int>& getLeadingDims() const { return ld; }

        uint64_t getSize() const { return size; }

        void div(const T alpha, bool conja, const Derived& A,
                                bool conjb, const Derived& B, const T beta)
        {
            using std::conj;

            assert(size == A.size && size == B.size);

            if (conja)
            {
                if (conjb)
                {
                    for (uint64_t i = 0;i < size;i++)
                    {
                        if (std::abs(B.data[i]) > DBL_MIN)
                        {
                            data[i] = beta*data[i] + alpha*conj(A.data[i])/conj(B.data[i]);
                        }
                    }
                }
                else
                {
                    for (uint64_t i = 0;i < size;i++)
                    {
                        if (std::abs(B.data[i]) > DBL_MIN)
                        {
                            data[i] = beta*data[i] + alpha*conj(A.data[i])/B.data[i];
                        }
                    }
                }
            }
            else
            {
                if (conjb)
                {
                    for (uint64_t i = 0;i < size;i++)
                    {
                        if (std::abs(B.data[i]) > DBL_MIN)
                        {
                            data[i] = beta*data[i] + alpha*A.data[i]/conj(B.data[i]);
                        }
                    }
                }
                else
                {
                    for (uint64_t i = 0;i < size;i++)
                    {
                        if (std::abs(B.data[i]) > DBL_MIN)
                        {
                            data[i] = beta*data[i] + alpha*A.data[i]/B.data[i];
                        }
                    }
                }
            }
        }

        void invert(const T alpha, bool conja, const Derived& A, const T beta)
        {
            using std::conj;

            assert(size == A.size);

            if (conja)
            {
                for (uint64_t i = 0;i < size;i++)
                {
                    if (std::abs(A.data[i]) > DBL_MIN)
                    {
                        data[i] = beta*data[i] + alpha/conj(A.data[i]);
                    }
                }
            }
            else
            {
                for (uint64_t i = 0;i < size;i++)
                {
                    if (std::abs(A.data[i]) > DBL_MIN)
                    {
                        data[i] = beta*data[i] + alpha/A.data[i];
                    }
                }
            }
        }

        virtual void print(FILE* fp) const = 0;

        virtual void print(std::ostream& stream) const = 0;

        T* getData() { return data; };

        const T* getData() const { return data; };

        T dot(bool conja, const Derived& A, const int* idx_A,
              bool conjb,                   const int* idx_B) const
        {
            Derived dt;
            dt.mult(1, conja,            A, idx_A,
                       conjb, getDerived(), idx_B,
                    0,                       NULL);
            return dt.getData()[0];
        }
};

}
}

#endif
