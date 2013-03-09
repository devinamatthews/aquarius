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

#ifndef _AQUARIUS_TENSOR_DENSE_TENSOR_HPP_
#define _AQUARIUS_TENSOR_DENSE_TENSOR_HPP_

#include <vector>
#include <algorithm>

#include "local_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <typename T> class PackedTensor;

template <typename T>
class DenseTensor : public LocalTensor< DenseTensor<T>,T >
{
    friend class PackedTensor<T>;
    friend class IndexedTensor< DenseTensor<T>,T >;

    INHERIT_FROM_LOCAL_TENSOR(DenseTensor<T>,T)

    public:
        DenseTensor(const DenseTensor& t, const T val)
        : LocalTensor< DenseTensor<T>,T >(0, (int*)NULL, (int*)NULL, 1, false)
        {
            data_[0] = val;
        }

        /*
        DenseTensor(const PackedTensor& A)
        : LocalTensor<DenseTensor>(A.ndim_, A.len_, (int*)NULL, A.size_, (bool)false)
        {
            CHECK_RETURN_VALUE(
            tensor_densify(&ndim_, len_, A.sym_));

            copy(A.data_, A.data_+size_, data_);
        }

        DenseTensor(const PackedTensor& A, const CopyType type)
        : LocalTensor<DenseTensor>(A.ndim_, A.len_, (int*)NULL, A.size_, (double*)NULL)
        {
            CHECK_RETURN_VALUE(
            tensor_densify(&ndim_, len_, A.sym_));

            switch (type)
            {
                case CLONE:
                    data_ = new double[size_];
                    copy(A.data_, A.data_+size_, data_);
                    isAlloced = true;
                    break;
                case REFERENCE:
                    data_ = A.data_;
                    isAlloced = false;
                    break;
                case REPLACE:
                    data_ = A.data_;
                    isAlloced = A.isAlloced;
                    const_cast<PackedTensor&>(A).isAlloced = false;
                    break;
            }
        }
        */

        DenseTensor(const DenseTensor<T>& A, const typename LocalTensor<DenseTensor<T>,T>::CopyType type=CLONE)
        : LocalTensor< DenseTensor<T>,T >(A, type) {}

        DenseTensor(const int ndim, const int *len, T* data, const bool zero=false)
        : LocalTensor< DenseTensor<T>,T >(ndim, len, NULL, getSize(ndim, len, NULL), data, zero) {}

        DenseTensor(const int ndim, const int *len, const bool zero=true)
        : LocalTensor< DenseTensor<T>,T >(ndim, len, NULL, getSize(ndim, len, NULL), zero) {}

        DenseTensor(const int ndim, const int *len, const int *ld, T* data, const bool zero=false)
        : LocalTensor< DenseTensor<T>,T >(ndim, len, ld, getSize(ndim, len, ld), data, zero) {}

        DenseTensor(const int ndim, const int *len, const int *ld, const bool zero=true)
        : LocalTensor< DenseTensor<T>,T >(ndim, len, ld, getSize(ndim, len, ld), zero) {}

        static uint64_t getSize(const int ndim, const int *len, const int *ld)
        {
            int64_t r = tensor_size_dense(ndim, len, ld);

            #ifdef VALIDATE_INPUTS
            CHECK_RETURN_VALUE(r);
            #endif //VALIDATE_INPUTS

            return r;
        }

        void print(FILE* fp) const
        {
            CHECK_RETURN_VALUE(
            tensor_print_dense(fp, data_, ndim_, len_, ld_));
        }

        void print(std::ostream& stream) const
        {
            #ifdef VALIDATE_INPUTS
            VALIDATE_TENSOR_THROW(ndim_, len_, ld_, NULL);
            #endif //VALIDATE_INPUTS

            std::vector<size_t> stride(ndim_);
            if (ld_ == NULL)
            {
                if (ndim_ > 0) stride[0] = 1;
                for (int i = 1;i < ndim_;i++) stride[i] = stride[i-1]*len_[i-1];
            }
            else
            {
                if (ndim_ > 0) stride[0] = ld_[0];
                for (int i = 1;i < ndim_;i++) stride[i] = stride[i-1]*ld_[i];
            }

            size_t size;
            if (ndim_ > 0)
            {
                size = stride[ndim_-1]*len_[ndim_-1];
            }
            else
            {
                size = 1;
            }

            size_t off = 0;

            /*
             * loop over elements in A
             */
            std::vector<int> pos(ndim_, 0);
            for (bool done = false;!done;)
            {
                #ifdef CHECK_BOUNDS
                if (off < 0 || off >= size) throw OutOfBoundsError();
                #endif //CHECK_BOUNDS

                for (int i = 0;i < ndim_;i++) stream << pos[i] << ' ';
                stream << std::scientific << std::setprecision(15) << data_[off] << '\n';

                for (int i = 0;i < ndim_;i++)
                {
                    if (pos[i] == len_[i]-1)
                    {
                        pos[i] = 0;
                        off -= stride[i]*(len_[i]-1);

                        if (i == ndim_-1)
                        {
                            done = true;
                            break;
                        }
                    }
                    else
                    {
                        pos[i]++;
                        off += stride[i];
                        break;
                    }
                }

                if (ndim_ == 0) done = true;
            }
            /*
             * end loop over A
             */
        }

        void mult(const T alpha, const DenseTensor<T>& A, const int* idx_A,
                                 const DenseTensor<T>& B, const int* idx_B,
                  const T beta,                           const int* idx_C)
        {
            CHECK_RETURN_VALUE(
            tensor_mult_dense_(alpha, A.data_, A.ndim_, A.len_, A.ld_, idx_A,
                                      B.data_, B.ndim_, B.len_, B.ld_, idx_B,
                               beta,    data_,   ndim_,   len_,   ld_, idx_C));
        }

        void sum(const T alpha, const DenseTensor<T>& A, const int* idx_A,
                 const T beta,                           const int* idx_B)
        {
            CHECK_RETURN_VALUE(
            tensor_sum_dense_(alpha, A.data_, A.ndim_, A.len_, A.ld_, idx_A,
                              beta,    data_,   ndim_,   len_,   ld_, idx_B));
        }

        void scale(const T alpha, const int* idx_A)
        {
            CHECK_RETURN_VALUE(
            tensor_scale_dense_(alpha, data_, ndim_, len_, ld_, idx_A));
        }

        /*
        void unpack(const PackedTensor& A)
        {
        #ifdef VALIDATE_INPUTS
            if (size_ != tensor_size_dense(A.ndim_, A.len_, NULL)) throw LengthMismatchError();
        #endif //VALIDATE_INPUTS

            CHECK_RETURN_VALUE(
            tensor_unpack(A.data_, data_, A.ndim_, A.len_, A.sym_));
        }
        */

        const DenseTensor slice(const int* start, const int* len) const
        {
            T* B;
            int ndim_B;
            std::vector<int> len_B(ndim_);
            std::vector<int> ldb(ndim_);

            CHECK_RETURN_VALUE(
            tensor_slice_dense(data_,   ndim_,         len_,        ld_,
                                  &B, &ndim_B, len_B.data(), ldb.data(),
                               start, len));

            return DenseTensor(ndim_B, len_B.data(), ldb.data(), B);
        }

        DenseTensor slice(const int* start, const int* len)
        {
            T* B;
            int ndim_B;
            std::vector<int> len_B(ndim_);
            std::vector<int> ldb(ndim_);

            CHECK_RETURN_VALUE(
            tensor_slice_dense(data_,   ndim_,         len_,        ld_,
                                  &B, &ndim_B, len_B.data(), ldb.data(),
                               start, len));

            return DenseTensor(ndim_B, len_B.data(), ldb.data(), B);
        }
};

template <typename T>
inline double scalar(const IndexedTensor< DenseTensor<T>,T >& other)
{
    T res;
    DenseTensor<T> dt(0, (int*)NULL, &res);
    dt[""] = other;
    return res;
}

template <typename T>
inline double scalar(const IndexedTensorMult< DenseTensor<T>,T >& other)
{
    T res;
    DenseTensor<T> dt(0, (int*)NULL, &res);
    dt[""] = other;
    return res;
}

}
}

#endif
