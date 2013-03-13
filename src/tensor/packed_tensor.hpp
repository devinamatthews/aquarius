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

#ifndef _AQUARIUS_TENSOR_PACKED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_PACKED_TENSOR_HPP_

#include <ostream>
#include <iostream>
#include <vector>
#include <cstdio>
#include <stdint.h>
#include <cstring>
#include <cassert>
#include <string>
#include <algorithm>
#include <iomanip>

#include "memory/memory.h"

#include "local_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <typename T> class DenseTensor;

template <class T>
class PackedTensor : public LocalTensor<PackedTensor<T>,T>
{
    friend class DenseTensor<T>;

    INHERIT_FROM_LOCAL_TENSOR(PackedTensor<T>,T)

    protected:
        int* sym_;

    public:
        PackedTensor(const PackedTensor& t, const T val)
        : Tensor<PackedTensor<T>,T>(*this),
          LocalTensor<PackedTensor,T>(t, val)
        {
            sym_ = SAFE_MALLOC(int, 1);
            sym_[0] = NS;
        }

        PackedTensor(const PackedTensor<T>& A, const typename LocalTensor<PackedTensor<T>,T>::CopyType type=CLONE)
        : Tensor<PackedTensor<T>,T>(*this),
          LocalTensor< PackedTensor<T>,T >(A, type)
        {
            sym_ = SAFE_MALLOC(int, ndim_);
            std::copy(A.sym_, A.sym_+ndim_, sym_);
        }

        PackedTensor(const int ndim, const int *len, const int *sym, T* data, const bool zero=false)
        : Tensor<PackedTensor<T>,T>(*this),
          LocalTensor< PackedTensor<T>,T >(ndim, len, NULL, getSize(ndim, len, NULL, sym), data, zero)
        {
            #ifdef VALIDATE_INPUTS
            VALIDATE_TENSOR_THROW(ndim, len, NULL, sym);
            #endif //VALIDATE_INPUTS

            sym_ = SAFE_MALLOC(int, ndim_);
            std::copy(sym, sym+ndim_, sym_);
        }

        PackedTensor(const int ndim, const int *len, const int *sym, const bool zero=true)
        : Tensor<PackedTensor<T>,T>(*this),
          LocalTensor< PackedTensor<T>,T >(ndim, len, NULL, getSize(ndim, len, NULL, sym), zero)
        {
            #ifdef VALIDATE_INPUTS
            VALIDATE_TENSOR_THROW(ndim, len, NULL, sym);
            #endif //VALIDATE_INPUTS

            sym_ = SAFE_MALLOC(int, ndim_);
            std::copy(sym, sym+ndim_, sym_);
        }

        PackedTensor(const int ndim, const int *len, const int *ld, const int *sym, T* data, const bool zero=false)
        : Tensor<PackedTensor<T>,T>(*this),
          LocalTensor< PackedTensor<T>,T >(ndim, len, ld, getSize(ndim, len, ld, sym), data, zero)
        {
            #ifdef VALIDATE_INPUTS
            VALIDATE_TENSOR_THROW(ndim, len, ld, sym);
            #endif //VALIDATE_INPUTS

            sym_ = SAFE_MALLOC(int, ndim_);
            std::copy(sym, sym+ndim_, sym_);
        }

        PackedTensor(const int ndim, const int *len, const int *ld, const int *sym, const bool zero=true)
        : Tensor<PackedTensor<T>,T>(*this),
          LocalTensor< PackedTensor<T>,T >(ndim, len, ld, getSize(ndim, len, ld, sym), zero)
        {
            #ifdef VALIDATE_INPUTS
            VALIDATE_TENSOR_THROW(ndim, len, ld, sym);
            #endif //VALIDATE_INPUTS

            sym_ = SAFE_MALLOC(int, ndim_);
            std::copy(sym, sym+ndim_, sym_);
        }

        ~PackedTensor()
        {
            FREE(sym_);
        }

        static uint64_t getSize(const int ndim, const int *len, const int *ld, const int *sym)
        {
            return tensor_size(ndim, len, ld, sym);
        }

        const int* getSymmetry() const
        {
            return sym_;
        }

        void print(FILE* fp) const
        {
            CHECK_RETURN_VALUE(
            tensor_print(fp, data_, ndim_, len_, ld_, sym_));
        }

        void print(std::ostream& stream) const
        {
            #ifdef VALIDATE_INPUTS
            VALIDATE_TENSOR_THROW(ndim_, len_, ld_, sym_);
            #endif //VALIDATE_INPUTS

            std::vector<int> idx(ndim_);
            first_packed_indices(ndim_, len_, sym_, idx.data());

            bool done = false;
            for (size_t k = 0;!done;k++)
            {
                #ifdef CHECK_BOUNDS
                if (k < 0 || k >= size_) throw OutOfBoundsError();
                #endif //CHECK_BOUNDS

                for (int i = 0;i < ndim_;i++) stream << idx[i];
                stream << std::scientific << std::setprecision(15) << data_[k] << '\n';

                done = !next_packed_indices(ndim_, len_, sym_, idx.data());
            }
        }

        void mult(const T alpha, bool conja, const IndexableTensor<PackedTensor<T>,T>& A, const int* idx_A,
                                 bool conjb, const IndexableTensor<PackedTensor<T>,T>& B, const int* idx_B,
                  const T beta,                                                           const int* idx_C)
        {
            CHECK_RETURN_VALUE(
            tensor_mult_(alpha,     A.data_,     A.ndim_,     A.len_,     A.ld_, A.sym_, idx_A,
                                    B.data_,     B.ndim_,     B.len_,     B.ld_, B.sym_, idx_B,
                         beta,  data_, ndim_, len_, ld_,   sym_, idx_C));
        }

        void sum(const T alpha, bool conja, const IndexableTensor<PackedTensor<T>,T>& A, const int* idx_A,
                 const T beta,                                                           const int* idx_B)
        {
            CHECK_RETURN_VALUE(
            tensor_sum_(alpha,     A.data_,     A.ndim_,     A.len_,     A.ld_, A.sym_, idx_A,
                        beta,  data_, ndim_, len_, ld_,   sym_, idx_B));
        }

        void scale(const T alpha, const int* idx_A)
        {
            CHECK_RETURN_VALUE(
            tensor_scale_(alpha, data_, ndim_, len_, ld_, sym_, idx_A));
        }

        /*
        void pack(const DenseTensor& A)
        {
            #ifdef VALIDATE_INPUTS
            if (size_ != tensor_size(A.ndim_, A.len_, NULL, sym_)) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            CHECK_RETURN_VALUE(
            tensor_pack(A.data_, data_, A.ndim_, A.len_, sym_));
        }

        void symmetrize(const DenseTensor& A)
        {
            #ifdef VALIDATE_INPUTS
            if (size_ != tensor_size(A.ndim_, A.len_, NULL, sym_)) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            CHECK_RETURN_VALUE(
            tensor_symmetrize(A.data_, data_, A.ndim_, A.len_, sym_));
        }
        */

        void resym(const T alpha, const IndexableTensor<PackedTensor<T>,T>& A, const char* idx_A,
                   const T beta,                                               const char* idx_B)
        {
            std::vector<int> idx_A_(A.ndim_);
            std::vector<int> idx_B_(ndim_);

            for (int i = 0;i <     A.ndim_;i++) idx_A_[i] = idx_A[i];
            for (int i = 0;i < ndim_;i++) idx_B_[i] = idx_B[i];

            resym(alpha, A, idx_A_.data(),
                   beta,    idx_B_.data());
        }

        void resym(const T alpha, const IndexableTensor<PackedTensor<T>,T>& A, const int* idx_A,
                   const T beta,                                               const int* idx_B)
        {
            CHECK_RETURN_VALUE(
            tensor_resym_(alpha,     A.data_,     A.ndim_,     A.len_,     A.ld_, A.sym_, idx_A,
                          beta,  data_, ndim_, len_, ld_,   sym_, idx_B));
        }
};

}
}

#endif
