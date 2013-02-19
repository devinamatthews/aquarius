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

#include <vector>
#include <algorithm>

#include "local_tensor.hpp"

namespace aquarius
{
namespace tensor
{

class PackedTensor : public LocalTensor<PackedTensor>
{
    friend class DenseTensor;
    friend class IndexedTensor<PackedTensor>;

    INHERIT_FROM_LOCAL_TENSOR(PackedTensor)

    protected:
        int* sym_;

    public:
        PackedTensor()
        : LocalTensor<PackedTensor>(0, NULL, NULL, 1)
        {
            sym_ = SAFE_MALLOC(int, 1);
            sym_[0] = NS;
        }

        PackedTensor(const PackedTensor& A, const CopyType type=CLONE)
        : LocalTensor<PackedTensor>(A, type)
        {
            sym_ = SAFE_MALLOC(int, ndim_);
            std::copy(A.sym_, A.sym_+ndim_, sym_);
        }

        PackedTensor(const int ndim, const int *len, const int *sym, double* data, const bool zero=false)
        : LocalTensor<PackedTensor>(ndim, len, NULL, getSize(ndim, len, NULL, sym), data, zero)
        {
            #ifdef VALIDATE_INPUTS
            VALIDATE_TENSOR_THROW(ndim, len, NULL, sym);
            #endif //VALIDATE_INPUTS

            sym_ = SAFE_MALLOC(int, ndim_);
            std::copy(sym, sym+ndim_, sym_);
        }

        PackedTensor(const int ndim, const int *len, const int *sym, const bool zero=true)
        : LocalTensor<PackedTensor>(ndim, len, NULL, getSize(ndim, len, NULL, sym), zero)
        {
            #ifdef VALIDATE_INPUTS
            VALIDATE_TENSOR_THROW(ndim, len, NULL, sym);
            #endif //VALIDATE_INPUTS

            sym_ = SAFE_MALLOC(int, ndim_);
            std::copy(sym, sym+ndim_, sym_);
        }

        PackedTensor(const int ndim, const int *len, const int *ld, const int *sym, double* data, const bool zero=false)
        : LocalTensor<PackedTensor>(ndim, len, ld, getSize(ndim, len, ld, sym), data, zero)
        {
            #ifdef VALIDATE_INPUTS
            VALIDATE_TENSOR_THROW(ndim, len, ld, sym);
            #endif //VALIDATE_INPUTS

            sym_ = SAFE_MALLOC(int, ndim_);
            std::copy(sym, sym+ndim_, sym_);
        }

        PackedTensor(const int ndim, const int *len, const int *ld, const int *sym, const bool zero=true)
        : LocalTensor<PackedTensor>(ndim, len, ld, getSize(ndim, len, ld, sym), zero)
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
                stream << scientific << setprecision(15) << data_[k] << '\n';

                done = !next_packed_indices(ndim_, len_, sym_, idx.data());
            }
        }

        void mult(const double alpha, const PackedTensor& A, const int* idx_A,
                                                    const PackedTensor& B, const int* idx_B,
                                const double beta,                         const int* idx_C)
        {
            CHECK_RETURN_VALUE(
            tensor_mult_(alpha, A.data_, A.ndim_, A.len_, A.ld_, A.sym_, idx_A,
                                B.data_, B.ndim_, B.len_, B.ld_, B.sym_, idx_B,
                         beta,    data_,   ndim_,   len_,   ld_,   sym_, idx_C));
        }

        void sum(const double alpha, const PackedTensor& A, const int* idx_A,
                               const double beta,                         const int* idx_B)
        {
            CHECK_RETURN_VALUE(
            tensor_sum_(alpha, A.data_, A.ndim_, A.len_, A.ld_, A.sym_, idx_A,
                        beta,    data_,   ndim_,   len_,   ld_,   sym_, idx_B));
        }

        void scale(const double alpha, const int* idx_A)
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

        void resym(const double alpha, const PackedTensor& A, const char* idx_A,
                   const double beta,                         const char* idx_B)
        {
            std::vector<int> idx_A_(A.ndim_);
            std::vector<int> idx_B_(ndim_);

            for (int i = 0;i < A.ndim_;i++) idx_A_[i] = idx_A[i];
            for (int i = 0;i < ndim_;i++) idx_B_[i] = idx_B[i];

            resym(alpha, A, idx_A_.data(),
                   beta,    idx_B_.data());
        }

        void resym(const double alpha, const PackedTensor& A, const int* idx_A,
                   const double beta,                         const int* idx_B)
        {
            CHECK_RETURN_VALUE(
            tensor_resym_(alpha, A.data_, A.ndim_, A.len_, A.ld_, A.sym_, idx_A,
                          beta,    data_,   ndim_,   len_,   ld_,   sym_, idx_B));
        }
};

template<>
double scalar(const IndexedTensor<PackedTensor>& other)
{
    double res;
    PackedTensor dt(0, (int*)NULL, (int*)NULL, &res);
    dt[""] = other;
    return res;
}

template<>
double scalar(const IndexedTensorMult<PackedTensor>& other)
{
    double res;
    PackedTensor dt(0, (int*)NULL, (int*)NULL, &res);
    dt[""] = other;
    return res;
}

}
}

#endif
