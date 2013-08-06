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

#include "packed_tensor.hpp"

using namespace std;
using namespace aquarius::tensor;

template <class T>
PackedTensor<T>::PackedTensor(T val)
: LocalTensor<PackedTensor,T>(val), sym(0) {}

template <class T>
PackedTensor<T>::PackedTensor(const PackedTensor<T>& A, T val)
: LocalTensor< PackedTensor<T>,T >(val), sym(0) {}

template <class T>
PackedTensor<T>::PackedTensor(const PackedTensor<T>& A)
: LocalTensor< PackedTensor<T>,T >(A), sym(A.sym) {}

template <class T>
PackedTensor<T>::PackedTensor(PackedTensor<T>& A, typename LocalTensor<PackedTensor<T>,T>::CopyType type)
: LocalTensor< PackedTensor<T>,T >(A, type), sym(A.sym) {}

template <class T>
PackedTensor<T>::PackedTensor(int ndim, const vector<int>& len, const vector<int>& sym, T* data, bool zero)
: LocalTensor< PackedTensor<T>,T >(ndim, len, vector<int>(), getSize(ndim, len, vector<int>(), sym), data, zero),
  sym(sym)
{
    assert(sym.size() == ndim);
    #ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR_THROW(ndim, len.data(), ld.data(), sym.data());
    #endif //VALIDATE_INPUTS
}

template <class T>
PackedTensor<T>::PackedTensor(int ndim, const vector<int>& len, const vector<int>& sym, bool zero)
: LocalTensor< PackedTensor<T>,T >(ndim, len, vector<int>(), getSize(ndim, len, vector<int>(), sym), zero),
  sym(sym)
{
    assert(sym.size() == ndim);
    #ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR_THROW(ndim, len.data(), ld.data(), sym.data());
    #endif //VALIDATE_INPUTS
}

template <class T>
PackedTensor<T>::PackedTensor(int ndim, const vector<int>& len, const vector<int>& ld, const vector<int>& sym, T* data, bool zero)
: LocalTensor< PackedTensor<T>,T >(ndim, len, ld, getSize(ndim, len, ld, sym), data, zero), sym(sym)
{
    assert(sym.size() == ndim);
    #ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR_THROW(ndim, len.data(), ld.data(), sym.data());
    #endif //VALIDATE_INPUTS
}

template <class T>
PackedTensor<T>::PackedTensor(int ndim, const vector<int>& len, const vector<int>& ld, const vector<int>& sym, bool zero)
: LocalTensor< PackedTensor<T>,T >(ndim, len, ld, getSize(ndim, len, ld, sym), zero), sym(sym)
{
    assert(sym.size() == ndim);
    #ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR_THROW(ndim, len.data(), ld.data(), sym.data());
    #endif //VALIDATE_INPUTS
}

template <class T>
uint64_t PackedTensor<T>::getSize(int ndim, const vector<int>& len, const vector<int>& ld, const vector<int>& sym)
{
    return tensor_size(ndim, len.data(), (ld.size() == 0 ? NULL : ld.data()), sym.data());
}

template <class T>
void PackedTensor<T>::print(FILE* fp) const
{
    CHECK_RETURN_VALUE(
    tensor_print(fp, data, ndim, len.data(), ld.data(), sym.data()));
}

template <class T>
void PackedTensor<T>::print(ostream& stream) const
{
    #ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR_THROW(ndim, len.data(), ld.data(), sym.data());
    #endif //VALIDATE_INPUTS

    vector<int> idx(ndim);
    first_packed_indices(ndim, len.data(), sym.data(), idx.data());

    bool done = false;
    for (size_t k = 0;!done;k++)
    {
        #ifdef CHECK_BOUNDS
        if (k < 0 || k >= size_) throw OutOfBoundsError();
        #endif //CHECK_BOUNDS

        for (int i = 0;i < ndim;i++) stream << idx[i];
        stream << scientific << setprecision(15) << data[k] << '\n';

        done = !next_packed_indices(ndim, len.data(), sym.data(), idx.data());
    }
}

template <class T>
void PackedTensor<T>::mult(const T alpha, bool conja, const PackedTensor<T>& A, const int* idx_A,
                                          bool conjb, const PackedTensor<T>& B, const int* idx_B,
                           const T  beta,                                       const int* idx_C)
{
    CHECK_RETURN_VALUE(
    tensor_mult_(alpha, A.data, A.ndim, A.len.data(), A.ld.data(), A.sym.data(), idx_A,
                        B.data, B.ndim, B.len.data(), B.ld.data(), B.sym.data(), idx_B,
                 beta,    data,   ndim,   len.data(),   ld.data(),   sym.data(), idx_C));
}

template <class T>
void PackedTensor<T>::sum(const T alpha, bool conja, const PackedTensor<T>& A, const int* idx_A,
                          const T  beta,                                       const int* idx_B)
{
    CHECK_RETURN_VALUE(
    tensor_sum_(alpha, A.data, A.ndim, A.len.data(), A.ld.data(), A.sym.data(), idx_A,
                beta,    data,   ndim,   len.data(),   ld.data(),   sym.data(), idx_B));
}

template <class T>
void PackedTensor<T>::scale(const T alpha, const int* idx_A)
{
    CHECK_RETURN_VALUE(
    tensor_scale_(alpha, data, ndim, len.data(), ld.data(), sym.data(), idx_A));
}

/*
void pack(const DenseTensor& A)
{
    #ifdef VALIDATE_INPUTS
    if (size_ != tensor_size(A.ndim, A.len, NULL, sym)) throw LengthMismatchError();
    #endif //VALIDATE_INPUTS

    CHECK_RETURN_VALUE(
    tensor_pack(A.data, data, A.ndim, A.len, sym));
}

void symmetrize(const DenseTensor& A)
{
    #ifdef VALIDATE_INPUTS
    if (size_ != tensor_size(A.ndim, A.len, NULL, sym)) throw LengthMismatchError();
    #endif //VALIDATE_INPUTS

    CHECK_RETURN_VALUE(
    tensor_symmetrize(A.data, data, A.ndim, A.len, sym));
}
*/

template <class T>
void PackedTensor<T>::resym(const T alpha, const PackedTensor<T>& A, const char* idx_A,
                            const T  beta,                           const char* idx_B)
{
    vector<int> idx_A_(A.ndim);
    vector<int> idx_B_(ndim);

    for (int i = 0;i < A.ndim;i++) idx_A_[i] = idx_A[i];
    for (int i = 0;i <   ndim;i++) idx_B_[i] = idx_B[i];

    resym(alpha, A, idx_A_.data(),
           beta,    idx_B_.data());
}

template <class T>
void PackedTensor<T>::resym(const T alpha, const PackedTensor<T>& A, const int* idx_A,
                            const T beta,                            const int* idx_B)
{
    CHECK_RETURN_VALUE(
    tensor_resym(alpha, A.data, A.ndim, A.len.data(), A.ld.data(), A.sym.data(), idx_A,
                  beta,   data,   ndim,   len.data(),   ld.data(),   sym.data(), idx_B));
}

//INSTANTIATE_SPECIALIZATIONS(PackedTensor);
