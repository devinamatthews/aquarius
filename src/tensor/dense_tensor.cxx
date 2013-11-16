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

#include "dense_tensor.hpp"

using namespace std;
using namespace aquarius::tensor;

template <typename T>
DenseTensor<T>::DenseTensor(const string& name, T val)
: LocalTensor<DenseTensor<T>,T>(name, val) {}

/*
DenseTensor(const PackedTensor& A)
: LocalTensor<DenseTensor>(A.ndim, A.len, (int*)NULL, A.size_, (bool)false)
{
    CHECK_RETURN_VALUE(
    tensor_densify(&ndim, len, A.sym_));

    copy(A.data, A.data+size_, data);
}

DenseTensor(const PackedTensor& A, const CopyType type)
: LocalTensor<DenseTensor>(A.ndim, A.len, (int*)NULL, A.size_, (double*)NULL)
{
    CHECK_RETURN_VALUE(
    tensor_densify(&ndim, len, A.sym_));

    switch (type)
    {
        case CLONE:
            data = new double[size_];
            copy(A.data, A.data+size_, data);
            isAlloced = true;
            break;
        case REFERENCE:
            data = A.data;
            isAlloced = false;
            break;
        case REPLACE:
            data = A.data;
            isAlloced = A.isAlloced;
            const_cast<PackedTensor&>(A).isAlloced = false;
            break;
    }
}
*/

template <typename T>
DenseTensor<T>::DenseTensor(const string& name, const DenseTensor<T>& A, T val)
: LocalTensor< DenseTensor<T>,T >(name, val) {}

template <typename T>
DenseTensor<T>::DenseTensor(const DenseTensor<T>& A)
: LocalTensor< DenseTensor<T>,T >(A) {}

template <typename T>
DenseTensor<T>::DenseTensor(const string& name, const DenseTensor<T>& A)
: LocalTensor< DenseTensor<T>,T >(name, A) {}

template <typename T>
DenseTensor<T>::DenseTensor(const string& name, DenseTensor<T>& A, typename LocalTensor<DenseTensor<T>,T>::CopyType type)
: LocalTensor< DenseTensor<T>,T >(name, A, type) {}

template <typename T>
DenseTensor<T>::DenseTensor(const string& name, int ndim, const vector<int>& len, T* data, bool zero)
: LocalTensor< DenseTensor<T>,T >(name, ndim, len, vector<int>(), getSize(ndim, len, vector<int>()), data, zero) {}

template <typename T>
DenseTensor<T>::DenseTensor(const string& name, int ndim, const vector<int>& len, bool zero)
: LocalTensor< DenseTensor<T>,T >(name, ndim, len, vector<int>(), getSize(ndim, len, vector<int>()), zero) {}

template <typename T>
DenseTensor<T>::DenseTensor(const string& name, int ndim, const vector<int>& len, const vector<int>& ld, T* data, bool zero)
: LocalTensor< DenseTensor<T>,T >(name, ndim, len, ld, getSize(ndim, len, ld), data, zero) {}

template <typename T>
DenseTensor<T>::DenseTensor(const string& name, int ndim, const vector<int>& len, const vector<int>& ld, bool zero)
: LocalTensor< DenseTensor<T>,T >(name, ndim, len, ld, getSize(ndim, len, ld), zero) {}

template <typename T>
uint64_t DenseTensor<T>::getSize(int ndim, const vector<int>& len, const vector<int>& ld)
{
    int64_t r = tensor_size_dense(ndim, len.data(), (ld.size() == 0 ? NULL : ld.data()));

    #ifdef VALIDATE_INPUTS
    CHECK_RETURN_VALUE(r);
    #endif //VALIDATE_INPUTS

    return r;
}

template <typename T>
void DenseTensor<T>::print(FILE* fp) const
{
    CHECK_RETURN_VALUE(
    tensor_print_dense(fp, data, ndim, len.data(), ld.data()));
}

template <typename T>
void DenseTensor<T>::print(ostream& stream) const
{
    #ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR_THROW(ndim, len.data(), ld.data(), NULL);
    #endif //VALIDATE_INPUTS

    vector<size_t> stride(ndim);
    if (ndim > 0) stride[0] = ld[0];
    for (int i = 1;i < ndim;i++) stride[i] = stride[i-1]*ld[i];

    size_t size;
    if (ndim > 0)
    {
        size = stride[ndim-1]*len[ndim-1];
    }
    else
    {
        size = 1;
    }

    size_t off = 0;

    /*
     * loop over elements in A
     */
    vector<int> pos(ndim, 0);
    for (bool done = false;!done;)
    {
        #ifdef CHECK_BOUNDS
        if (off < 0 || off >= size) throw OutOfBoundsError();
        #endif //CHECK_BOUNDS

        for (int i = 0;i < ndim;i++) stream << pos[i] << ' ';
        stream << scientific << setprecision(15) << data[off] << '\n';

        for (int i = 0;i < ndim;i++)
        {
            if (pos[i] == len[i]-1)
            {
                pos[i] = 0;
                off -= stride[i]*(len[i]-1);

                if (i == ndim-1)
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

        if (ndim == 0) done = true;
    }
    /*
     * end loop over A
     */
}

template <typename T>
void DenseTensor<T>::mult(const T alpha, bool conja, const DenseTensor<T>& A, const string& idx_A,
                                         bool conjb, const DenseTensor<T>& B, const string& idx_B,
                          const T beta,                                       const string& idx_C)
{
    vector<int> idx_A_(    A.ndim);
    vector<int> idx_B_(    B.ndim);
    vector<int> idx_C_(this->ndim);

    for (int i = 0;i <     A.ndim;i++) idx_A_[i] = idx_A[i];
    for (int i = 0;i <     B.ndim;i++) idx_B_[i] = idx_B[i];
    for (int i = 0;i < this->ndim;i++) idx_C_[i] = idx_C[i];

    CHECK_RETURN_VALUE(
    tensor_mult_dense_(alpha, A.data, A.ndim, A.len.data(), A.ld.data(), idx_A_.data(),
                              B.data, B.ndim, B.len.data(), B.ld.data(), idx_B_.data(),
                       beta,    data,   ndim,   len.data(),   ld.data(), idx_C_.data()));
}

template <typename T>
void DenseTensor<T>::sum(const T alpha, bool conja, const DenseTensor<T>& A, const string& idx_A,
                         const T beta,                                       const string& idx_B)
{
    vector<int> idx_A_(    A.ndim);
    vector<int> idx_B_(this->ndim);

    for (int i = 0;i <     A.ndim;i++) idx_A_[i] = idx_A[i];
    for (int i = 0;i < this->ndim;i++) idx_B_[i] = idx_B[i];

    CHECK_RETURN_VALUE(
    tensor_sum_dense_(alpha, A.data, A.ndim, A.len.data(), A.ld.data(), idx_A_.data(),
                      beta,    data,   ndim,   len.data(),   ld.data(), idx_B_.data()));
}

template <typename T>
void DenseTensor<T>::scale(const T alpha, const string& idx_A)
{
    vector<int> idx_A_(this->ndim);

    for (int i = 0;i < this->ndim;i++) idx_A_[i] = idx_A[i];

    CHECK_RETURN_VALUE(
    tensor_scale_dense_(alpha, data, ndim, len.data(), ld.data(), idx_A_.data()));
}

/*
void unpack(const PackedTensor& A)
{
#ifdef VALIDATE_INPUTS
    if (size_ != tensor_size_dense(A.ndim, A.len, NULL)) throw LengthMismatchError();
#endif //VALIDATE_INPUTS

    CHECK_RETURN_VALUE(
    tensor_unpack(A.data, data, A.ndim, A.len, A.sym_));
}

template <typename T>
DenseTensor<T> DenseTensor<T>::slice(const vector<int>& start, const vector<int>& len)
{
    T* B;
    int ndimB;
    vector<int> lenB(ndim);
    vector<int> ldB(ndim);

    CHECK_RETURN_VALUE(
    tensor_slice_dense(data,   ndim,  len.data(),  ld.data(),
                         &B, &ndimB, lenB.data(), ldB.data(),
                       start.data(), len.data()));

    lenB.resize(ndimB);
    ldB.resize(ndimB);

    return DenseTensor<T>(ndimB, lenB, ldB, B);
}
*/

INSTANTIATE_SPECIALIZATIONS(DenseTensor);
