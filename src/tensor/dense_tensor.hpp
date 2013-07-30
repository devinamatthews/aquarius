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

template <typename T> class PackedTensor;

template <typename T>
class DenseTensor : public LocalTensor< DenseTensor<T>,T >
{
    friend class PackedTensor<T>;

    INHERIT_FROM_LOCAL_TENSOR(DenseTensor<T>,T)

    public:
        DenseTensor(const DenseTensor& t, const T val);

        DenseTensor(const DenseTensor<T>& A,
                const typename LocalTensor<DenseTensor<T>,T>::CopyType type=LocalTensor<DenseTensor<T>,T>::CLONE);

        DenseTensor(const int ndim, const int *len, T* data, const bool zero=false);

        DenseTensor(const int ndim, const int *len, const bool zero=true);

        DenseTensor(const int ndim, const int *len, const int *ld, T* data, const bool zero=false);

        DenseTensor(const int ndim, const int *len, const int *ld, const bool zero=true);

        static uint64_t getSize(const int ndim, const int *len, const int *ld);

        void print(FILE* fp) const;

        void print(std::ostream& stream) const;

        void mult(const T alpha, bool conja, const DenseTensor<T>& A, const int* idx_A,
                                 bool conjb, const DenseTensor<T>& B, const int* idx_B,
                  const T beta,                                       const int* idx_C);

        void sum(const T alpha, bool conja, const DenseTensor<T>& A, const int* idx_A,
                 const T beta,                                       const int* idx_B);

        void scale(const T alpha, const int* idx_A);

        DenseTensor slice(const int* start, const int* len);
};

}
}

#endif
