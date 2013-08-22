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

#ifndef _AQUARIUS_TENSOR_DIST_TENSOR_HPP_
#define _AQUARIUS_TENSOR_DIST_TENSOR_HPP_

#include <ostream>
#include <iostream>
#include <vector>
#include <cstdio>
#include <stdint.h>
#include <cstring>
#include <cassert>
#include <string>
#include <algorithm>
#include <cfloat>

#include "ctf.hpp"

#include "memory/memory.h"
#include "util/distributed.hpp"

#include "util.h"
#include "indexable_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template<typename T>
int conv_idx(const int ndim, const T* cidx, int*& iidx);

template<typename T>
int conv_idx(const int ndim_A, const T* cidx_A, int*& iidx_A,
             const int ndim_B, const T* cidx_B, int*& iidx_B);

template<typename T>
int conv_idx(const int ndim_A, const T* cidx_A, int*& iidx_A,
             const int ndim_B, const T* cidx_B, int*& iidx_B,
             const int ndim_C, const T* cidx_C, int*& iidx_C);

template <typename T>
class DistTensor : public IndexableTensor< DistTensor<T>,T >, public Distributed<T>
{
    INHERIT_FROM_INDEXABLE_TENSOR(DistTensor<T>,T)

    protected:
        int tid;
        int *len_;
        int *sym_;

    public:
        using Distributed<T>::ctf;

        DistTensor(const DistTensor& t, const T val);

        DistTensor(tCTF_World<T>& ctf);

        DistTensor(const DistTensor<T>& A,
                   const bool copy=true,
                   const bool zero=false);

        DistTensor(tCTF_World<T>& ctf,
                   const int ndim, const int *len, const int *sym,
                   const bool zero=true);

        ~DistTensor();

        void set_name(char const * name_);

        const int* getLengths() const { return len_; }

        const int* getSymmetry() const { return sym_; }

        T* getRawData(int64_t& size);

        const T* getRawData(int64_t& size) const;

        void getLocalData(int64_t& npair, tkv_pair<T>*& pairs) const;

        void getRemoteData(int64_t npair, tkv_pair<T>* pairs) const;

        void writeRemoteData(int64_t npair, tkv_pair<T>* pairs);

        void writeRemoteData(double alpha, double beta, int64_t npair, tkv_pair<T>* pairs);

        void getAllData(int64_t& npair, T*& vals) const;

        void getAllData(int64_t& npair, T*& vals, const int rank) const;

        void div(const T alpha, bool conja, const DistTensor<T>& A,
                                bool conjb, const DistTensor<T>& B, const T beta);

        void invert(const T alpha, bool conja, const DistTensor<T>& A, const T beta);

        void print(FILE* fp, double cutoff = 0.0) const;

        void compare(FILE* fp, const DistTensor<T>& other, double cutoff = 0.0) const;

        double reduce(CTF_OP op) const;

        void mult(const T alpha, bool conja, const DistTensor<T>& A, const int *idx_A,
                                 bool conjb, const DistTensor<T>& B, const int *idx_B,
                  const T  beta,                                     const int *idx_C);

        void sum(const T alpha, bool conja, const DistTensor<T>& A, const int *idx_A,
                 const T  beta,                                     const int *idx_B);

        void scale(const T alpha, const int* idx_A);

        T dot(bool conja, const DistTensor<T>& A, const int* idx_A,
              bool conjb,                         const int* idx_B) const;
};

}
}

#endif
