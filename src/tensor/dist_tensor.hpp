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
#include "indexabletensor.hpp"

namespace aquarius
{
namespace tensor
{

template<typename T>
int conv_idx(const int ndim, const T* cidx, int*& iidx)
{
    iidx = SAFE_MALLOC(int, ndim);

    int n = 0;
    for (int i = 0;i < ndim;i++)
    {
        int j;
        for (j = 0;j < i;j++)
        {
            if (cidx[i] == cidx[j])
            {
                iidx[i] = iidx[j];
                break;
            }
        }
        if (j == i)
        {
            iidx[i] = n++;
        }
    }

    return n;
}

template<typename T>
int conv_idx(const int ndim_A, const T* cidx_A, int*& iidx_A,
             const int ndim_B, const T* cidx_B, int*& iidx_B)
{
    iidx_B = SAFE_MALLOC(int, ndim_B);

    int n = conv_idx(ndim_A, cidx_A, iidx_A);
    for (int i = 0;i < ndim_B;i++)
    {
        int j;
        for (j = 0;j < ndim_A;j++)
        {
            if (cidx_B[i] == cidx_A[j])
            {
                iidx_B[i] = iidx_A[j];
                break;
            }
        }
        if (j == ndim_A)
        {
            for (j = 0;j < i;j++)
            {
                if (cidx_B[i] == cidx_B[j])
                {
                    iidx_B[i] = iidx_B[j];
                    break;
                }
            }
            if (j == i)
            {
                iidx_B[i] = n++;
            }
        }
    }

    return n;
}

template<typename T>
int conv_idx(const int ndim_A, const T* cidx_A, int*& iidx_A,
             const int ndim_B, const T* cidx_B, int*& iidx_B,
             const int ndim_C, const T* cidx_C, int*& iidx_C)
{
    iidx_C = SAFE_MALLOC(int, ndim_C);

    int n = conv_idx(ndim_A, cidx_A, iidx_A,
                     ndim_B, cidx_B, iidx_B);
    for (int i = 0;i < ndim_C;i++)
    {
        int j;
        for (j = 0;j < ndim_B;j++)
        {
            if (cidx_C[i] == cidx_B[j])
            {
                iidx_C[i] = iidx_B[j];
                break;
            }
        }
        if (j == ndim_B)
        {
            for (j = 0;j < ndim_A;j++)
            {
                if (cidx_C[i] == cidx_A[j])
                {
                    iidx_C[i] = iidx_A[j];
                    break;
                }
            }
            if (j == ndim_A)
            {
                for (j = 0;j < i;j++)
                {
                    if (cidx_C[i] == cidx_C[j])
                    {
                        iidx_C[i] = iidx_C[j];
                        break;
                    }
                }
                if (j == i)
                {
                    iidx_C[i] = n++;
                }
            }
        }
    }

    return n;
}

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

        DistTensor(const DistTensor& t, const T val)
        : IndexableTensor< DistTensor<T>,T >(),
          Distributed<T>(t.ctf)
        {
            len_ = SAFE_MALLOC(int, ndim_);
            sym_ = SAFE_MALLOC(int, ndim_);

            int ret = ctf.ctf->define_tensor(ndim_, len_, sym_, &tid);
            assert(ret == DIST_TENSOR_SUCCESS);

            if (this->rank == 0)
            {
                tkv_pair<T> p(0, val);
                writeRemoteData(1, &p);
            }
            else
            {
                writeRemoteData(0, NULL);
            }
        }

        DistTensor(tCTF_World<T>& ctf)
        : IndexableTensor< DistTensor<T>,T >(), Distributed<T>(ctf)
        {
            len_ = SAFE_MALLOC(int, ndim_);
            sym_ = SAFE_MALLOC(int, ndim_);

            int ret = ctf.ctf->define_tensor(ndim_, len_, sym_, &tid);
            assert(ret == DIST_TENSOR_SUCCESS);
        }

        DistTensor(const DistTensor<T>& A,
                   const bool copy=true,
                   const bool zero=false)
        : IndexableTensor< DistTensor<T>,T >(A.ndim_), Distributed<T>(A.ctf)
        {
            len_ = SAFE_MALLOC(int, ndim_);
            sym_ = SAFE_MALLOC(int, ndim_);

            std::copy(A.len_, A.len_+ndim_, len_);
            std::copy(A.sym_, A.sym_+ndim_, sym_);

            int ret = ctf.ctf->define_tensor(ndim_, len_, sym_, &tid);
            assert(ret == DIST_TENSOR_SUCCESS);

            if (copy)
            {
                ret = ctf.ctf->copy_tensor(A.tid, tid);
                assert(ret == DIST_TENSOR_SUCCESS);
            }
            else if (zero)
            {
                int64_t size;
                T* raw_data = getRawData(size);
                std::fill(raw_data, raw_data+size, (T)0);
            }
        }

        DistTensor(tCTF_World<T>& ctf,
                   const int ndim, const int *len, const int *sym,
                   const bool zero=true)
        : IndexableTensor< DistTensor<T>,T >(ndim), Distributed<T>(ctf)
        {
            len_ = SAFE_MALLOC(int, ndim_);
            sym_ = SAFE_MALLOC(int, ndim_);

            std::copy(len, len+ndim_, len_);
            std::copy(sym, sym+ndim_, sym_);

            #ifdef VALIDATE_INPUTS
            validate_tensor(ndim, len, NULL, sym);
            #endif //VALIDATE_INPUTS

            int ret = ctf.ctf->define_tensor(ndim, len, sym, &tid);
            assert(ret == DIST_TENSOR_SUCCESS);

            if (zero)
            {
                int64_t size;
                T* raw_data = getRawData(size);
                std::fill(raw_data, raw_data+size, (T)0);
            }
        }

        virtual ~DistTensor()
        {
            FREE(len_);
            FREE(sym_);
            ctf.ctf->clean_tensor(tid);
        }

        const int* getLengths() const
        {
            return len_;
        }

        const int* getSymmetry() const
        {
            return sym_;
        }

        T* getRawData(int64_t& size)
        {
            T *data;
            int ret = ctf.ctf->get_raw_data(tid, &data, &size);
            assert(ret == DIST_TENSOR_SUCCESS);
            return data;
        }

        const T* getRawData(int64_t& size) const
        {
            T *data;
            int ret = ctf.ctf->get_raw_data(tid, &data, &size);
            assert(ret == DIST_TENSOR_SUCCESS);
            return data;
        }

        void getLocalData(int64_t& npair, tkv_pair<T>*& pairs) const
        {
            int ret = ctf.ctf->read_local_tensor(tid, &npair, &pairs);
            assert(ret == DIST_TENSOR_SUCCESS);
        }

        void getRemoteData(int64_t npair, tkv_pair<T>* pairs) const
        {
            int ret = ctf.ctf->read_tensor(tid, npair, pairs);
            assert(ret == DIST_TENSOR_SUCCESS);
        }

        void writeRemoteData(int64_t npair, tkv_pair<T>* pairs)
        {
            int ret = ctf.ctf->write_tensor(tid, npair, pairs);
            assert(ret == DIST_TENSOR_SUCCESS);
        }

        void writeRemoteData(double alpha, double beta, int64_t npair, tkv_pair<T>* pairs)
        {
            int ret = ctf.ctf->write_tensor(tid, npair, alpha, beta, pairs);
            assert(ret == DIST_TENSOR_SUCCESS);
        }

        void getAllData(int64_t& npair, T*& vals) const
        {
            int ret = ctf.ctf->allread_tensor(tid, &npair, &vals);
            assert(ret == DIST_TENSOR_SUCCESS);
        }

        void getAllData(int64_t& npair, T*& vals, const int rank) const
        {
            if (this->rank == rank)
            {
                std::vector<tkv_pair<T> > pairs;
                std::vector<int> idx(ndim_, 0);

                first_packed_indices(ndim_, len_, sym_, idx.data());

                do
                {
                    int64_t key = 0, stride = 1;
                    for (int i = 0;i < ndim_;i++)
                    {
                        key += idx[i]*stride;
                        stride *= len_[i];
                    }
                    pairs.push_back(tkv_pair<T>(key, (T)0));
                }
                while (next_packed_indices(ndim_, len_, sym_, idx.data()));

                int ret = ctf.ctf->read_tensor(tid, pairs.size(), pairs.data());
                assert(ret == DIST_TENSOR_SUCCESS);

                npair = pairs.size();
                vals = (T*)malloc(sizeof(T)*npair);
                assert(vals != NULL);

                for (size_t i = 0;i < npair;i++)
                {
                    vals[i] = pairs[i].d;
                }
            }
            else
            {
                int ret = ctf.ctf->read_tensor(tid, 0, NULL);
                assert(ret == DIST_TENSOR_SUCCESS);
            }
        }

        void div(const T alpha, bool conja, const DistTensor<T>& A,
                                bool conjb, const DistTensor<T>& B, const T beta)
        {
            ctf.ctf->align(A.tid, tid);
            ctf.ctf->align(B.tid, tid);
            int64_t size, size_A, size_B;
            T* raw_data = getRawData(size);
            const T* raw_data_A = A.getRawData(size_A);
            const T* raw_data_B = B.getRawData(size_B);
            assert(size == size_A);
            assert(size == size_B);
            if (conja)
            {
                if (conjb)
                {
                    for (int64_t i = 0;i < size;i++)
                    {
                        if (std::abs(raw_data_B[i]) > DBL_MIN)
                        {
                            raw_data[i] = beta*raw_data[i] + alpha*std::conj(raw_data_A[i])/std::conj(raw_data_B[i]);
                        }
                    }
                }
                else
                {
                    for (int64_t i = 0;i < size;i++)
                    {
                        if (std::abs(raw_data_B[i]) > DBL_MIN)
                        {
                            raw_data[i] = beta*raw_data[i] + alpha*std::conj(raw_data_A[i])/raw_data_B[i];
                        }
                    }
                }
            }
            else
            {
                if (conjb)
                {
                    for (int64_t i = 0;i < size;i++)
                    {
                        if (std::abs(raw_data_B[i]) > DBL_MIN)
                        {
                            raw_data[i] = beta*raw_data[i] + alpha*raw_data_A[i]/std::conj(raw_data_B[i]);
                        }
                    }
                }
                else
                {
                    for (int64_t i = 0;i < size;i++)
                    {
                        if (std::abs(raw_data_B[i]) > DBL_MIN)
                        {
                            raw_data[i] = beta*raw_data[i] + alpha*raw_data_A[i]/raw_data_B[i];
                        }
                    }
                }
            }
        }

        void invert(const T alpha, bool conja, const DistTensor<T>& A, const T beta)
        {
            ctf.ctf->align(tid, A.tid);
            int64_t size, size_A;
            T* raw_data = getRawData(size);
            const T* raw_data_A = A.getRawData(size_A);
            assert(size == size_A);
            if (conja)
            {
                for (int64_t i = 0;i < size;i++)
                {
                    if (std::abs(raw_data_A[i]) > DBL_MIN)
                    {
                        raw_data[i] = beta*raw_data[i] + alpha/std::conj(raw_data_A[i]);
                    }
                }
            }
            else
            {
                for (int64_t i = 0;i < size;i++)
                {
                    if (std::abs(raw_data_A[i]) > DBL_MIN)
                    {
                        raw_data[i] = beta*raw_data[i] + alpha/raw_data_A[i];
                    }
                }
            }
        }

        void print(FILE* fp) const
        {
            ctf.ctf->print_tensor(fp, tid);
        }

        void compare(FILE* fp, const DistTensor<T>& other) const
        {
            ctf.ctf->compare_tensor(fp, tid, other.tid);
        }

        T reduce(CTF_OP op) const
        {
            int ret;
            T ans = 0.0;
            ret = ctf.ctf->reduce_tensor(tid, op, &ans);
            assert(ret == DIST_TENSOR_SUCCESS);
            return ans;
        }

        void mult(const T alpha, bool conja, const DistTensor<T>& A, const int *idx_A,
                                 bool conjb, const DistTensor<T>& B, const int *idx_B,
                  const T  beta,                                     const int *idx_C)
        {
            int ret;
            CTF_ctr_type_t tp;

            tp.tid_A = A.tid;
            tp.tid_B = B.tid;
            tp.tid_C =   tid;
            conv_idx(    A.ndim_, idx_A, tp.idx_map_A,
                         B.ndim_, idx_B, tp.idx_map_B,
                     ndim_, idx_C, tp.idx_map_C);

            ret = ctf.ctf->contract(&tp, alpha, beta);

            FREE(tp.idx_map_A);
            FREE(tp.idx_map_B);
            FREE(tp.idx_map_C);
            assert(ret == DIST_TENSOR_SUCCESS);
        }

        void sum(const T alpha, bool conja, const DistTensor<T>& A, const int *idx_A,
                 const T  beta,                                     const int *idx_B)
        {
            int ret;
            CTF_sum_type_t st;

            st.tid_A = A.tid;
            st.tid_B =   tid;
            conv_idx(    A.ndim_, idx_A, st.idx_map_A,
                     ndim_, idx_B, st.idx_map_B);

            ret = ctf.ctf->sum_tensors(&st, alpha, beta);

            FREE(st.idx_map_A);
            FREE(st.idx_map_B);
            assert(ret == DIST_TENSOR_SUCCESS);
        }

        void scale(const T alpha, const int* idx_A)
        {
            int ret;
            int * idx_map_A;

            conv_idx(ndim_, idx_A, idx_map_A);

            ret = ctf.ctf->scale_tensor(alpha, tid, idx_map_A);

            FREE(idx_map_A);
            assert(ret == DIST_TENSOR_SUCCESS);
        }

        T dot(bool conja, const DistTensor<T>& A, const int* idx_A,
              bool conjb,                         const int* idx_B) const
        {
            DistTensor<T> dt(A, (T)0);
            int64_t n;
            T ret, *val;
            dt.mult(1, conja,     A, idx_A,
                       conjb, *this, idx_B,
                    0,                NULL);
            dt.getAllData(n, val);
            assert(n==1);
            ret = val[0];
            free(val);
            return ret;
        }
};

}
}

#endif
