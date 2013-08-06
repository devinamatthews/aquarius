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

#include "dist_tensor.hpp"

using namespace std;
using namespace aquarius::tensor;

template<typename T>
int aquarius::tensor::conv_idx(int ndim, const T* cidx, int*& iidx)
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

template int aquarius::tensor::conv_idx(int, const  int*, int*&);
template int aquarius::tensor::conv_idx(int, const char*, int*&);

template<typename T>
int aquarius::tensor::conv_idx(int ndimA, const T* cidx_A, int*& iidx_A,
             int ndimB, const T* cidx_B, int*& iidx_B)
{
    iidx_B = SAFE_MALLOC(int, ndimB);

    int n = conv_idx(ndimA, cidx_A, iidx_A);
    for (int i = 0;i < ndimB;i++)
    {
        int j;
        for (j = 0;j < ndimA;j++)
        {
            if (cidx_B[i] == cidx_A[j])
            {
                iidx_B[i] = iidx_A[j];
                break;
            }
        }
        if (j == ndimA)
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

template int aquarius::tensor::conv_idx(int, const  int*, int*&,
                      int, const  int*, int*&);
template int aquarius::tensor::conv_idx(int, const char*, int*&,
                      int, const char*, int*&);

template<typename T>
int aquarius::tensor::conv_idx(int ndimA, const T* cidx_A, int*& iidx_A,
             int ndimB, const T* cidx_B, int*& iidx_B,
             int ndimC, const T* cidx_C, int*& iidx_C)
{
    iidx_C = SAFE_MALLOC(int, ndimC);

    int n = conv_idx(ndimA, cidx_A, iidx_A,
                     ndimB, cidx_B, iidx_B);
    for (int i = 0;i < ndimC;i++)
    {
        int j;
        for (j = 0;j < ndimB;j++)
        {
            if (cidx_C[i] == cidx_B[j])
            {
                iidx_C[i] = iidx_B[j];
                break;
            }
        }
        if (j == ndimB)
        {
            for (j = 0;j < ndimA;j++)
            {
                if (cidx_C[i] == cidx_A[j])
                {
                    iidx_C[i] = iidx_A[j];
                    break;
                }
            }
            if (j == ndimA)
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

template int aquarius::tensor::conv_idx(int, const  int*, int*&,
                      int, const  int*, int*&,
                      int, const  int*, int*&);
template int aquarius::tensor::conv_idx(int, const char*, int*&,
                      int, const char*, int*&,
                      int, const char*, int*&);

/*
 * Create a scalar (0-dimensional tensor)
 */
template <typename T>
DistTensor<T>::DistTensor(Arena<T>& arena, T scalar)
: IndexableTensor< DistTensor<T>,T >(), Distributed<T>(arena), len(0), sym(0)
{
    allocate();

    int64_t size;
    T* raw_data = getRawData(size);
    fill(raw_data, raw_data+size, scalar);
}

/*
 * Create a scalar (0-dimensional tensor) on the same arena as A
 */
template <typename T>
DistTensor<T>::DistTensor(const DistTensor<T>& A, T scalar)
: IndexableTensor< DistTensor<T>,T >(), Distributed<T>(A.arena), len(0), sym(0)
{
    allocate();

    int64_t size;
    T* raw_data = getRawData(size);
    fill(raw_data, raw_data+size, scalar);
}

/*
 * Create a tensor of the same size and shape as A, optionally copying or zeroing the data
 */
template <typename T>
DistTensor<T>::DistTensor(const DistTensor<T>& A, bool copy, bool zero)
: IndexableTensor< DistTensor<T>,T >(A.ndim), Distributed<T>(A.arena), len(A.len), sym(A.sym)
{
    allocate();

    if (copy)
    {
        int ret = arena.ctf.ctf->copy_tensor(A.tid, tid);
        assert(ret == DIST_TENSOR_SUCCESS);
    }
    else if (zero)
    {
        int64_t size;
        T* raw_data = getRawData(size);
        fill(raw_data, raw_data+size, (T)0);
    }
}

/*
 * Create a tensor of the specified size and shape, optionally zeroing the data
 */
template <typename T>
DistTensor<T>::DistTensor(Arena<T>& arena, int ndim, const vector<int>& len, const vector<int>& sym,
                          bool zero)
: IndexableTensor< DistTensor<T>,T >(ndim), Distributed<T>(arena), len(len), sym(sym)
{
    assert(len.size() == ndim);
    assert(sym.size() == ndim);

    #ifdef VALIDATE_INPUTS
    validate_tensor(ndim, len.data(), NULL, sym.data());
    #endif //VALIDATE_INPUTS

    allocate();

    if (zero)
    {
        int64_t size;
        T* raw_data = getRawData(size);
        fill(raw_data, raw_data+size, (T)0);
    }
}

template <typename T>
DistTensor<T>::~DistTensor()
{
    free();
}

template <typename T>
void DistTensor<T>::allocate()
{
    int ret = arena.ctf.ctf->define_tensor(ndim, len.data(), sym.data(), &tid);
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::free()
{
    int ret = arena.ctf.ctf->clean_tensor(tid);
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::resize(int ndim, const vector<int>& len, const vector<int>& sym, bool zero)
{
    assert(len.size() == ndim);
    assert(sym.size() == ndim);

    free();
    this->ndim = ndim;
    this->len = len;
    this->sym = sym;
    allocate();

    if (zero)
    {
        int64_t size;
        T* raw_data = getRawData(size);
        fill(raw_data, raw_data+size, (T)0);
    }
}

template <typename T>
T* DistTensor<T>::getRawData(int64_t& size)
{
    return const_cast<T*>(const_cast<const DistTensor<T>&>(*this).getRawData(size));
}

template <typename T>
const T* DistTensor<T>::getRawData(int64_t& size) const
{
    T *data;
    int ret = arena.ctf.ctf->get_raw_data(tid, &data, &size);
    assert(ret == DIST_TENSOR_SUCCESS);
    return data;
}

template <typename T>
void DistTensor<T>::getLocalData(vector<tkv_pair<T> >& pairs) const
{
    int64_t npair;
    tkv_pair<T> *data;
    int ret = arena.ctf.ctf->read_local_tensor(tid, &npair, &data);
    assert(ret == DIST_TENSOR_SUCCESS);
    pairs.assign(data, data+npair);
    ::free(data);
}

template <typename T>
void DistTensor<T>::getRemoteData(vector<tkv_pair<T> >& pairs) const
{
    int ret = arena.ctf.ctf->read_tensor(tid, pairs.size(), pairs.data());
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::getRemoteData() const
{
    int ret = arena.ctf.ctf->read_tensor(tid, 0, NULL);
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::writeRemoteData(const vector<tkv_pair<T> >& pairs)
{
    int ret = arena.ctf.ctf->write_tensor(tid, pairs.size(), pairs.data());
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::writeRemoteData()
{
    int ret = arena.ctf.ctf->write_tensor(tid, 0, NULL);
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::writeRemoteData(double alpha, double beta, const vector<tkv_pair<T> >& pairs)
{
    int ret = arena.ctf.ctf->write_tensor(tid, pairs.size(), alpha, beta, pairs.data());
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::writeRemoteData(double alpha, double beta)
{
    int ret = arena.ctf.ctf->write_tensor(tid, 0, alpha, beta, NULL);
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::getAllData(vector<T>& vals) const
{

    getAllData(vals, 0);
    int64_t npair = vals.size();
    this->comm.Bcast(&npair, 1, MPI_LONG_LONG, 0);
    if (this->rank != 0) vals.resize(npair);
    this->comm.Bcast(vals.data(), npair, MPI_TYPE_<T>::value(), 0);
}

template <typename T>
void DistTensor<T>::getAllData(vector<T>& vals, int rank) const
{
    if (this->rank == rank)
    {
        vector<tkv_pair<T> > pairs;
        vector<int> idx(ndim, 0);

        first_packed_indices(ndim, len.data(), sym.data(), idx.data());

        do
        {
            int64_t key = 0, stride = 1;
            for (int i = 0;i < ndim;i++)
            {
                key += idx[i]*stride;
                stride *= len[i];
            }
            pairs.push_back(tkv_pair<T>(key, (T)0));
        }
        while (next_packed_indices(ndim, len.data(), sym.data(), idx.data()));

        int ret = arena.ctf.ctf->read_tensor(tid, pairs.size(), pairs.data());
        assert(ret == DIST_TENSOR_SUCCESS);

        sort(pairs.begin(), pairs.end());
        size_t npair = pairs.size();
        vals.resize(npair);

        for (size_t i = 0;i < npair;i++)
        {
            vals[i] = pairs[i].d;
        }
    }
    else
    {
        int ret = arena.ctf.ctf->read_tensor(tid, 0, NULL);
        assert(ret == DIST_TENSOR_SUCCESS);
    }
}

template <typename T>
void DistTensor<T>::div(T alpha, bool conja, const DistTensor<T>& A,
                                 bool conjb, const DistTensor<T>& B, T beta)
{
    arena.ctf.ctf->align(A.tid, tid);
    arena.ctf.ctf->align(B.tid, tid);
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
                if (abs(raw_data_B[i]) > DBL_MIN)
                {
                    raw_data[i] = beta*raw_data[i] + alpha*conj(raw_data_A[i])/conj(raw_data_B[i]);
                }
            }
        }
        else
        {
            for (int64_t i = 0;i < size;i++)
            {
                if (abs(raw_data_B[i]) > DBL_MIN)
                {
                    raw_data[i] = beta*raw_data[i] + alpha*conj(raw_data_A[i])/raw_data_B[i];
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
                if (abs(raw_data_B[i]) > DBL_MIN)
                {
                    raw_data[i] = beta*raw_data[i] + alpha*raw_data_A[i]/conj(raw_data_B[i]);
                }
            }
        }
        else
        {
            for (int64_t i = 0;i < size;i++)
            {
                if (abs(raw_data_B[i]) > DBL_MIN)
                {
                    raw_data[i] = beta*raw_data[i] + alpha*raw_data_A[i]/raw_data_B[i];
                }
            }
        }
    }
}

template <typename T>
void DistTensor<T>::invert(T alpha, bool conja, const DistTensor<T>& A, T beta)
{
    arena.ctf.ctf->align(tid, A.tid);
    int64_t size, size_A;
    T* raw_data = getRawData(size);
    const T* raw_data_A = A.getRawData(size_A);
    assert(size == size_A);
    if (conja)
    {
        for (int64_t i = 0;i < size;i++)
        {
            if (abs(raw_data_A[i]) > DBL_MIN)
            {
                raw_data[i] = beta*raw_data[i] + alpha/conj(raw_data_A[i]);
            }
        }
    }
    else
    {
        for (int64_t i = 0;i < size;i++)
        {
            if (abs(raw_data_A[i]) > DBL_MIN)
            {
                raw_data[i] = beta*raw_data[i] + alpha/raw_data_A[i];
            }
        }
    }
}

template <typename T>
void DistTensor<T>::print(FILE* fp, double cutoff) const
{
    arena.ctf.ctf->print_tensor(fp, tid, cutoff);
}

template <typename T>
void DistTensor<T>::compare(FILE* fp, const DistTensor<T>& other, double cutoff) const
{
    arena.ctf.ctf->compare_tensor(fp, tid, other.tid, cutoff);
}

template <typename T>
typename real_type<T>::type DistTensor<T>::norm(int p) const
{
    int ret;
    T ans = (T)0;
    if (p == 00)
    {
        ret = arena.ctf.ctf->reduce_tensor(tid, CTF_OP_MAXABS, &ans);
    }
    else if (p == 1)
    {
        ret = arena.ctf.ctf->reduce_tensor(tid, CTF_OP_SUMABS, &ans);
    }
    else if (p == 2)
    {
        ret = arena.ctf.ctf->reduce_tensor(tid, CTF_OP_SQNRM2, &ans);
    }
    assert(ret == DIST_TENSOR_SUCCESS);
    return abs(ans);
}

template <typename T>
void DistTensor<T>::mult(T alpha, bool conja, const DistTensor<T>& A, const int *idx_A,
                                  bool conjb, const DistTensor<T>& B, const int *idx_B,
                         T  beta,                                     const int *idx_C)
{
    int ret;
    CTF_ctr_type_t tp;

    tp.tid_A = A.tid;
    tp.tid_B = B.tid;
    tp.tid_C =   tid;
    conv_idx(    A.ndim, idx_A, tp.idx_map_A,
                 B.ndim, idx_B, tp.idx_map_B,
                   ndim, idx_C, tp.idx_map_C);

    ret = arena.ctf.ctf->contract(&tp, alpha, beta);

    FREE(tp.idx_map_A);
    FREE(tp.idx_map_B);
    FREE(tp.idx_map_C);
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::sum(T alpha, bool conja, const DistTensor<T>& A, const int *idx_A,
                        T  beta,                                     const int *idx_B)
{
    int ret;
    CTF_sum_type_t st;

    st.tid_A = A.tid;
    st.tid_B =   tid;
    conv_idx(    A.ndim, idx_A, st.idx_map_A,
                   ndim, idx_B, st.idx_map_B);

    ret = arena.ctf.ctf->sum_tensors(&st, alpha, beta);

    FREE(st.idx_map_A);
    FREE(st.idx_map_B);
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
void DistTensor<T>::scale(T alpha, const int* idx_A)
{
    int ret;
    int * idx_map_A;

    conv_idx(ndim, idx_A, idx_map_A);

    ret = arena.ctf.ctf->scale_tensor(alpha, tid, idx_map_A);

    FREE(idx_map_A);
    assert(ret == DIST_TENSOR_SUCCESS);
}

template <typename T>
T DistTensor<T>::dot(bool conja, const DistTensor<T>& A, const int* idx_A,
                     bool conjb,                         const int* idx_B) const
{
    DistTensor<T> dt(A.arena);
    int64_t n;
    vector<T> val;
    dt.mult(1, conja,     A, idx_A,
               conjb, *this, idx_B,
            0,                NULL);
    dt.getAllData(val);
    assert(val.size()==1);
    return val[0];
}

INSTANTIATE_SPECIALIZATIONS(DistTensor);
