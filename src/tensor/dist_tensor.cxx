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

/*
 * Create a scalar (0-dimensional tensor)
 */
template <typename T>
DistTensor<T>::DistTensor(const Arena& arena, T scalar)
: IndexableTensor< DistTensor<T>,T >(), Resource(arena), len(0), sym(0)
{
    allocate();
    *dt = scalar;
}

/*
 * Create a scalar (0-dimensional tensor) on the same arena as A
 */
template <typename T>
DistTensor<T>::DistTensor(const DistTensor<T>& A, T scalar)
: IndexableTensor< DistTensor<T>,T >(), Resource(A.arena), len(0), sym(0)
{
    allocate();
    *dt = scalar;
}

/*
 * Create a tensor of the same size and shape as A, optionally copying or zeroing the data
 */
template <typename T>
DistTensor<T>::DistTensor(const DistTensor<T>& A, bool copy, bool zero)
: IndexableTensor< DistTensor<T>,T >(A.ndim), Resource(A.arena), len(A.len), sym(A.sym)
{
    allocate();

    if (copy)
    {
        *this = A;
    }
    else if (zero)
    {
        *dt = (T)0;
    }
}

/*
 * Create a tensor of the specified size and shape, optionally zeroing the data
 */
template <typename T>
DistTensor<T>::DistTensor(const Arena& arena, int ndim, const vector<int>& len, const vector<int>& sym,
                          bool zero)
: IndexableTensor< DistTensor<T>,T >(ndim), Resource(arena), len(len), sym(sym)
{
    assert(len.size() == ndim);
    assert(sym.size() == ndim);

    #ifdef VALIDATE_INPUTS
    validate_tensor(ndim, len.data(), NULL, sym.data());
    #endif //VALIDATE_INPUTS

    allocate();
    if (zero) *dt = (T)0;
}

template <typename T>
DistTensor<T>::~DistTensor()
{
    free();
}

template <typename T>
void DistTensor<T>::allocate()
{
    dt = new tCTF_Tensor<T>(ndim, len.data(), sym.data(), arena.ctf<T>(), "FIXME", 1);
}

template <typename T>
void DistTensor<T>::free()
{
    delete dt;
}

template <typename T>
void DistTensor<T>::resize(int ndim, const vector<int>& len, const vector<int>& sym, bool zero)
{
    assert(len.size() == ndim);
    assert(sym.size() == ndim);

    this->ndim = ndim;
    this->len = len;
    this->sym = sym;

    free();
    allocate();
    if (zero) *dt = (T)0;
}

template <typename T>
T* DistTensor<T>::getRawData(int64_t& size)
{
    return const_cast<T*>(const_cast<const DistTensor<T>&>(*this).getRawData(size));
}

template <typename T>
const T* DistTensor<T>::getRawData(int64_t& size) const
{
    long_int size_;
    T* data = dt->get_raw_data(&size_);
    size = size_;
    return data;
}

template <typename T>
void DistTensor<T>::getLocalData(vector<tkv_pair<T> >& pairs) const
{
    int64_t npair;
    tkv_pair<T> *data;
    dt->get_local_data(&npair, &data);
    pairs.assign(data, data+npair);
    if (npair > 0) ::free(data);
}

template <typename T>
void DistTensor<T>::getRemoteData(vector<tkv_pair<T> >& pairs) const
{
    dt->get_remote_data(pairs.size(), pairs.data());
}

template <typename T>
void DistTensor<T>::getRemoteData() const
{
    dt->get_remote_data(0, NULL);
}

template <typename T>
void DistTensor<T>::writeRemoteData(const vector<tkv_pair<T> >& pairs)
{
    dt->write_remote_data(pairs.size(), pairs.data());
}

template <typename T>
void DistTensor<T>::writeRemoteData()
{
    dt->write_remote_data(0, NULL);
}

template <typename T>
void DistTensor<T>::writeRemoteData(double alpha, double beta, const vector<tkv_pair<T> >& pairs)
{
    dt->add_remote_data(pairs.size(), alpha, beta, pairs.data());
}

template <typename T>
void DistTensor<T>::writeRemoteData(double alpha, double beta)
{
    dt->add_remote_data(0, alpha, beta, NULL);
}

template <typename T>
void DistTensor<T>::getAllData(vector<T>& vals) const
{

    getAllData(vals, 0);
    int64_t npair = vals.size();
    arena.Bcast(&npair, 1, 0);
    if (rank != 0) vals.resize(npair);
    arena.Bcast(vals, 0);
}

template <typename T>
void DistTensor<T>::getAllData(vector<T>& vals, int rank) const
{
    if (arena.rank == rank)
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

        dt->get_remote_data(pairs.size(), pairs.data());

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
        dt->get_remote_data(0, NULL);
    }
}

template <typename T>
void DistTensor<T>::div(T alpha, bool conja, const DistTensor<T>& A,
                                 bool conjb, const DistTensor<T>& B, T beta)
{
    dt->align(*A.dt);
    dt->align(*B.dt);
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
    dt->align(*A.dt);
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
    dt->print(fp, cutoff);
}

template <typename T>
void DistTensor<T>::compare(FILE* fp, const DistTensor<T>& other, double cutoff) const
{
    dt->compare(*other.dt, fp, cutoff);
}

template <typename T>
typename real_type<T>::type DistTensor<T>::norm(int p) const
{
    T ans = (T)0;
    if (p == 00)
    {
        ans = dt->reduce(CTF_OP_MAXABS);
    }
    else if (p == 1)
    {
        ans = dt->reduce(CTF_OP_SUMABS);
    }
    else if (p == 2)
    {
        ans = sqrt(dt->reduce(CTF_OP_SQNRM2));
    }
    return abs(ans);
}

template <typename T>
void DistTensor<T>::mult(T alpha, bool conja, const DistTensor<T>& A, const string& idx_A,
                                  bool conjb, const DistTensor<T>& B, const string& idx_B,
                         T  beta,                                     const string& idx_C)
{
    dt->contract(alpha, *A.dt, idx_A.c_str(),
                        *B.dt, idx_B.c_str(),
                  beta,        idx_C.c_str());
}

template <typename T>
void DistTensor<T>::sum(T alpha, bool conja, const DistTensor<T>& A, const string& idx_A,
                        T  beta,                                     const string& idx_B)
{
    dt->sum(alpha, *A.dt, idx_A.c_str(),
             beta,        idx_B.c_str());
}

template <typename T>
void DistTensor<T>::scale(T alpha, const string& idx_A)
{
    dt->scale(alpha, idx_A.c_str());
}

template <typename T>
T DistTensor<T>::dot(bool conja, const DistTensor<T>& A, const string& idx_A,
                     bool conjb,                         const string& idx_B) const
{
    DistTensor<T> dt(A.arena);
    vector<T> val;
    dt.mult(1, conja,     A, idx_A,
               conjb, *this, idx_B,
            0,                  "");
    dt.getAllData(val);
    assert(val.size()==1);
    return val[0];
}

INSTANTIATE_SPECIALIZATIONS(DistTensor);
