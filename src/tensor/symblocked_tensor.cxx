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

#include "symblocked_tensor.hpp"

using namespace std;
using namespace aquarius::tensor;
using namespace aquarius::symmetry;
using namespace aquarius::task;

template <class T>
SymmetryBlockedTensor<T>::SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& other)
: IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>(other), Resource(other.arena),
  group(other.group), len(other.len), sym(other.sym) {}

template <class T>
SymmetryBlockedTensor<T>::SymmetryBlockedTensor(SymmetryBlockedTensor<T>* other)
: IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>(other->ndim, 0), Resource(other->arena),
  group(other->group), len(other->len), sym(other->sym)
{
    tensors.swap(other->tensors);
    delete other;
}

template <class T>
SymmetryBlockedTensor<T>::SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& other, T scalar)
: IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>(0, 0), Resource(other.arena),
  group(other.group), len(0), sym(0)
{
    tensors.resize(1, NULL);
    tensors[0].tensor = new DistTensor<T>(other.arena, scalar);
}

template <class T>
SymmetryBlockedTensor<T>::SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& A,
                                                const vector<vector<int> >& start_A,
                                                const vector<vector<int> >& len_A)
: IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>(A.ndim, 0), Resource(A.arena),
  group(A.group), len(len_A), sym(A.sym)
{
    allocate(false);
    slice((T)1, false, A, start_A, (T)0);
}

template <class T>
SymmetryBlockedTensor<T>::SymmetryBlockedTensor(const Arena& arena, const PointGroup& group,
                                                int ndim, const vector<vector<int> >& len,
                                                const vector<int>& sym, bool zero)
: IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>(ndim, 0), Resource(arena),
  group(group), len(len), sym(sym)
{
    assert(sym.size() == ndim);
    assert(len.size() == ndim);
    allocate(zero);
}

template <class T>
void SymmetryBlockedTensor<T>::allocate(bool zero)
{
    int n = group.getNumIrreps();
    vector<Representation> irreps;
    for (int i = 0;i < n;i++) irreps.push_back(group.getIrrep(i));

    int ntensors = 1;
    vector<int> sublen(ndim);

    for (int i = 0;i < ndim;i++)
    {
        assert(len[i].size() == n);
        ntensors *= n;
        sublen[i] = len[i][0];
    }

    tensors.resize(ntensors);

    int t = 0;
    vector<int> idx(ndim, 0);
    vector<Representation> prod(ndim+1, group.totallySymmetricIrrep());
    for (bool done = false;!done;t++)
    {
        if (prod[0].isTotallySymmetric())
        {
            vector<int> subsym(sym);

            bool ok = true;
            for (int i = 0;i < ndim-1;i++)
            {
                if (sym[i] != NS)
                {
                    if (idx[i] < idx[i+1])
                    {
                        subsym[i] = NS;
                    }
                    else if (idx[i] > idx[i+1])
                    {
                        ok = false;
                    }
                }
            }

            assert(t < ntensors);
            if (ok) tensors[t].tensor = new DistTensor<T>(arena, ndim, sublen, subsym, zero);
        }

        for (int i = 0;i < ndim;i++)
        {
            idx[i] = (idx[i] == n-1 ? 0 : idx[i]+1);
            sublen[i] = len[i][idx[i]];

            if (idx[i] != 0)
            {
                for (int j = i;j >= 0;j--)
                {
                    prod[j] = prod[j+1];
                    prod[j] *= irreps[idx[j]];
                }
                break;
            }
            else if (i == ndim-1)
            {
                done = true;
            }
        }

        if (ndim == 0) done = true;
    }
}

template <class T>
DistTensor<T>& SymmetryBlockedTensor<T>::operator()(const vector<int>& irreps)
{
    return const_cast<DistTensor<T>&>(const_cast<const SymmetryBlockedTensor<T>&>(*this)(irreps));
}

template <class T>
const DistTensor<T>& SymmetryBlockedTensor<T>::operator()(const vector<int>& irreps) const
{
    assert(irreps.size() == this->ndim);

    int n = group.getNumIrreps();

    int off = 0;
    int stride = 1;
    for (int i = 0;i < ndim;i++)
    {
        off += stride*irreps[i];
        stride *= n;
    }

    assert(tensors[off].tensor != NULL);

    return *tensors[off].tensor;
}

template <class T>
bool SymmetryBlockedTensor<T>::exists(const vector<int>& irreps) const
{
    assert(irreps.size() == this->ndim);

    int n = group.getNumIrreps();

    int off = 0;
    int stride = 1;
    for (int i = 0;i < ndim;i++)
    {
        off += stride*irreps[i];
        stride *= n;
    }

    return tensors[off].tensor != NULL;
}

template <class T>
void SymmetryBlockedTensor<T>::slice(T alpha, bool conja, const SymmetryBlockedTensor<T>& A,
                                     const vector<vector<int> >& start_A, T beta)
{
    int n = group.getNumIrreps();
    slice(alpha, conja, A, start_A, beta, vector<vector<int> >(this->ndim,vector<int>(n,0)), len);
}

template <class T>
void SymmetryBlockedTensor<T>::slice(T alpha, bool conja, const SymmetryBlockedTensor<T>& A,
                                     T beta, const vector<vector<int> >& start_B)
{
    int n = group.getNumIrreps();
    slice(alpha, conja, A, vector<vector<int> >(this->ndim,vector<int>(n,0)), beta, start_B, A.len);
}

template <class T>
void SymmetryBlockedTensor<T>::slice(T alpha, bool conja, const SymmetryBlockedTensor<T>& A,
                                                          const vector<vector<int> >& start_A,
                                     T  beta,             const vector<vector<int> >& start_B,
                                                          const vector<vector<int> >& len)
{
    assert(this->ndim == A.ndim);

    int n = group.getNumIrreps();

    vector<vector<int> > end_A(ndim);
    vector<vector<int> > end_B(ndim);

    for (int i = 0;i < ndim;i++)
    {
        end_A[i].resize(n);
        end_B[i].resize(n);

        assert(sym[i] == A.sym[i] && sym[i] == NS);
        for (int j = 0;j < n;j++)
        {
            end_A[i][j] = start_A[i][j]+len[i][j];
            end_B[i][j] = start_B[i][j]+len[i][j];
            assert(start_A[i][j] >= 0);
            assert(start_B[i][j] >= 0);
            assert(end_A[i][j] <= A.len[i][j]);
            assert(end_B[i][j] <= this->len[i][j]);
        }
    }

    vector<int> stride(ndim,1);
    for (int i = 1;i < ndim;i++) stride[i] = stride[i-1]*n;

    int off_A = 0;
    vector<int> iA(ndim, 0);
    vector<int> start_A_sub(ndim);
    vector<int> start_B_sub(ndim);
    vector<int> len_sub(ndim);
    for (bool doneA = false;!doneA;)
    {
        if (tensors[off_A] != NULL && A.tensors[off_A] != NULL)
        {
            for (int i = 0;i < ndim;i++)
            {
                start_A_sub[i] = start_A[i][iA[i]];
                start_B_sub[i] = start_B[i][iA[i]];
                len_sub[i] = len[i][iA[i]];
            }

            tensors[off_A].tensor->slice(alpha, conja, *A.tensors[off_A].tensor, start_A_sub,
                                          beta,                                  start_B_sub, len_sub);
        }

        for (int i = 0;i < ndim;i++)
        {
            iA[i]++;
            off_A += stride[i];

            if (iA[i] == n)
            {
                iA[i] = 0;
                off_A -= stride[i]*n;
                if (i == ndim-1) doneA = true;
            }
            else
            {
                break;
            }
        }

        if (this->ndim == 0) doneA = true;
    }
}

template <class T>
vector<int> SymmetryBlockedTensor<T>::getStrides(const string& indices, const int ndim,
                                                 const int len, const string& idx_A)
{
    vector<int> strides(indices.size(), 0);
    vector<int> stride_A(ndim);

    if (ndim > 0) stride_A[0] = 1;
    for (int i = 1;i < ndim;i++)
    {
        stride_A[i] = stride_A[i-1]*len;
    }

    for (int i = 0;i < indices.size();i++)
    {
        for (int j = 0;j < ndim;j++)
        {
            if (indices[i] == idx_A[j])
            {
                strides[i] += stride_A[j];
            }
        }
    }

    return strides;
}

template <class T>
void SymmetryBlockedTensor<T>::mult(T alpha, bool conja, const SymmetryBlockedTensor<T>& A, const string& idx_A,
                                             bool conjb, const SymmetryBlockedTensor<T>& B, const string& idx_B,
                                    T  beta,                                                const string& idx_C)
{
    assert(group == A.group);
    assert(group == B.group);

    int n = group.getNumIrreps();

    string idx_A_ = idx_A;
    string idx_B_ = idx_B;
    string idx_C_ = idx_C;

    double f1 = align_symmetric_indices(A.ndim, idx_A_, A.sym.data(),
                                        B.ndim, idx_B_, B.sym.data(),
                                          ndim, idx_C_,   sym.data());
    f1 *= overcounting_factor(A.ndim, idx_A_, A.sym.data(),
                              B.ndim, idx_B_, B.sym.data(),
                                ndim, idx_C_,   sym.data());

    string inds_AB = idx_A_+idx_B_;
    string inds_C = idx_C_;

    uniq(inds_AB);
    uniq(inds_C);
    exclude(inds_AB, inds_C);

    int nAB = inds_AB.size();
    int nC = inds_C.size();

    vector<int> stride_A_AB = getStrides(inds_AB, A.ndim, n, idx_A_);
    vector<int> stride_B_AB = getStrides(inds_AB, B.ndim, n, idx_B_);
    vector<int> stride_A_C = getStrides(inds_C, A.ndim, n, idx_A_);
    vector<int> stride_B_C = getStrides(inds_C, B.ndim, n, idx_B_);
    vector<int> stride_C = getStrides(inds_C, ndim, n, idx_C_);

    int off_A = 0;
    int off_B = 0;
    int off_C = 0;
    vector<int> iC(nC, 0);
    for (bool doneC = false;!doneC;)
    {
        double beta_ = beta;

        vector<int> iAB(nAB, 0);
        for (bool doneAB = false;!doneAB;)
        {
            if (A.tensors[off_A] != NULL &&
                B.tensors[off_B] != NULL &&
                  tensors[off_C] != NULL)
            {
                double f2 = overcounting_factor(A.ndim, idx_A_, A.tensors[off_A].tensor->getSymmetry().data(),
                                                B.ndim, idx_B_, B.tensors[off_B].tensor->getSymmetry().data(),
                                                  ndim, idx_C_,   tensors[off_C].tensor->getSymmetry().data());

                tensors[off_C].tensor->mult(alpha*f1/f2, conja, *A.tensors[off_A].tensor, idx_A_,
                                                         conjb, *B.tensors[off_B].tensor, idx_B_,
                                            beta_      ,                                  idx_C_);

                beta_ = 1.0;
            }

            for (int i = 0;i < nAB;i++)
            {
                iAB[i]++;
                off_A += stride_A_AB[i];
                off_B += stride_B_AB[i];

                if (iAB[i] >= n)
                {
                    off_A -= stride_A_AB[i]*iAB[i];
                    off_B -= stride_B_AB[i]*iAB[i];
                    iAB[i] = 0;
                    if (i == nAB-1) doneAB = true;
                }
                else
                {
                    break;
                }
            }

            if (nAB == 0) doneAB = true;
        }

        for (int i = 0;i < nC;i++)
        {
            iC[i]++;
            off_A += stride_A_C[i];
            off_B += stride_B_C[i];
            off_C += stride_C[i];

            if (iC[i] >= n)
            {
                off_A -= stride_A_C[i]*iC[i];
                off_B -= stride_B_C[i]*iC[i];
                off_C -= stride_C[i]*iC[i];
                iC[i] = 0;
                if (i == nC-1) doneC = true;
            }
            else
            {
                break;
            }
        }

        if (nC == 0) doneC = true;
    }
}

template <class T>
void SymmetryBlockedTensor<T>::sum(T alpha, bool conja, const SymmetryBlockedTensor<T>& A, const string& idx_A,
                                   T  beta,                                                const string& idx_B)
{
    assert(group == A.group);

    int n = group.getNumIrreps();

    string idx_A_ = idx_A;
    string idx_B_ = idx_B;

    double f = align_symmetric_indices(A.ndim, idx_A_, A.sym.data(),
                                         ndim, idx_B_,   sym.data());

    string inds_A = idx_A_;
    string inds_B = idx_B_;

    uniq(inds_A);
    uniq(inds_B);
    exclude(inds_A, inds_B);

    int nA = inds_A.size();
    int nB = inds_B.size();

    vector<int> stride_A_A = getStrides(inds_A, A.ndim, n, idx_A_);
    vector<int> stride_A_B = getStrides(inds_B, A.ndim, n, idx_A_);
    vector<int> stride_B = getStrides(inds_B, ndim, n, idx_B_);

    int off_A = 0;
    int off_B = 0;
    vector<int> iB(nB, 0);
    for (bool doneB = false;!doneB;)
    {
        double beta_ = beta;

        vector<int> iA(nA, 0);
        for (bool doneA = false;!doneA;)
        {
            if (A.tensors[off_A] != NULL &&
                  tensors[off_B] != NULL)
            {
                tensors[off_B].tensor->sum(alpha*f, conja, *A.tensors[off_A].tensor, idx_A_,
                                             beta_,                                  idx_B_);
            }

            for (int i = 0;i < nA;i++)
            {
                iA[i]++;
                off_A += stride_A_A[i];

                if (iA[i] == n)
                {
                    iA[i] = 0;
                    off_A -= stride_A_A[i]*n;
                    if (i == nA-1) doneA = true;
                }
                else
                {
                    break;
                }
            }

            if (nA == 0) doneA = true;

            beta_ = 1.0;
        }

        for (int i = 0;i < nB;i++)
        {
            iB[i]++;
            off_A += stride_A_B[i];
            off_B += stride_B[i];

            if (iB[i] == n)
            {
                iB[i] = 0;
                off_A -= stride_A_B[i]*n;
                off_B -= stride_B[i]*n;
                if (i == nB-1) doneB = true;
            }
            else
            {
                break;
            }
        }

        if (nB == 0) doneB = true;
    }
}

template <class T>
void SymmetryBlockedTensor<T>::scale(T alpha, const string& idx_A)
{
    int n = group.getNumIrreps();

    string inds_A = idx_A;

    uniq(inds_A);

    int nA = inds_A.size();

    vector<int> stride_A = getStrides(inds_A, ndim, n, idx_A);

    int off_A = 0;
    vector<int> iA(nA, 0);
    for (bool doneA = false;!doneA;)
    {
        if (tensors[off_A] != NULL)
        {
            tensors[off_A].tensor->scale(alpha, idx_A);
        }

        for (int i = 0;i < nA;i++)
        {
            iA[i]++;
            off_A += stride_A[i];

            if (iA[i] == n)
            {
                iA[i] = 0;
                off_A -= stride_A[i]*n;
                if (i == nA-1) doneA = true;
            }
            else
            {
                break;
            }
        }

        if (nA == 0) doneA = true;
    }
}

template <class T>
T SymmetryBlockedTensor<T>::dot(bool conja, const SymmetryBlockedTensor<T>& A, const string& idx_A,
                                bool conjb,                                    const string& idx_B) const
{
    assert(group == A.group);

    T sum = (T)0;
    int n = group.getNumIrreps();

    string idx_A_ = idx_A;
    string idx_B_ = idx_B;

    string nll;
    double f1 = align_symmetric_indices(A.ndim, idx_A_, A.sym.data(),
                                          ndim, idx_B_,   sym.data(),
                                              0,   nll,         NULL);
    f1 *= overcounting_factor(A.ndim, idx_A_, A.sym.data(),
                                ndim, idx_B_,   sym.data(),
                                    0,   nll,         NULL);

    string inds_AB = idx_A_+idx_B_;

    uniq(inds_AB);
    int nAB = inds_AB.size();

    vector<int> stride_A_AB = getStrides(inds_AB, A.ndim, n, idx_A_);
    vector<int> stride_B_AB = getStrides(inds_AB,   ndim, n, idx_B_);

    int off_A = 0;
    int off_B = 0;
    vector<int> iAB(nAB, 0);
    for (bool doneAB = false;!doneAB;)
    {
        if (A.tensors[off_A] != NULL &&
              tensors[off_B] != NULL)
        {
            double f2 = overcounting_factor(A.ndim, idx_A_, A.tensors[off_A].tensor->getSymmetry().data(),
                                              ndim, idx_B_,   tensors[off_B].tensor->getSymmetry().data(),
                                                 0,    nll,                                          NULL);

            sum += (f1/f2)*tensors[off_B].tensor->dot(conja, *A.tensors[off_A].tensor, idx_A_,
                                                      conjb,                           idx_B_);
        }

        for (int i = 0;i < nAB;i++)
        {
            iAB[i]++;
            off_A += stride_A_AB[i];
            off_B += stride_B_AB[i];

            if (iAB[i] == n)
            {
                iAB[i] = 0;
                off_A -= stride_A_AB[i]*n;
                off_B -= stride_B_AB[i]*n;
                if (i == nAB-1) doneAB = true;
            }
            else
            {
                break;
            }
        }

        if (nAB == 0) doneAB = true;
    }

    return sum;
}

template <class T>
void SymmetryBlockedTensor<T>::weight(const vector<const vector<vector<T> >*>& d)
{
    int n = group.getNumIrreps();

    vector<int> stride(ndim,1);
    for (int i = 1;i < ndim;i++) stride[i] = stride[i-1]*n;

    int off_A = 0;
    vector<int> iA(this->ndim, 0);
    vector<const vector<T>*> dsub(this->ndim);
    for (bool doneA = false;!doneA;)
    {
        if (tensors[off_A] != NULL)
        {
            for (int i = 0;i < this->ndim;i++) dsub[i] = &((*d[i])[iA[i]]);
            tensors[off_A].tensor->weight(dsub);
        }

        for (int i = 0;i < this->ndim;i++)
        {
            iA[i]++;
            off_A += stride[i];

            if (iA[i] == n)
            {
                iA[i] = 0;
                off_A -= stride[i]*n;
                if (i == this->ndim-1) doneA = true;
            }
            else
            {
                break;
            }
        }

        if (ndim == 0) doneA = true;
    }
}

template <class T>
typename std::real_type<T>::type SymmetryBlockedTensor<T>::norm(int p) const
{
    typename std::real_type<T>::type nrm = 0;

    int n = group.getNumIrreps();

    vector<int> stride(ndim,1);
    for (int i = 1;i < ndim;i++) stride[i] = stride[i-1]*n;

    int off_A = 0;
    vector<int> iA(this->ndim, 0);
    vector<const vector<T>*> dsub(this->ndim);
    for (bool doneA = false;!doneA;)
    {
        if (tensors[off_A] != NULL)
        {
            double factor = 1;
            const vector<int>& subsym = tensors[off_A].tensor->getSymmetry();
            for (int i = 0;i < ndim;)
            {
                int j; for (j = i;sym[j] != NS;j++); j++;

                int m = j-i;
                for (int k = i;k < j;)
                {
                    int l; for (l = k;subsym[l] != NS;l++); l++;
                    int o = l-k;
                    factor *= binom(m,o);
                    m -= o;
                    k = l;
                }

                i = j;
            }

            if      (p == 2) factor = sqrt(factor);
            else if (p == 0) factor = 1;

            nrm += factor*tensors[off_A].tensor->norm(p);
        }

        for (int i = 0;i < this->ndim;i++)
        {
            iA[i]++;
            off_A += stride[i];

            if (iA[i] == n)
            {
                iA[i] = 0;
                off_A -= stride[i]*n;
                if (i == this->ndim-1) doneA = true;
            }
            else
            {
                break;
            }
        }

        if (ndim == 0) doneA = true;
    }

    return nrm;
}

INSTANTIATE_SPECIALIZATIONS(SymmetryBlockedTensor);
