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

template <class T>
SymmetryBlockedTensor<T>::SymmetryBlockedTensor(tCTF_World<T>& ctf, const PointGroup& group)
: IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>(0, 0), Distributed<T>(ctf),
  group_(group), len_(), sym_()
{
    tensors_.resize(1, NULL);
    tensors_[0] = new DistTensor<T>(ctf);
}

template <class T>
SymmetryBlockedTensor<T>::SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& other)
: IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>(other), Distributed<T>(other.ctf),
  group_(other.group_), len_(other.len_), sym_(other.sym_) {}

template <class T>
SymmetryBlockedTensor<T>::SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& other, T scalar)
: IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>(0, 1), Distributed<T>(other.ctf),
  group_(other.group_), len_(0), sym_(0)
{
    tensors_[0].tensor_ = new DistTensor<T>(other.ctf, 0, NULL, NULL);
    *tensors_[0].tensor_ = scalar;
}

template <class T>
SymmetryBlockedTensor<T>::SymmetryBlockedTensor(tCTF_World<T>& ctf, const PointGroup& group,
                                                const int ndim, const int **len, const int *sym, const bool zero)
: IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>(ndim, 0), Distributed<T>(ctf),
  group_(group), len_(ndim), sym_(ndim)
{
    int n = group.getNumIrreps();
    const Representation* irreps = group.getIrreps();

    sym_.assign(sym, sym+ndim_);

    int ntensors = 1;
    vector<int> sublen(ndim_);

    for (int i = 0;i < ndim_;i++)
    {
        ntensors *= n;
        len_[i].assign(len[i], len[i]+n);
        sublen[i] = len_[i][0];
    }

    tensors_.resize(ntensors);

    int t = 0;
    vector<int> idx(ndim_, 0);
    vector<Representation> prod(ndim_+1, group.totallySymmetricIrrep());
    for (bool done = false;!done;t++)
    {
        cout << idx << prod[0].isTotallySymmetric() << endl;

        if (prod[0].isTotallySymmetric())
        {
            vector<int> subsym(sym, sym+n);

            bool ok = true;
            for (int i = 0;i < ndim_-1;i++)
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

            if (ok) tensors_[t].tensor_ = new DistTensor<T>(ctf, ndim, sublen.data(), subsym.data(), zero);
        }

        for (int i = 0;i < ndim_;i++)
        {
            idx[i]++;

            if (idx[i] < n)
            {
                idx[i] = 0;
                sublen[i] = len_[i][idx[i]];

                for (int j = i;j >= 0;j--)
                {
                    prod[j] = irreps[idx[j]]*prod[j+1];
                }

                break;
            }
            else if (i == ndim_-1)
            {
                done = true;
            }
            else
            {
                sublen[i] = len_[i][idx[i]];
            }
        }

        if (ndim_ == 0) done = true;
    }
}

template <class T>
vector<int> SymmetryBlockedTensor<T>::getStrides(const vector<int> indices, const int ndim,
                                                 const int len, const int* idx_A)
{
    vector<int> strides(indices.size(), 0);
    vector<int> stride_A(ndim);

    stride_A[0] = 1;
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
void SymmetryBlockedTensor<T>::mult(const double alpha, bool conja, const SymmetryBlockedTensor<T>& A, const int* idx_A_,
                                                        bool conjb, const SymmetryBlockedTensor<T>& B, const int* idx_B_,
                                     const double beta,                                                const int* idx_C_)
{
    assert(group_ == A.group_);
    assert(group_ == B.group_);

    int n = group_.getNumIrreps();

    vector<int> idx_A(idx_A_, idx_A_+A.ndim_);
    vector<int> idx_B(idx_B_, idx_B_+B.ndim_);
    vector<int> idx_C(idx_C_, idx_C_+  ndim_);

    double f1 = align_symmetric_indices(A.ndim_, idx_A, A.sym_.data(),
                                        B.ndim_, idx_B, B.sym_.data(),
                                          ndim_, idx_C,   sym_.data());
    f1 *= overcounting_factor(A.ndim_, idx_A, A.sym_.data(),
                              B.ndim_, idx_B, B.sym_.data(),
                                ndim_, idx_C,   sym_.data());

    vector<int> inds_AB = idx_A+idx_B;
    vector<int> inds_C = idx_C;

    uniq(inds_AB);
    uniq(inds_C);
    exclude(inds_AB, inds_C);

    int nAB = inds_AB.size();
    int nC = inds_C.size();

    vector<int> stride_A_AB = getStrides(inds_AB, A.ndim_, n, idx_A.data());
    vector<int> stride_B_AB = getStrides(inds_AB, B.ndim_, n, idx_B.data());
    vector<int> stride_A_C = getStrides(inds_C, A.ndim_, n, idx_A.data());
    vector<int> stride_B_C = getStrides(inds_C, B.ndim_, n, idx_B.data());
    vector<int> stride_C = getStrides(inds_C, ndim_, n, idx_C.data());

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
            if (A.tensors_[off_A] != NULL &&
                B.tensors_[off_B] != NULL &&
                  tensors_[off_C] != NULL)
            {
                double f2 = overcounting_factor(A.ndim_, idx_A, A.tensors_[off_A].tensor_->getSymmetry(),
                                                B.ndim_, idx_B, B.tensors_[off_B].tensor_->getSymmetry(),
                                                  ndim_, idx_C,   tensors_[off_C].tensor_->getSymmetry());

                tensors_[off_C].tensor_->mult(alpha*f1/f2, conja, *A.tensors_[off_A].tensor_, idx_A.data(),
                                                           conjb, *B.tensors_[off_B].tensor_, idx_B.data(),
                                              beta_      ,                                    idx_C.data());
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

            beta_ = 1.0;
        }

        for (int i = 0;i < nC;i++)
        {
            iC[i]++;
            off_A += stride_A_C[i];
            off_B += stride_B_C[i];
            off_C += stride_C[i];

            if (iC[i] == n)
            {
                iC[i] = 0;
                off_A -= stride_A_C[i]*n;
                off_B -= stride_B_C[i]*n;
                off_C -= stride_C[i]*n;
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
void SymmetryBlockedTensor<T>::sum(const double alpha, bool conja, const SymmetryBlockedTensor<T>& A, const int* idx_A_,
                                    const double beta,                                                const int* idx_B_)
{
    assert(group_ == A.group_);

    int n = group_.getNumIrreps();

    vector<int> idx_A(idx_A_, idx_A_+A.ndim_);
    vector<int> idx_B(idx_B_, idx_B_+  ndim_);

    double f = align_symmetric_indices(A.ndim_, idx_A, A.sym_.data(),
                                         ndim_, idx_B,   sym_.data());

    vector<int> inds_A = idx_A;
    vector<int> inds_B = idx_B;

    uniq(inds_A);
    uniq(inds_B);
    exclude(inds_A, inds_B);

    int nA = inds_A.size();
    int nB = inds_B.size();

    vector<int> stride_A_A = getStrides(inds_A, A.ndim_, n, idx_A.data());
    vector<int> stride_A_B = getStrides(inds_B, A.ndim_, n, idx_A.data());
    vector<int> stride_B = getStrides(inds_B, ndim_, n, idx_B.data());

    int off_A = 0;
    int off_B = 0;
    vector<int> iB(nB, 0);
    for (bool doneB = false;!doneB;)
    {
        double beta_ = beta;

        vector<int> iA(nA, 0);
        for (bool doneA = false;!doneA;)
        {
            if (A.tensors_[off_A] != NULL &&
                  tensors_[off_B] != NULL)
            {
                tensors_[off_B].tensor_->sum(alpha*f, conja, *A.tensors_[off_A].tensor_, idx_A.data(),
                                             beta_  ,                                    idx_B.data());
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
void SymmetryBlockedTensor<T>::scale(const double alpha, const int* idx_A)
{
    int n = group_.getNumIrreps();

    vector<int> inds_A(idx_A, idx_A+ndim_);

    uniq(inds_A);

    int nA = inds_A.size();

    vector<int> stride_A = getStrides(inds_A, ndim_, n, idx_A);

    int off_A = 0;
    vector<int> iA(nA, 0);
    for (bool doneA = false;!doneA;)
    {
        if (tensors_[off_A] != NULL)
        {
            tensors_[off_A].tensor_->scale(alpha, idx_A);
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
T SymmetryBlockedTensor<T>::dot(bool conja, const SymmetryBlockedTensor<T>& A, const int* idx_A_,
                                bool conjb,                                    const int* idx_B_) const
{
    assert(group_ == A.group_);

    T sum = (T)0;
    int n = group_.getNumIrreps();

    vector<int> idx_A(idx_A_, idx_A_+A.ndim_);
    vector<int> idx_B(idx_B_, idx_B_+  ndim_);
    vector<int> nll;

    double f1 = align_symmetric_indices(A.ndim_, idx_A, A.sym_.data(),
                                          ndim_, idx_B,   sym_.data(),
                                              0,   nll,          NULL);
    f1 *= overcounting_factor(A.ndim_, idx_A, A.sym_.data(),
                                ndim_, idx_B,   sym_.data(),
                                    0,   nll,          NULL);

    vector<int> inds_AB = idx_A+idx_B;

    uniq(inds_AB);
    int nAB = inds_AB.size();

    vector<int> stride_A_AB = getStrides(inds_AB, A.ndim_, n, idx_A.data());
    vector<int> stride_B_AB = getStrides(inds_AB,   ndim_, n, idx_B.data());

    int off_A = 0;
    int off_B = 0;
    vector<int> iAB(nAB, 0);
    for (bool doneAB = false;!doneAB;)
    {
        if (A.tensors_[off_A] != NULL &&
              tensors_[off_B] != NULL)
        {
            double f2 = overcounting_factor(A.ndim_, idx_A, A.tensors_[off_A].tensor_->getSymmetry(),
                                              ndim_, idx_B,   tensors_[off_B].tensor_->getSymmetry(),
                                                  0,   nll,                                     NULL);

            sum += (f1/f2)*tensors_[off_B].tensor_->dot(conja, *A.tensors_[off_A].tensor_, idx_A.data(),
                                                        conjb,                             idx_B.data());
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

INSTANTIATE_SPECIALIZATIONS(SymmetryBlockedTensor);
