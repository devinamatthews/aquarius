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

#ifndef _AQUARIUS_TENSOR_SYMBLOCKED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_SYMBLOCKED_TENSOR_HPP_

#include <vector>
#include <iostream>
#include <cassert>

#include "../src/dist_tensor/sym_indices.hxx"

#include "stl_ext/stl_ext.hpp"

#include "dist_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <class T>
class SymmetryBlockedTensor : public IndexableTensor< SymmetryBlockedTensor<T> >, public Distributed<T>
{
    INHERIT_FROM_INDEXABLE_TENSOR(SymmetryBlockedTensor<T>, T);

    protected:
        const symmetry::PointGroup& group_;
        std::vector< std::vector<int> > len_;
        std::vector<int> sym_;
        std::vector< DistTensor<T>* > tensors;

    public:
        SymmetryBlockedTensor(tCTF_World<T>& ctf, const symmetry::PointGroup& group)
        : IndexableTensor< SymmetryBlockedTensor<T> >(), Distributed<T>(ctf),
          group_(group), len_(), sym_(), tensors()
        {
            tensors.resize(1, NULL);
            tensors[0] = new DistTensor<T>(ctf);
        }

        SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& other)
        : IndexableTensor< SymmetryBlockedTensor<T> >(other.ndim_), Distributed(other.ctf),
          group_(other.group_), len_(other.len_), sym_(other.sym_), tensors(other.tensors)
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL) tensors[i] = new DistTensor<T>(*tensors[i]);
            }
        }

        SymmetryBlockedTensor(tCTF_World<T>& ctf, const symmetry::PointGroup& group,
                              const int ndim, const int **len, const int *sym, const bool zero=true)
        : IndexableTensor< SymmetryBlockedTensor<T> >(ndim), Distributed<T>(ctf),
          group_(group), len_(ndim), sym_(ndim), tensors()
        {
            int n = group.getNumIrreps();
            const symmetry::Representation* irreps = group.getIrreps();

            sym_.assign(sym, sym+ndim_);

            int ntensors = 1;
            std::vector<int> sublen(ndim_);

            for (int i = 0;i < ndim_;i++)
            {
                ntensors *= n;
                len_[i].assign(len[i], len[i]+n);
                sublen[i] = len_[i][0];
            }

            tensors.resize(ntensors, NULL);

            std::vector<int> idx(ndim_, 0);
            std::vector<symmetry::Representation> prod(ndim_+1, group.totallySymmetricIrrep());
            for (bool done = false;!done;)
            {
                std::cout << idx << prod[0].isTotallySymmetric() << std::endl;

                if (prod[0].isTotallySymmetric())
                {
                    std::vector<int> subsym(sym, sym+n);

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

                    if (ok) tensors[i] = new DistTensor<T>(ctf, ndim, sublen, subsym, zero);
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

        ~SymmetryBlockedTensor()
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL) delete tensors[i];
            }
        }

        virtual SymmetryBlockedTensor<T>& operator=(const double val)
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL) tensors[i] = val;
            }

            return *this;
        }

        static std::vector<int> getStrides(const std::vector<int> indices, const int ndim,
                                           const int len, const int* idx_A)
        {
            std::vector<int> strides(indices.size(), 0);
            std::vector<int> stride_A(ndim);

            stride_A[0] = 1;
            for (int i = 1;i < ndim;i++)
            {
                stride_A[i] = stride_A[i-1]*n;
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

        virtual void mult(const double alpha, const SymmetryBlockedTensor<T>& A, const int* idx_A,
                                              const SymmetryBlockedTensor<T>& B, const int* idx_B,
                          const double beta,                                     const int* idx_C)
        {
            assert(group_ == A.group_);
            assert(group_ == B.group_);

            int n = group_.getNumIrreps();

            double f1 = align_symmetric_indices(A.ndim_, idx_A, A.sym_,
                                               B.ndim_, idx_B, B.sym_,
                                                 ndim_, idx_C,   sym_);
            f1 *= overcounting_factor(A.ndim_, idx_A, A.sym_,
                                      B.ndim_, idx_B, B.sym_,
                                        ndim_, idx_C,   sym_);

            std::vector<int> inds_AB(idx_A, idx_A+A.ndim_);
            idx_AB.insert(idx_AB.end(), idx_B, idx_B+B.ndim_);

            std::vector<int> inds_C(idx_C, idx_C+ndim_);

            uniq(inds_AB);
            uniq(inds_C);
            exclude(inds_AB, inds_C);

            int nAB = inds_AB.size();
            int nC = inds_C.size();

            std::vector<int> stride_A_AB = getStrides(inds_AB, A.ndim_, n, idx_A);
            std::vector<int> stride_B_AB = getStrides(inds_AB, B.ndim_, n, idx_B);
            std::vector<int> stride_A_C = getStrides(inds_C, A.ndim_, n, idx_A);
            std::vector<int> stride_B_C = getStrides(inds_C, B.ndim_, n, idx_B);
            std::vector<int> stride_C = getStrides(inds_C, ndim_, n, idx_C);

            int off_A = 0;
            int off_B = 0;
            int off_C = 0;
            std::vector<int> iC(nC, 0);
            for (bool doneC = false;!doneC;)
            {
                double beta_ = beta;

                std::vector<int> iAB(nAB, 0);
                for (bool doneAB = false;!doneAB;)
                {
                    if (A.tensors[off_A] != NULL &&
                        B.tensors[off_B] != NULL &&
                          tensors[off_C] != NULL)
                    {
                        double f2 = overcounting_factor(A.ndim_, idx_A, A.tensors[off_A]->sym_,
                                                        B.ndim_, idx_B, B.tensors[off_B]->sym_,
                                                          ndim_, idx_C,   tensors[off_C]->sym_);

                        tensors[off_C]->mult(alpha*f1/f2, *A.tensors[off_A], idx_A,
                                                          *B.tensors[off_B], idx_B,
                                             beta_      ,                    idx_C);
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

        virtual void sum(const double alpha, const SymmetryBlockedTensor<T>& A, const int* idx_A,
                         const double beta,                                     const int* idx_B)
        {
            assert(group_ == A.group_);

            int n = group_.getNumIrreps();

            double f = align_symmetric_indices(A.ndim_, idx_A, A.sym_,
                                                 ndim_, idx_B,   sym_);

            std::vector<int> inds_A(idx_A, idx_A+A.ndim_);
            std::vector<int> inds_B(idx_B, idx_B+ndim_);

            uniq(inds_A);
            uniq(inds_B);
            exclude(inds_A, inds_B);

            int nA = inds_A.size();
            int nB = inds_B.size();

            std::vector<int> stride_A_A = getStrides(inds_A, A.ndim_, n, idx_A);
            std::vector<int> stride_A_B = getStrides(inds_B, A.ndim_, n, idx_A);
            std::vector<int> stride_B = getStrides(inds_B, ndim_, n, idx_B);

            int off_A = 0;
            int off_B = 0;
            std::vector<int> iB(nB, 0);
            for (bool doneB = false;!doneB;)
            {
                double beta_ = beta;

                std::vector<int> iA(nA, 0);
                for (bool doneA = false;!doneA;)
                {
                    if (A.tensors[off_A] != NULL &&
                          tensors[off_B] != NULL)
                    {
                        tensors[off_B]->sum(alpha*f, *A.tensors[off_A], idx_A,
                                            beta_  ,                    idx_B);
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

        virtual void scale(const double alpha, const int* idx_A)
        {
            int n = group_.getNumIrreps();

            std::vector<int> inds_A(idx_A, idx_A+ndim_);

            uniq(inds_A);

            int nA = inds_A.size();

            std::vector<int> stride_A = getStrides(inds_A, ndim_, n, idx_A);

            int off_A = 0;
            std::vector<int> iA(nA, 0);
            for (bool doneA = false;!doneA;)
            {
                if (tensors[off_A] != NULL)
                {
                    tensors[off_A]->scale(alpha, idx_A);
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
};

}
}

#endif
