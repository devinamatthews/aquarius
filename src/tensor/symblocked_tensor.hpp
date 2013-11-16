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
#include <cassert>
#include <string>
#include <algorithm>
#include <map>

#include "ctf.hpp"
#include "../src/dist_tensor/sym_indices.hxx"

#include "symmetry/symmetry.hpp"
#include "task/task.hpp"

#include "composite_tensor.hpp"
#include "ctf_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <class T>
class SymmetryBlockedTensor : public IndexableCompositeTensor<SymmetryBlockedTensor<T>,CTFTensor<T>,T>,
                              public task::Resource
{
    INHERIT_FROM_INDEXABLE_COMPOSITE_TENSOR(SymmetryBlockedTensor<T>,CTFTensor<T>,T)

    protected:
        const symmetry::PointGroup& group;
        std::vector< std::vector<int> > len;
        std::vector<int> sym;
        std::vector<double> factor;
        std::vector<std::vector<int> > reorder;
        static std::map<const tCTF_World<T>*,std::map<const symmetry::PointGroup*,std::pair<int,SymmetryBlockedTensor<T>*> > > scalars;

        static std::vector<int> getStrides(const std::string& indices, const int ndim,
                                           const int len, const std::string& idx_A);

        void allocate(bool zero);

        CTFTensor<T>& operator()(const std::vector<int>& irreps);

        const CTFTensor<T>& operator()(const std::vector<int>& irreps) const;

        void register_scalar();

        void unregister_scalar();

        SymmetryBlockedTensor<T>& scalar() const;

    public:
        SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& other);

        SymmetryBlockedTensor(SymmetryBlockedTensor<T>* other);

        SymmetryBlockedTensor(const std::string& name, const SymmetryBlockedTensor<T>& other);

        SymmetryBlockedTensor(const std::string& name, SymmetryBlockedTensor<T>* other);

        SymmetryBlockedTensor(const std::string& name, const SymmetryBlockedTensor<T>& other, T scalar);

        SymmetryBlockedTensor(const std::string& name, const SymmetryBlockedTensor<T>& A,
                              const std::vector<std::vector<int> >& start_A,
                              const std::vector<std::vector<int> >& len_A);

        SymmetryBlockedTensor(const std::string& name, const Arena& arena, const symmetry::PointGroup& group,
                              int ndim, const std::vector<std::vector<int> >& len,
                              const std::vector<int>& sym, bool zero=true);

        ~SymmetryBlockedTensor();

        const symmetry::PointGroup& getGroup() const { return group; }

        const std::vector<std::vector<int> >& getLengths() const { return len; }

        const std::vector<int>& getSymmetry() const { return sym; }

        bool exists(const std::vector<int>& irreps) const;

        T* getRawData(const std::vector<int>& irreps, int64_t& size)
        {
            return (*this)(irreps).getRawData(size);
        }

        const T* getRawData(const std::vector<int>& irreps, int64_t& size) const
        {
            return (*this)(irreps).getRawData(size);
        }

        template <typename Container>
        void getLocalData(const std::vector<int>& irreps, Container& pairs) const
        {
            (*this)(irreps).getLocalData(pairs);
        }

        template <typename Container>
        void getRemoteData(const std::vector<int>& irreps, Container& pairs) const
        {
            (*this)(irreps).getRemoteData(pairs);
        }

        void getRemoteData(const std::vector<int>& irreps) const
        {
            (*this)(irreps).getRemoteData();
        }

        template <typename Container>
        void writeRemoteData(const std::vector<int>& irreps, const Container& pairs)
        {
            (*this)(irreps).writeRemoteData(pairs);
        }

        void writeRemoteData(const std::vector<int>& irreps)
        {
            (*this)(irreps).writeRemoteData();
        }

        template <typename Container>
        void writeRemoteData(const std::vector<int>& irreps, double alpha, double beta, const Container& pairs)
        {
            (*this)(irreps).writeRemoteData(alpha, beta, pairs);
        }

        void writeRemoteData(const std::vector<int>& irreps, double alpha, double beta)
        {
            (*this)(irreps).writeRemoteData(alpha, beta);
        }

        template <typename Container>
        void getAllData(const std::vector<int>& irreps, Container& vals) const
        {
            (*this)(irreps).getAllData(vals);
        }

        template <typename Container>
        void getAllData(const std::vector<int>& irreps, Container& vals, int rank) const
        {
            (*this)(irreps).getAllData(vals, rank);
        }

        void getAllData(const std::vector<int>& irreps, int rank) const
        {
            (*this)(irreps).getAllData(rank);
        }

        void slice(T alpha, bool conja, const SymmetryBlockedTensor<T>& A,
                   const std::vector<std::vector<int> >& start_A, T beta);

        void slice(T alpha, bool conja, const SymmetryBlockedTensor<T>& A,
                   T beta, const std::vector<std::vector<int> >& start_B);

        void slice(T alpha, bool conja, const SymmetryBlockedTensor<T>& A,
                                        const std::vector<std::vector<int> >& start_A,
                   T  beta,             const std::vector<std::vector<int> >& start_B,
                                        const std::vector<std::vector<int> >& len);

        virtual void mult(T alpha, bool conja, const SymmetryBlockedTensor<T>& A, const std::string& idx_A,
                                   bool conjb, const SymmetryBlockedTensor<T>& B, const std::string& idx_B,
                          T beta,                                                 const std::string& idx_C);

        virtual void sum(T alpha, bool conja, const SymmetryBlockedTensor<T>& A, const std::string& idx_A,
                         T beta,                                                 const std::string& idx_B);

        virtual void scale(T alpha, const std::string& idx_A);

        virtual T dot(bool conja, const SymmetryBlockedTensor<T>& A, const std::string& idx_A,
                      bool conjb,                                    const std::string& idx_B) const;

        void weight(const std::vector<const std::vector<std::vector<T> >*>& d);

        typename std::real_type<T>::type norm(int p) const;
};

}
}

#endif
