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

#ifndef _AQUARIUS_TENSOR_CTF_TENSOR_HPP_
#define _AQUARIUS_TENSOR_CTF_TENSOR_HPP_

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
#include "util/stl_ext.hpp"
#include "task/task.hpp"

#include "util.h"
#include "indexable_tensor.hpp"

namespace std
{
    template <typename T> struct is_pod<tkv_pair<T> > { static const bool value = true; };
}

namespace aquarius
{
namespace tensor
{

template <typename T>
class CTFTensor : public IndexableTensor< CTFTensor<T>,T >, public task::Resource
{
    INHERIT_FROM_INDEXABLE_TENSOR(CTFTensor<T>,T)

    protected:
        tCTF_Tensor<T>* dt;
        std::vector<int> len;
        std::vector<int> sym;
        static std::map<const tCTF_World<T>*,std::pair<int,CTFTensor<T>*> > scalars;

        void allocate();

        void free();

        void register_scalar();

        void unregister_scalar();

        CTFTensor<T>& scalar() const;

    public:
        CTFTensor(const Arena& arena, T scalar = (T)0);

        CTFTensor(const CTFTensor<T>& A, T scalar);

        CTFTensor(const CTFTensor<T>& A, bool copy=true, bool zero=false);

        CTFTensor(CTFTensor<T>* A);

        CTFTensor(const CTFTensor<T>& A, const std::vector<int>& start_A, const std::vector<int>& len_A);

        CTFTensor(const Arena& arena, int ndim, const std::vector<int>& len, const std::vector<int>& sym,
                   bool zero=true);

        ~CTFTensor();

        void resize(int ndim, const std::vector<int>& len, const std::vector<int>& sym, bool zero);

        const std::vector<int>& getLengths() const { return len; }

        const std::vector<int>& getSymmetry() const { return sym; }

        T* getRawData(int64_t& size);

        const T* getRawData(int64_t& size) const;

        template <typename Container>
        void getLocalData(Container& pairs) const
        {
            int64_t npair;
            tkv_pair<T> *data;
            dt->read_local(&npair, &data);
            pairs.assign(data, data+npair);
            if (npair > 0) ::free(data);
        }

        template <typename Container>
        void getRemoteData(Container& pairs) const
        {
            dt->read(pairs.size(), pairs.data());
        }

        void getRemoteData() const
        {
            dt->read(0, NULL);
        }

        template <typename Container>
        void writeRemoteData(const Container& pairs)
        {
            dt->write(pairs.size(), pairs.data());
        }

        void writeRemoteData()
        {
            dt->write(0, NULL);
        }

        template <typename Container>
        void writeRemoteData(double alpha, double beta, const Container& pairs)
        {
            dt->write(pairs.size(), alpha, beta, pairs.data());
        }

        void writeRemoteData(double alpha, double beta)
        {
            dt->write(0, alpha, beta, NULL);
        }

        template <typename Container>
        void getAllData(Container& vals) const
        {

            getAllData(vals, 0);
            int64_t npair = vals.size();
            this->arena.Bcast(&npair, 1, 0);
            if (this->rank != 0) vals.resize(npair);
            this->arena.Bcast(vals, 0);
        }

        template <typename Container>
        void getAllData(Container& vals, int rank) const
        {
            assert(this->rank == rank);

            for (int i = 0;i < ndim;i++)
            {
                if (len[i] == 0)
                {
                    vals.clear();
                    return;
                }
            }

            std::vector<tkv_pair<T> > pairs;
            std::vector<int> idx(ndim, 0);

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

            dt->read(pairs.size(), pairs.data());

            sort(pairs.begin(), pairs.end());
            size_t npair = pairs.size();
            vals.resize(npair);

            for (size_t i = 0;i < npair;i++)
            {
                vals[i] = pairs[i].d;
            }
        }

        void getAllData(int rank) const
        {
            assert(this->arena.rank != rank);
            dt->read(0, NULL);
        }

        void slice(T alpha, bool conja, const CTFTensor<T>& A,
                   const std::vector<int>& start_A, T beta);

        void slice(T alpha, bool conja, const CTFTensor<T>& A,
                   T beta, const std::vector<int>& start_B);

        void slice(T alpha, bool conja, const CTFTensor<T>& A, const std::vector<int>& start_A,
                   T  beta,                                     const std::vector<int>& start_B,
                                                                const std::vector<int>& len);

        void div(T alpha, bool conja, const CTFTensor<T>& A,
                          bool conjb, const CTFTensor<T>& B, T beta);

        void invert(T alpha, bool conja, const CTFTensor<T>& A, T beta);

        void weight(const std::vector<const std::vector<T>*>& d);

        void print(FILE* fp, double cutoff = -1.0) const;

        void compare(FILE* fp, const CTFTensor<T>& other, double cutoff = 0.0) const;

        typename std::real_type<T>::type norm(int p) const;

        void mult(T alpha, bool conja, const CTFTensor<T>& A, const std::string& idx_A,
                           bool conjb, const CTFTensor<T>& B, const std::string& idx_B,
                  T  beta,                                     const std::string& idx_C);

        void sum(T alpha, bool conja, const CTFTensor<T>& A, const std::string& idx_A,
                 T  beta,                                     const std::string& idx_B);

        void scale(T alpha, const std::string& idx_A);

        T dot(bool conja, const CTFTensor<T>& A, const std::string& idx_A,
              bool conjb,                         const std::string& idx_B) const;
};

}
}

#endif
