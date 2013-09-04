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

#include "ctf.hpp"
#include "../src/dist_tensor/sym_indices.hxx"

#include "symmetry/symmetry.hpp"
#include "util/distributed.hpp"

#include "composite_tensor.hpp"
#include "dist_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <class T>
class SymmetryBlockedTensor : public IndexableCompositeTensor<SymmetryBlockedTensor<T>,DistTensor<T>,T>, public Distributed
{
    INHERIT_FROM_INDEXABLE_COMPOSITE_TENSOR(SymmetryBlockedTensor<T>,DistTensor<T>,T)

    protected:
        const symmetry::PointGroup& group;
        std::vector< std::vector<int> > len;
        std::vector<int> sym;

    public:
        SymmetryBlockedTensor(const Arena& arena, const symmetry::PointGroup& group, T scalar = (T)0);

        SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& other);

        SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& other, T scalar);

        SymmetryBlockedTensor(const Arena& arena, const symmetry::PointGroup& group,
                              int ndim, const std::vector<std::vector<int> >& len,
                              const std::vector<int>& sym, bool zero=true);

        virtual ~SymmetryBlockedTensor() {}

        static std::vector<int> getStrides(const std::string& indices, const int ndim,
                                           const int len, const std::string& idx_A);

        virtual void mult(const double alpha, bool conja, const SymmetryBlockedTensor<T>& A, const std::string& idx_A,
                                              bool conjb, const SymmetryBlockedTensor<T>& B, const std::string& idx_B,
                          const double beta,                                                 const std::string& idx_C);

        virtual void sum(const double alpha, bool conja, const SymmetryBlockedTensor<T>& A, const std::string& idx_A,
                         const double beta,                                                 const std::string& idx_B);

        virtual void scale(const double alpha, const std::string& idx_A);

        virtual T dot(bool conja, const SymmetryBlockedTensor<T>& A, const std::string& idx_A,
                      bool conjb,                                    const std::string& idx_B) const;
};

}
}

#endif
