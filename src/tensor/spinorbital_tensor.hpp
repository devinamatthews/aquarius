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

#ifndef _AQUARIUS_TENSOR_SPINORBITAL_TENSOR_HPP_
#define _AQUARIUS_TENSOR_SPINORBITAL_TENSOR_HPP_

#include <iostream>
#include <vector>
#include <cstring>
#include <cassert>
#include <string>
#include <algorithm>

#include "autocc/autocc.hpp"
#include "memory/memory.h"

#include "dist_tensor.hpp"
#include "composite_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template<class T>
class SpinorbitalTensor : public IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T>
{
    INHERIT_FROM_INDEXABLE_COMPOSITE_TENSOR(SpinorbitalTensor<T>,DistTensor<T>,T)

    protected:
        struct SpinCase
        {
            DistTensor<T>* tensor;
            std::vector<autocc::Line> logical;
            int nA, nE, nM, nI;
            std::vector<int> log_to_phys;
            double permFactor;

            SpinCase(DistTensor<T>& tensor) : tensor(&tensor) {}
        };

        std::vector<autocc::Line> logical;
        int nA, nE, nM, nI;
        int spin;
        //double permFactor;
        std::vector<SpinCase> cases;

    public:
        void set_name(char const * name_);

        SpinorbitalTensor(const SpinorbitalTensor<T>& t, const T val);

        SpinorbitalTensor(const SpinorbitalTensor<T>& other);

        SpinorbitalTensor(const autocc::Manifold& left, const autocc::Manifold& right, const int spin=0);

        SpinorbitalTensor(const std::string& logical, const int spin=0);

        virtual ~SpinorbitalTensor() {}

        void addSpinCase(DistTensor<T>* tensor, std::string logical, std::string physical, double factor = 1.0, bool isAlloced = true);

        void addSpinCase(DistTensor<T>& tensor, std::string logical, std::string physical, double factor = 1.0, bool isAlloced = false);

        void addSpinCase(DistTensor<T>* tensor, const autocc::Manifold& alpha_left,
                         const autocc::Manifold& alpha_right, double factor = 1.0, bool isAlloced = true);

        void addSpinCase(DistTensor<T>& tensor, const autocc::Manifold& alpha_left,
                         const autocc::Manifold& alpha_right, double factor = 1.0, bool isAlloced = false);

        DistTensor<T>& operator()(int nA, int nM, int nE, int nI);

        const DistTensor<T>& operator()(int nA, int nM, int nE, int nI) const;

        void mult(const T alpha, bool conja, const SpinorbitalTensor<T>& A_, const int* idx_A,
                                 bool conjb, const SpinorbitalTensor<T>& B_, const int* idx_B,
                  const T beta_,                                             const int* idx_C);

        void sum(const T alpha, bool conja, const SpinorbitalTensor<T>& A_, const int* idx_A,
                 const T beta_,                                             const int* idx_B);

        void scale(const T alpha, const int* idx_A);

        void div(const T alpha, bool conja, const SpinorbitalTensor<T>& A,
                                bool conjb, const SpinorbitalTensor<T>& B, const T beta);

        void invert(const T alpha, bool conja, const SpinorbitalTensor<T>& A, const T beta);

        T dot(bool conja, const SpinorbitalTensor<T>& A, const int* idx_A,
              bool conjb,                                const int* idx_B) const;

    protected:
        void addSpinCase(DistTensor<T>& tensor, std::vector<autocc::Line> logical, std::vector<autocc::Line> physical, double factor, bool isAlloced);

        static void matchTypes(const int nin_A, const int nout_A, const std::vector<autocc::Line>& log_A, const int* idx_A,
                               const int nin_B, const int nout_B, const std::vector<autocc::Line>& log_B, const int* idx_B,
                               std::vector<autocc::Line>& out_A, std::vector<autocc::Line>& in_A,
                               std::vector<autocc::Line>& pout_A, std::vector<autocc::Line>& hout_A,
                               std::vector<autocc::Line>& pin_A, std::vector<autocc::Line>& hin_A,
                               std::vector<autocc::Line>& sum_A);
};

}
}

#endif
