/* Copyin (c) 2013, Devin Matthews
 * All ins reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyin
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyin
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
#include "operator/space.hpp"
#include "util/stl_ext.hpp"
#include "task/task.hpp"

#include "dist_tensor.hpp"
#include "composite_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template<class T>
class SpinorbitalTensor : public IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T>, public task::Resource
{
    INHERIT_FROM_INDEXABLE_COMPOSITE_TENSOR(SpinorbitalTensor<T>,DistTensor<T>,T)

    public:
        SpinorbitalTensor(const SpinorbitalTensor<T>& t, const T val);

        SpinorbitalTensor(const SpinorbitalTensor<T>& other);

        SpinorbitalTensor(const Arena& arena,
                          const std::vector<op::Space>& spaces,
                          const std::vector<int>& nout,
                          const std::vector<int>& nin, int spin=0);

        DistTensor<T>& operator()(const std::vector<int>& alpha_out,
                                  const std::vector<int>& alpha_in);

        const DistTensor<T>& operator()(const std::vector<int>& alpha_out,
                                        const std::vector<int>& alpha_in) const;

        void mult(const T alpha, bool conja, const SpinorbitalTensor<T>& A_, const std::string& idx_A,
                                 bool conjb, const SpinorbitalTensor<T>& B_, const std::string& idx_B,
                  const T beta_,                                             const std::string& idx_C);

        void sum(const T alpha, bool conja, const SpinorbitalTensor<T>& A_, const std::string& idx_A,
                 const T beta_,                                             const std::string& idx_B);

        void scale(const T alpha, const std::string& idx_A);

        void weight(const std::vector<const std::vector<T>*>& da,
                    const std::vector<const std::vector<T>*>& db);

        T dot(bool conja, const SpinorbitalTensor<T>& A, const std::string& idx_A,
              bool conjb,                                const std::string& idx_B) const;

    protected:
        struct SpinCase
        {
            DistTensor<T> *tensor;
            std::vector<int> alpha_out, alpha_in;

            void construct(SpinorbitalTensor<T>& t,
                           const std::vector<int>& alpha_out,
                           const std::vector<int>& alpha_in);
        };

        std::vector<op::Space> spaces;
        std::vector<int> nout, nin;
        int spin;
        std::vector<SpinCase> cases;

        /*
        static void matchTypes(const int nin_A, const int nout_A, const std::vector<autocc::Line>& log_A, const std::string& idx_A,
                               const int nin_B, const int nout_B, const std::vector<autocc::Line>& log_B, const std::string& idx_B,
                               std::vector<autocc::Line>& out_A, std::vector<autocc::Line>& in_A,
                               std::vector<autocc::Line>& pout_A, std::vector<autocc::Line>& hout_A,
                               std::vector<autocc::Line>& pin_A, std::vector<autocc::Line>& hin_A,
                               std::vector<autocc::Line>& sum_A);
         */
};

}
}

#endif
