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

#ifndef _AQUARIUS_SCF_CHOLESKYSCF_HPP_
#define _AQUARIUS_SCF_CHOLESKYSCF_HPP_

#include "operator/cholesky.hpp"

#include "scf.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class CholeskyUHF : public UHF<T>
{
    protected:
        const op::CholeskyIntegrals<T>& chol;
        tensor::DistTensor<T> *J;
        tensor::DistTensor<T> *JD;
        tensor::DistTensor<T> *La_occ, *Lb_occ;
        tensor::DistTensor<T> *LDa_occ, *LDb_occ;

    public:
        CholeskyUHF(const input::Config& config, const op::CholeskyIntegrals<T>& chol);

        ~CholeskyUHF();

        const op::CholeskyIntegrals<T>& getCholesky() const { return chol; }

    protected:
        void buildFock();
};

}
}

#endif
