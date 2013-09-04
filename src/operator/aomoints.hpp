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

#ifndef _AQUARIUS_OPERATOR_AOMOINTS_HPP_
#define _AQUARIUS_OPERATOR_AOMOINTS_HPP_

#include "scf/aoscf.hpp"
#include "integrals/2eints.hpp"

#include "moints.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class AOMOIntegrals : public MOIntegrals<T>
{
    public:
        AOMOIntegrals(const std::string& name, const input::Config& config);

    private:
        enum Side {PQ, RS};
        enum Index {A, B};

        struct abrs_integrals;

        struct pqrs_integrals : Distributed
        {
            int np, nq, nr, ns;
            size_t nints;
            T *ints;
            idx4_t *idxs;

            /*
             * Read integrals in and break (pq|rs)=(rs|pq) symmetry
             */
            pqrs_integrals(int norb, const integrals::ERI& aoints);

            pqrs_integrals(abrs_integrals& abrs);

            void free();

            void sortInts(bool rles, size_t& nrs, size_t*& rscount);

            /*
             * Redistribute integrals such that each node has all pq for each rs pair
             */
            void collect(bool rles);
        };

        struct abrs_integrals : Distributed
        {
            int na, nb, nr, ns;
            size_t nints;
            T *ints;
            size_t nrs;
            idx2_t *rs;

            /*
             * Expand (p_i q_j|r_k s_l) into (pq|r_k s_l), where pleq = true indicates
             * that (pq|rs) = (qp|rs) and only p_i <= q_j is stored in the input
             */
            abrs_integrals(pqrs_integrals& pqrs, const bool pleq);

            /*
             * Transform (ab|rs) -> (cb|rs) (index = A) or (ab|rs) -> (ac|rs) (index = B)
             *
             * C is ldc*nc if trans = 'N' and ldc*[na|nb] if trans = 'T'
             */
            abrs_integrals transform(Index index, const char trans, const int nc, const double* C, const int ldc);

            void transcribe(tensor::DistTensor<T>& tensor, bool assymij, bool assymkl, bool reverse);

            void free();
        };

    protected:
        void run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif
