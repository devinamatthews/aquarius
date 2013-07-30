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

#ifndef _AQUARIUS_SCF_AOMOINTS_HPP_
#define _AQUARIUS_SCF_AOMOINTS_HPP_

#include "aoints.hpp"
#include "aoscf.hpp"
#include "moints.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class AOMOIntegrals : public MOIntegrals<T>
{
    public:
        AOMOIntegrals(const AOUHF<T>& uhf);

    private:
        using Distributed<T>::nproc;
        using Distributed<T>::rank;
        enum Side {PQ, RS};
        enum Index {A, B};

        static bool sortIntsByRS(const idx4_t& i1, const idx4_t& i2);

        struct abrs_integrals;

        struct pqrs_integrals : Distributed<T>
        {
            int np, nq, nr, ns;
            size_t nints;
            T *ints;
            idx4_t *idxs;

            using Distributed<T>::nproc;
            using Distributed<T>::rank;

            /*
             * Read integrals in and break (pq|rs)=(rs|pq) symmetry
             */
            pqrs_integrals(const AOIntegrals<T>& aoints);

            pqrs_integrals(abrs_integrals& abrs);

            void free();

            /*
             * Redistribute integrals such that each node has all pq for each rs pair
             */
            void collect(const bool rles);
        };

        struct abrs_integrals : Distributed<T>
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
        void doTransformation(const AOIntegrals<T>& ints);
};

}
}

#endif
