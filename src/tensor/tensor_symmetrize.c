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

#include "tensor.h"
#include "util.h"

int tensor_symmetrize(const double* restrict A, double* restrict B, const int ndim, const int* restrict len, const int* restrict sym)
{
    int i, j, sgn;
    int idx[ndim];
    size_t stride[ndim];
    int fac;
    size_t off_A, off_B;
    size_t size_A, size_B;

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(ndim, len, NULL, sym);
#endif //VALIDATE_INPUTS

    stride[0] = 1;
    for (i = 0;i < ndim - 1;i++)
        stride[i + 1] = stride[i] * len[i];

    size_A = tensor_size_dense(ndim, len, NULL);
    size_B = tensor_size(ndim, len, NULL, sym);

    /*
     * initialize B with zeroes
     */
    dzero(size_B, B, 1);

    /*
     * initialize index string to lexicographically minimum values
     */
    first_packed_indices(ndim, len, sym, idx);

    /*
     * loop over index combinations which are valid for packed storage (i.e. i <(<=) j <(<=) k ...)
     * these correspond to unique, contiguous elements of A
     */
    for (off_B = 0;;off_B++)
    {
#ifdef CHECK_BOUNDS
        if (off_B < 0 || off_B > size_B) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

        /*
         * the algorithm in next_permutation generates only distinct permutations for repeated indices
         * if there are symmetry-related indices which are equal, we must include an additional factor
         * for identical permutations
         */
        fac = 1;
        j = 1;
        for (i = 1;i < ndim;i++)
        {
            if (sym[i - 1] == SY && idx[i - 1] == idx[i])
            {
                j++;
            }
            else
            {
                if (j > 1) fac *= factorial(j);
                j = 1;
            }
        }
        fac *= factorial(j);

        /*
         * for each of these combinations, loop over all permutations in each symmetry-related index group
         * to inizialize all redundant elements of B
         */
        sgn = 1;
        while (sgn != 0)
        {
            off_A = 0;
            for (i = 0;i < ndim;i++) off_A += stride[i]*idx[i];

#ifdef CHECK_BOUNDS
            if (off_A < 0 || off_A > size_A) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

            B[off_B] += fac * sgn * A[off_A];

            /*
             * loop over symmetry-related index groups
             * increment each group to next permutation
             */
            sgn *= next_permutation(ndim, sym, idx);
        }

        /*
         * increment the index string to the next in lexicographic order
         */
        if (!next_packed_indices(ndim, len, sym, idx)) break;
    }

    return TENSOR_SUCCESS;
}
