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

/*
 * transform an unpacked (redundant) tensor into packed form, only the ordered (i <(<=) j <(<=) k ...) elements of A are accessed
 *
 * A(input)             - tensor to pack
 * B(output)            - packed tensor
 * ndim(input)          - number of dimensions in tensor
 * len(input)           - length of each dimension is tensor, must be the same for symmetry-related indices
 * sym(input)           - sym[i] = 1 iff ind[i] <= ind[i+1] in A,
 *                        sym[i] = -1 iff ind[i] < ind[i+1] in A,
 *                        sym[i] = 0 otherwise
 *
 * performs NO input checking
 */
int tensor_pack(const double* restrict A, double* restrict B, const int ndim, const int* restrict len, const int* restrict sym)
{
    int i, j;
    int idx[ndim];
    size_t stride[ndim];
    bool done;
    int fac;
    size_t off_A, off_B;
    size_t size_A, size_B;

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(ndim, len, NULL, sym);
#endif //VALIDATE_INPUTS

    fac = 1;

//    /*
//     * multiply by a factor of m! for each group of m symmetry-related indices
//     *
//     * (not necessary unless we really need a normalized transformation)
//     */
//    j = 1;
//    for (i = 1;i < ndim;i++)
//    {
//        if (sym[i-1] != 0)
//        {
//            j++;
//        }
//        else
//        {
//            if (j > 1) fac *= fact(j);
//            j = 1;
//        }
//    }
//    fac *= fact(j);

    size_A = tensor_size_dense(ndim, len, NULL);
    size_B = tensor_size(ndim, len, NULL, sym);

    stride[0] = 1;
    for (i = 0;i < ndim - 1;i++)
        stride[i + 1] = stride[i] * len[i];

    first_packed_indices(ndim, len, sym, idx);

    done = false;
    for (off_B = 0;!done;off_B++)
    {
        off_A = 0;
        for (i = 0;i < ndim;i++) off_A += idx[i]*stride[i];

#ifdef CHECK_BOUNDS
        if (off_A < 0 || off_A >= size_A) return TENSOR_OUT_OF_BOUNDS;
        if (off_B < 0 || off_B >= size_B) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

        B[off_B] = fac * A[off_A];

        done = !next_packed_indices(ndim, len, sym, idx);
    }

    return TENSOR_SUCCESS;
}
