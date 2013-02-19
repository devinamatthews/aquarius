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
#include <string.h>

int tensor_scale_dense_(const double alpha, double* restrict A, const int ndim_A, const int* restrict len_A, const int* restrict lda, const int* restrict idx_A)
{
    int i, j;
    bool found, done;
    int ndim_uniq;
    int len_uniq[ndim_A];
    size_t stride[ndim_A];
    size_t off;
    size_t inc[ndim_A];
    size_t size;
    int pos[ndim_A];

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(ndim_A, len_A, lda, NULL);
#endif //VALIDATE_INPUTS

    if (lda == NULL)
    {
        if (ndim_A > 0) stride[0] = 1;
        for (i = 1;i < ndim_A;i++) stride[i] = stride[i-1]*len_A[i-1];
    }
    else
    {
        if (ndim_A > 0) stride[0] = lda[0];
        for (i = 1;i < ndim_A;i++) stride[i] = stride[i-1]*lda[i];
    }

    if (ndim_A > 0)
    {
        size = stride[ndim_A-1]*len_A[ndim_A-1];
    }
    else
    {
        size = 1;
    }

    ndim_uniq = 0;
    memset(inc, 0, ndim_A*sizeof(size_t));

    for (i = 0;i < ndim_A;i++)
    {
        found = false;

        for (j = i+1;j < ndim_A;j++)
        {
            if (idx_A[i] == idx_A[j])
            {
#ifdef VALIDATE_INPUTS
                if (len_A[i] != len_A[j]) return TENSOR_LENGTH_MISMATCH;
#endif //VALIDATE_INPUTS
                found = true;
#ifndef VALIDATE_INPUTS
                break;
#endif //VALIDATE_INPUTS
            }
        }

        if (!found)
        {
            len_uniq[ndim_uniq] = len_A[i];

            for (j = 0;j < ndim_A;j++)
            {
                if (idx_A[j] == idx_A[i]) inc[ndim_uniq] += stride[j];
            }

            ndim_uniq++;
        }
    }

    off = 0;

    /*
     * loop over elements in A
     */
    memset(pos, 0, ndim_uniq*sizeof(int));
    for (done = false;!done;)
    {
#ifdef CHECK_BOUNDS
        if (off < 0 || off >= size) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

        if (alpha == 0.0)
        {
            A[off] = 0.0;
        }
        else
        {
            A[off] *= alpha;
        }

        for (i = 0;i < ndim_uniq;i++)
        {
            if (pos[i] == len_uniq[i] - 1)
            {
                pos[i] = 0;
                off -= inc[i]*(len_uniq[i]-1);

                if (i == ndim_uniq - 1)
                {
                    done = true;
                    break;
                }
            }
            else
            {
                pos[i]++;
                off += inc[i];
                break;
            }
        }

        if (ndim_uniq == 0) done = true;
    }
    /*
     * end loop over A
     */

    return TENSOR_SUCCESS;
}
