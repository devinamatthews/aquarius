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

int tensor_print_dense(FILE* fp, const double* A, const int ndim_A, const int* len_A, const int* lda)
{
    size_t off, size;
    int i;
    size_t stride[ndim_A];
    int pos[ndim_A];
    bool done;

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

    off = 0;

    /*
     * loop over elements in A
     */
    memset(pos, 0, ndim_A*sizeof(int));
    for (done = false;!done;)
    {
#ifdef CHECK_BOUNDS
        if (off < 0 || off >= size) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

        for (i = 0;i < ndim_A;i++) fprintf(fp, "%d ", pos[i]);
        fprintf(fp, "%.15e\n", A[off]);

        for (i = 0;i < ndim_A;i++)
        {
            if (pos[i] == len_A[i] - 1)
            {
                pos[i] = 0;
                off -= stride[i]*(len_A[i]-1);

                if (i == ndim_A - 1)
                {
                    done = true;
                    break;
                }
            }
            else
            {
                off += stride[i];
                pos[i]++;
                break;
            }
        }

        if (ndim_A == 0) done = true;
    }
    /*
     * end loop over A
     */

    return TENSOR_SUCCESS;
}
