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

int tensor_slice_dense(const double*  A, const int  ndim_A, const int* len_A, const int* lda,
                             double** B,       int* ndim_B,       int* len_B,       int* ldb,
                       const int* start, const int* len)
{
    int i;
    size_t ld, stride;

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(ndim_A, len_A, lda, NULL);
#endif //VALIDATE_INPUTS

    *B = (double*)A;
    *ndim_B = 0;

    if (ndim_A > 0 && lda != NULL)
    {
        stride = lda[0];
        ld = lda[0];
    }
    else
    {
        stride = 1;
        ld = 1;
    }

    for (i = 0;i < ndim_A;i++)
    {
#ifdef VALIDATE_INPUTS
        if (start[i] < 0 || start[i] >= len_A[i]) return TENSOR_INVALID_START;
        if (len[i] < 0 || len[i] > len_A[i]) return TENSOR_INVALID_LENGTH;
        if (start[i]+MAX(1,len[i]) > len_A[i]) return TENSOR_LENGTH_MISMATCH;
#endif //VALIDATE_INPUTS

        *B += stride*start[i];

        if (len[i] == 0)
        {
            if (lda != NULL && i < ndim_A-1)
            {
                ld *= lda[i+1];
            }
            else
            {
                ld *= len_A[i];
            }
        }
        else
        {
            ldb[*ndim_B] = ld;
            len_B[*ndim_B] = len[i];
            (*ndim_B)++;

            if (lda != NULL && i < ndim_A-1)
            {
                ld = lda[i+1];
            }
            else
            {
                ld = len_A[i];
            }
        }

        if (lda == NULL || i == ndim_A-1)
        {
            stride *= len_A[i];
        }
        else
        {
            stride *= lda[i+1];
        }
    }

    return TENSOR_SUCCESS;
}
