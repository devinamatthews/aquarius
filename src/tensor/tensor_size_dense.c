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

int64_t tensor_size_dense(const int ndim, const int* len, const int* ld)
{
    int i;
    int64_t size;

#ifdef VALIDATE_INPUTS
    if (ndim < 0) return TENSOR_INVALID_NDIM;
#endif //VALIDATE_INPUTS

    if (ld == NULL || ndim == 0)
    {
        size = 1;
    }
    else
    {
#ifdef VALIDATE_INPUTS
        if (ld[0] <= 0) return TENSOR_INVALID_LD;
#endif //VALIDATE_INPUTS
        size = ld[0];
    }

    for (i = 0;i < ndim;i++)
    {
#ifdef VALIDATE_INPUTS
        if (len[i] <= 0) return TENSOR_INVALID_LENGTH;
#endif //VALIDATE_INPUTS

        if (ld == NULL || i == ndim-1)
        {
            size *= len[i];
        }
        else
        {
#ifdef VALIDATE_INPUTS
            if (ld[i+1] <= 0) return TENSOR_INVALID_LD;
#endif //VALIDATE_INPUTS
            size *= ld[i+1];
        }
    }

    return size;
}
