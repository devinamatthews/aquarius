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

int64_t tensor_size(const int ndim, const int* len, const int* ld, const int* sym)
{
    int i, k, mp;
    int64_t size, tmp;

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(ndim, len, ld, sym);
#endif //VALIDATE_INPUTS

    k = 1;
    tmp = 1;
    size = (ld == NULL || ndim == 0 ? 1 : ld[0]);
    if (ndim > 0) mp = len[0];
    for (i = 0;i < ndim;i++)
    {
        tmp = (tmp*mp)/k;
        k++;
        switch (sym[i])
        {
            case AS:
            case SH:
                mp--;
                break;
            case SY:
            case NS:
                mp++;
                break;
        }

        if (sym[i] == NS)
        {
            if (ld == NULL || i == ndim-1)
            {
                size *= ld[i+1];
            }
            else
            {
                size *= tmp;
            }
            k = 1;
            tmp = 1;
            if (i < ndim-1) mp = len[i+1];
        }
    }

    return size;
}
