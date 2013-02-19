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

int tensor_densify(int* ndim, int* len, const int* sym)
{
    int i, k, mp;
    int nd;
    size_t tmp;
    int len_dense[*ndim];

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(*ndim, len, NULL, sym);
#endif //VALIDATE_INPUTS

    k = 1;
    tmp = 1;
    mp = len[0];
    nd = 0;
    for (i = 0;i < *ndim;i++)
    {
        tmp = (tmp * mp) / k;
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
            len_dense[nd] = tmp;
            nd++;
            k = 1;
            tmp = 1;
            if (i < *ndim - 1) mp = len[i + 1];
        }
    }

    *ndim = nd;
    memcpy(len, len_dense, nd*sizeof(int));

    return TENSOR_SUCCESS;
}
