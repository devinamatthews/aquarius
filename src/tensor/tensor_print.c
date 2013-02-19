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

int tensor_print(FILE* fp, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A)
{
    size_t k, size_A;
    int i;
    int idx[ndim_A];
    bool done;

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(ndim_A, len_A, NULL, sym_A);
#endif //VALIDATE_INPUTS

    size_A = tensor_size(ndim_A, len_A, NULL, sym_A);

    first_packed_indices(ndim_A, len_A, sym_A, idx);

    done = false;
    for (k = 0;!done;k++)
    {
#ifdef CHECK_BOUNDS
        if (k < 0 || k >= size_A) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

        for (i = 0;i < ndim_A;i++) fprintf(fp, "%d ", idx[i]);
        fprintf(fp, "%.15e\n", A[k]);

        done = !next_packed_indices(ndim_A, len_A, sym_A, idx);
    }

    return TENSOR_SUCCESS;
}
