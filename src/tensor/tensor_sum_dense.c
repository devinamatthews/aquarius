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

/**
 * Sum a tensor (presumably operated on in one or more ways) onto a second.
 *
 * This form generalizes all of the unary operations trace, transpose, diagonal, and replicate, which may be performed
 * in any combination. Even in the case that only a subset of the elements of B are written to by the operation, all elements
 * of B are first scaled by beta. Replication is performed in-place.
 *
 * \author Devin Matthews
 * \date Sep. 28 2011
 */

#include "tensor.h"
#include "util.h"
#include <string.h>

int tensor_sum_dense_(const double alpha, const double* restrict A, const int ndim_A, const int* restrict len_A, const int* restrict lda, const int* restrict idx_A,
                      const double beta,        double* restrict B, const int ndim_B, const int* restrict len_B, const int* restrict ldb, const int* restrict idx_B)
{
    int i, j;
    bool found, done_A, done_B, done_AB;
    double temp;
    int ndim_uniq_A;
    int ndim_uniq_B;
    int ndim_uniq_AB;
    int len_uniq_A[ndim_A];
    int len_uniq_B[ndim_B];
    int len_uniq_AB[ndim_A+ndim_B];
    size_t stride_A[ndim_A];
    size_t stride_B[ndim_B];
    size_t off_A, off_B;
    size_t inc_A_A[ndim_A];
    size_t inc_B_B[ndim_B];
    size_t inc_A_AB[ndim_A+ndim_B];
    size_t inc_B_AB[ndim_A+ndim_B];
    size_t size_A, size_B;
    int pos_A[ndim_A];
    int pos_B[ndim_B];
    int pos_AB[ndim_A+ndim_B];

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(ndim_A, len_A, lda, NULL);
    VALIDATE_TENSOR(ndim_B, len_B, ldb, NULL);
#endif //VALIDATE_INPUTS

    if (lda == NULL)
    {
        if (ndim_A > 0) stride_A[0] = 1;
        for (i = 1;i < ndim_A;i++) stride_A[i] = stride_A[i-1]*len_A[i-1];
    }
    else
    {
        if (ndim_A > 0) stride_A[0] = lda[0];
        for (i = 1;i < ndim_A;i++) stride_A[i] = stride_A[i-1]*lda[i];
    }

    if (ldb == NULL)
    {
        if (ndim_B > 0) stride_B[0] = 1;
        for (i = 1;i < ndim_B;i++) stride_B[i] = stride_B[i-1]*len_B[i-1];
    }
    else
    {
        if (ndim_B > 0) stride_B[0] = ldb[0];
        for (i = 1;i < ndim_B;i++) stride_B[i] = stride_B[i-1]*ldb[i];
    }

    if (ndim_A > 0)
    {
        size_A = stride_A[ndim_A-1]*len_A[ndim_A-1];
    }
    else
    {
        size_A = 1;
    }

    if (ndim_B > 0)
    {
        size_B = stride_B[ndim_B-1]*len_B[ndim_B-1];
    }
    else
    {
        size_B = 1;
    }

    ndim_uniq_A = 0;
    ndim_uniq_B = 0;
    ndim_uniq_AB = 0;
    memset(inc_A_A, 0, ndim_A*sizeof(size_t));
    memset(inc_B_B, 0, ndim_B*sizeof(size_t));
    memset(inc_A_AB, 0, (ndim_A+ndim_B)*sizeof(size_t));
    memset(inc_B_AB, 0, (ndim_A+ndim_B)*sizeof(size_t));

    /*
     * sort the indices into three groups: B only, A and B, and A only
     */
    for (i = 0;i < ndim_B;i++)
    {
        found = false;

        for (j = i+1;j < ndim_B;j++)
        {
            if (idx_B[i] == idx_B[j])
            {
#ifdef VALIDATE_INPUTS
                if (len_B[i] != len_B[j]) return TENSOR_LENGTH_MISMATCH;
#endif //VALIDATE_INPUTS
                found = true;
#ifndef VALIDATE_INPUTS
                break;
#endif //VALIDATE_INPUTS
            }
        }

#ifndef VALIDATE_INPUTS
        if (!found) {
#endif //VALIDATE_INPUTS
        for (j = 0;j < ndim_A;j++)
        {
            if (idx_B[i] == idx_A[j])
            {
#ifdef VALIDATE_INPUTS
                if (len_B[i] != len_A[j]) return TENSOR_LENGTH_MISMATCH;
#endif //VALIDATE_INPUTS
                found = true;
#ifndef VALIDATE_INPUTS
                break;
#endif //VALIDATE_INPUTS
            }
        }
#ifndef VALIDATE_INPUTS
        }
#endif //VALIDATE_INPUTS

        /*
         * idx_B[i] is in B but not A
         */
        if (!found)
        {
            len_uniq_B[ndim_uniq_B] = len_B[i];

            for (j = 0;j <= i;j++)
            {
                if (idx_B[j] == idx_B[i]) inc_B_B[ndim_uniq_B] += stride_B[j];
            }

            ndim_uniq_B++;
        }
    }

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
            for (j = 0;j < ndim_B;j++)
            {
                if (idx_B[j] == idx_A[i])
                {
                    found = true;
                    break;
                }
            }

            /*
             * idx_A[i] is in A and B
             */
            if (found)
            {
                len_uniq_AB[ndim_uniq_AB] = len_A[i];

                for (j = 0;j < ndim_A;j++)
                {
                    if (idx_A[j] == idx_A[i]) inc_A_AB[ndim_uniq_AB] += stride_A[j];
                }

                for (j = 0;j < ndim_B;j++)
                {
                    if (idx_B[j] == idx_A[i]) inc_B_AB[ndim_uniq_AB] += stride_B[j];
                }

                ndim_uniq_AB++;
            }
            /*
             * idx_A[i] is in A but not B
             */
            else
            {
                len_uniq_A[ndim_uniq_A] = len_A[i];

                for (j = 0;j < ndim_A;j++)
                {
                    if (idx_A[j] == idx_A[i]) inc_A_A[ndim_uniq_A] += stride_A[j];
                }

                ndim_uniq_A++;
            }
        }
    }

    off_A = 0;
    off_B = 0;

    /*
     * loop over elements in the first replicate of B (will also change off_A)
     */
    memset(pos_AB, 0, ndim_uniq_AB*sizeof(int));
    for (done_AB = false;!done_AB;)
    {
        temp = 0.0;

        /*
         * loop over elements in A to be summed onto this element of B
         */
        memset(pos_A, 0, ndim_uniq_A*sizeof(int));
        for (done_A = false;!done_A;)
        {
#ifdef CHECK_BOUNDS
            if (off_A < 0 || off_A >= size_A) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

            temp += A[off_A];

            for (i = 0;i < ndim_uniq_A;i++)
            {
                if (pos_A[i] == len_uniq_A[i] - 1)
                {
                    pos_A[i] = 0;
                    off_A -= inc_A_A[i]*(len_uniq_A[i]-1);

                    if (i == ndim_uniq_A - 1)
                    {
                        done_A = true;
                        break;
                    }
                }
                else
                {
                    pos_A[i]++;
                    off_A += inc_A_A[i];
                    break;
                }
            }

            if (ndim_uniq_A == 0) done_A = true;
        }
        /*
         * end loop over A
         */

        temp *= alpha;

        /*
         * loop over replicates of B
         */
        memset(pos_B, 0, ndim_uniq_B*sizeof(int));
        for (done_B = false;!done_B;)
        {
#ifdef CHECK_BOUNDS
            if (off_B < 0 || off_B >= size_B) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

            if (beta == 0.0)
            {
                B[off_B] = temp;
            }
            else
            {
                B[off_B] = temp + beta*B[off_B];
            }

            for (i = 0;i < ndim_uniq_B;i++)
            {
                if (pos_B[i] == len_uniq_B[i] - 1)
                {
                    pos_B[i] = 0;
                    off_B -= inc_B_B[i]*(len_uniq_B[i]-1);

                    if (i == ndim_uniq_B - 1)
                    {
                        done_B = true;
                        break;
                    }
                }
                else
                {
                    pos_B[i]++;
                    off_B += inc_B_B[i];
                    break;
                }
            }

            if (ndim_uniq_B == 0) done_B = true;
        }
        /*
         * end loop over B
         */

        for (i = 0;i < ndim_uniq_AB;i++)
        {
            if (pos_AB[i] == len_uniq_AB[i] - 1)
            {
                pos_AB[i] = 0;
                off_A -= inc_A_AB[i]*(len_uniq_AB[i]-1);
                off_B -= inc_B_AB[i]*(len_uniq_AB[i]-1);

                if (i == ndim_uniq_AB - 1)
                {
                    done_AB = true;
                    break;
                }
            }
            else
            {
                pos_AB[i]++;
                off_A += inc_A_AB[i];
                off_B += inc_B_AB[i];
                break;
            }
        }

        if (ndim_uniq_AB == 0) done_AB = true;
    }
    /*
     * end loop over AB
     */

    return TENSOR_SUCCESS;
}
