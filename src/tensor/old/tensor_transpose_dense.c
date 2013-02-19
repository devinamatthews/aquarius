/**
 * Transpose a tensor and sum onto a second.
 *
 * The general form for a transposition operation is \f$ab\ldots \rightarrow P(ab\ldots)\f$
 * where P is some permutation.
 *
 * \author Devin Matthews
 * \date Sep. 28 2011
 */

#include "tensor.h"
#include "util.h"
#include <string.h>

int tensor_transpose_dense_(const double alpha, const double* restrict A, const int ndim_A, const int* restrict len_A, const int* restrict lda, const int* restrict idx_A,
                            const double beta,        double* restrict B, const int ndim_B, const int* restrict len_B, const int* restrict ldb, const int* restrict idx_B)
{
    int i, j;
    bool found, done;
    size_t stride_A[ndim_A];
    size_t stride_B[ndim_B];
    size_t reorder[ndim_B];
    size_t off_A, off_B;
    size_t size_A, size_B;
    int idx[ndim_A];

#ifdef VALIDATE_INPUTS
    if (ndim_A != ndim_B) return TENSOR_INDEX_MISMATCH;
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

    for (i = 0;i < ndim_A;i++)
    {
#ifdef VALIDATE_INPUTS
        /*
         * make sure that each index appears only once in A...
         */
        for (j = i+1;j < ndim_A;j++)
        {
           if (idx_A[i] == idx_A[j]) return TENSOR_INDEX_MISMATCH;
        }

        /*
         * ...and once in B
         */
        found = false;
#endif //VALIDATE_INPUTS
        for (j = 0;j < ndim_B;j++)
        {
           if (idx_A[i] == idx_B[j])
           {
#ifdef VALIDATE_INPUTS
               if (found) return TENSOR_INDEX_MISMATCH;
               if (len_A[i] != len_B[j]) return TENSOR_LENGTH_MISMATCH;
               found = true;
#endif //VALIDATE_INPUTS
               reorder[i] = j;
#ifndef VALIDATE_INPUTS
               break;
#endif //VALIDATE_INPUTS
           }
        }
#ifdef VALIDATE_INPUTS
        if (!found) return TENSOR_INDEX_MISMATCH;
#endif //VALIDATE_INPUTS
    }

    memset(idx, 0, ndim_A*sizeof(int));
    off_A = 0;
    off_B = 0;

    for (done = false;!done;)
    {
#ifdef CHECK_BOUNDS
        if (off_A < 0 || off_A >= size_A) return TENSOR_OUT_OF_BOUNDS;
        if (off_B < 0 || off_B >= size_B) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

        if (beta == 0.0)
        {
            B[off_B] = alpha*A[off_A];
        }
        else
        {
            B[off_B] = alpha*A[off_A] + beta*B[off_B];
        }

        for (i = 0;i < ndim_A;i++)
        {
            if (idx[i] == len_A[i] - 1)
            {
                idx[i] = 0;
                off_A -= stride_A[i]*(len_A[i]-1);
                off_B -= stride_B[reorder[i]]*(len_A[i]-1);

                if (i == ndim_A - 1)
                {
                    done = true;
                    break;
                }
            }
            else
            {
                idx[i]++;
                off_A += stride_A[i];
                off_B += stride_B[reorder[i]];
                break;
            }
        }

        if (ndim_A == 0) done = true;
    }

    return TENSOR_SUCCESS;
}
