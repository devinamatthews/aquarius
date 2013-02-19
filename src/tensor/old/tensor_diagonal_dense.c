/**
 * Extract the (semi)diagonal of a tensor and sum onto a second
 *
 * The general form for a diagonal extraction is \f$ab\ldots k*l*\ldots \rightarrow ab\ldots k*l*\ldots where k*
 * denotes the index k appearing one or more times (which need not be the same number of times in A and B).
 * Indices may be transposed.
 *
 * \author Devin Matthews
 * \date Oct. 1 2011
 */

#include "tensor.h"
#include "util.h"
#include <string.h>

int tensor_diagonal_dense_(const double alpha, const double* restrict A, const int ndim_A, const int* restrict len_A, const int* restrict lda, const int* restrict idx_A,
                           const double beta,        double* restrict B, const int ndim_B, const int* restrict len_B, const int* restrict ldb, const int* restrict idx_B)
{
    int i, j;
    bool found, done_AB;
    double temp;
    int ndim_uniq_AB;
    int len_uniq_AB[ndim_A+ndim_B];
    size_t stride_A[ndim_A];
    size_t stride_B[ndim_B];
    size_t off_A, off_B;
    size_t inc_A_AB[ndim_A+ndim_B];
    size_t inc_B_AB[ndim_A+ndim_B];
    size_t size_A, size_B;
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

    ndim_uniq_AB = 0;
    memset(inc_A_AB, 0, (ndim_A+ndim_B)*sizeof(size_t));
    memset(inc_B_AB, 0, (ndim_A+ndim_B)*sizeof(size_t));

    /*
     * all indices must appear in both A and B, but any number of times
     */
    for (i = 0;i < ndim_A;i++)
    {
#ifdef VALIDATE_INPUTS
        found = false;
        for (j = 0;j < ndim_B;j++)
        {
            if (idx_A[i] == idx_B[j])
            {
                if (len_A[i] != len_B[j]) return TENSOR_LENGTH_MISMATCH;
                found = true;
            }
        }
        if (!found) return TENSOR_INDEX_MISMATCH;
#endif //VALIDATE_INPUTS

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
            len_uniq_AB[ndim_uniq_AB] = len_A[i];

            for (j = 0;j <= i;j++)
            {
                if (idx_A[j] == idx_A[i]) inc_A_AB[ndim_uniq_AB] += stride_A[j];
            }

            for (j = 0;j < ndim_B;j++)
            {
                if (idx_B[j] == idx_A[i]) inc_B_AB[ndim_uniq_AB] += stride_B[j];
            }

            ndim_uniq_AB++;
        }
    }

#ifdef VALIDATE_INPUTS
    for (i = 0;i < ndim_B;i++)
    {
        found = false;
        for (j = 0;j < ndim_A;j++)
        {
            if (idx_B[i] == idx_A[j]) found = true;
        }
        if (!found) return TENSOR_INDEX_MISMATCH;
    }
#endif //VALIDATE_INPUTS

    off_A = 0;
    off_B = 0;

    /*
     * loop over elements in A and B
     */
    memset(pos_AB, 0, ndim_uniq_AB*sizeof(int));
    for (done_AB = false;!done_AB;)
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
