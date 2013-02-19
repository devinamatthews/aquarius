/**
 * Replicate a tensor and sum onto a second
 *
 * The general form for a replication operation is \f$ab\ldots \rightarrow ab\ldots c*d*\ldots\f$ where c*
 * denotes the index c appearing one or more times. Any indices may be transposed. Replication is
 * performed in-place (meaning that no additional scratch space is required). Even in the case that only a
 * subset of the elements of B are written two, all elements of B are first scaled by beta.
 *
 * \author Devin Matthews
 * \date Sep. 30 2011
 */

#include "tensor.h"
#include "util.h"
#include <string.h>

int tensor_replicate_dense_(const double alpha, const double* restrict A, const int ndim_A, const int* restrict len_A, const int* restrict lda, const int* restrict idx_A,
                            const double beta,        double* restrict B, const int ndim_B, const int* restrict len_B, const int* restrict ldb, const int* restrict idx_B)
{
    int i, j;
    double temp;
    bool found, done_B, done_AB;
    int ndim_uniq_B;
    int ndim_uniq_AB;
    int len_uniq_B[ndim_B];
    int len_uniq_AB[ndim_A+ndim_B];
    size_t stride_A[ndim_A];
    size_t stride_B[ndim_B];
    size_t off_A, off_B;
    size_t inc_B_B[ndim_B];
    size_t inc_A_AB[ndim_A+ndim_B];
    size_t inc_B_AB[ndim_A+ndim_B];
    size_t size_A, size_B;
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

    ndim_uniq_B = 0;
    ndim_uniq_AB = 0;
    memset(inc_B_B, 0, ndim_B*sizeof(size_t));
    memset(inc_A_AB, 0, (ndim_A+ndim_B)*sizeof(size_t));
    memset(inc_B_AB, 0, (ndim_A+ndim_B)*sizeof(size_t));

    /*
     * sort the indices into two groups: B only and A and B
     * indices that appear in A may only appear once in each tensor
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
            }
        }

        for (j = 0;j < ndim_A;j++)
        {
            if (idx_B[i] == idx_A[j])
            {
#ifdef VALIDATE_INPUTS
                if (len_B[i] != len_A[j]) return TENSOR_LENGTH_MISMATCH;
#endif //VALIDATE_INPUTS
                found = true;
            }
        }

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
#ifdef VALIDATE_INPUTS
        for (j = i+1;j < ndim_A;j++)
        {
            if (idx_A[i] == idx_A[j]) return TENSOR_INDEX_MISMATCH;
        }

        found = false;

        for (j = 0;j < ndim_B;j++)
        {
            if (idx_A[i] == idx_B[j])
            {
                if (found) return TENSOR_INDEX_MISMATCH;
                found = true;
            }
        }

        if (!found) return TENSOR_INDEX_MISMATCH;
#endif //VALIDATE_INPUTS

        /*
         * idx_A[i] is in A and B
         */
        len_uniq_AB[ndim_uniq_AB] = len_A[i];

        inc_A_AB[ndim_uniq_AB] = stride_A[i];

        for (j = 0;j < ndim_B;j++)
        {
            if (idx_B[j] == idx_A[i]) inc_B_AB[ndim_uniq_AB] += stride_B[j];
        }

        ndim_uniq_AB++;
    }

    off_A = 0;
    off_B = 0;

    /*
     * loop over elements in this replicate of B and in A
     */
    memset(pos_AB, 0, ndim_uniq_AB*sizeof(int));
    for (done_AB = false;!done_AB;)
    {
#ifdef CHECK_BOUNDS
        if (off_A < 0 || off_A >= size_A) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

        temp = alpha*A[off_A];

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
