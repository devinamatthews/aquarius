/**
 * Contract two tensors into a third.
 *
 * The general form for a contraction is \f$ab\ldots ef\ldots \times ef\ldots cd\ldots \rightarrow ab\ldots cd\ldots\f$
 * where the indices ef... will be summed over. Indices may be transposed in any tensor. Any index group may be empty
 * (in the case that ef... is empty, this reduces to an outer product).
 *
 * \author Devin Matthews
 * \date Oct. 3 2011
 */

#include "tensor.h"
#include "util.h"
#include <string.h>

int tensor_contract_dense_(const double alpha, const double* restrict A, const int ndim_A, const int* restrict len_A, const int* restrict lda, const int* restrict idx_A,
                                               const double* restrict B, const int ndim_B, const int* restrict len_B, const int* restrict ldb, const int* restrict idx_B,
                           const double beta,        double* restrict C, const int ndim_C, const int* restrict len_C, const int* restrict ldc, const int* restrict idx_C)
{
    int i, j;
    bool found, done_AB, done_ABC;
    double temp;
    int ndim_uniq_AB;
    int ndim_uniq_ABC;
    int len_uniq_AB[ndim_A+ndim_B];
    int len_uniq_ABC[ndim_A+ndim_B+ndim_C];
    size_t stride_A[ndim_A];
    size_t stride_B[ndim_B];
    size_t stride_C[ndim_C];
    size_t off_A, off_B, off_C;
    size_t inc_A_AB[ndim_A+ndim_B];
    size_t inc_B_AB[ndim_A+ndim_B];
    size_t inc_A_ABC[ndim_A+ndim_B+ndim_C];
    size_t inc_B_ABC[ndim_A+ndim_B+ndim_C];
    size_t inc_C_ABC[ndim_A+ndim_B+ndim_C];
    size_t size_A, size_B, size_C;
    int pos_AB[ndim_A+ndim_B];
    int pos_ABC[ndim_A+ndim_B+ndim_C];

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(ndim_A, len_A, lda, NULL);
    VALIDATE_TENSOR(ndim_B, len_B, ldb, NULL);
    VALIDATE_TENSOR(ndim_C, len_C, ldc, NULL);
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

    if (ldc == NULL)
    {
        if (ndim_C > 0) stride_C[0] = 1;
        for (i = 1;i < ndim_C;i++) stride_C[i] = stride_C[i-1]*len_C[i-1];
    }
    else
    {
        if (ndim_C > 0) stride_C[0] = ldc[0];
        for (i = 1;i < ndim_C;i++) stride_C[i] = stride_C[i-1]*ldc[i];
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

    if (ndim_C > 0)
    {
        size_C = stride_C[ndim_C-1]*len_C[ndim_C-1];
    }
    else
    {
        size_C = 1;
    }

    ndim_uniq_AB = 0;
    ndim_uniq_ABC = 0;
    memset(inc_A_AB, 0, (ndim_A+ndim_B)*sizeof(size_t));
    memset(inc_B_AB, 0, (ndim_A+ndim_B)*sizeof(size_t));
    memset(inc_A_ABC, 0, (ndim_A+ndim_B+ndim_C)*sizeof(size_t));
    memset(inc_B_ABC, 0, (ndim_A+ndim_B+ndim_C)*sizeof(size_t));
    memset(inc_C_ABC, 0, (ndim_A+ndim_B+ndim_C)*sizeof(size_t));

    for (i = 0;i < ndim_C;i++)
    {
#ifdef VALIDATE_INPUTS
        for (j = i+1;j < ndim_C;j++)
        {
            if (idx_C[i] == idx_C[j]) return TENSOR_INDEX_MISMATCH;
        }

        found = false;

        for (j = 0;j < ndim_A;j++)
        {
            if (idx_C[i] == idx_A[j])
            {
                if (found) return TENSOR_INDEX_MISMATCH;
                if (len_C[i] != len_A[j]) return TENSOR_LENGTH_MISMATCH;
                found = true;
            }
        }

        for (j = 0;j < ndim_B;j++)
        {
            if (idx_C[i] == idx_B[j])
            {
                if (found) return TENSOR_INDEX_MISMATCH;
                if (len_C[i] != len_B[j]) return TENSOR_LENGTH_MISMATCH;
                found = true;
            }
        }

        if (!found) return TENSOR_INDEX_MISMATCH;
#endif //VALIDATE_INPUTS

        len_uniq_ABC[ndim_uniq_ABC] = len_C[i];

        for (j = 0;j < ndim_A;j++)
        {
            if (idx_A[j] == idx_C[i]) inc_A_ABC[ndim_uniq_ABC] += stride_A[j];
        }

        for (j = 0;j < ndim_B;j++)
        {
            if (idx_B[j] == idx_C[i]) inc_B_ABC[ndim_uniq_ABC] += stride_B[j];
        }

        for (j = 0;j < ndim_C;j++)
        {
            if (idx_C[j] == idx_C[i]) inc_C_ABC[ndim_uniq_ABC] += stride_C[j];
        }

        ndim_uniq_ABC++;
    }

    for (i = 0;i < ndim_A;i++)
    {
#ifdef VALIDATE_INPUTS
        for (j = i+1;j < ndim_A;j++)
        {
            if (idx_A[i] == idx_A[j]) return TENSOR_INDEX_MISMATCH;
        }
#endif //VALIDATE_INPUTS

        found = false;
        for (j = 0;j < ndim_C;j++)
        {
            if (idx_A[i] == idx_C[j])
            {
                found = true;
                break;
            }
        }
        if (found) continue;

#ifdef VALIDATE_INPUTS
        found = false;
        for (j = 0;j < ndim_B;j++)
        {
            if (idx_A[i] == idx_B[j])
            {
                if (found) return TENSOR_INDEX_MISMATCH;
                if (len_A[i] != len_B[j]) return TENSOR_LENGTH_MISMATCH;
                found = true;
            }
        }
        if (!found) return TENSOR_INDEX_MISMATCH;
#endif //VALIDATE_INPUTS

        /*
         * idx_A[i] is in both A and B
         */
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

#ifdef VALIDATE_INPUTS
    for (i = 0;i < ndim_B;i++)
    {
        found = false;
        for (j = 0;j < ndim_C;j++)
        {
            if (idx_B[i] == idx_C[j])
            {
                found = true;
                break;
            }
        }
        if (found) continue;

        found = false;
        for (j = 0;j < ndim_A;j++)
        {
            if (idx_B[i] == idx_A[j])
            {
                if (found) return TENSOR_INDEX_MISMATCH;
                if (len_B[i] != len_A[j]) return TENSOR_LENGTH_MISMATCH;
                found = true;
            }
        }
        if (!found) return TENSOR_INDEX_MISMATCH;
    }
#endif //VALIDATE_INPUTS

    off_A = 0;
    off_B = 0;
    off_C = 0;

    /*
     * loop over elements in C (will also change off_A and off_B)
     */
    memset(pos_ABC, 0, ndim_uniq_ABC*sizeof(int));
    for (done_ABC = false;!done_ABC;)
    {
#ifdef CHECK_BOUNDS
        if (off_C < 0 || off_C >= size_C) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

        temp = 0.0;

        /*
         * loop over elements in A an B to be summed onto this element of C
         */
        memset(pos_AB, 0, ndim_uniq_AB*sizeof(int));
        for (done_AB = false;!done_AB;)
        {
#ifdef CHECK_BOUNDS
            if (off_A < 0 || off_A >= size_A) return TENSOR_OUT_OF_BOUNDS;
            if (off_B < 0 || off_B >= size_B) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

            temp += A[off_A]*B[off_B];

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

        temp *= alpha;

        if (beta == 0.0)
        {
            C[off_C] = temp;
        }
        else
        {
            C[off_C] = temp + beta*C[off_C];
        }

        for (i = 0;i < ndim_uniq_ABC;i++)
        {
            if (pos_ABC[i] == len_uniq_ABC[i] - 1)
            {
                pos_ABC[i] = 0;
                off_A -= inc_A_ABC[i]*(len_uniq_ABC[i]-1);
                off_B -= inc_B_ABC[i]*(len_uniq_ABC[i]-1);
                off_C -= inc_C_ABC[i]*(len_uniq_ABC[i]-1);

                if (i == ndim_uniq_ABC - 1)
                {
                    done_ABC = true;
                    break;
                }
            }
            else
            {
                pos_ABC[i]++;
                off_A += inc_A_ABC[i];
                off_B += inc_B_ABC[i];
                off_C += inc_C_ABC[i];
                break;
            }
        }

        if (ndim_uniq_ABC == 0) done_ABC = true;
    }
    /*
     * end loop over ABC
     */

    return TENSOR_SUCCESS;
}
