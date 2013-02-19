/**
 * Weight a tensor by a second and sum onto a third
 *
 * The general form for a weighting is \f$ab\ldots ef\ldots \times ef\ldots cd\ldots \rightarrow ab\ldots cd\ldots ef\ldots\f$
 * with no indices being summed over. Indices may be transposed in any tensor. Any index group may be empty
 * (in the case that ef... is empty, this reduces to an outer product).
 *
 * \author Devin Matthews
 * \date Oct. 4 2011
 */

#include "tensor.h"
#include "util.h"
#include <string.h>

int tensor_weight_dense_(const double alpha, const double* restrict A, const int ndim_A, const int* restrict len_A, const int* restrict lda, const int* restrict idx_A,
                                             const double* restrict B, const int ndim_B, const int* restrict len_B, const int* restrict ldb, const int* restrict idx_B,
                         const double beta,        double* restrict C, const int ndim_C, const int* restrict len_C, const int* restrict ldc, const int* restrict idx_C)
{
    int i, j;
    bool found, done_AC, done_BC, done_ABC;
    double temp;
    int ndim_uniq_AC;
    int ndim_uniq_BC;
    int ndim_uniq_ABC;
    int len_uniq_AC[ndim_A+ndim_C];
    int len_uniq_BC[ndim_B+ndim_C];
    int len_uniq_ABC[ndim_A+ndim_B+ndim_C];
    size_t stride_A[ndim_A];
    size_t stride_B[ndim_B];
    size_t stride_C[ndim_C];
    size_t off_A, off_B, off_C;
    size_t inc_A_AC[ndim_A+ndim_C];
    size_t inc_C_AC[ndim_A+ndim_C];
    size_t inc_B_BC[ndim_B+ndim_C];
    size_t inc_C_BC[ndim_B+ndim_C];
    size_t inc_A_ABC[ndim_A+ndim_B+ndim_C];
    size_t inc_B_ABC[ndim_A+ndim_B+ndim_C];
    size_t inc_C_ABC[ndim_A+ndim_B+ndim_C];
    size_t size_A, size_B, size_C;
    int pos_AC[ndim_A+ndim_C];
    int pos_BC[ndim_B+ndim_C];
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

    ndim_uniq_AC = 0;
    ndim_uniq_BC = 0;
    ndim_uniq_ABC = 0;
    memset(inc_A_AC, 0, (ndim_A+ndim_C)*sizeof(size_t));
    memset(inc_C_AC, 0, (ndim_A+ndim_C)*sizeof(size_t));
    memset(inc_B_BC, 0, (ndim_B+ndim_C)*sizeof(size_t));
    memset(inc_C_BC, 0, (ndim_B+ndim_C)*sizeof(size_t));
    memset(inc_A_ABC, 0, (ndim_A+ndim_B+ndim_C)*sizeof(size_t));
    memset(inc_B_ABC, 0, (ndim_A+ndim_B+ndim_C)*sizeof(size_t));
    memset(inc_C_ABC, 0, (ndim_A+ndim_B+ndim_C)*sizeof(size_t));

    /*
     * sort the indices into three groups: A and C, B and C, and A, B, and C
     */
    for (i = 0;i < ndim_A;i++)
    {
#ifdef VALIDATE_INPUTS
        for (j = i+1;j < ndim_A;j++)
        {
            if (idx_A[i] == idx_A[j]) return TENSOR_INDEX_MISMATCH;
        }

        found = false;

        for (j = 0;j < ndim_C;j++)
        {
            if (idx_C[j] == idx_A[i])
            {
                if (found) return TENSOR_INDEX_MISMATCH;
                if (len_A[i] != len_C[j]) return TENSOR_LENGTH_MISMATCH;
                found = true;
            }
        }

        if (!found) return TENSOR_INDEX_MISMATCH;
#endif //VALIDATE_INPUTS

        found = false;

        for (j = 0;j < ndim_B;j++)
        {
            if (idx_B[j] == idx_A[i])
            {
#ifdef VALIDATE_INPUTS
                if (found) return TENSOR_INDEX_MISMATCH;
                if (len_A[i] != len_B[j]) return TENSOR_LENGTH_MISMATCH;
#endif //VALIDATE_INPUTS
                found = true;
#ifndef VALIDATE_INPUTS
                break;
#endif //VALIDATE_INPUTS
            }
        }

        if (found)
        {
            len_uniq_ABC[ndim_uniq_ABC] = len_A[i];

            inc_A_ABC[ndim_uniq_ABC] = stride_A[i];

            for (j = 0;j < ndim_B;j++)
            {
                if (idx_B[j] == idx_A[i])
                {
                    inc_B_ABC[ndim_uniq_ABC] = stride_B[j];
                    break;
                }
            }

            for (j = 0;j < ndim_C;j++)
            {
                if (idx_C[j] == idx_A[i])
                {
                    inc_C_ABC[ndim_uniq_ABC] = stride_C[j];
                    break;
                }
            }

            ndim_uniq_ABC++;
        }
        else
        {
            len_uniq_AC[ndim_uniq_AC] = len_A[i];

            inc_A_AC[ndim_uniq_AC] = stride_A[i];

            for (j = 0;j < ndim_C;j++)
            {
                if (idx_C[j] == idx_A[i])
                {
                    inc_C_AC[ndim_uniq_AC] = stride_C[j];
                    break;
                }
            }

            ndim_uniq_AC++;
        }
    }

    for (i = 0;i < ndim_B;i++)
    {
#ifdef VALIDATE_INPUTS
        for (j = i+1;j < ndim_B;j++)
        {
            if (idx_B[i] == idx_B[j]) return TENSOR_INDEX_MISMATCH;
        }

        found = false;

        for (j = 0;j < ndim_C;j++)
        {
            if (idx_C[j] == idx_B[i])
            {
                if (found) return TENSOR_INDEX_MISMATCH;
                if (len_B[i] != len_C[j]) return TENSOR_LENGTH_MISMATCH;
                found = true;
            }
        }

        if (!found) return TENSOR_INDEX_MISMATCH;
#endif //VALIDATE_INPUTS

        found = false;

        for (j = 0;j < ndim_A;j++)
        {
            if (idx_A[j] == idx_B[i])
            {
                found = true;
                break;
            }
        }

        if (found) continue;

        len_uniq_BC[ndim_uniq_BC] = len_B[i];

        inc_B_BC[ndim_uniq_BC] = stride_B[i];

        for (j = 0;j < ndim_C;j++)
        {
            if (idx_C[j] == idx_B[i])
            {
                inc_C_BC[ndim_uniq_BC] = stride_C[j];
                break;
            }
        }

        ndim_uniq_BC++;
    }

#ifdef VALIDATE_INPUTS
    for (i = 0;i < ndim_C;i++)
    {
        found = false;

        for (j = 0;j < ndim_A;j++)
        {
            if (idx_A[j] == idx_C[i])
            {
                found = true;
                break;
            }
        }

        for (j = 0;j < ndim_B;j++)
        {
            if (idx_B[j] == idx_C[i])
            {
                found = true;
                break;
            }
        }

        if (!found) return TENSOR_INDEX_MISMATCH;
    }
#endif //VALIDATE_INPUTS

    off_A = 0;
    off_B = 0;
    off_C = 0;

    /*
     * loop over weighted elements
     */
    memset(pos_ABC, 0, ndim_uniq_ABC*sizeof(int));
    for (done_ABC = false;!done_ABC;)
    {
        if (tensor_size_dense(ndim_uniq_AC, len_uniq_AC, NULL) < tensor_size_dense(ndim_uniq_BC, len_uniq_BC, NULL))
        {
            /*
             * loop over elements in A (will also change off_C)
             */
            memset(pos_AC, 0, ndim_uniq_AC*sizeof(int));
            for (done_AC = false;!done_AC;)
            {
#ifdef CHECK_BOUNDS
                if (off_A < 0 || off_A >= size_A) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

                temp = alpha*A[off_A];

                /*
                 * loop over elements in B (will also change off_C)
                 */
                memset(pos_BC, 0, ndim_uniq_BC*sizeof(int));
                for (done_BC = false;!done_BC;)
                {
#ifdef CHECK_BOUNDS
                    if (off_B < 0 || off_B >= size_B) return TENSOR_OUT_OF_BOUNDS;
                    if (off_C < 0 || off_C >= size_C) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

                    if (beta == 0.0)
                    {
                        C[off_C] = temp*B[off_B];
                    }
                    else
                    {
                        C[off_C] = temp*B[off_B] + beta*C[off_C];
                    }

                    for (i = 0;i < ndim_uniq_BC;i++)
                    {
                        if (pos_BC[i] == len_uniq_BC[i] - 1)
                        {
                            pos_BC[i] = 0;
                            off_B -= inc_B_BC[i]*(len_uniq_BC[i]-1);
                            off_C -= inc_C_BC[i]*(len_uniq_BC[i]-1);

                            if (i == ndim_uniq_BC - 1)
                            {
                                done_BC = true;
                                break;
                            }
                        }
                        else
                        {
                            pos_BC[i]++;
                            off_B += inc_B_BC[i];
                            off_C += inc_C_BC[i];
                            break;
                        }
                    }

                    if (ndim_uniq_BC == 0) done_BC = true;
                }
                /*
                 * end loop over BC
                 */

                for (i = 0;i < ndim_uniq_AC;i++)
                {
                    if (pos_AC[i] == len_uniq_AC[i] - 1)
                    {
                        pos_AC[i] = 0;
                        off_A -= inc_A_AC[i]*(len_uniq_AC[i]-1);
                        off_C -= inc_C_AC[i]*(len_uniq_AC[i]-1);

                        if (i == ndim_uniq_AC - 1)
                        {
                            done_AC = true;
                            break;
                        }
                    }
                    else
                    {
                        pos_AC[i]++;
                        off_A += inc_A_AC[i];
                        off_C += inc_C_AC[i];
                        break;
                    }
                }

                if (ndim_uniq_AC == 0) done_AC = true;
            }
            /*
             * end loop over AC
             */
        }
        else
        {
            /*
             * loop over elements in B (will also change off_C)
             */
            memset(pos_BC, 0, ndim_uniq_BC*sizeof(int));
            for (done_BC = false;!done_BC;)
            {
#ifdef CHECK_BOUNDS
                if (off_B < 0 || off_B >= size_B) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

                temp = alpha*B[off_B];

                /*
                 * loop over elements in A (will also change off_C)
                 */
                memset(pos_AC, 0, ndim_uniq_AC*sizeof(int));
                for (done_AC = false;!done_AC;)
                {
#ifdef CHECK_BOUNDS
                    if (off_A < 0 || off_A >= size_A) return TENSOR_OUT_OF_BOUNDS;
                    if (off_C < 0 || off_C >= size_C) return TENSOR_OUT_OF_BOUNDS;
#endif //CHECK_BOUNDS

                    if (beta == 0.0)
                    {
                        C[off_C] = temp*A[off_A];
                    }
                    else
                    {
                        C[off_C] = temp*A[off_A] + beta*C[off_C];
                    }

                    for (i = 0;i < ndim_uniq_AC;i++)
                    {
                        if (pos_AC[i] == len_uniq_AC[i] - 1)
                        {
                            pos_AC[i] = 0;
                            off_A -= inc_A_AC[i]*(len_uniq_AC[i]-1);
                            off_C -= inc_C_AC[i]*(len_uniq_AC[i]-1);

                            if (i == ndim_uniq_AC - 1)
                            {
                                done_AC = true;
                                break;
                            }
                        }
                        else
                        {
                            pos_AC[i]++;
                            off_A += inc_A_AC[i];
                            off_C += inc_C_AC[i];
                            break;
                        }
                    }

                    if (ndim_uniq_AC == 0) done_AC = true;
                }
                /*
                 * end loop over AC
                 */

                for (i = 0;i < ndim_uniq_BC;i++)
                {
                    if (pos_BC[i] == len_uniq_BC[i] - 1)
                    {
                        pos_BC[i] = 0;
                        off_B -= inc_B_BC[i]*(len_uniq_BC[i]-1);
                        off_C -= inc_C_BC[i]*(len_uniq_BC[i]-1);

                        if (i == ndim_uniq_BC - 1)
                        {
                            done_BC = true;
                            break;
                        }
                    }
                    else
                    {
                        pos_BC[i]++;
                        off_B += inc_B_BC[i];
                        off_C += inc_C_BC[i];
                        break;
                    }
                }

                if (ndim_uniq_BC == 0) done_BC = true;
            }
            /*
             * end loop over BC
             */
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
