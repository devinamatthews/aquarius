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
#include <stdlib.h>
#include <string.h>

int sort_ints(const void *a, const void *b)
{
    return *(const int*)a-*(const int*)b;
}

int tensor_sum_(const double alpha, const double* restrict A, const int ndim_A, const int* restrict len_A, const int* restrict ld_A, const int* restrict sym_A, const int* restrict idx_A,
                const double beta,        double* restrict B, const int ndim_B, const int* restrict len_B, const int* restrict ld_B, const int* restrict sym_B, const int* restrict idx_B)
{
    int i, j, k, l, ncc;
    int ngroup;
    int idx[ndim_A+ndim_B];
    int idx_perm[ndim_A+ndim_B];
    int cc_len[ndim_A+ndim_B];
    int len[ndim_A+ndim_B];
    int sym[ndim_A+ndim_B];
    int sym_perm[ndim_A+ndim_B];
    int which[ndim_A+ndim_B];
    int stride[ndim_A+ndim_B];
    int pos[ndim_A+ndim_B];
    int skip[ndim_A+ndim_B];
    int map[ndim_A+ndim_B];
    int idx_A_sorted[ndim_A];
    int idx_B_sorted[ndim_B];
    int idx_A_perm[ndim_A];
    int idx_B_perm[ndim_B];
    int stride_A[ndim_A];
    int stride_B[ndim_B];
    int pos_A, pos_B;
    int skip_A[ndim_A];
    int skip_B[ndim_B];
    int group_A[ndim_A];
    int group_B[ndim_B];
    size_t size_A, size_B, tmp;
    bool G[ndim_A+ndim_B][ndim_A+ndim_B];
    int permutation[ndim_A+ndim_B];
    int unit_start[ndim_A+ndim_B];
    int unit_len[ndim_A+ndim_B];
    int nunit;
    int sgn;

#ifdef VALIDATE_INPUTS
    VALIDATE_TENSOR(ndim_A, len_A, ld_A, sym_A);
    VALIDATE_TENSOR(ndim_B, len_B, ld_B, sym_B);

    for (i = 0;i < ndim_A;i++)
    {
        for (j = i+1;j < ndim_A;j++)
        {
            if (idx_A[i] == idx_A[j])
            {
                if (len_A[i] != len_A[j]) return TENSOR_LENGTH_MISMATCH;
                if ((sym_A[i] == SY && (sym_A[j] == AS || sym_A[j] == SH)) ||
                    (sym_A[j] == SY && (sym_A[i] == AS || sym_A[i] == SH))) return TENSOR_SYMMETRY_MISMATCH;
            }
        }

        for (j = 0;j < ndim_B;j++)
        {
            if (idx_A[i] == idx_B[j])
            {
                if (len_A[i] != len_B[j]) return TENSOR_LENGTH_MISMATCH;
                if ((sym_A[i] == SY && (sym_B[j] == AS || sym_B[j] == SH)) ||
                    (sym_B[j] == SY && (sym_A[i] == AS || sym_A[i] == SH))) return TENSOR_SYMMETRY_MISMATCH;
            }
        }
    }

    for (i = 0;i < ndim_B;i++)
    {
        for (j = i+1;j < ndim_B;j++)
        {
            if (idx_B[i] == idx_B[j])
            {
                if (len_B[i] != len_B[j]) return TENSOR_LENGTH_MISMATCH;
                if ((sym_B[i] == SY && (sym_B[j] == AS || sym_B[j] == SH)) ||
                    (sym_B[j] == SY && (sym_B[i] == AS || sym_B[i] == SH))) return TENSOR_SYMMETRY_MISMATCH;
            }
        }
    }
#endif //VALIDATE_INPUTS

    /*
    for (i = 0;i < ndim_A;i++) printf("%c ", idx_A[i]);
    printf("\n");
    for (i = 0;i < ndim_B;i++) printf("%c ", idx_B[i]);
    printf("\n");
    */

    memcpy(idx_A_sorted, idx_A, sizeof(int)*ndim_A);
    memcpy(idx_B_sorted, idx_B, sizeof(int)*ndim_B);

    /*
     * Sort indices of A within groups
     */
    for (i = 0;i < ndim_A;i = j)
    {
        for (j = i+1;sym_A[j-1] != NS;j++);
        qsort(idx_A_sorted+i, j-i, sizeof(int), sort_ints);
    }

    /*
     * Sort indices of B within groups
     */
    for (i = 0;i < ndim_B;i = j)
    {
        for (j = i+1;sym_B[j-1] != NS;j++);
        qsort(idx_B_sorted+i, j-i, sizeof(int), sort_ints);
    }

    /*
     * Set up stride, skip, group, and pos information
     */
    tensor_info(ndim_A, len_A, ld_A, sym_A,
                group_A, skip_A, stride_A, &size_A);
    tensor_info(ndim_B, len_B, ld_B, sym_B,
                group_B, skip_B, stride_B, &size_B);

    /*
    for (i = 0;i < ndim_A+ndim_B;i++)
    {
        for (j = 0;j < ndim_A+ndim_B;j++) printf("%d ", G[i][j]);
        printf("\n");
    }
    */

    /*
     * Generate the adjacency matrix for the indices, where adjacency
     * is determined by either relation by symmetry or a repeated label
     */
    index_connectivity(ndim_A, sym_A, idx_A, ndim_B, sym_B, idx_B, G);
    /*
     * The connected components of this graph are the clusters, within which
     * we must permute indices to achieve the proper symmetry
     */
    ncc = connected_components(ndim_A+ndim_B, G, idx, cc_len);

    memset(map, 0, sizeof(int)*(ndim_A+ndim_B));

    /*
    j = 0;
    for (i = 0;i < ncc;i++)
    {
        printf("%d: %d\n", i, cc_len[i]);
        for (k = j;k < j+cc_len[i];k++) printf("%d ", idx[k]);
        printf("\n");
        j += cc_len[i];
    }
    */

    /*
     * Sort global index list within clusters
     * This ensures that if indices always appear together,
     * they are always seen in the same order and can be simply tested
     */
    j = 0;
    for (i = 0;i < ncc;i++)
    {
        map[j] = CLUSTER_|UNIT_;
        for (k = j;k < j+cc_len[i];k++)
        {
            for (l = j+cc_len[i]-1;l > k;l--)
            {
                if (LOOKUP2(idx,idx[l]) < LOOKUP2(idx,idx[l-1]))
                {
                    tmp = idx[l];
                    idx[l] = idx[l-1];
                    idx[l-1] = tmp;
                }
            }
        }
        j += cc_len[i];
    }

    /*
     * Find units, groups of indices which always appear together
     * Also set up permutation-invariant global arrays (sym and len)
     */
    nunit = 0;
    for (i = 0;i < ndim_A+ndim_B;i++)
    {
        //printf("%d = %d\n", i, idx[i]);
        sym[i] = LOOKUP2(sym,idx[i]);
        if (sym[i] == SH) sym[i] = AS;
        len[i] = LOOKUP2(len,idx[i]);

        if (i > 0 && LOOKUP2(idx,idx[i]) != LOOKUP2(idx,idx[i-1]))
        {
            for (j = 0;j < ndim_A;j++)
            {
                if ((idx_A_sorted[j] == LOOKUP2(idx,idx[i]) &&
                     (j == 0 ||
                      idx_A_sorted[j-1] != LOOKUP2(idx,idx[i-1]) ||
                      group_A[j] != group_A[j-1])) ||
                    (idx_A_sorted[j] == LOOKUP2(idx,idx[i-1]) &&
                     (j == ndim_A-1 ||
                      idx_A_sorted[j+1] != LOOKUP2(idx,idx[i]) ||
                      group_A[j+1] != group_A[j]))) map[i] |= UNIT_;
            }

            for (j = 0;j < ndim_B;j++)
            {
                if ((idx_B_sorted[j] == LOOKUP2(idx,idx[i]) &&
                     (j == 0 ||
                      idx_B_sorted[j-1] != LOOKUP2(idx,idx[i-1]) ||
                      group_B[j] != group_B[j-1])) ||
                    (idx_B_sorted[j] == LOOKUP2(idx,idx[i-1]) &&
                     (j == ndim_B-1 ||
                      idx_B_sorted[j+1] != LOOKUP2(idx,idx[i]) ||
                      group_B[j+1] != group_B[j]))) map[i] |= UNIT_;
            }
        }

        if (map[i]&UNIT_)
        {
            nunit++;
            unit_start[nunit-1] = i;
            unit_len[nunit-1] = 0;
        }
        unit_len[nunit-1]++;
        permutation[i] = nunit-1;
    }

    /*
     * Fix symmetry:
     * - Give repeated indices EQ symmetry
     * - Make sure that all indices in the cluster have some
     *   non-NS symmetry except the last (which must be NS)
     */
    j = 0;
    for (i = 0;i < ncc;i++)
    {
        tmp = NS;
        for (k = j;k < j+cc_len[i];k++) if (sym[k] != NS) tmp = sym[k];

        for (k = j;k < j+cc_len[i]-1;k++)
        {
            if (LOOKUP2(idx,idx[k]) == LOOKUP2(idx,idx[k+1]))
            {
                sym[k] = EQ;
            }
            else
            {
                sym[k] = tmp;
            }
        }

        sym[j+cc_len[i]-1] = NS;
        j += cc_len[i];
    }

    do
    {
        k = 0;
        l = 0;
        pos_A = 0;
        pos_B = 0;
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            idx_perm[i] = idx[unit_start[permutation[i]]++];

            if (idx_perm[i] > ndim_A)
            {
                pos[i] = pos_B++;
                stride[i] = stride_B[idx_perm[i]-ndim_A];
                skip[i] = skip_B[idx_perm[i]-ndim_A];
                which[i] = 1;
                if (sym_B[l] == NS) pos_B = 0;
                idx_B_perm[l++] = idx_B_sorted[idx_perm[i]-ndim_A];
            }
            else
            {
                pos[i] = pos_A++;
                stride[i] = stride_A[idx_perm[i]];
                skip[i] = skip_A[idx_perm[i]];
                which[i] = 0;
                if (sym_A[k] == NS) pos_A = 0;
                idx_A_perm[l++] = idx_A_sorted[idx_perm[i]];
            }

            sym_perm[i] = sym[i];
            if (i < ndim_A+ndim_B-1 && sym_perm[i] != NS &&
                idx_perm[i] > idx_perm[i+1] &&
                !G[idx_perm[i]][idx_perm[i+1]])
            {
                sym_perm[i] = (sym_perm[i] == AS ? SY : AS);
            }
        }

        for (i = 0;i < nunit;i++) unit_start[i] -= unit_len[i];

        sgn = tensor_sign(ndim_A, sym_A, idx_A_sorted, idx_A_perm) *
              tensor_sign(ndim_B, sym_B, idx_B_sorted, idx_B_perm);

        char* sym_names[4] = {"  ","<="," <"," ="};

        printf("permutation (%+d):", sgn);
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            printf("%3d  ", permutation[i]);
        }
        printf("\n");

        printf("indices: ");
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            printf("  %c%s", LOOKUP2(idx,idx_perm[i]), sym_names[sym_perm[i]]);
        }
        printf("\n");

        printf("len    : ");
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            printf("%3d  ", len[i]);
        }
        printf("\n");

        printf("group  : ");
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            printf("%3d  ", LOOKUP2(group,idx_perm[i]));
        }
        printf("\n");

        printf("pos    : ");
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            printf("%3d  ", pos[i]);
        }
        printf("\n");

        printf("stride : ");
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            printf("%3d  ", stride[i]);
        }
        printf("\n");

        printf("skip   : ");
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            printf("%3d  ", skip[i]);
        }
        printf("\n");

        printf("which  : ");
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            printf("%3d  ", which[i]);
        }
        printf("\n");

        printf("map    : ");
        for (i = 0;i < ndim_A+ndim_B;i++)
        {
            printf("%3d  ", map[i]);
        }
        printf("\n");

        //return tensor_iterate(alpha*sgn, 2, (double*const[]){A, B}, (const size_t[]){size_A, size_B},
        //                      ndim_A+ndim_B, len, sym, which, stride, pos, skip);
    } while (next_permutation(ndim_A+ndim_B, sym, permutation) != 0);

    return TENSOR_SUCCESS;
}
