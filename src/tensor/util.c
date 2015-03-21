#include "tensor.h"
#include "util.h"

void tensor_info(const int ndim, const int* len, const int* ld, const int* sym,
                 int* group, int* skip, int* stride, size_t* size)
{
    int i;
    int g = 0;
    int p = 0;
    size_t tmp = (ndim == 0 ? 1 : len[0]);
    size_t size_;

    size_ = (ld == NULL || ndim == 0 ? 1 : ld[0]);
    for (i = 0;i < ndim;i++)
    {
        if (group != NULL) group[i] = g;
        if (stride != NULL) stride[i] = size_;
        if (skip != NULL)
            skip[i] = (i == 0 || sym == NULL || sym[i-1] == NS || ld == NULL ? 0 : ld[i]);

        if (sym == NULL || sym[i] == NS)
        {
            size_ *= (ld == NULL || i == ndim-1 ? tmp : ld[i+1]);
            tmp = (i == ndim-1 ? 1 : len[i+1]);
            g++;
            p = 0;
        }
        else if (sym[i] == SY)
        {
            p--;
            tmp = (tmp*(len[i]-p))/(1-p);
        }
        else
        {
            p++;
            tmp = (tmp*(len[i]-p))/(1+p);
        }
    }

    if (size != NULL) *size = size_;
}

int tensor_iterate(const double alpha, const int ntensor, double * const restrict * restrict data, const size_t* restrict size,
                   const int ndim, const int* restrict len, const int* restrict sym,
                   const int* restrict which, const int* restrict stride, const int* restrict pos,
                   const int* restrict skip)
{
    int i, j, minidx, lsym;
    bool done, atmax;
    size_t off[ntensor];
    int idx[ndim];
    double tmp;

    assert(ntensor >= 2);

    if (ndim > 0) idx[0] = 0;
    for (i = 1;i < ndim;i++)
    {
        switch (sym[i-1])
        {
            case NS: idx[i] = 0; break;
            case SY: idx[i] = idx[i-1]; break;
            case AS: idx[i] = idx[i-1]+1; break;
            case EQ: idx[i] = idx[i-1]; break;
        }
        if (idx[i] >= len[i]) return TENSOR_SUCCESS;
    }

    memset(off, 0, sizeof(size_t)*ntensor);

    for (done = false;!done;)
    {
#ifdef CHECK_BOUNDS
        for (i = 0;i < ntensor;i++)
        {
            if (off[i] < 0 || off[i] >= size[i]) return TENSOR_OUT_OF_BOUNDS;
        }
#endif //CHECK_BOUNDS

        tmp = data[0][off[0]];
        for (i = 1;i < ntensor-1;i++) tmp *= data[i][off[i]];
        data[ntensor-1][off[ntensor-1]] += alpha*tmp;

        if (ndim == 0) done = true;
        for (i = 0;i < ndim;)
        {
            lsym = (i == 0 ? NS : sym[i-1]);

            for (j = i;sym[j] == EQ;j++);
            switch (sym[j])
            {
                case NS: atmax = (idx[j] == len[j]-1); break;
                case SY: atmax = (idx[j] == idx[j+1]); break;
                case AS: atmax = (idx[j] == idx[j+1]-1); break;
            }

            if (atmax)
            {
                if (j == ndim) done = true;

                switch (lsym)
                {
                    case NS: minidx = 0; break;
                    case SY: minidx = idx[i-1]; break;
                    case AS: minidx = idx[i-1]+1; break;
                }

                for (;i <= j;i++)
                {
                    off[which[i]] -= offset_between_indices(skip[i]+minidx, skip[i]+i, pos[i], lsym, stride[i]);
                    idx[i] = minidx;
                }
            }
            else
            {
                for (;i <= j;i++)
                {
                    off[which[i]] += offset_to_next_index(skip[i]+i, pos[i], lsym, stride[i]);
                    idx[i]++;
                }
                break;
            }
        }
    }

    return TENSOR_SUCCESS;
}

void index_connectivity(const int ndim_A, const int* sym_A, const int* idx_A,
                        const int ndim_B, const int* sym_B, const int* idx_B,
                        bool G[ndim_A+ndim_B][ndim_A+ndim_B])
{
    int i, j;

    memset(G, 0, sizeof(bool)*(ndim_A+ndim_B)*(ndim_A+ndim_B));

    for (i = 0;i < ndim_A;i++)
    {
        G[i][i] = true;
        for (j = i+1;sym_A[j-1] != NS;j++) G[i][j] = true;

        for (j = i+1;j < ndim_A;j++)
        {
            if (idx_A[i] == idx_A[j]) G[i][j] = true;
        }

        for (j = 0;j < ndim_B;j++)
        {
            if (idx_A[i] == idx_B[j]) G[i][ndim_A+j] = true;
        }
    }

    for (i = 0;i < ndim_B;i++)
    {
        G[ndim_A+i][ndim_A+i] = true;
        for (j = i+1;sym_B[j-1] != NS;j++) G[ndim_A+i][ndim_A+j] = true;

        for (j = i+1;j < ndim_B;j++)
        {
            if (idx_B[i] == idx_B[j]) G[ndim_A+i][ndim_A+j] = true;
        }
    }
}

void dzero(const int n, double* a, const int inca)
{
    int i;
    int ia;

    if (inca == 1)
    {
#ifdef NO_IEEE754
        m = n%5;
        for (i = 0;i < m;i++)
        {
            a[i] = 0.0;
        }
        for (i = m;i < n;i += 5)
        {
            a[i  ] = 0.0;
            a[i+1] = 0.0;
            a[i+2] = 0.0;
            a[i+3] = 0.0;
            a[i+4] = 0.0;
        }
#else
        /*
         * for IEE754, we use the special value +0 = 0x0000000000000000
         */
        memset(a, 0, n*sizeof(double));
#endif
    }
    else
    {
        ia = 0;
        for (i = 0;i < n;i++)
        {
            a[ia] = 0.0;
            ia += inca;
        }
    }
}

void izero(const int n, int* a, const int inca)
{
    int i;
    int ia;

    if (inca == 1)
    {
        memset(a, 0, n*sizeof(int));
    }
    else
    {
        ia = 0;
        for (i = 0;i < n;i++)
        {
            a[ia] = 0;
            ia += inca;
        }
    }
}

int idot(const int n, const int* restrict x, const int incx, const int* restrict y, const int incy)
{
    int i, m;
    int ix, iy;
    int sum;

    sum = 0;

    if (incx == 1 && incy == 1)
    {
        m = n%5;
        for (i = 0;i < m;i++)
        {
            sum += x[i]*y[i];
        }
        for (i = m;i < n;i += 5)
        {
            sum += x[i  ]*y[i  ];
            sum += x[i+1]*y[i+1];
            sum += x[i+2]*y[i+2];
            sum += x[i+3]*y[i+3];
            sum += x[i+4]*y[i+4];
        }
    }
    else
    {
        ix = 0;
        iy = 0;
        for (i = 0;i < n;i++)
        {
            sum += x[ix] * y[iy];
            ix += incx;
            iy += incy;
        }
    }

    return sum;
}

/*
 * n!
 */
uint64_t factorial(const int n)
{
    int i;
    uint64_t j;

    j = 1;

    for (i = n;i > 1;i--)
    {
        j = j * i;
    }

    return j;
}

/**
 * Binomial coefficient.
 *
 * \f[
 * \left(\begin{array}{c}a\\b\end{array}\right) = \frac{a!}{b!(a-b)!},\;0 \le b \le a
 * \f]
 */
uint64_t binomial(const int a, const int b)
{
    int i;
    uint64_t j;

    if (b < 0 || b > a) return 0;

    j = 1;
    for (i = 1;i <= (b < a-b ? b : a-b);i++)
    {
        j = (j * (a - i + 1)) / i;
    }

    return j;
}

size_t offset_to_next_index(const int current_index, const int group_pos, const int sym, const size_t stride)
{
    if (group_pos > 0)
    {
        switch (sym)
        {
            case AS:
            case SH:
                return stride * binomial(current_index, group_pos);
            case SY:
                return stride * binomial(current_index + group_pos, group_pos);
            default:
                return 0;
        }
    }
    else
    {
        return stride;
    }
}

size_t offset_between_indices(const int index1, const int index2, const int group_pos, const int sym, const size_t stride)
{
    if (group_pos > 0)
    {
        switch (sym)
        {
            case AS:
            case SH:
                return stride * (binomial(index2, group_pos+1) - binomial(index1, group_pos+1));
            case SY:
                return stride * (binomial(index2 + group_pos, group_pos+1) -
                                 binomial(index1 + group_pos, group_pos+1));
            default:
                return 0;
        }
    }
    else
    {
        return stride * (index2 - index1);
    }
}

size_t offset_to_index(const int index, const int group_pos, const int sym, const size_t stride)
{
    if (group_pos > 0)
    {
        switch (sym)
        {
            case AS:
            case SH:
                return stride * binomial(index, group_pos+1);
            case SY:
                return stride * binomial(index + group_pos, group_pos+1);
            default:
                return 0;
        }
    }
    else
    {
        return stride * index;
    }
}

/*
 * reverse the sequence in s of length n, returning the number of swaps performed
 */
int reverse_sequence(const int n, int* s)
{
    int tmp, i;

    for (i = 0;i < n / 2;i++)
    {
        tmp = s[i];
        s[i] = s[n-1-i];
        s[n-1-i] = tmp;
    }

    return n / 2;
}

int tensor_sign(const int ndim, const int* sym, const int* idx1, const int* idx2)
{
    int i, j;
    int sign = 1;

    for (i = 0;i < ndim;i = j)
    {
        for (j = i+1;sym[j-1] != NS;j++);
        if (sym[i] == AS) sign *= relative_sign(j-i, idx1+i, idx2+i);
    }

    return sign;
}

int relative_sign(const int n, const int* s1, const int* s2)
{
    int i, j, k;
    int sign = 1;
    bool seen[n];

    memset(seen, 0, sizeof(bool)*n);

    for (i = 0;i < n;i++)
    {
        if (seen[i]) continue;

        j = i;
        while (true)
        {
            for (k = 0;k < n && (!(s1[k] == s2[j]) || seen[k]);k++);
            j = k;
            seen[j] = true;
            if (j == i) break;
            sign = -sign;
        }
    }

    return sign;
}

int next_permutation(const int ndim, const int* sym, int* idx)
{
    int i, i0, l, m, sgn, tmp;

    sgn = 1;

    i0 = 0;
    for (i = 0;i < ndim;i++)
    {
        if (sym[i] == NS)
        {
            /*
             * we now have a symmetry-related group idx[i0]...idx[i]
             * increment to the next permutation
             */

            /*
             * find the largest l s.t. idx[l] < idx[l+1]
             * thus, idx[l+1]...idx[i] is non-increasing
             */
            for (l = i - 1;l >= i0;l--)
            {
                if (idx[l] < idx[l + 1]) break;
            }

            /*
             * if no such l is found, then this is the last permutation
             */
            if (l < i0)
            {
                if (i > 0)
                {
                    /*
                     * the last permutation is in reverse order, fix it and adjust sign
                     * sym[i-1] may not be in the current group, but only when this is a group of size 1
                     * in which case reverse_sequence will return 0
                     */
                    if (reverse_sequence(i + 1 - i0, &idx[i0]) % 2 == 1 && sym[i-1] == AS) sgn = -sgn;
                }

                if (i == ndim - 1) return 0;
            }
            /*
             * otherwise, increment to the next permutation
             */
            else
            {
                /*
                 * find the largest index m s.t. idx[l] < idx[m]
                 */
                for (m = i;m > l;m--)
                {
                    if (idx[l] < idx[m]) break;
                }

                /*
                 * swap idx[l] with idx[m]
                 * sym[i-1] must be in the current group since there is more than
                 * one permutation and hence more than one index
                 */
                tmp = idx[l];
                idx[l] = idx[m];
                idx[m] = tmp;
                if (sym[i-1] == AS) sgn = -sgn;

                /*
                 * and reverse idx[l+1]...idx[i]
                 */
                if (reverse_sequence(i - l, &idx[l + 1]) % 2 == 1 && sym[i-1] == AS) sgn = -sgn;

                break;
            }

            i0 = i + 1;
        }
    }

    return (ndim > 0 ? sgn : 0);
}


void indices_from_labels(const int n, const int* labels, int* used_labels, int* indices)
{
    int i, j;
    bool found;

    for (i = 0;i < n;i++)
    {
        found = false;
        for (j = 0;used_labels[j] != '\0';j++)
        {
            if (used_labels[j] == labels[i])
            {
                indices[i] = j;
                found = true;
                break;
            }
        }

        if (!found)
        {
            used_labels[j] = labels[i];
            indices[i] = j;
        }
    }
}

void reorder_from_indices(const int n, const int* idx1, const int* idx2, int* reorder)
{
    int i, j;

    for (i = 0;i < n;i++)
    {
        for (j = 0;j < n;j++)
        {
            if (idx2[j] == idx1[i])
            {
                reorder[i] = j;
                break;
            }
        }
    }
}

int validate_tensor(const int ndim, const int* len, const int* ld, const int* sym)
{
    int i, k;
    uint64_t size;

    if (ndim < 0) return TENSOR_INVALID_NDIM;

    k = 0;
    size = 1;
    for (i = 0;i < ndim;i++)
    {
        if (len[i] <= 0) return TENSOR_INVALID_LENGTH;

        if (ld != NULL)
        {
            if (i == 0 || !sym || sym[i-1] == NS)
            {
                if (ld[i] <= 0) return TENSOR_INVALID_LD;
                if (ld[i] < size) return TENSOR_LD_TOO_SMALL;
            }
            else
            {
                if (sym[i] != NS && ld[i] != ld[i+1]) return TENSOR_INVALID_LD;
                if (ld[i] < 0) return TENSOR_INVALID_LD;
            }
        }

        if (sym != NULL)
        {
            if (sym[i] != NS && sym[i] != SY &&
                sym[i] != AS && sym[i] != SH) return TENSOR_INVALID_SYMMETRY;

            // the last index can't be < or <= the next one!
            if (i == ndim-1 && sym[i] != NS) return TENSOR_SYMMETRY_MISMATCH;

            // make sure that the edge length is long enough for antisymmetric indices
            if (i == 0 || sym[i-1] == NS)
            {
                size = len[i];
                k = 0;
            }
            else if (sym[i-1] == AS || sym[i-1] == SH)
            {
                size = (size*(len[i]-k))/(1+k);
                k++;
            }
            else
            {
                size = (size*(len[i]-k))/(1-k);
                k--;
            }

            if (len[i] < MAX(k,0)) return TENSOR_INVALID_LENGTH;

            /*
             * check the Golden Rules (1-3):
             *
             * 1) Thou shalt have symmetry-related indices arranged in ascending order
             */
            // enforced by the input format
            /*
             * 2) Thou shalt have the same edge length for all symmetry-related indices.
             */
            if (i < ndim-1 && sym[i] != NS && len[i] != len[i+1]) return TENSOR_LENGTH_MISMATCH;
            /*
             * 3) Thou shalt have the same symmetry relation (AS, SY, or SH) for all symmetry-related indices.
             */
            if (i < ndim-1 && sym[i] != sym[i+1] && sym[i] != NS && sym[i+1] != NS) return TENSOR_SYMMETRY_MISMATCH;
        }
        else
        {
            size = len[i];
        }
    }

    return 0;
}

int connected_components(const int n, const bool G[n][n], int cc[n], int len[n])
{
    int i;
    int ncc = 0;
    bool seen[n];
    int *curcc = cc;

    memset(seen, 0, sizeof(bool)*n);

    for (i = 0;i < n;i++)
    {
        if (seen[i]) continue;
        len[ncc++] = connected_component(n, G, i, seen, &curcc);
    }

    return ncc;
}

int connected_component(const int n, const bool G[n][n], const int i, bool seen[n], int** cc)
{
    int len = 1;
    int j;

    if (seen[i]) return 0;

    seen[i] = true;
    *((*cc)++) = i;
    for (j = i+1;j < n;j++)
    {
        if (!G[i][j]) continue;
        len += connected_component(n, G, j, seen, cc);
    }

    return len;
}
