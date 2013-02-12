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
 * ARE DISCLAIMED. IN NO EVENT SHALL EDGAR SOLOMONIK BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "internal.h"

shell_t* SLIDE_new_shell(const center_t* pos, const int L, const int nprim, const int ncontr, const bool spherical, const bool contaminants,
                         const double* exponents, const double* coefficients, const int idx[8])
{
    shell_t* shell;
    int i;

    shell = MALLOC(shell_t, 1);
    if (shell == NULL) return NULL;

    shell->pos = (center_t*)pos;
    shell->L = L;
    shell->nprim = nprim;
    shell->ncontr = ncontr;
    shell->spherical = spherical;
    shell->contaminants = contaminants;
    shell->exponents = MALLOC(double, shell->nprim);
    shell->coefficients = MALLOC(double, shell->nprim*shell->ncontr);
    if (shell->exponents == NULL ||
        shell->coefficients == NULL)
    {
        if (shell->exponents != NULL) FREE(shell->exponents);
        if (shell->coefficients != NULL) FREE(shell->coefficients);
        FREE(shell);
        return NULL;
    }

    memcpy(shell->exponents, exponents, sizeof(double)*shell->nprim);
    memcpy(shell->coefficients, coefficients, sizeof(double)*shell->nprim*shell->ncontr);

    memcpy(shell->idx, idx, sizeof(int)*8);

    if (spherical && !contaminants)
    {
        shell->nfunc = angmom[L]->nspher;
    }
    else
    {
        shell->nfunc = angmom[L]->ncart;
    }

    if (SLIDE_shell_sym(shell) == -1)
    {
        SLIDE_free_shell(shell);
        return NULL;
    }

    SLIDE_shell_norm(shell);

    if (ordering == SLIDE_ORDER_SICF || ordering == SLIDE_ORDER_SIFC)
    {
        for (i = 0;i < order-1;i++)
        {
            shell->idx[i+1] = shell->idx[i] + shell->nfunc_irrep[i]*shell->ncontr;
        }
    }
    else if (ordering == SLIDE_ORDER_SFIC || ordering == SLIDE_ORDER_SFCI || ordering == SLIDE_ORDER_SCFI)
    {
        for (i = 0;i < order-1;i++)
        {
            shell->idx[i+1] = shell->idx[i];
        }
    }
    else if (ordering == SLIDE_ORDER_SCIF)
    {
        for (i = 0;i < order-1;i++)
        {
            shell->idx[i+1] = shell->idx[i] + shell->nfunc_irrep[i];
        }
    }

    return shell;
}

shell_t* SLIDE_copy_shell(const shell_t* shell, const bool spherical, const bool contaminants)
{
    shell_t* s;

    s = MALLOC(shell_t, 1);
    if (s == NULL) return NULL;

    memcpy(s, shell, sizeof(shell_t));

    s->exponents = MALLOC(double, shell->nprim);
    s->coefficients = MALLOC(double, shell->nprim*shell->ncontr);
    s->func_irrep = (int(*)[8])MALLOC(int, shell->nfunc*8);
    s->irrep_pos = (int(*)[8])MALLOC(int, shell->nfunc*8);
    s->irreps = (int(*)[8])MALLOC(int, shell->nfunc*8);
    if (s->exponents == NULL ||
        s->coefficients == NULL ||
        s->func_irrep == NULL ||
        s->irrep_pos == NULL ||
        s->irreps == NULL)
    {
        if (s->exponents != NULL) FREE(s->exponents);
        if (s->coefficients != NULL) FREE(s->coefficients);
        if (s->func_irrep != NULL) FREE(s->func_irrep);
        if (s->irrep_pos != NULL) FREE(s->irrep_pos);
        if (s->irreps != NULL) FREE(s->irreps);
        FREE(s);
        return NULL;
    }

    memcpy(s->exponents, shell->exponents, sizeof(double)*shell->nprim);
    memcpy(s->coefficients, shell->coefficients, sizeof(double)*shell->nprim*shell->ncontr);
    memcpy(s->func_irrep, shell->func_irrep, sizeof(int)*shell->nfunc*8);
    memcpy(s->irrep_pos, shell->irrep_pos, sizeof(int)*shell->nfunc*8);
    memcpy(s->irreps, shell->irreps, sizeof(int)*shell->nfunc*8);

    return s;
}

void SLIDE_free_shell(shell_t* shell)
{
    if (shell != NULL)
    {
        if (shell->func_irrep != NULL) FREE(shell->func_irrep);
        if (shell->irrep_pos != NULL) FREE(shell->irrep_pos);
        if (shell->irreps != NULL) FREE(shell->irreps);
        if (shell->exponents != NULL) FREE(shell->exponents);
        if (shell->coefficients != NULL) FREE(shell->coefficients);

        FREE(shell);
    }
}

int SLIDE_shell_sym(shell_t* s)
{
    int irrep, func, i, j, nonzero;
    int proj[8];
    center_t* c;

    c = s->pos;

    s->irreps = (int(*)[8])MALLOC(int, s->nfunc * 8);
    if (s->irreps == NULL) return -1;

    s->func_irrep = (int(*)[8])MALLOC(int, s->nfunc * 8);
    if (s->func_irrep == NULL) return -1;

    s->irrep_pos = (int(*)[8])MALLOC(int, s->nfunc * 8);
    if (s->irrep_pos == NULL) return -1;

    memset(s->nfunc_irrep, 0, sizeof(int)*8);

    /*
     * the tricky part: determine the irrep of each final SO function
     * each AO function will give n SO functions, where n is the number of symmetry-equivalent atoms associated to this shell
     *
     * loop through the functions and for each irrep, determine the projection of this function onto the degenerate centers
     * if this is non-zero, then this is one of the irreps for this function
     *
     * this could potentially be simplified by using the DCR of the atom's stabilizer, but this works
     */
    for (func = 0;func < s->nfunc;func++)
    {
        i = 0;

        for (irrep = 0;irrep < order;irrep++)
        {
            SLIDE_shell_func_centers(s, func, irrep, proj);

            /*
             * check if the projection is non-zero
             */
            nonzero = 0;
            for (j = 0;j < c->degeneracy;j++)
                nonzero += abs(proj[j]);
            if (nonzero > 0)
            {
                s->irrep_pos[func][irrep] = i;
                s->func_irrep[func][i] = s->nfunc_irrep[irrep];
                s->irreps[func][i] = irrep;
                s->nfunc_irrep[irrep]++;
                i++;
            }
        }
    }

    return 0;
}

void SLIDE_shell_func_centers(const shell_t* shell, const int func, const int irrep, int* proj)
{
    int j, op;

    for (j = 0;j < shell->pos->degeneracy;j++)
        proj[j] = 0;

    /*
     * do the projection, using the characters of the irrep and the parity of either the cartesian or spherical functions
     */
    for (op = 0;op < order;op++)
    {
        if (shell->spherical)
        {
            proj[shell->pos->centermap[op]] += chars[irrep][op] * angmom[shell->L]->spherparity[ops[op]][spherical_order[shell->L][func]];
        }
        else
        {
            proj[shell->pos->centermap[op]] += chars[irrep][op] * angmom[shell->L]->cartparity[ops[op]][cartesian_order[shell->L][func]];
        }
    }

    for (j = 0;j < shell->pos->degeneracy;j++)
        proj[j] /= (order/shell->pos->degeneracy);
}

void SLIDE_shell_norm(shell_t* s)
{
    double PI2_N34 = 0.25197943553838073034791409490358;

    double norm;
    int i, j, k, m, n;

    m = s->ncontr;
    n = s->nprim;

    for (i = 0;i < m;i++)
    {
        norm = 0.0;
        for (j = 0;j < n;j++)
        {
            for (k = 0;k < n;k++)
            {
                norm += s->coefficients[i * n + j] * s->coefficients[i * n + k] * pow(2.0 * sqrt(s->exponents[j] * s->exponents[k]) / (s->exponents[j] + s->exponents[k]), (double)s->L + 1.5);
            }
        }

        for (j = 0;j < n;j++)
        {
            s->coefficients[i * n + j] *= PI2_N34 * pow(4 * s->exponents[j], ((double)s->L + 1.5) / 2) / sqrt(norm);
        }
    }
}

//TODO: incorporate cartesian_order etc.
void SLIDE_shell_ao_to_so(const shell_t* shell, const int ao_ordering, double* aoso, const int ld)
{
    int i, j, k, l, m, nc;
    int so, irrep;
    int ao, ao_stride, nao;
    int proj[8];
    int ncontr, nfunc, degen;

    ncontr = shell->ncontr;
    nfunc = shell->nfunc;
    degen = shell->pos->degeneracy;
    nc = angmom[shell->L]->ncart;
    nao = nc*ncontr*degen;

    for (i = 0;i < nfunc;i++)
    {
        for (j = 0;j < ncontr;j++)
        {
            for (k = 0;k < degen;k++)
            {
                irrep = shell->irreps[i][k];

                switch (ordering)
                {
                    case SLIDE_ORDER_ISCF:
                    case SLIDE_ORDER_SICF:
                        so = j*shell->nfunc_irrep[irrep] + shell->func_irrep[i][k];
                        for (l = 0;l < irrep;l++) so += shell->nfunc_irrep[l]*shell->ncontr;
                        break;
                    case SLIDE_ORDER_ISFC:
                    case SLIDE_ORDER_SIFC:
                        so = shell->func_irrep[i][k]*shell->ncontr + j;
                        for (l = 0;l < irrep;l++) so += shell->nfunc_irrep[l]*shell->ncontr;
                        break;
                    case SLIDE_ORDER_SCIF:
                        so = j*shell->nfunc*shell->pos->degeneracy + shell->func_irrep[i][k];
                        for (l = 0;l < irrep;l++) so += shell->nfunc_irrep[l];
                        break;
                    case SLIDE_ORDER_SFIC:
                        so = i*shell->ncontr*shell->pos->degeneracy + shell->irrep_pos[i][irrep]*shell->ncontr + j;
                        break;
                    case SLIDE_ORDER_SCFI:
                        so = j*shell->nfunc*shell->pos->degeneracy + i*shell->pos->degeneracy + shell->irrep_pos[i][irrep];
                        break;
                    case SLIDE_ORDER_SFCI:
                        so = i*shell->ncontr*shell->pos->degeneracy + j*shell->pos->degeneracy + shell->irrep_pos[i][irrep];
                        break;
                }

                zero(nao, aoso+ld*so, 1);

                SLIDE_shell_func_centers(shell, i, irrep, proj);

                for (l = j;l < j+1;l++)
                {
                    for (m = 0;m < degen;m++)
                    {
                        switch (ao_ordering)
                        {
                            case SLIDE_SHELL_ORDER_FCD:
                                ao = m+l*degen;
                                ao_stride = ncontr*degen;
                                break;
                            case SLIDE_SHELL_ORDER_FDC:
                                ao = l+m*ncontr;
                                ao_stride = ncontr*degen;
                                break;
                            case SLIDE_SHELL_ORDER_CFD:
                                ao = m+l*nc*degen;
                                ao_stride = degen;
                                break;
                            case SLIDE_SHELL_ORDER_CDF:
                                ao = (m+l*degen)*nc;
                                ao_stride = 1;
                                break;
                            case SLIDE_SHELL_ORDER_DFC:
                                ao = l+m*nc*ncontr;
                                ao_stride = ncontr;
                                break;
                            case SLIDE_SHELL_ORDER_DCF:
                                ao = (l+m*ncontr)*nc;
                                ao_stride = 1;
                                break;
                        }

                        daxpy(nc, proj[m], angmom[shell->L]->cart2spher+i*nc, 1, aoso+so*ld+ao, ao_stride);
                    }
                }
            }
        }
    }
}
