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

void transpose(const int m, const int n, const double alpha, const double* A, const int lda,
		                                 const double  beta,       double* B, const int ldb)
{
    size_t i, j;

    if (alpha == 1.0)
    {
    	if (beta == 0.0)
    	{
			for (i = 0;i < m;i++)
			{
				for (j = 0;j < n;j++)
				{
					B[i*ldb + j] = A[j*lda + i];
				}
			}
    	}
    	else if (beta == 1.0)
    	{
			for (i = 0;i < m;i++)
			{
				for (j = 0;j < n;j++)
				{
					B[i*ldb + j] += A[j*lda + i];
				}
			}
    	}
    	else
    	{
			for (i = 0;i < m;i++)
			{
				for (j = 0;j < n;j++)
				{
					B[i*ldb + j] = beta*B[i*ldb + j] + A[j*lda + i];
				}
			}
    	}
    }
    else
    {
    	if (beta == 0.0)
    	{
			for (i = 0;i < m;i++)
			{
				for (j = 0;j < n;j++)
				{
					B[i*ldb + j] = alpha*A[j*lda + i];
				}
			}
    	}
    	else if (beta == 1.0)
    	{
			for (i = 0;i < m;i++)
			{
				for (j = 0;j < n;j++)
				{
					B[i*ldb + j] += alpha*A[j*lda + i];
				}
			}
    	}
    	else
    	{
			for (i = 0;i < m;i++)
			{
				for (j = 0;j < n;j++)
				{
					B[i*ldb + j] = beta*B[i*ldb + j] + alpha*A[j*lda + i];
				}
			}
    	}
    }
}

void ao2so2(const shell_t* a, const shell_t* b, const int r, double* restrict aointegrals, double* restrict sointegrals)
{
    int i, j, e, f, length, irrepe, irrepf;
    double fac;

    length = a->ncontr*b->ncontr;

    for (j = 0;j < b->nfunc;j++)
    {
        for (i = 0;i < a->nfunc;i++)
        {
            for (f = 0;f < b->pos->degeneracy;f++)
            {
                for (e = 0;e < a->pos->degeneracy;e++)
                {
                    irrepe = a->irreps[i][e];
                    irrepf = b->irreps[j][f];
                    if (dirprd[irrepe][irrepf] != 0) continue;

                    if (b->spherical)
                    {
                        fac = angmom[b->L]->spherparity[r][j]*chars[irrepf][r];
                    }
                    else
                    {
                        fac = angmom[b->L]->cartparity[r][j]*chars[irrepf][r];
                    }

                    daxpy(length, fac, aointegrals, 1, sointegrals, 1);

                    sointegrals += length;
                }
            }
            aointegrals += length;
        }
    }
}

void ao2so4(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int r, const int t, const int st,
            double* restrict aointegrals, double* restrict sointegrals)
{
    int i, j, k, l, e, f, g, h;
    double fac;
    int length, irrepe, irrepf, irrepg, irreph;

    length = a->ncontr*b->ncontr*c->ncontr*d->ncontr;

    for (l = 0;l < d->nfunc;l++)
    {
        for (k = 0;k < c->nfunc;k++)
        {
            for (j = 0;j < b->nfunc;j++)
            {
                for (i = 0;i < a->nfunc;i++)
                {
                    for (h = 0;h < d->pos->degeneracy;h++)
                    {
                        for (g = 0;g < c->pos->degeneracy;g++)
                        {
                            for (f = 0;f < b->pos->degeneracy;f++)
                            {
                                for (e = 0;e < a->pos->degeneracy;e++)
                                {
                                    irrepe = a->irreps[i][e];
                                    irrepf = b->irreps[j][f];
                                    irrepg = c->irreps[k][g];
                                    irreph = d->irreps[l][h];
                                    if (dirprd[dirprd[irrepe][irrepf]][dirprd[irrepg][irreph]] != 0) continue;

                                    if (b->spherical)
                                    {
                                        fac = angmom[b->L]->spherparity[r][j]*chars[irrepf][r];
                                    }
                                    else
                                    {
                                        fac = angmom[b->L]->cartparity[r][j]*chars[irrepf][r];
                                    }

                                    if (c->spherical)
                                    {
                                        fac *= angmom[c->L]->spherparity[t][k]*chars[irrepg][t];
                                    }
                                    else
                                    {
                                        fac *= angmom[c->L]->cartparity[t][k]*chars[irrepg][t];
                                    }

                                    if (d->spherical)
                                    {
                                        fac *= angmom[d->L]->spherparity[st][l]*chars[irreph][st];
                                    }
                                    else
                                    {
                                        fac *= angmom[d->L]->cartparity[st][l]*chars[irreph][st];
                                    }

                                    daxpy(length, fac, aointegrals, 1, sointegrals, 1);

                                    sointegrals += length;
                                }
                            }
                        }
                    }
                    aointegrals += length;
                }
            }
        }
    }
}

void cart2spher2r(const shell_t* a, const shell_t* b, const int nprim, double* restrict buf1, double* restrict buf2)
{
    size_t m, n, k;

    m = b->nfunc;
    n = angmom[a->L]->ncart*nprim;
    k = angmom[b->L]->ncart;

    if (b->spherical)
    {
        // [b,j]' x [xa,b]' = [j,xa]
        dgemm('T', 'T', m, n, k, 1.0, angmom[b->L]->cart2spher, k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [xa,j]' = [j,xa]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = a->nfunc;
    n = b->nfunc*nprim;
    k = angmom[a->L]->ncart;

    if (a->spherical)
    {
        // [a,i]' x [jx,a]' = [i,jx]
        dgemm('T', 'T', m, n, k, 1.0, angmom[a->L]->cart2spher, k, buf2, n, 0.0, buf1, m);
    }
    else
    {
        // [jx,i]' = [i,jx]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    dcopy(m*n, buf1, 1, buf2, 1);
}

void cart2spher2l(const shell_t* a, const shell_t* b, const int nprim, double* restrict buf1, double* restrict buf2)
{
    size_t m, n, k;

    n = a->nfunc;
    m = angmom[b->L]->ncart*nprim;
    k = angmom[a->L]->ncart;

    if (a->spherical)
    {
        // [a,bx]' x [a,i]  = [bx,i]
        dgemm('T', 'N', m, n, k, 1.0, buf1, k, angmom[b->L]->cart2spher, k, 0.0, buf2, m);
    }
    else
    {
        // [i,bx]' = [bx,i]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    n = b->nfunc;
    m = a->nfunc*nprim;
    k = angmom[b->L]->ncart;

    if (b->spherical)
    {
        // [b,xi]' x [b,j] = [xi,j]
        dgemm('T', 'N', m, n, k, 1.0, buf2, k, angmom[a->L]->cart2spher, k, 0.0, buf1, m);
    }
    else
    {
        // [j,xi]' = [xi,j]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    dcopy(m*n, buf1, 1, buf2, 1);
}

void cart2spher4r(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int nprim, double* restrict buf1, double* restrict buf2)
{
    size_t m, n, k;

    m = d->nfunc;
    n = angmom[a->L]->ncart*angmom[b->L]->ncart*angmom[c->L]->ncart*nprim;
    k = angmom[d->L]->ncart;

    if (d->spherical)
    {
        // [d,l]' x [xabc,d]' = [l,xabc]
        dgemm('T', 'T', m, n, k, 1.0, angmom[d->L]->cart2spher, k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [xabc,l]' = [l,xabc]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = c->nfunc;
    n = d->nfunc*angmom[a->L]->ncart*angmom[b->L]->ncart*nprim;
    k = angmom[c->L]->ncart;

    if (c->spherical)
    {
        // [c,k]' x [lxab,c]' = [k,lxab]
        dgemm('T', 'T', m, n, k, 1.0, angmom[c->L]->cart2spher, k, buf2, n, 0.0, buf1, m);
    }
    else
    {
        // [lxab,k]' = [k,lxab]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    m = b->nfunc;
    n = c->nfunc*d->nfunc*angmom[a->L]->ncart*nprim;
    k = angmom[b->L]->ncart;

    if (b->spherical)
    {
        // [b,j]' x [klxa,b]' = [j,klxa]
        dgemm('T', 'T', m, n, k, 1.0, angmom[b->L]->cart2spher, k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [klxa,j]' = [j,klxa]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = a->nfunc;
    n = b->nfunc*c->nfunc*d->nfunc*nprim;
    k = angmom[a->L]->ncart;

    if (a->spherical)
    {
        // [a,i]' x [jklx,a]' = [i,jklx]
        dgemm('T', 'T', m, n, k, 1.0, angmom[a->L]->cart2spher, k, buf2, n, 0.0, buf1, m);
    }
    else
    {
        // [jklx,i]' = [i,jklx]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    dcopy(m*n, buf1, 1, buf2, 1);
}

void cart2spher4l(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int nprim, double* restrict buf1, double* restrict buf2)
{
    size_t m, n, k;

    m = a->nfunc;
    n = angmom[b->L]->ncart*angmom[c->L]->ncart*angmom[d->L]->ncart*nprim;
    k = angmom[a->L]->ncart;

    if (a->spherical)
    {
        // [a,bcdx]' x [a,i]' = [bcdx,i]
        dgemm('T', 'N', m, n, k, 1.0, buf1, k, angmom[d->L]->cart2spher, k, 0.0, buf2, m);
    }
    else
    {
        // [i,bcdx]' = [bcdx,i]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = b->nfunc;
    n = a->nfunc*angmom[c->L]->ncart*angmom[d->L]->ncart*nprim;
    k = angmom[b->L]->ncart;

    if (b->spherical)
    {
        // [b,cdxi]' x [b,j]' = [cdxi,j]
        dgemm('T', 'N', m, n, k, 1.0, buf2, k, angmom[c->L]->cart2spher, k, 0.0, buf1, m);
    }
    else
    {
        // [j,cdxi]' = [cdxi,j]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    m = c->nfunc;
    n = a->nfunc*b->nfunc*angmom[d->L]->ncart*nprim;
    k = angmom[c->L]->ncart;

    if (c->spherical)
    {
        // [c,dxij]' x [c,k]' = [dxij,k]
        dgemm('T', 'N', m, n, k, 1.0, buf1, k, angmom[b->L]->cart2spher, k, 0.0, buf2, m);
    }
    else
    {
        // [k,dxij]' = [dxij,k]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = d->nfunc;
    n = a->nfunc*b->nfunc*c->nfunc*nprim;
    k = angmom[d->L]->ncart;

    if (d->spherical)
    {
        // [d,xijk]' x [d,l]' = [xijk,l]
        dgemm('T', 'N', m, n, k, 1.0, buf2, k, angmom[a->L]->cart2spher, k, 0.0, buf1, m);
    }
    else
    {
        // [l,xijk]' = [xijk,l]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    dcopy(m*n, buf1, 1, buf2, 1);
}

void prim2contr2r(const shell_t* a, const shell_t* b, const int nfunc, double* restrict buf1, double* restrict buf2)
{
    size_t m, n, k;

    // [b,j]' x [xa,b]' = [j,xa]
    m = b->ncontr;
    n = a->nprim*nfunc;
    k = b->nprim;
    dgemm('T', 'T', m, n, k, 1.0, b->coefficients, k, buf1, n, 0.0, buf2, m);

    // [a,i]' x [jx,a]' = [i,jx]
    m = a->ncontr;
    n = b->ncontr*nfunc;
    k = a->nprim;
    dgemm('T', 'T', m, n, k, 1.0, a->coefficients, k, buf2, n, 0.0, buf1, m);

    dcopy(m*n, buf1, 1, buf2, 1);
}

void prim2contr2l(const shell_t* a, const shell_t* b, const int nfunc, double* restrict buf1, double* restrict buf2)
{
    size_t m, n, k;

    // [a,bx]' x [a,i] = [bx,i]
    m = a->ncontr;
    n = b->nprim*nfunc;
    k = a->nprim;
    dgemm('T', 'N', m, n, k, 1.0, buf1, k, b->coefficients, k, 0.0, buf2, m);

    // [b,xi]' x [b,j] = [xi,j]
    m = b->ncontr;
    n = a->ncontr*nfunc;
    k = b->nprim;
    dgemm('T', 'N', m, n, k, 1.0, buf2, k, a->coefficients, k, 0.0, buf1, m);

    dcopy(m*n, buf1, 1, buf2, 1);
}

void prim2contr4r(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int nfunc, double* restrict buf1, double* restrict buf2)
{
    size_t m, n, k;

    // [d,l]' x [xabc,d]' = [l,xabc]
    m = d->ncontr;
    n = a->nprim*b->nprim*c->nprim*nfunc;
    k = d->nprim;
    dgemm('T', 'T', m, n, k, 1.0, d->coefficients, k, buf1, n, 0.0, buf2, m);

    // [c,k]' x [lxab,c]' = [k,lxab]
    m = c->ncontr;
    n = d->ncontr*a->nprim*b->nprim*nfunc;
    k = c->nprim;
    dgemm('T', 'T', m, n, k, 1.0, c->coefficients, k, buf2, n, 0.0, buf1, m);

    // [b,j]' x [klxa,b]' = [j,klxa]
    m = b->ncontr;
    n = c->ncontr*d->ncontr*a->nprim*nfunc;
    k = b->nprim;
    dgemm('T', 'T', m, n, k, 1.0, b->coefficients, k, buf1, n, 0.0, buf2, m);

    // [a,i]' x [jkl,xa]' = [i,jklx]
    m = a->ncontr;
    n = b->ncontr*c->ncontr*d->ncontr*nfunc;
    k = a->nprim;
    dgemm('T', 'T', m, n, k, 1.0, a->coefficients, k, buf2, n, 0.0, buf1, m);

    dcopy(m*n, buf1, 1, buf2, 1);
}

void prim2contr4l(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int nfunc, double* restrict buf1, double* restrict buf2)
{
    size_t m, n, k;

    // [a,bcdx]' x [a,i] = [bcdx,i]
    m = a->ncontr;
    n = b->nprim*c->nprim*d->nprim*nfunc;
    k = a->nprim;
    dgemm('T', 'N', m, n, k, 1.0, buf1, k, a->coefficients, k, 0.0, buf2, m);

    // [b,cdxi]' x [b,j] = [cdxi,j]
    m = b->ncontr;
    n = a->ncontr*c->nprim*d->nprim*nfunc;
    k = b->nprim;
    dgemm('T', 'N', m, n, k, 1.0, buf2, k, b->coefficients, k, 0.0, buf1, m);

    // [c,dxij]' x [c,k] = [dxij,k]
    m = c->ncontr;
    n = a->ncontr*b->ncontr*d->nprim*nfunc;
    k = c->nprim;
    dgemm('T', 'N', m, n, k, 1.0, buf1, k, c->coefficients, k, 0.0, buf2, m);

    // [d,xijk]' x [d,l] = [xijk,l]
    m = d->ncontr;
    n = a->ncontr*b->ncontr*c->ncontr*nfunc;
    k = d->nprim;
    dgemm('T', 'N', m, n, k, 1.0, buf2, k, d->coefficients, k, 0.0, buf1, m);

    dcopy(m*n, buf1, 1, buf2, 1);
}
