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

#include "internal.h"

/**
 * generate the roots and weights of the Rys quadrature
 *
 * see Golub, G. H.; Welsch, J. H. Math. Comput. 23, 221-230 (1969)
 *     K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 */
void rysquad(const double T, const int n, double* restrict rt, double* restrict wt)
{
    double a[n], b[n - 1], tmp;
    int i, j, k;
    double R[n + 1][n + 1];
    double fm0;

    for (i = 0;i < n + 1;i++)
    {
        for (j = 0;j < n + 1;j++)
        {
            R[i][j] = fm(T, i + j);
        }
    }

    for (i = 0;i < n + 1;i++)
    {
        for (j = i;j < n + 1;j++)
        {
            for (k = 0;k <= i - 1;k++)
            {
                R[i][j] -= R[k][i] * R[k][j] / R[k][k];
            }
        }
    }

    for (i = 0;i < n + 1;i++)
    {
        tmp = sqrt(fabs(1 / R[i][i]));

        for (j = 0;j < n + 1;j++)
        {
            R[i][j] *= tmp;
        }
    }

    a[0] = R[0][1] / R[0][0];

    for (i = 0;i < n - 1;i++)
    {
        a[i + 1] = R[i + 1][i + 2] / R[i + 1][i + 1] - R[i][i + 1] / R[i][i];
        b[i] = R[i + 1][i + 1] / R[i][i];
    }

    fm0 = fm(T, 0);
    quadrature(n, a, b, fm0, rt, wt);
}

/**
 * generate the roots and weights of the quadrature from the eigenvalues and first element of each eigenvector
 */
void quadrature(const int n, double* restrict a, double* restrict b, const double mu0, double* restrict rt, double* restrict wt)
{
    int i, info;
    double Z[n * n], work[n * 2 - 2];

    info = dstev('V', n, a, b, Z, n, work);
    if (info != 0)
    {
        fprintf(stderr, "rys.c: quadrature: info in dstev = %d\n", info);
        exit(1);
    }

    for (i = 0;i < n;i++)
    {
        rt[i] = a[i];
        wt[i] = Z[i * n] * Z[i * n] * mu0;
    }
}
