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

#include <stdbool.h>

#define TAYLOR_N 8

const static int TMAX[40] = { 33, 37, 40, 43, 46, 49, 51, 53, 56, 58,
                              60, 62, 64, 66, 68, 70, 72, 74, 76, 78,
                              80, 82, 83, 85, 87, 89, 90, 92, 94, 96,
                              97, 99, 101, 102, 104, 106, 108, 110, 112, 114 };

static double FMTABLE[2281][40];

/**
 * Compute \f$F_m(T)\f$
 *
 * by evaluating \f$ \frac{T^{-m-1/2}}{2} \gamma(m+\frac{1}{2},T) \f$
 * where \f$\gamma(x,z)\f$ is the lower incomplete gamma function.
 */
double fm(double T, int m)
{
    const double epsilon = 1e-15;

    double sum, delt, ap;

    PROFILE_FLOPS(1);
    ap = m + 0.5;

    /* Jeff: logical comparison with FP is problematic. */
    assert( !(T<0.0) );
    if (T < 1.e-13)
    {
        PROFILE_FLOPS(DIV_FLOPS);
        return 0.5 / ap;
    }

    if (T > TMAX[m]) return fmasymptotic(T, m);

    PROFILE_FLOPS(2+EXP_FLOPS+DIV_FLOPS);
    sum = 1.0 / (2 * ap * exp(T));
    delt = sum;
    ap += 1.0;

    while (true)
    {
        PROFILE_FLOPS(3+2*DIV_FLOPS);
        delt *= T / ap;
        ap += 1.0;
        if (fabs(delt) / sum < epsilon) break;
        sum += delt;
    }

    return sum;
}

double fmtaylor(double T, int n)
{
    double ans, fac;
    int i;
    int tidx;
    double tr, trmt;

    PROFILE_FLOPS((TAYLOR_N+1)*(3+DIV_FLOPS));

    tidx = (int)round(T * 20.0);
    tr = (double)tidx / 20.0;
    trmt = tr - T;

    ans = 0.0;
    fac = 1.0;

    for (i = 0;i < TAYLOR_N;i++)
    {
        ans += FMTABLE[tidx][n + i] * fac;
        fac *= trmt / (double)(i + 1);
    }

    return ans;
}

double fmasymptotic(double T, int m)
{
    const double PI = 3.1415926535897932384626433832795;

    double sum;
    int i;

    PROFILE_FLOPS(EXP_FLOPS+SQRT_FLOPS+2*DIV_FLOPS+(m-1)*LOG_FLOPS+m+1);

    sum = -m * log(2 * T);
    for (i = 3;i < 2 * m;i += 2)
        sum += log((double)i);

    return exp(sum) * sqrt(PI / T) / 2;
}

void fmrecursive(double T, int n, double* array)
{
    int i;
    double emt;

    if (T > TMAX[n])
    {
        array[n] = fmasymptotic(T, n);
    }
    else
    {
        array[n] = fmtaylor(T, n);
    }

    if (n == 0) return;

    emt = exp(-T);

    for (i = n;i > 0;i--)
    {
        array[i - 1] = (2 * T * array[i] + emt) / (2 * i - 1);
    }

    PROFILE_FLOPS((5+DIV_FLOPS)*(n-1)+EXP_FLOPS);
}

void calcfmtable()
{
    // store tabulated points for T=0,TMAX[m],0.05, m=0,39
    #pragma omp parallel for
    for (int m = 0;m < 40;m++)
    {
        for (int Tidx = 0;Tidx < TMAX[m] * 20 + 1;Tidx++)
        {
            FMTABLE[Tidx][m] = fm(((double)Tidx) / 20.0, m);
        }
    }
}
