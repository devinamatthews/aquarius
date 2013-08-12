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

#ifndef _AQUARIUS_UTIL_H_
#define _AQUARIUS_UTIL_H_

#include "mpi.h"

#ifdef __cplusplus
#ifdef USE_ELEMENTAL
#include "elemental.hpp"
#endif
#endif

#include <math.h>
#include <stdio.h>
#include <stdint.h>

/*
#define INSTANTIATE_SPECIALIZATIONS(name) \
template class name<double>; \
template class name<float>; \
template class name<std::complex<double> >; \
template class name<std::complex<float> >;
*/

#define INSTANTIATE_SPECIALIZATIONS(name) \
template class name<double>;

#define INSTANTIATE_SPECIALIZATIONS_2(name,extra1) \
template class name<double,extra1>;

#define INSTANTIATE_SPECIALIZATIONS_3(name,extra1,extra2) \
template class name<double,extra1,extra2>;

#define CONCAT(...) __VA_ARGS__

#define ASSERT(x,...)               \
if (!(x))                           \
{                                   \
    fprintf(stderr, "%s: ", __func__); \
    fprintf(stderr,__VA_ARGS__);    \
    fprintf(stderr,"\n");           \
    abort();                        \
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define ABS(a) ((a) < 0 ? -(a) : (a))

#define SIZE(m,n,lda) ((m)*(n) == 0 ? 0 : (lda)*((n)-1)+(m))

#define MASTER if (gpid() == 0)

#define NOT_MASTER if (gpid() != 0)

#define RANK(r) if (gpid() == (r))

#define NOT_RANK(r) if (gpid() != (r))

#define PRINT(...) MASTER printf(__VA_ARGS__)

#define ERROR(...) ERROR_AT(__FILE__, __LINE__, __VA_ARGS__)

#define ERROR_AT(file, line, ...) \
do \
{ \
        fprintf(stderr, "%s(%d): ", file, line); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n"); \
    abort(); \
} while (0)

#ifdef DEBUG

#define DPRINTF(...) \
MASTER \
{ \
    printf("%s(%d): ", __FILE__, __LINE__); \
    printf(__VA_ARGS__); \
}

#define DPRINTFC(...) \
MASTER \
{ \
    printf(__VA_ARGS__); \
}

#define ADPRINTF(...) \
do \
{ \
    printf("%s(%d): ", __FILE__, __LINE__); \
    printf(__VA_ARGS__); \
} while (0)

#define ADPRINTFC(...) \
do \
{ \
    printf(__VA_ARGS__); \
} while (0)

#else

#define DPRINTF(...)

#define DPRINTFC(...)

#define ADPRINTF(...)

#define ADPRINTFC(...)

#endif

static inline int gpid()
{
    static int pid = -1;
    if (pid == -1) MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    return pid;
}

static inline int alland(const int what)
{
    int ret = what;
    MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    return ret;
}

static inline int allor(const int what)
{
    int ret = what;
    MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    return ret;
}

static inline double allmax(const double what)
{
    double ret = what;
    MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return ret;
}

static inline double allsum(const double what)
{
    double ret = what;
    MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return ret;
}

static inline long allsuml(const long what)
{
    long ret = what;
    MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    return ret;
}

#ifdef __cplusplus

static inline long allsum(const long what)
{
    long ret = what;
    MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    return ret;
}

static inline long long allsum(const long long what)
{
    long ret = what;
    MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return ret;
}

static inline unsigned long allsum(const unsigned long what)
{
    long ret = what;
    MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    return ret;
}

static inline unsigned long long allsum(const unsigned long long what)
{
    long ret = what;
    MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return ret;
}

template <typename T>
void transpose(const size_t m, const size_t n, const T alpha, const T* A, const size_t lda,
                                               const T  beta,       T* B, const size_t ldb)
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

#endif

static inline int roundup(int x, int y)
{
    return ((x+y-1)/y)*y;
}

static inline double dist2(const double* a, const double* b)
{
    return (a[0]-b[0])*(a[0]-b[0]) +
           (a[1]-b[1])*(a[1]-b[1]) +
           (a[2]-b[2])*(a[2]-b[2]);
}

static inline double dist(const double* a, const double* b)
{
    return sqrt(dist2(a, b));
}

static inline int binom(int a, int b)
{
    int i, j;

    if (b < 0 || b > a) return 0;

    j = 1;
    for (i = 1;i <= MIN(b,a-b);i++)
    {
        j = (j*(a-i+1))/i;
    }

    return j;
}

static inline int64_t fact(int n)
{
    int i;
    int64_t j;

    j = 1;

    for (i = n;i > 1;i--)
    {
        j = j*i;
    }

    return j;
}

static inline int64_t dfact(int n)
{
    int i;
    int64_t j;

    j = 1;

    for (i = n;i > 1;i -= 2)
    {
        j = j*i;
    }

    return j;
}

#endif
