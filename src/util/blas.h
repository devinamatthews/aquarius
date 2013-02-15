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

#ifndef _AQUARIUS_UTIL_BLAS_H_
#define _AQUARIUS_UTIL_BLAS_H_

#include "fortran.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * Level 1 BLAS, FORTRAN prototypes
 */
void    FFUNC(drotg,DROTG)  (double* a, double* b, double* c, double* s);
void    FFUNC(drotmg,DROTMG)(double* d1, double* d2, double* a, const double* b, double* param);
void    FFUNC(drot,DROT)    (const integer* n,                            double* x, const integer* incx,       double* y, const integer* incy, const double* c, const double* s);
void    FFUNC(drotm,DROTM)  (const integer* n,                            double* x, const integer* incx,       double* y, const integer* incy, double* param);
void    FFUNC(dswap,DSWAP)  (const integer* n,                            double* x, const integer* incx,       double* y, const integer* incy);
void    FFUNC(dscal,DSCAL)  (const integer* n, const double* alpha,       double* x, const integer* incx);
void    FFUNC(dcopy,DCOPY)  (const integer* n,                      const double* x, const integer* incx,       double* y, const integer* incy);
void    FFUNC(daxpy,DAXPY)  (const integer* n, const double* alpha, const double* x, const integer* incx,       double* y, const integer* incy);
double  FFUNC(ddot,DDOT)    (const integer* n,                      const double* x, const integer* incx, const double* y, const integer* incy);
double  FFUNC(dnrm2,DNRM2)  (const integer* n,                      const double* x, const integer* incx);
double  FFUNC(dasum,DASUM)  (const integer* n,                      const double* x, const integer* incx);
integer FFUNC(idamax,IDAMAX)(const integer* n,                      const double* x, const integer* incx);

/*
 * Level 2 BLAS, FORTRAN prototypes
 */
void FFUNC(dgemv,DGEMV)(                  const char* trans,                   const integer* m, const integer* n,                                       const double* alpha, const double* a, const integer* lda, const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void FFUNC(dgbmv,DGBMV)(                  const char* trans,                   const integer* m, const integer* n, const integer* kl, const integer* ku, const double* alpha, const double* a, const integer* lda, const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void FFUNC(dsymv,DSYMV)(const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* a, const integer* lda, const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void FFUNC(dsbmv,DSBMV)(const char* uplo,                                                        const integer* n, const integer* k,                     const double* alpha, const double* a, const integer* lda, const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void FFUNC(dspmv,DSPMV)(const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* ap,                    const double* x, const integer* incx, const double* beta, double* y, const integer* incy);
void FFUNC(dtrmv,DTRMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                            const double* a, const integer* lda,       double* x, const integer* incx);
void FFUNC(dtbmv,DTBMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                          const double* a, const integer* lda,       double* x, const integer* incx);
void FFUNC(dtpmv,DTPMV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                            const double* ap,                          double* x, const integer* incx);
void FFUNC(dtrsv,DTRSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                            const double* a, const integer* lda,       double* x, const integer* incx);
void FFUNC(dtbsv,DTBSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n, const integer* k,                                          const double* a, const integer* lda,       double* x, const integer* incx);
void FFUNC(dtpsv,DTPSV)(const char* uplo, const char* trans, const char* diag,                   const integer* n,                                                            const double* ap,                          double* x, const integer* incx);
void FFUNC(dger,DGER)  (                                                       const integer* m, const integer* n,                                       const double* alpha, const double* x, const integer* incx, const double* y, const integer* incy, double* a, const integer* lda);
void FFUNC(dsyr,DSYR)  (const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* x, const integer* incx,                                       double* a, const integer* lda);
void FFUNC(dspr,DSPR)  (const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* x, const integer* incx,                                       double* ap);
void FFUNC(dsyr2,DSYR2)(const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* x, const integer* incx, const double* y, const integer* incy, double* a, const integer* lda);
void FFUNC(dspr2,DSPR2)(const char* uplo,                                                        const integer* n,                                       const double* alpha, const double* x, const integer* incx, const double* y, const integer* incy, double* ap);

/*
 * Level 3 BLAS, FORTRAN prototypes
 */
void FFUNC(dgemm,DGEMM)  (                                    const char* transa, const char* transb,                   const integer* m, const integer* n, const integer* k, const double* alpha, const double* a, const integer* lda, const double* b, const integer* ldb, const double* beta, double* c, const integer* ldc);
void FFUNC(dsymm,DSYMM)  (const char* side, const char* uplo,                                                           const integer* m, const integer* n,                   const double* alpha, const double* a, const integer* lda, const double* b, const integer* ldb, const double* beta, double* c, const integer* ldc);
void FFUNC(dsyrk,DSYRK)  (                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const double* alpha, const double* a, const integer* lda,                                      const double* beta, double* c, const integer* ldc);
void FFUNC(dsyr2k,DSYR2K)(                  const char* uplo, const char* trans,                                                          const integer* n, const integer* k, const double* alpha, const double* a, const integer* lda, const double* b, const integer* ldb, const double* beta, double* c, const integer* ldc);
void FFUNC(dtrmm,DTRMM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const double* alpha, const double* a, const integer* lda,       double* b, const integer* ldb);
void FFUNC(dtrsm,DTRSM)  (const char* side, const char* uplo, const char* transa,                     const char* diag, const integer* m, const integer* n,                   const double* alpha, const double* a, const integer* lda,       double* b, const integer* ldb);

/*
 * Level 1 BLAS, C wrappers
 */
static inline void   c_drotg (double* a, double* b, double* c, double* s)
{
    FFUNC(drotg,DROTG)(a, b, c, s);
}

static inline void   c_drotmg(double* d1, double* d2, double* a, const double b, double* param)
{
    FFUNC(drotmg,DROTMG)(d1, d2, a, &b, param);
}

static inline void   c_drot  (const integer n,                           double* x, const integer incx,       double* y, const integer incy, const double c, const double s)
{
    FFUNC(drot,DROT)(&n, x, &incx, y, &incy, &c, &s);
}

static inline void   c_drotm (const integer n,                           double* x, const integer incx,       double* y, const integer incy, double* param)
{
    FFUNC(drotm,DROTM)(&n, x, &incx, y, &incy, param);
}

static inline void   c_dswap (const integer n,                           double* x, const integer incx,       double* y, const integer incy)
{
    FFUNC(dswap,DSWAP)(&n, x, &incx, y, &incy);
}

static inline void   c_dscal (const integer n, const double alpha,       double* x, const integer incx)
{
    FFUNC(dscal,DSCAL)(&n, &alpha, x, &incx);
}

static inline void   c_dcopy (const integer n,                     const double* x, const integer incx,       double* y, const integer incy)
{
    FFUNC(dcopy,DCOPY)(&n, x, &incx, y, &incy);
}

static inline void   c_daxpy (const integer n, const double alpha, const double* x, const integer incx,       double* y, const integer incy)
{
    FFUNC(daxpy,DAXPY)(&n, &alpha, x, &incx, y, &incy);
}

static inline double c_ddot  (const integer n,                     const double* x, const integer incx, const double* y, const integer incy)
{
    return FFUNC(ddot,DDOT)(&n, x, &incx, y, &incy);
}

static inline double c_dnrm2 (const integer n,                     const double* x, const integer incx)
{
    return FFUNC(dnrm2,DNRM2)(&n, x, &incx);
}

static inline double c_dasum (const integer n,                     const double* x, const integer incx)
{
    return FFUNC(dasum,DASUM)(&n, x, &incx);
}

static inline integer    c_idamax(const integer n,                     const double* x, const integer incx)
{
    return FFUNC(idamax,IDAMAX)(&n, x, &incx)-1;
}

/*
 * Level 2 BLAS, C wrappers
 */
static inline void c_dgemv(                 const char trans,                  const integer m, const integer n,                             const double alpha, const double* a, const integer lda, const double* x, const integer incx, const double beta, double* y, const integer incy)
{
    FFUNC(dgemv,DGEMV)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_dgbmv(                 const char trans,                  const integer m, const integer n, const integer kl, const integer ku, const double alpha, const double* a, const integer lda, const double* x, const integer incx, const double beta, double* y, const integer incy)
{
    FFUNC(dgbmv,DGBMV)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_dsymv(const char uplo,                                                 const integer n,                             const double alpha, const double* a, const integer lda, const double* x, const integer incx, const double beta, double* y, const integer incy)
{
    FFUNC(dsymv,DSYMV)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_dsbmv(const char uplo,                                                 const integer n, const integer k,                const double alpha, const double* a, const integer lda, const double* x, const integer incx, const double beta, double* y, const integer incy)
{
    FFUNC(dsbmv,DSBMV)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

static inline void c_dspmv(const char uplo,                                                 const integer n,                             const double alpha, const double* ap,               const double* x, const integer incx, const double beta, double* y, const integer incy)
{
    FFUNC(dspmv,DSPMV)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
}

static inline void c_dtrmv(const char uplo, const char trans, const char diag,              const integer n,                                                 const double* a, const integer lda,       double* x, const integer incx)
{
    FFUNC(dtrmv,DTRMV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_dtbmv(const char uplo, const char trans, const char diag,              const integer n, const integer k,                                    const double* a, const integer lda,       double* x, const integer incx)
{
    FFUNC(dtbmv,DTBMV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_dtpmv(const char uplo, const char trans, const char diag,              const integer n,                                                 const double* ap,                     double* x, const integer incx)
{
    FFUNC(dtpmv,DTPMV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_dtrsv(const char uplo, const char trans, const char diag,              const integer n,                                                 const double* a, const integer lda,       double* x, const integer incx)
{
    FFUNC(dtrsv,DTRSV)(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

static inline void c_dtbsv(const char uplo, const char trans, const char diag,              const integer n, const integer k,                                    const double* a, const integer lda,       double* x, const integer incx)
{
    FFUNC(dtbsv,DTBSV)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

static inline void c_dtpsv(const char uplo, const char trans, const char diag,              const integer n,                                                 const double* ap,                     double* x, const integer incx)
{
    FFUNC(dtpsv,DTPSV)(&uplo, &trans, &diag, &n, ap, x, &incx);
}

static inline void c_dger (                                                    const integer m, const integer n,                             const double alpha, const double* x, const integer incx, const double* y, const integer incy, double* a, const integer lda)
{
    FFUNC(dger,DGER)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_dsyr (const char uplo,                                                 const integer n,                             const double alpha, const double* x, const integer incx,                                  double* a, const integer lda)
{
    FFUNC(dsyr,DSYR)(&uplo, &n, &alpha, x, &incx, a, &lda);
}

static inline void c_dspr (const char uplo,                                                 const integer n,                             const double alpha, const double* x, const integer incx,                                  double* ap)
{
    FFUNC(dspr,DSPR)(&uplo, &n, &alpha, x, &incx, ap);
}

static inline void c_dsyr2(const char uplo,                                                 const integer n,                             const double alpha, const double* x, const integer incx, const double* y, const integer incy, double* a, const integer lda)
{
    FFUNC(dsyr2,DSYR2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void c_dspr2(const char uplo,                                                 const integer n,                             const double alpha, const double* x, const integer incx, const double* y, const integer incy, double* ap)
{
    FFUNC(dspr2,DSPR2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
}

/*
 * Level 3 BLAS, C wrappers
 */
static inline void c_dgemm (                                  const char transa, const char transb,                  const integer m, const integer n, const integer k, const double alpha, const double* a, const integer lda, const double* b, const integer ldb, const double beta, double* c, const integer ldc)
{
    FFUNC(dgemm,DGEMM)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_dsymm (const char side, const char uplo,                                                        const integer m, const integer n,              const double alpha, const double* a, const integer lda, const double* b, const integer ldb, const double beta, double* c, const integer ldc)
{
    FFUNC(dsymm,DSYMM)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_dsyrk (                 const char uplo, const char trans,                                                   const integer n, const integer k, const double alpha, const double* a, const integer lda,                                 const double beta, double* c, const integer ldc)
{
    FFUNC(dsyrk,DSYRK)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

static inline void c_dsyr2k(                 const char uplo, const char trans,                                                   const integer n, const integer k, const double alpha, const double* a, const integer lda, const double* b, const integer ldb, const double beta, double* c, const integer ldc)
{
    FFUNC(dsyr2k,DSYR2K)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

static inline void c_dtrmm (const char side, const char uplo, const char transa,                    const char diag, const integer m, const integer n,              const double alpha, const double* a, const integer lda,       double* b, const integer ldb)
{
    FFUNC(dtrmm,DTRMM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

static inline void c_dtrsm (const char side, const char uplo, const char transa,                    const char diag, const integer m, const integer n,              const double alpha, const double* a, const integer lda,       double* b, const integer ldb)
{
    FFUNC(dtrsm,DTRSM)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

#ifdef __cplusplus
}
#endif

/*
 * #define more familiar names for the C versions
 */
#define drotg  c_drotg
#define drotmg c_drotmg
#define drot   c_drot
#define drotm  c_drotm
#define dswap  c_dswap
#define dscal  c_dscal
#define dcopy  c_dcopy
#define daxpy  c_daxpy
#define ddot   c_ddot
#define dnrm2  c_dnrm2
#define dasum  c_dasum
#define idamax c_idamax
#define dgemv  c_dgemv
#define dsymv  c_dsymv
#define dsbmv  c_dsbmv
#define dspmv  c_dspmv
#define dtrmv  c_dtrmv
#define dtbmv  c_dtbmv
#define dtpmv  c_dtpmv
#define dtrsv  c_dtrsv
#define dtbsv  c_dtbsv
#define dtpsv  c_dtpsv
#define dger   c_dger
#define dsyr   c_dsyr
#define dspr   c_dspr
#define dsyr2  c_dsyr2
#define dspr2  c_dspr2
#define dgemm  c_dgemm
#define dsymm  c_dsymm
#define dsyrk  c_dsyrk
#define ssyr2k c_ssyr2k
#define dsyr2k c_dsyr2k
#define dtrmm  c_dtrmm
#define dtrsm  c_dtrsm

#endif
