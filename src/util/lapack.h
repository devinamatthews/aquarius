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

#ifndef _AQUARIUS_UTIL_LAPACK_H_
#define _AQUARIUS_UTIL_LAPACK_H_

#include "memory/memory.h"

#include "fortran.h"
#include "blas.h"

enum {AXBX=1,
       ABX=2,
       BAX=3};

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * FORTRAN prototypes
 */
void FFUNC(sgetrf,SGETRF)( integer* m, integer* n, float* a, integer* lda,
                    integer* ipiv, integer *info );

void FFUNC(dgetrf,DGETRF)( integer* m, integer* n, double* a, integer* lda,
                    integer* ipiv, integer *info );

void FFUNC(cgetrf,CGETRF)( integer* m, integer* n, scomplex* a, integer* lda,
                    integer* ipiv, integer *info );

void FFUNC(zgetrf,ZGETRF)( integer* m, integer* n, dcomplex* a, integer* lda,
                    integer* ipiv, integer *info );

void FFUNC(sgbtrf,SGBTRF)( integer* m, integer* n, integer* kl,
                    integer* ku, float* ab, integer* ldab,
                    integer* ipiv, integer *info );

void FFUNC(dgbtrf,DGBTRF)( integer* m, integer* n, integer* kl,
                    integer* ku, double* ab, integer* ldab,
                    integer* ipiv, integer *info );

void FFUNC(cgbtrf,CGBTRF)( integer* m, integer* n, integer* kl,
                    integer* ku, scomplex* ab, integer* ldab,
                    integer* ipiv, integer *info );

void FFUNC(zgbtrf,ZGBTRF)( integer* m, integer* n, integer* kl,
                    integer* ku, dcomplex* ab, integer* ldab,
                    integer* ipiv, integer *info );

void FFUNC(sgttrf,SGTTRF)( integer* n, float* dl, float* d, float* du,
                    float* du2, integer* ipiv, integer *info );

void FFUNC(dgttrf,DGTTRF)( integer* n, double* dl, double* d, double* du,
                    double* du2, integer* ipiv, integer *info );

void FFUNC(cgttrf,CGTTRF)( integer* n, scomplex* dl, scomplex* d, scomplex* du,
                    scomplex* du2, integer* ipiv, integer *info );

void FFUNC(zgttrf,ZGTTRF)( integer* n, dcomplex* dl, dcomplex* d, dcomplex* du,
                    dcomplex* du2, integer* ipiv, integer *info );

void FFUNC(spotrf,SPOTRF)( char* uplo, integer* n, float* a, integer* lda,
                    integer *info );

void FFUNC(dpotrf,DPOTRF)( char* uplo, integer* n, double* a, integer* lda,
                    integer *info );

void FFUNC(cpotrf,CPOTRF)( char* uplo, integer* n, scomplex* a, integer* lda,
                    integer *info );

void FFUNC(zpotrf,ZPOTRF)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    integer *info );

void FFUNC(spstrf,SPSTRF)( char* uplo, integer* n, float* a, integer* lda,
                    integer* piv, integer* rank, float* tol,
                    float* work, integer *info );

void FFUNC(dpstrf,DPSTRF)( char* uplo, integer* n, double* a, integer* lda,
                    integer* piv, integer* rank, double* tol,
                    double* work, integer *info );

void FFUNC(cpstrf,CPSTRF)( char* uplo, integer* n, scomplex* a, integer* lda,
                    integer* piv, integer* rank, float* tol,
                    float* work, integer *info );

void FFUNC(zpstrf,ZPSTRF)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    integer* piv, integer* rank, double* tol,
                    double* work, integer *info );

void FFUNC(spftrf,SPFTRF)( char* transr, char* uplo, integer* n, float* a,
                    integer *info );

void FFUNC(dpftrf,DPFTRF)( char* transr, char* uplo, integer* n, double* a,
                    integer *info );

void FFUNC(cpftrf,CPFTRF)( char* transr, char* uplo, integer* n, scomplex* a,
                    integer *info );

void FFUNC(zpftrf,ZPFTRF)( char* transr, char* uplo, integer* n, dcomplex* a,
                    integer *info );

void FFUNC(spptrf,SPPTRF)( char* uplo, integer* n, float* ap, integer *info );

void FFUNC(dpptrf,DPPTRF)( char* uplo, integer* n, double* ap, integer *info );

void FFUNC(cpptrf,CPPTRF)( char* uplo, integer* n, scomplex* ap, integer *info );

void FFUNC(zpptrf,ZPPTRF)( char* uplo, integer* n, dcomplex* ap, integer *info );

void FFUNC(spbtrf,SPBTRF)( char* uplo, integer* n, integer* kd, float* ab,
                    integer* ldab, integer *info );

void FFUNC(dpbtrf,DPBTRF)( char* uplo, integer* n, integer* kd, double* ab,
                    integer* ldab, integer *info );

void FFUNC(cpbtrf,CPBTRF)( char* uplo, integer* n, integer* kd, scomplex* ab,
                    integer* ldab, integer *info );

void FFUNC(zpbtrf,ZPBTRF)( char* uplo, integer* n, integer* kd, dcomplex* ab,
                    integer* ldab, integer *info );

void FFUNC(spttrf,SPTTRF)( integer* n, float* d, float* e, integer *info );

void FFUNC(dpttrf,DPTTRF)( integer* n, double* d, double* e, integer *info );

void FFUNC(cpttrf,CPTTRF)( integer* n, float* d, scomplex* e, integer *info );

void FFUNC(zpttrf,ZPTTRF)( integer* n, double* d, dcomplex* e, integer *info );

void FFUNC(ssytrf,SSYTRF)( char* uplo, integer* n, float* a, integer* lda,
                    integer* ipiv, float* work, integer* lwork,
                    integer *info );

void FFUNC(dsytrf,DSYTRF)( char* uplo, integer* n, double* a, integer* lda,
                    integer* ipiv, double* work, integer* lwork,
                    integer *info );

void FFUNC(csytrf,CSYTRF)( char* uplo, integer* n, scomplex* a, integer* lda,
                    integer* ipiv, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zsytrf,ZSYTRF)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    integer* ipiv, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(chetrf,CHETRF)( char* uplo, integer* n, scomplex* a, integer* lda,
                    integer* ipiv, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zhetrf,ZHETRF)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    integer* ipiv, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(ssptrf,SSPTRF)( char* uplo, integer* n, float* ap, integer* ipiv,
                    integer *info );

void FFUNC(dsptrf,DSPTRF)( char* uplo, integer* n, double* ap, integer* ipiv,
                    integer *info );

void FFUNC(csptrf,CSPTRF)( char* uplo, integer* n, scomplex* ap, integer* ipiv,
                    integer *info );

void FFUNC(zsptrf,ZSPTRF)( char* uplo, integer* n, dcomplex* ap, integer* ipiv,
                    integer *info );

void FFUNC(chptrf,CHPTRF)( char* uplo, integer* n, scomplex* ap, integer* ipiv,
                    integer *info );

void FFUNC(zhptrf,ZHPTRF)( char* uplo, integer* n, dcomplex* ap, integer* ipiv,
                    integer *info );

void FFUNC(sgetrs,SGETRS)( char* trans, integer* n, integer* nrhs,
                    const float* a, integer* lda, const integer* ipiv,
                    float* b, integer* ldb, integer *info );

void FFUNC(dgetrs,DGETRS)( char* trans, integer* n, integer* nrhs,
                    const double* a, integer* lda, const integer* ipiv,
                    double* b, integer* ldb, integer *info );

void FFUNC(cgetrs,CGETRS)( char* trans, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, const integer* ipiv,
                    scomplex* b, integer* ldb, integer *info );

void FFUNC(zgetrs,ZGETRS)( char* trans, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, const integer* ipiv,
                    dcomplex* b, integer* ldb, integer *info );

void FFUNC(sgbtrs,SGBTRS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const float* ab, integer* ldab,
                    const integer* ipiv, float* b, integer* ldb,
                    integer *info );

void FFUNC(dgbtrs,DGBTRS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const double* ab, integer* ldab,
                    const integer* ipiv, double* b, integer* ldb,
                    integer *info );

void FFUNC(cgbtrs,CGBTRS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const scomplex* ab, integer* ldab,
                    const integer* ipiv, scomplex* b, integer* ldb,
                    integer *info );

void FFUNC(zgbtrs,ZGBTRS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const dcomplex* ab, integer* ldab,
                    const integer* ipiv, dcomplex* b, integer* ldb,
                    integer *info );

void FFUNC(sgttrs,SGTTRS)( char* trans, integer* n, integer* nrhs,
                    const float* dl, const float* d, const float* du,
                    const float* du2, const integer* ipiv, float* b,
                    integer* ldb, integer *info );

void FFUNC(dgttrs,DGTTRS)( char* trans, integer* n, integer* nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* du2, const integer* ipiv, double* b,
                    integer* ldb, integer *info );

void FFUNC(cgttrs,CGTTRS)( char* trans, integer* n, integer* nrhs,
                    const scomplex* dl, const scomplex* d, const scomplex* du,
                    const scomplex* du2, const integer* ipiv, scomplex* b,
                    integer* ldb, integer *info );

void FFUNC(zgttrs,ZGTTRS)( char* trans, integer* n, integer* nrhs,
                    const dcomplex* dl, const dcomplex* d, const dcomplex* du,
                    const dcomplex* du2, const integer* ipiv, dcomplex* b,
                    integer* ldb, integer *info );

void FFUNC(spotrs,SPOTRS)( char* uplo, integer* n, integer* nrhs,
                    const float* a, integer* lda, float* b,
                    integer* ldb, integer *info );

void FFUNC(dpotrs,DPOTRS)( char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, double* b,
                    integer* ldb, integer *info );

void FFUNC(cpotrs,CPOTRS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, scomplex* b,
                    integer* ldb, integer *info );

void FFUNC(zpotrs,ZPOTRS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, dcomplex* b,
                    integer* ldb, integer *info );

void FFUNC(spftrs,SPFTRS)( char* transr, char* uplo, integer* n, integer* nrhs,
                    const float* a, float* b, integer* ldb,
                    integer *info );

void FFUNC(dpftrs,DPFTRS)( char* transr, char* uplo, integer* n, integer* nrhs,
                    const double* a, double* b, integer* ldb,
                    integer *info );

void FFUNC(cpftrs,CPFTRS)( char* transr, char* uplo, integer* n, integer* nrhs,
                    const scomplex* a, scomplex* b, integer* ldb,
                    integer *info );

void FFUNC(zpftrs,ZPFTRS)( char* transr, char* uplo, integer* n, integer* nrhs,
                    const dcomplex* a, dcomplex* b, integer* ldb,
                    integer *info );

void FFUNC(spptrs,SPPTRS)( char* uplo, integer* n, integer* nrhs,
                    const float* ap, float* b, integer* ldb,
                    integer *info );

void FFUNC(dpptrs,DPPTRS)( char* uplo, integer* n, integer* nrhs,
                    const double* ap, double* b, integer* ldb,
                    integer *info );

void FFUNC(cpptrs,CPPTRS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* ap, scomplex* b, integer* ldb,
                    integer *info );

void FFUNC(zpptrs,ZPPTRS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* ap, dcomplex* b, integer* ldb,
                    integer *info );

void FFUNC(spbtrs,SPBTRS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const float* ab, integer* ldab, float* b,
                    integer* ldb, integer *info );

void FFUNC(dpbtrs,DPBTRS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const double* ab, integer* ldab, double* b,
                    integer* ldb, integer *info );

void FFUNC(cpbtrs,CPBTRS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const scomplex* ab, integer* ldab, scomplex* b,
                    integer* ldb, integer *info );

void FFUNC(zpbtrs,ZPBTRS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const dcomplex* ab, integer* ldab, dcomplex* b,
                    integer* ldb, integer *info );

void FFUNC(spttrs,SPTTRS)( integer* n, integer* nrhs, const float* d,
                    const float* e, float* b, integer* ldb,
                    integer *info );

void FFUNC(dpttrs,DPTTRS)( integer* n, integer* nrhs, const double* d,
                    const double* e, double* b, integer* ldb,
                    integer *info );

void FFUNC(cpttrs,CPTTRS)( integer* n, integer* nrhs, const float* d,
                    const scomplex* e, scomplex* b, integer* ldb,
                    integer *info );

void FFUNC(zpttrs,ZPTTRS)( integer* n, integer* nrhs, const double* d,
                    const dcomplex* e, dcomplex* b, integer* ldb,
                    integer *info );

void FFUNC(ssytrs,SSYTRS)( char* uplo, integer* n, integer* nrhs,
                    const float* a, integer* lda, const integer* ipiv,
                    float* b, integer* ldb, integer *info );

void FFUNC(dsytrs,DSYTRS)( char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, const integer* ipiv,
                    double* b, integer* ldb, integer *info );

void FFUNC(csytrs,CSYTRS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, const integer* ipiv,
                    scomplex* b, integer* ldb, integer *info );

void FFUNC(zsytrs,ZSYTRS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, const integer* ipiv,
                    dcomplex* b, integer* ldb, integer *info );

void FFUNC(chetrs,CHETRS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, const integer* ipiv,
                    scomplex* b, integer* ldb, integer *info );

void FFUNC(zhetrs,ZHETRS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, const integer* ipiv,
                    dcomplex* b, integer* ldb, integer *info );

void FFUNC(ssptrs,SSPTRS)( char* uplo, integer* n, integer* nrhs,
                    const float* ap, const integer* ipiv, float* b,
                    integer* ldb, integer *info );

void FFUNC(dsptrs,DSPTRS)( char* uplo, integer* n, integer* nrhs,
                    const double* ap, const integer* ipiv, double* b,
                    integer* ldb, integer *info );

void FFUNC(csptrs,CSPTRS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* ap, const integer* ipiv, scomplex* b,
                    integer* ldb, integer *info );

void FFUNC(zsptrs,ZSPTRS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* ap, const integer* ipiv, dcomplex* b,
                    integer* ldb, integer *info );

void FFUNC(chptrs,CHPTRS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* ap, const integer* ipiv, scomplex* b,
                    integer* ldb, integer *info );

void FFUNC(zhptrs,ZHPTRS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* ap, const integer* ipiv, dcomplex* b,
                    integer* ldb, integer *info );

void FFUNC(strtrs,STRTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const float* a, integer* lda,
                    float* b, integer* ldb, integer *info );

void FFUNC(dtrtrs,DTRTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const double* a, integer* lda,
                    double* b, integer* ldb, integer *info );

void FFUNC(ctrtrs,CTRTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const scomplex* a, integer* lda,
                    scomplex* b, integer* ldb, integer *info );

void FFUNC(ztrtrs,ZTRTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const dcomplex* a, integer* lda,
                    dcomplex* b, integer* ldb, integer *info );

void FFUNC(stptrs,STPTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const float* ap, float* b,
                    integer* ldb, integer *info );

void FFUNC(dtptrs,DTPTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const double* ap, double* b,
                    integer* ldb, integer *info );

void FFUNC(ctptrs,CTPTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const scomplex* ap, scomplex* b,
                    integer* ldb, integer *info );

void FFUNC(ztptrs,ZTPTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const dcomplex* ap, dcomplex* b,
                    integer* ldb, integer *info );

void FFUNC(stbtrs,STBTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const float* ab,
                    integer* ldab, float* b, integer* ldb,
                    integer *info );

void FFUNC(dtbtrs,DTBTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const double* ab,
                    integer* ldab, double* b, integer* ldb,
                    integer *info );

void FFUNC(ctbtrs,CTBTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const scomplex* ab,
                    integer* ldab, scomplex* b, integer* ldb,
                    integer *info );

void FFUNC(ztbtrs,ZTBTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const dcomplex* ab,
                    integer* ldab, dcomplex* b, integer* ldb,
                    integer *info );

void FFUNC(sgecon,SGECON)( char* norm, integer* n, const float* a, integer* lda,
                    float* anorm, float* rcond, float* work,
                    integer* iwork, integer *info );

void FFUNC(dgecon,DGECON)( char* norm, integer* n, const double* a, integer* lda,
                    double* anorm, double* rcond, double* work,
                    integer* iwork, integer *info );

void FFUNC(cgecon,CGECON)( char* norm, integer* n, const scomplex* a, integer* lda,
                    float* anorm, float* rcond, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(zgecon,ZGECON)( char* norm, integer* n, const dcomplex* a, integer* lda,
                    double* anorm, double* rcond, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(sgbcon,SGBCON)( char* norm, integer* n, integer* kl, integer* ku,
                    const float* ab, integer* ldab, const integer* ipiv,
                    float* anorm, float* rcond, float* work,
                    integer* iwork, integer *info );

void FFUNC(dgbcon,DGBCON)( char* norm, integer* n, integer* kl, integer* ku,
                    const double* ab, integer* ldab, const integer* ipiv,
                    double* anorm, double* rcond, double* work,
                    integer* iwork, integer *info );

void FFUNC(cgbcon,CGBCON)( char* norm, integer* n, integer* kl, integer* ku,
                    const scomplex* ab, integer* ldab, const integer* ipiv,
                    float* anorm, float* rcond, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(zgbcon,ZGBCON)( char* norm, integer* n, integer* kl, integer* ku,
                    const dcomplex* ab, integer* ldab, const integer* ipiv,
                    double* anorm, double* rcond, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(sgtcon,SGTCON)( char* norm, integer* n, const float* dl,
                    const float* d, const float* du, const float* du2,
                    const integer* ipiv, float* anorm, float* rcond,
                    float* work, integer* iwork, integer *info );

void FFUNC(dgtcon,DGTCON)( char* norm, integer* n, const double* dl,
                    const double* d, const double* du, const double* du2,
                    const integer* ipiv, double* anorm, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(cgtcon,CGTCON)( char* norm, integer* n, const scomplex* dl,
                    const scomplex* d, const scomplex* du, const scomplex* du2,
                    const integer* ipiv, float* anorm, float* rcond,
                    scomplex* work, integer *info );

void FFUNC(zgtcon,ZGTCON)( char* norm, integer* n, const dcomplex* dl,
                    const dcomplex* d, const dcomplex* du, const dcomplex* du2,
                    const integer* ipiv, double* anorm, double* rcond,
                    dcomplex* work, integer *info );

void FFUNC(spocon,SPOCON)( char* uplo, integer* n, const float* a, integer* lda,
                    float* anorm, float* rcond, float* work,
                    integer* iwork, integer *info );

void FFUNC(dpocon,DPOCON)( char* uplo, integer* n, const double* a, integer* lda,
                    double* anorm, double* rcond, double* work,
                    integer* iwork, integer *info );

void FFUNC(cpocon,CPOCON)( char* uplo, integer* n, const scomplex* a, integer* lda,
                    float* anorm, float* rcond, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(zpocon,ZPOCON)( char* uplo, integer* n, const dcomplex* a, integer* lda,
                    double* anorm, double* rcond, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(sppcon,SPPCON)( char* uplo, integer* n, const float* ap, float* anorm,
                    float* rcond, float* work, integer* iwork,
                    integer *info );

void FFUNC(dppcon,DPPCON)( char* uplo, integer* n, const double* ap, double* anorm,
                    double* rcond, double* work, integer* iwork,
                    integer *info );

void FFUNC(cppcon,CPPCON)( char* uplo, integer* n, const scomplex* ap, float* anorm,
                    float* rcond, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zppcon,ZPPCON)( char* uplo, integer* n, const dcomplex* ap, double* anorm,
                    double* rcond, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(spbcon,SPBCON)( char* uplo, integer* n, integer* kd, const float* ab,
                    integer* ldab, float* anorm, float* rcond,
                    float* work, integer* iwork, integer *info );

void FFUNC(dpbcon,DPBCON)( char* uplo, integer* n, integer* kd, const double* ab,
                    integer* ldab, double* anorm, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(cpbcon,CPBCON)( char* uplo, integer* n, integer* kd, const scomplex* ab,
                    integer* ldab, float* anorm, float* rcond,
                    scomplex* work, float* iwork, integer *info );

void FFUNC(zpbcon,ZPBCON)( char* uplo, integer* n, integer* kd, const dcomplex* ab,
                    integer* ldab, double* anorm, double* rcond,
                    dcomplex* work, double* iwork, integer *info );

void FFUNC(sptcon,SPTCON)( integer* n, const float* d, const float* e,
                    float* anorm, float* rcond, float* work,
                    integer *info );

void FFUNC(dptcon,DPTCON)( integer* n, const double* d, const double* e,
                    double* anorm, double* rcond, double* work,
                    integer *info );

void FFUNC(cptcon,CPTCON)( integer* n, const float* d, const scomplex* e,
                    float* anorm, float* rcond, float* work,
                    integer *info );

void FFUNC(zptcon,ZPTCON)( integer* n, const double* d, const dcomplex* e,
                    double* anorm, double* rcond, double* work,
                    integer *info );

void FFUNC(ssycon,SSYCON)( char* uplo, integer* n, const float* a, integer* lda,
                    const integer* ipiv, float* anorm, float* rcond,
                    float* work, integer* iwork, integer *info );

void FFUNC(dsycon,DSYCON)( char* uplo, integer* n, const double* a, integer* lda,
                    const integer* ipiv, double* anorm, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(csycon,CSYCON)( char* uplo, integer* n, const scomplex* a, integer* lda,
                    const integer* ipiv, float* anorm, float* rcond,
                    scomplex* work, integer *info );

void FFUNC(zsycon,ZSYCON)( char* uplo, integer* n, const dcomplex* a, integer* lda,
                    const integer* ipiv, double* anorm, double* rcond,
                    dcomplex* work, integer *info );

void FFUNC(checon,CHECON)( char* uplo, integer* n, const scomplex* a, integer* lda,
                    const integer* ipiv, float* anorm, float* rcond,
                    scomplex* work, integer *info );

void FFUNC(zhecon,ZHECON)( char* uplo, integer* n, const dcomplex* a, integer* lda,
                    const integer* ipiv, double* anorm, double* rcond,
                    dcomplex* work, integer *info );

void FFUNC(sspcon,SSPCON)( char* uplo, integer* n, const float* ap,
                    const integer* ipiv, float* anorm, float* rcond,
                    float* work, integer* iwork, integer *info );

void FFUNC(dspcon,DSPCON)( char* uplo, integer* n, const double* ap,
                    const integer* ipiv, double* anorm, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(cspcon,CSPCON)( char* uplo, integer* n, const scomplex* ap,
                    const integer* ipiv, float* anorm, float* rcond,
                    scomplex* work, integer *info );

void FFUNC(zspcon,ZSPCON)( char* uplo, integer* n, const dcomplex* ap,
                    const integer* ipiv, double* anorm, double* rcond,
                    dcomplex* work, integer *info );

void FFUNC(chpcon,CHPCON)( char* uplo, integer* n, const scomplex* ap,
                    const integer* ipiv, float* anorm, float* rcond,
                    scomplex* work, integer *info );

void FFUNC(zhpcon,ZHPCON)( char* uplo, integer* n, const dcomplex* ap,
                    const integer* ipiv, double* anorm, double* rcond,
                    dcomplex* work, integer *info );

void FFUNC(strcon,STRCON)( char* norm, char* uplo, char* diag, integer* n,
                    const float* a, integer* lda, float* rcond,
                    float* work, integer* iwork, integer *info );

void FFUNC(dtrcon,DTRCON)( char* norm, char* uplo, char* diag, integer* n,
                    const double* a, integer* lda, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(ctrcon,CTRCON)( char* norm, char* uplo, char* diag, integer* n,
                    const scomplex* a, integer* lda, float* rcond,
                    scomplex* work, float* rwork, integer *info );

void FFUNC(ztrcon,ZTRCON)( char* norm, char* uplo, char* diag, integer* n,
                    const dcomplex* a, integer* lda, double* rcond,
                    dcomplex* work, double* rwork, integer *info );

void FFUNC(stpcon,STPCON)( char* norm, char* uplo, char* diag, integer* n,
                    const float* ap, float* rcond, float* work,
                    integer* iwork, integer *info );

void FFUNC(dtpcon,DTPCON)( char* norm, char* uplo, char* diag, integer* n,
                    const double* ap, double* rcond, double* work,
                    integer* iwork, integer *info );

void FFUNC(ctpcon,CTPCON)( char* norm, char* uplo, char* diag, integer* n,
                    const scomplex* ap, float* rcond, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(ztpcon,ZTPCON)( char* norm, char* uplo, char* diag, integer* n,
                    const dcomplex* ap, double* rcond, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(stbcon,STBCON)( char* norm, char* uplo, char* diag, integer* n,
                    integer* kd, const float* ab, integer* ldab,
                    float* rcond, float* work, integer* iwork,
                    integer *info );

void FFUNC(dtbcon,DTBCON)( char* norm, char* uplo, char* diag, integer* n,
                    integer* kd, const double* ab, integer* ldab,
                    double* rcond, double* work, integer* iwork,
                    integer *info );

void FFUNC(ctbcon,CTBCON)( char* norm, char* uplo, char* diag, integer* n,
                    integer* kd, const scomplex* ab, integer* ldab,
                    float* rcond, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(ztbcon,ZTBCON)( char* norm, char* uplo, char* diag, integer* n,
                    integer* kd, const dcomplex* ab, integer* ldab,
                    double* rcond, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(sgerfs,SGERFS)( char* trans, integer* n, integer* nrhs,
                    const float* a, integer* lda, const float* af,
                    integer* ldaf, const integer* ipiv, const float* b,
                    integer* ldb, float* x, integer* ldx, float* ferr,
                    float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dgerfs,DGERFS)( char* trans, integer* n, integer* nrhs,
                    const double* a, integer* lda, const double* af,
                    integer* ldaf, const integer* ipiv, const double* b,
                    integer* ldb, double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(cgerfs,CGERFS)( char* trans, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, const scomplex* af,
                    integer* ldaf, const integer* ipiv, const scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, float* ferr,
                    float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zgerfs,ZGERFS)( char* trans, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, const dcomplex* af,
                    integer* ldaf, const integer* ipiv, const dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, double* ferr,
                    double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(sgerfsx,SGERFSX)( char* trans, char* equed, integer* n, integer* nrhs,
                     const float* a, integer* lda, const float* af,
                     integer* ldaf, const integer* ipiv, const float* r,
                     const float* c, const float* b, integer* ldb,
                     float* x, integer* ldx, float* rcond, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     float* work, integer* iwork, integer *info );

void FFUNC(dgerfsx,DGERFSX)( char* trans, char* equed, integer* n, integer* nrhs,
                     const double* a, integer* lda, const double* af,
                     integer* ldaf, const integer* ipiv, const double* r,
                     const double* c, const double* b, integer* ldb,
                     double* x, integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(cgerfsx,CGERFSX)( char* trans, char* equed, integer* n, integer* nrhs,
                     const scomplex* a, integer* lda, const scomplex* af,
                     integer* ldaf, const integer* ipiv, const float* r,
                     const float* c, const scomplex* b, integer* ldb,
                     scomplex* x, integer* ldx, float* rcond, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     scomplex* work, float* rwork, integer *info );

void FFUNC(zgerfsx,ZGERFSX)( char* trans, char* equed, integer* n, integer* nrhs,
                     const dcomplex* a, integer* lda, const dcomplex* af,
                     integer* ldaf, const integer* ipiv, const double* r,
                     const double* c, const dcomplex* b, integer* ldb,
                     dcomplex* x, integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     dcomplex* work, double* rwork, integer *info );

void FFUNC(sgbrfs,SGBRFS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const float* ab, integer* ldab,
                    const float* afb, integer* ldafb,
                    const integer* ipiv, const float* b, integer* ldb,
                    float* x, integer* ldx, float* ferr, float* berr,
                    float* work, integer* iwork, integer *info );

void FFUNC(dgbrfs,DGBRFS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const double* ab, integer* ldab,
                    const double* afb, integer* ldafb,
                    const integer* ipiv, const double* b, integer* ldb,
                    double* x, integer* ldx, double* ferr, double* berr,
                    double* work, integer* iwork, integer *info );

void FFUNC(cgbrfs,CGBRFS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const scomplex* ab, integer* ldab,
                    const scomplex* afb, integer* ldafb,
                    const integer* ipiv, const scomplex* b, integer* ldb,
                    scomplex* x, integer* ldx, float* ferr, float* berr,
                    scomplex* work, float* rwork, integer *info );

void FFUNC(zgbrfs,ZGBRFS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const dcomplex* ab, integer* ldab,
                    const dcomplex* afb, integer* ldafb,
                    const integer* ipiv, const dcomplex* b, integer* ldb,
                    dcomplex* x, integer* ldx, double* ferr, double* berr,
                    dcomplex* work, double* rwork, integer *info );

void FFUNC(sgbrfsx,SGBRFSX)( char* trans, char* equed, integer* n, integer* kl,
                     integer* ku, integer* nrhs, const float* ab,
                     integer* ldab, const float* afb, integer* ldafb,
                     const integer* ipiv, const float* r, const float* c,
                     const float* b, integer* ldb, float* x,
                     integer* ldx, float* rcond, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     float* work, integer* iwork, integer *info );

void FFUNC(dgbrfsx,DGBRFSX)( char* trans, char* equed, integer* n, integer* kl,
                     integer* ku, integer* nrhs, const double* ab,
                     integer* ldab, const double* afb, integer* ldafb,
                     const integer* ipiv, const double* r, const double* c,
                     const double* b, integer* ldb, double* x,
                     integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(cgbrfsx,CGBRFSX)( char* trans, char* equed, integer* n, integer* kl,
                     integer* ku, integer* nrhs, const scomplex* ab,
                     integer* ldab, const scomplex* afb, integer* ldafb,
                     const integer* ipiv, const float* r, const float* c,
                     const scomplex* b, integer* ldb, scomplex* x,
                     integer* ldx, float* rcond, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     scomplex* work, float* rwork, integer *info );

void FFUNC(zgbrfsx,ZGBRFSX)( char* trans, char* equed, integer* n, integer* kl,
                     integer* ku, integer* nrhs, const dcomplex* ab,
                     integer* ldab, const dcomplex* afb, integer* ldafb,
                     const integer* ipiv, const double* r, const double* c,
                     const dcomplex* b, integer* ldb, dcomplex* x,
                     integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     dcomplex* work, double* rwork, integer *info );

void FFUNC(sgtrfs,SGTRFS)( char* trans, integer* n, integer* nrhs,
                    const float* dl, const float* d, const float* du,
                    const float* dlf, const float* df, const float* duf,
                    const float* du2, const integer* ipiv, const float* b,
                    integer* ldb, float* x, integer* ldx, float* ferr,
                    float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dgtrfs,DGTRFS)( char* trans, integer* n, integer* nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* dlf, const double* df, const double* duf,
                    const double* du2, const integer* ipiv, const double* b,
                    integer* ldb, double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(cgtrfs,CGTRFS)( char* trans, integer* n, integer* nrhs,
                    const scomplex* dl, const scomplex* d, const scomplex* du,
                    const scomplex* dlf, const scomplex* df, const scomplex* duf,
                    const scomplex* du2, const integer* ipiv, const scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, float* ferr,
                    float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zgtrfs,ZGTRFS)( char* trans, integer* n, integer* nrhs,
                    const dcomplex* dl, const dcomplex* d, const dcomplex* du,
                    const dcomplex* dlf, const dcomplex* df, const dcomplex* duf,
                    const dcomplex* du2, const integer* ipiv, const dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, double* ferr,
                    double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(sporfs,SPORFS)( char* uplo, integer* n, integer* nrhs,
                    const float* a, integer* lda, const float* af,
                    integer* ldaf, const float* b, integer* ldb,
                    float* x, integer* ldx, float* ferr, float* berr,
                    float* work, integer* iwork, integer *info );

void FFUNC(dporfs,DPORFS)( char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, const double* af,
                    integer* ldaf, const double* b, integer* ldb,
                    double* x, integer* ldx, double* ferr, double* berr,
                    double* work, integer* iwork, integer *info );

void FFUNC(cporfs,CPORFS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, const scomplex* af,
                    integer* ldaf, const scomplex* b, integer* ldb,
                    scomplex* x, integer* ldx, float* ferr, float* berr,
                    scomplex* work, float* rwork, integer *info );

void FFUNC(zporfs,ZPORFS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, const dcomplex* af,
                    integer* ldaf, const dcomplex* b, integer* ldb,
                    dcomplex* x, integer* ldx, double* ferr, double* berr,
                    dcomplex* work, double* rwork, integer *info );

void FFUNC(sporfsx,SPORFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const float* a, integer* lda, const float* af,
                     integer* ldaf, const float* s, const float* b,
                     integer* ldb, float* x, integer* ldx, float* rcond,
                     float* berr, integer* n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer* nparams, float* params, float* work,
                     integer* iwork, integer *info );

void FFUNC(dporfsx,DPORFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const double* a, integer* lda, const double* af,
                     integer* ldaf, const double* s, const double* b,
                     integer* ldb, double* x, integer* ldx, double* rcond,
                     double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, double* work,
                     integer* iwork, integer *info );

void FFUNC(cporfsx,CPORFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const scomplex* a, integer* lda, const scomplex* af,
                     integer* ldaf, const float* s, const scomplex* b,
                     integer* ldb, scomplex* x, integer* ldx, float* rcond,
                     float* berr, integer* n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer* nparams, float* params, scomplex* work,
                     float* rwork, integer *info );

void FFUNC(zporfsx,ZPORFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const dcomplex* a, integer* lda, const dcomplex* af,
                     integer* ldaf, const double* s, const dcomplex* b,
                     integer* ldb, dcomplex* x, integer* ldx, double* rcond,
                     double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, dcomplex* work,
                     double* rwork, integer *info );

void FFUNC(spprfs,SPPRFS)( char* uplo, integer* n, integer* nrhs,
                    const float* ap, const float* afp, const float* b,
                    integer* ldb, float* x, integer* ldx, float* ferr,
                    float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dpprfs,DPPRFS)( char* uplo, integer* n, integer* nrhs,
                    const double* ap, const double* afp, const double* b,
                    integer* ldb, double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(cpprfs,CPPRFS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* ap, const scomplex* afp, const scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, float* ferr,
                    float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zpprfs,ZPPRFS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* ap, const dcomplex* afp, const dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, double* ferr,
                    double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(spbrfs,SPBRFS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const float* ab, integer* ldab, const float* afb,
                    integer* ldafb, const float* b, integer* ldb,
                    float* x, integer* ldx, float* ferr, float* berr,
                    float* work, integer* iwork, integer *info );

void FFUNC(dpbrfs,DPBRFS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const double* ab, integer* ldab, const double* afb,
                    integer* ldafb, const double* b, integer* ldb,
                    double* x, integer* ldx, double* ferr, double* berr,
                    double* work, integer* iwork, integer *info );

void FFUNC(cpbrfs,CPBRFS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const scomplex* ab, integer* ldab, const scomplex* afb,
                    integer* ldafb, const scomplex* b, integer* ldb,
                    scomplex* x, integer* ldx, float* ferr, float* berr,
                    scomplex* work, float* rwork, integer *info );

void FFUNC(zpbrfs,ZPBRFS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const dcomplex* ab, integer* ldab, const dcomplex* afb,
                    integer* ldafb, const dcomplex* b, integer* ldb,
                    dcomplex* x, integer* ldx, double* ferr, double* berr,
                    dcomplex* work, double* rwork, integer *info );

void FFUNC(sptrfs,SPTRFS)( integer* n, integer* nrhs, const float* d,
                    const float* e, const float* df, const float* ef,
                    const float* b, integer* ldb, float* x,
                    integer* ldx, float* ferr, float* berr, float* work,
                    integer *info );

void FFUNC(dptrfs,DPTRFS)( integer* n, integer* nrhs, const double* d,
                    const double* e, const double* df, const double* ef,
                    const double* b, integer* ldb, double* x,
                    integer* ldx, double* ferr, double* berr, double* work,
                    integer *info );

void FFUNC(cptrfs,CPTRFS)( integer* n, integer* nrhs, const float* d,
                    const scomplex* e, const float* df, const scomplex* ef,
                    const scomplex* b, integer* ldb, scomplex* x,
                    integer* ldx, float* ferr, float* berr, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(zptrfs,ZPTRFS)( integer* n, integer* nrhs, const double* d,
                    const dcomplex* e, const double* df, const dcomplex* ef,
                    const dcomplex* b, integer* ldb, dcomplex* x,
                    integer* ldx, double* ferr, double* berr, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(ssyrfs,SSYRFS)( char* uplo, integer* n, integer* nrhs,
                    const float* a, integer* lda, const float* af,
                    integer* ldaf, const integer* ipiv, const float* b,
                    integer* ldb, float* x, integer* ldx, float* ferr,
                    float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dsyrfs,DSYRFS)( char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, const double* af,
                    integer* ldaf, const integer* ipiv, const double* b,
                    integer* ldb, double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(csyrfs,CSYRFS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, const scomplex* af,
                    integer* ldaf, const integer* ipiv, const scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, float* ferr,
                    float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zsyrfs,ZSYRFS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, const dcomplex* af,
                    integer* ldaf, const integer* ipiv, const dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, double* ferr,
                    double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(cherfs,CHERFS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, const scomplex* af,
                    integer* ldaf, const integer* ipiv, const scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, float* ferr,
                    float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zherfs,ZHERFS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, const dcomplex* af,
                    integer* ldaf, const integer* ipiv, const dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, double* ferr,
                    double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(ssyrfsx,SSYRFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const float* a, integer* lda, const float* af,
                     integer* ldaf, const integer* ipiv, const float* s,
                     const float* b, integer* ldb, float* x,
                     integer* ldx, float* rcond, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     float* work, integer* iwork, integer *info );

void FFUNC(dsyrfsx,DSYRFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const double* a, integer* lda, const double* af,
                     integer* ldaf, const integer* ipiv, const double* s,
                     const double* b, integer* ldb, double* x,
                     integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(csyrfsx,CSYRFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const scomplex* a, integer* lda, const scomplex* af,
                     integer* ldaf, const integer* ipiv, const float* s,
                     const scomplex* b, integer* ldb, scomplex* x,
                     integer* ldx, float* rcond, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     scomplex* work, float* rwork, integer *info );

void FFUNC(zsyrfsx,ZSYRFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const dcomplex* a, integer* lda, const dcomplex* af,
                     integer* ldaf, const integer* ipiv, const double* s,
                     const dcomplex* b, integer* ldb, dcomplex* x,
                     integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     dcomplex* work, double* rwork, integer *info );

void FFUNC(cherfsx,CHERFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const scomplex* a, integer* lda, const scomplex* af,
                     integer* ldaf, const integer* ipiv, const float* s,
                     const scomplex* b, integer* ldb, scomplex* x,
                     integer* ldx, float* rcond, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     scomplex* work, float* rwork, integer *info );

void FFUNC(zherfsx,ZHERFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const dcomplex* a, integer* lda, const dcomplex* af,
                     integer* ldaf, const integer* ipiv, const double* s,
                     const dcomplex* b, integer* ldb, dcomplex* x,
                     integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     dcomplex* work, double* rwork, integer *info );

void FFUNC(ssprfs,SSPRFS)( char* uplo, integer* n, integer* nrhs,
                    const float* ap, const float* afp, const integer* ipiv,
                    const float* b, integer* ldb, float* x,
                    integer* ldx, float* ferr, float* berr, float* work,
                    integer* iwork, integer *info );

void FFUNC(dsprfs,DSPRFS)( char* uplo, integer* n, integer* nrhs,
                    const double* ap, const double* afp, const integer* ipiv,
                    const double* b, integer* ldb, double* x,
                    integer* ldx, double* ferr, double* berr, double* work,
                    integer* iwork, integer *info );

void FFUNC(csprfs,CSPRFS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* ap, const scomplex* afp, const integer* ipiv,
                    const scomplex* b, integer* ldb, scomplex* x,
                    integer* ldx, float* ferr, float* berr, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(zsprfs,ZSPRFS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* ap, const dcomplex* afp, const integer* ipiv,
                    const dcomplex* b, integer* ldb, dcomplex* x,
                    integer* ldx, double* ferr, double* berr, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(chprfs,CHPRFS)( char* uplo, integer* n, integer* nrhs,
                    const scomplex* ap, const scomplex* afp, const integer* ipiv,
                    const scomplex* b, integer* ldb, scomplex* x,
                    integer* ldx, float* ferr, float* berr, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(zhprfs,ZHPRFS)( char* uplo, integer* n, integer* nrhs,
                    const dcomplex* ap, const dcomplex* afp, const integer* ipiv,
                    const dcomplex* b, integer* ldb, dcomplex* x,
                    integer* ldx, double* ferr, double* berr, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(strrfs,STRRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const float* a, integer* lda,
                    const float* b, integer* ldb, const float* x,
                    integer* ldx, float* ferr, float* berr, float* work,
                    integer* iwork, integer *info );

void FFUNC(dtrrfs,DTRRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const double* a, integer* lda,
                    const double* b, integer* ldb, const double* x,
                    integer* ldx, double* ferr, double* berr, double* work,
                    integer* iwork, integer *info );

void FFUNC(ctrrfs,CTRRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const scomplex* a, integer* lda,
                    const scomplex* b, integer* ldb, const scomplex* x,
                    integer* ldx, float* ferr, float* berr, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(ztrrfs,ZTRRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const dcomplex* a, integer* lda,
                    const dcomplex* b, integer* ldb, const dcomplex* x,
                    integer* ldx, double* ferr, double* berr, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(stprfs,STPRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const float* ap, const float* b,
                    integer* ldb, const float* x, integer* ldx,
                    float* ferr, float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dtprfs,DTPRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const double* ap, const double* b,
                    integer* ldb, const double* x, integer* ldx,
                    double* ferr, double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(ctprfs,CTPRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const scomplex* ap, const scomplex* b,
                    integer* ldb, const scomplex* x, integer* ldx,
                    float* ferr, float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(ztprfs,ZTPRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const dcomplex* ap, const dcomplex* b,
                    integer* ldb, const dcomplex* x, integer* ldx,
                    double* ferr, double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(stbrfs,STBRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const float* ab,
                    integer* ldab, const float* b, integer* ldb,
                    const float* x, integer* ldx, float* ferr,
                    float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dtbrfs,DTBRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const double* ab,
                    integer* ldab, const double* b, integer* ldb,
                    const double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(ctbrfs,CTBRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const scomplex* ab,
                    integer* ldab, const scomplex* b, integer* ldb,
                    const scomplex* x, integer* ldx, float* ferr,
                    float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(ztbrfs,ZTBRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const dcomplex* ab,
                    integer* ldab, const dcomplex* b, integer* ldb,
                    const dcomplex* x, integer* ldx, double* ferr,
                    double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(sgetri,SGETRI)( integer* n, float* a, integer* lda,
                    const integer* ipiv, float* work, integer* lwork,
                    integer *info );

void FFUNC(dgetri,DGETRI)( integer* n, double* a, integer* lda,
                    const integer* ipiv, double* work, integer* lwork,
                    integer *info );

void FFUNC(cgetri,CGETRI)( integer* n, scomplex* a, integer* lda,
                    const integer* ipiv, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zgetri,ZGETRI)( integer* n, dcomplex* a, integer* lda,
                    const integer* ipiv, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(spotri,SPOTRI)( char* uplo, integer* n, float* a, integer* lda,
                    integer *info );

void FFUNC(dpotri,DPOTRI)( char* uplo, integer* n, double* a, integer* lda,
                    integer *info );

void FFUNC(cpotri,CPOTRI)( char* uplo, integer* n, scomplex* a, integer* lda,
                    integer *info );

void FFUNC(zpotri,ZPOTRI)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    integer *info );

void FFUNC(spftri,SPFTRI)( char* transr, char* uplo, integer* n, float* a,
                    integer *info );

void FFUNC(dpftri,DPFTRI)( char* transr, char* uplo, integer* n, double* a,
                    integer *info );

void FFUNC(cpftri,CPFTRI)( char* transr, char* uplo, integer* n, scomplex* a,
                    integer *info );

void FFUNC(zpftri,ZPFTRI)( char* transr, char* uplo, integer* n, dcomplex* a,
                    integer *info );

void FFUNC(spptri,SPPTRI)( char* uplo, integer* n, float* ap, integer *info );

void FFUNC(dpptri,DPPTRI)( char* uplo, integer* n, double* ap, integer *info );

void FFUNC(cpptri,CPPTRI)( char* uplo, integer* n, scomplex* ap, integer *info );

void FFUNC(zpptri,ZPPTRI)( char* uplo, integer* n, dcomplex* ap, integer *info );

void FFUNC(ssytri,SSYTRI)( char* uplo, integer* n, float* a, integer* lda,
                    const integer* ipiv, float* work, integer *info );

void FFUNC(dsytri,DSYTRI)( char* uplo, integer* n, double* a, integer* lda,
                    const integer* ipiv, double* work, integer *info );

void FFUNC(csytri,CSYTRI)( char* uplo, integer* n, scomplex* a, integer* lda,
                    const integer* ipiv, scomplex* work, integer *info );

void FFUNC(zsytri,ZSYTRI)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    const integer* ipiv, dcomplex* work, integer *info );

void FFUNC(chetri,CHETRI)( char* uplo, integer* n, scomplex* a, integer* lda,
                    const integer* ipiv, scomplex* work, integer *info );

void FFUNC(zhetri,ZHETRI)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    const integer* ipiv, dcomplex* work, integer *info );

void FFUNC(ssptri,SSPTRI)( char* uplo, integer* n, float* ap,
                    const integer* ipiv, float* work, integer *info );

void FFUNC(dsptri,DSPTRI)( char* uplo, integer* n, double* ap,
                    const integer* ipiv, double* work, integer *info );

void FFUNC(csptri,CSPTRI)( char* uplo, integer* n, scomplex* ap,
                    const integer* ipiv, scomplex* work, integer *info );

void FFUNC(zsptri,ZSPTRI)( char* uplo, integer* n, dcomplex* ap,
                    const integer* ipiv, dcomplex* work, integer *info );

void FFUNC(chptri,CHPTRI)( char* uplo, integer* n, scomplex* ap,
                    const integer* ipiv, scomplex* work, integer *info );

void FFUNC(zhptri,ZHPTRI)( char* uplo, integer* n, dcomplex* ap,
                    const integer* ipiv, dcomplex* work, integer *info );

void FFUNC(strtri,STRTRI)( char* uplo, char* diag, integer* n, float* a,
                    integer* lda, integer *info );

void FFUNC(dtrtri,DTRTRI)( char* uplo, char* diag, integer* n, double* a,
                    integer* lda, integer *info );

void FFUNC(ctrtri,CTRTRI)( char* uplo, char* diag, integer* n, scomplex* a,
                    integer* lda, integer *info );

void FFUNC(ztrtri,ZTRTRI)( char* uplo, char* diag, integer* n, dcomplex* a,
                    integer* lda, integer *info );

void FFUNC(stftri,STFTRI)( char* transr, char* uplo, char* diag, integer* n,
                    float* a, integer *info );

void FFUNC(dtftri,DTFTRI)( char* transr, char* uplo, char* diag, integer* n,
                    double* a, integer *info );

void FFUNC(ctftri,CTFTRI)( char* transr, char* uplo, char* diag, integer* n,
                    scomplex* a, integer *info );

void FFUNC(ztftri,ZTFTRI)( char* transr, char* uplo, char* diag, integer* n,
                    dcomplex* a, integer *info );

void FFUNC(stptri,STPTRI)( char* uplo, char* diag, integer* n, float* ap,
                    integer *info );

void FFUNC(dtptri,DTPTRI)( char* uplo, char* diag, integer* n, double* ap,
                    integer *info );

void FFUNC(ctptri,CTPTRI)( char* uplo, char* diag, integer* n, scomplex* ap,
                    integer *info );

void FFUNC(ztptri,ZTPTRI)( char* uplo, char* diag, integer* n, dcomplex* ap,
                    integer *info );

void FFUNC(sgeequ,SGEEQU)( integer* m, integer* n, const float* a,
                    integer* lda, float* r, float* c, float* rowcnd,
                    float* colcnd, float* amax, integer *info );

void FFUNC(dgeequ,DGEEQU)( integer* m, integer* n, const double* a,
                    integer* lda, double* r, double* c, double* rowcnd,
                    double* colcnd, double* amax, integer *info );

void FFUNC(cgeequ,CGEEQU)( integer* m, integer* n, const scomplex* a,
                    integer* lda, float* r, float* c, float* rowcnd,
                    float* colcnd, float* amax, integer *info );

void FFUNC(zgeequ,ZGEEQU)( integer* m, integer* n, const dcomplex* a,
                    integer* lda, double* r, double* c, double* rowcnd,
                    double* colcnd, double* amax, integer *info );

void FFUNC(sgeequb,SGEEQUB)( integer* m, integer* n, const float* a,
                    integer* lda, float* r, float* c, float* rowcnd,
                    float* colcnd, float* amax, integer *info );

void FFUNC(dgeequb,DGEEQUB)( integer* m, integer* n, const double* a,
                    integer* lda, double* r, double* c, double* rowcnd,
                    double* colcnd, double* amax, integer *info );

void FFUNC(cgeequb,CGEEQUB)( integer* m, integer* n, const scomplex* a,
                    integer* lda, float* r, float* c, float* rowcnd,
                    float* colcnd, float* amax, integer *info );

void FFUNC(zgeequb,ZGEEQUB)( integer* m, integer* n, const dcomplex* a,
                    integer* lda, double* r, double* c, double* rowcnd,
                    double* colcnd, double* amax, integer *info );

void FFUNC(sgbequ,SGBEQU)( integer* m, integer* n, integer* kl,
                    integer* ku, const float* ab, integer* ldab,
                    float* r, float* c, float* rowcnd, float* colcnd,
                    float* amax, integer *info );

void FFUNC(dgbequ,DGBEQU)( integer* m, integer* n, integer* kl,
                    integer* ku, const double* ab, integer* ldab,
                    double* r, double* c, double* rowcnd, double* colcnd,
                    double* amax, integer *info );

void FFUNC(cgbequ,CGBEQU)( integer* m, integer* n, integer* kl,
                    integer* ku, const scomplex* ab, integer* ldab,
                    float* r, float* c, float* rowcnd, float* colcnd,
                    float* amax, integer *info );

void FFUNC(zgbequ,ZGBEQU)( integer* m, integer* n, integer* kl,
                    integer* ku, const dcomplex* ab, integer* ldab,
                    double* r, double* c, double* rowcnd, double* colcnd,
                    double* amax, integer *info );

void FFUNC(sgbequb,SGBEQUB)( integer* m, integer* n, integer* kl,
                    integer* ku, const float* ab, integer* ldab,
                    float* r, float* c, float* rowcnd, float* colcnd,
                    float* amax, integer *info );

void FFUNC(dgbequb,DGBEQUB)( integer* m, integer* n, integer* kl,
                    integer* ku, const double* ab, integer* ldab,
                    double* r, double* c, double* rowcnd, double* colcnd,
                    double* amax, integer *info );

void FFUNC(cgbequb,CGBEQUB)( integer* m, integer* n, integer* kl,
                    integer* ku, const scomplex* ab, integer* ldab,
                    float* r, float* c, float* rowcnd, float* colcnd,
                    float* amax, integer *info );

void FFUNC(zgbequb,ZGBEQUB)( integer* m, integer* n, integer* kl,
                    integer* ku, const dcomplex* ab, integer* ldab,
                    double* r, double* c, double* rowcnd, double* colcnd,
                    double* amax, integer *info );

void FFUNC(spoequ,SPOEQU)( integer* n, const float* a, integer* lda, float* s,
                    float* scond, float* amax, integer *info );

void FFUNC(dpoequ,DPOEQU)( integer* n, const double* a, integer* lda, double* s,
                    double* scond, double* amax, integer *info );

void FFUNC(cpoequ,CPOEQU)( integer* n, const scomplex* a, integer* lda, float* s,
                    float* scond, float* amax, integer *info );

void FFUNC(zpoequ,ZPOEQU)( integer* n, const dcomplex* a, integer* lda, double* s,
                    double* scond, double* amax, integer *info );

void FFUNC(spoequb,SPOEQUB)( integer* n, const float* a, integer* lda, float* s,
                    float* scond, float* amax, integer *info );

void FFUNC(dpoequb,DPOEQUB)( integer* n, const double* a, integer* lda, double* s,
                    double* scond, double* amax, integer *info );

void FFUNC(cpoequb,CPOEQUB)( integer* n, const scomplex* a, integer* lda, float* s,
                    float* scond, float* amax, integer *info );

void FFUNC(zpoequb,ZPOEQUB)( integer* n, const dcomplex* a, integer* lda, double* s,
                    double* scond, double* amax, integer *info );

void FFUNC(sppequ,SPPEQU)( char* uplo, integer* n, const float* ap, float* s,
                    float* scond, float* amax, integer *info );

void FFUNC(dppequ,DPPEQU)( char* uplo, integer* n, const double* ap, double* s,
                    double* scond, double* amax, integer *info );

void FFUNC(cppequ,CPPEQU)( char* uplo, integer* n, const scomplex* ap, float* s,
                    float* scond, float* amax, integer *info );

void FFUNC(zppequ,ZPPEQU)( char* uplo, integer* n, const dcomplex* ap, double* s,
                    double* scond, double* amax, integer *info );

void FFUNC(spbequ,SPBEQU)( char* uplo, integer* n, integer* kd, const float* ab,
                    integer* ldab, float* s, float* scond, float* amax,
                    integer *info );

void FFUNC(dpbequ,DPBEQU)( char* uplo, integer* n, integer* kd, const double* ab,
                    integer* ldab, double* s, double* scond, double* amax,
                    integer *info );

void FFUNC(cpbequ,CPBEQU)( char* uplo, integer* n, integer* kd, const scomplex* ab,
                    integer* ldab, float* s, float* scond, float* amax,
                    integer *info );

void FFUNC(zpbequ,ZPBEQU)( char* uplo, integer* n, integer* kd, const dcomplex* ab,
                    integer* ldab, double* s, double* scond, double* amax,
                    integer *info );

void FFUNC(ssyequb,SSYEQUB)( char* uplo, integer* n, const float* a,
                     integer* lda, float* s, float* scond, float* amax,
                     float* work, integer *info );

void FFUNC(dsyequb,DSYEQUB)( char* uplo, integer* n, const double* a,
                     integer* lda, double* s, double* scond, double* amax,
                     double* work, integer *info );

void FFUNC(csyequb,CSYEQUB)( char* uplo, integer* n, const scomplex* a,
                     integer* lda, float* s, float* scond, float* amax,
                     scomplex* work, integer *info );

void FFUNC(zsyequb,ZSYEQUB)( char* uplo, integer* n, const dcomplex* a,
                     integer* lda, double* s, double* scond, double* amax,
                     dcomplex* work, integer *info );

void FFUNC(cheequb,CHEEQUB)( char* uplo, integer* n, const scomplex* a,
                     integer* lda, float* s, float* scond, float* amax,
                     scomplex* work, integer *info );

void FFUNC(zheequb,ZHEEQUB)( char* uplo, integer* n, const dcomplex* a,
                     integer* lda, double* s, double* scond, double* amax,
                     dcomplex* work, integer *info );

void FFUNC(sgesv,SGESV)( integer* n, integer* nrhs, float* a, integer* lda,
                   integer* ipiv, float* b, integer* ldb,
                   integer *info );

void FFUNC(dgesv,DGESV)( integer* n, integer* nrhs, double* a, integer* lda,
                   integer* ipiv, double* b, integer* ldb,
                   integer *info );

void FFUNC(cgesv,CGESV)( integer* n, integer* nrhs, scomplex* a, integer* lda,
                   integer* ipiv, scomplex* b, integer* ldb,
                   integer *info );

void FFUNC(zgesv,ZGESV)( integer* n, integer* nrhs, dcomplex* a, integer* lda,
                   integer* ipiv, dcomplex* b, integer* ldb,
                   integer *info );

void FFUNC(dsgesv,DSGESV)( integer* n, integer* nrhs, double* a, integer* lda,
                    integer* ipiv, double* b, integer* ldb, double* x,
                    integer* ldx, double* work, float* swork,
                    integer* iter, integer *info );

void FFUNC(zcgesv,ZCGESV)( integer* n, integer* nrhs, dcomplex* a, integer* lda,
                    integer* ipiv, dcomplex* b, integer* ldb, dcomplex* x,
                    integer* ldx, dcomplex* work, scomplex* swork, double* rwork,
                    integer* iter, integer *info );

void FFUNC(sgesvx,SGESVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    float* a, integer* lda, float* af, integer* ldaf,
                    integer* ipiv, char* equed, float* r, float* c,
                    float* b, integer* ldb, float* x, integer* ldx,
                    float* rcond, float* ferr, float* berr, float* work,
                    integer* iwork, integer *info );

void FFUNC(dgesvx,DGESVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    double* a, integer* lda, double* af, integer* ldaf,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer* ldb, double* x, integer* ldx,
                    double* rcond, double* ferr, double* berr, double* work,
                    integer* iwork, integer *info );

void FFUNC(cgesvx,CGESVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    scomplex* a, integer* lda, scomplex* af, integer* ldaf,
                    integer* ipiv, char* equed, float* r, float* c,
                    scomplex* b, integer* ldb, scomplex* x, integer* ldx,
                    float* rcond, float* ferr, float* berr, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(zgesvx,ZGESVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    dcomplex* a, integer* lda, dcomplex* af, integer* ldaf,
                    integer* ipiv, char* equed, double* r, double* c,
                    dcomplex* b, integer* ldb, dcomplex* x, integer* ldx,
                    double* rcond, double* ferr, double* berr, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(sgesvxx,SGESVXX)( char* fact, char* trans, integer* n, integer* nrhs,
                     float* a, integer* lda, float* af, integer* ldaf,
                     integer* ipiv, char* equed, float* r, float* c,
                     float* b, integer* ldb, float* x, integer* ldx,
                     float* rcond, float* rpvgrw, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     float* work, integer* iwork, integer *info );

void FFUNC(dgesvxx,DGESVXX)( char* fact, char* trans, integer* n, integer* nrhs,
                     double* a, integer* lda, double* af, integer* ldaf,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer* ldb, double* x, integer* ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(cgesvxx,CGESVXX)( char* fact, char* trans, integer* n, integer* nrhs,
                     scomplex* a, integer* lda, scomplex* af, integer* ldaf,
                     integer* ipiv, char* equed, float* r, float* c,
                     scomplex* b, integer* ldb, scomplex* x, integer* ldx,
                     float* rcond, float* rpvgrw, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     scomplex* work, float* rwork, integer *info );

void FFUNC(zgesvxx,ZGESVXX)( char* fact, char* trans, integer* n, integer* nrhs,
                     dcomplex* a, integer* lda, dcomplex* af, integer* ldaf,
                     integer* ipiv, char* equed, double* r, double* c,
                     dcomplex* b, integer* ldb, dcomplex* x, integer* ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     dcomplex* work, double* rwork, integer *info );

void FFUNC(sgbsv,SGBSV)( integer* n, integer* kl, integer* ku,
                   integer* nrhs, float* ab, integer* ldab,
                   integer* ipiv, float* b, integer* ldb,
                   integer *info );

void FFUNC(dgbsv,DGBSV)( integer* n, integer* kl, integer* ku,
                   integer* nrhs, double* ab, integer* ldab,
                   integer* ipiv, double* b, integer* ldb,
                   integer *info );

void FFUNC(cgbsv,CGBSV)( integer* n, integer* kl, integer* ku,
                   integer* nrhs, scomplex* ab, integer* ldab,
                   integer* ipiv, scomplex* b, integer* ldb,
                   integer *info );

void FFUNC(zgbsv,ZGBSV)( integer* n, integer* kl, integer* ku,
                   integer* nrhs, dcomplex* ab, integer* ldab,
                   integer* ipiv, dcomplex* b, integer* ldb,
                   integer *info );

void FFUNC(sgbsvx,SGBSVX)( char* fact, char* trans, integer* n, integer* kl,
                    integer* ku, integer* nrhs, float* ab,
                    integer* ldab, float* afb, integer* ldafb,
                    integer* ipiv, char* equed, float* r, float* c,
                    float* b, integer* ldb, float* x, integer* ldx,
                    float* rcond, float* ferr, float* berr, float* work,
                    integer* iwork, integer *info );

void FFUNC(dgbsvx,DGBSVX)( char* fact, char* trans, integer* n, integer* kl,
                    integer* ku, integer* nrhs, double* ab,
                    integer* ldab, double* afb, integer* ldafb,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer* ldb, double* x, integer* ldx,
                    double* rcond, double* ferr, double* berr, double* work,
                    integer* iwork, integer *info );

void FFUNC(cgbsvx,CGBSVX)( char* fact, char* trans, integer* n, integer* kl,
                    integer* ku, integer* nrhs, scomplex* ab,
                    integer* ldab, scomplex* afb, integer* ldafb,
                    integer* ipiv, char* equed, float* r, float* c,
                    scomplex* b, integer* ldb, scomplex* x, integer* ldx,
                    float* rcond, float* ferr, float* berr, scomplex* work,
                    float* rwork, integer *info );

void FFUNC(zgbsvx,ZGBSVX)( char* fact, char* trans, integer* n, integer* kl,
                    integer* ku, integer* nrhs, dcomplex* ab,
                    integer* ldab, dcomplex* afb, integer* ldafb,
                    integer* ipiv, char* equed, double* r, double* c,
                    dcomplex* b, integer* ldb, dcomplex* x, integer* ldx,
                    double* rcond, double* ferr, double* berr, dcomplex* work,
                    double* rwork, integer *info );

void FFUNC(sgbsvxx,SGBSVXX)( char* fact, char* trans, integer* n, integer* kl,
                     integer* ku, integer* nrhs, float* ab,
                     integer* ldab, float* afb, integer* ldafb,
                     integer* ipiv, char* equed, float* r, float* c,
                     float* b, integer* ldb, float* x, integer* ldx,
                     float* rcond, float* rpvgrw, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     float* work, integer* iwork, integer *info );

void FFUNC(dgbsvxx,DGBSVXX)( char* fact, char* trans, integer* n, integer* kl,
                     integer* ku, integer* nrhs, double* ab,
                     integer* ldab, double* afb, integer* ldafb,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer* ldb, double* x, integer* ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(cgbsvxx,CGBSVXX)( char* fact, char* trans, integer* n, integer* kl,
                     integer* ku, integer* nrhs, scomplex* ab,
                     integer* ldab, scomplex* afb, integer* ldafb,
                     integer* ipiv, char* equed, float* r, float* c,
                     scomplex* b, integer* ldb, scomplex* x, integer* ldx,
                     float* rcond, float* rpvgrw, float* berr,
                     integer* n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer* nparams, float* params,
                     scomplex* work, float* rwork, integer *info );

void FFUNC(zgbsvxx,ZGBSVXX)( char* fact, char* trans, integer* n, integer* kl,
                     integer* ku, integer* nrhs, dcomplex* ab,
                     integer* ldab, dcomplex* afb, integer* ldafb,
                     integer* ipiv, char* equed, double* r, double* c,
                     dcomplex* b, integer* ldb, dcomplex* x, integer* ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     dcomplex* work, double* rwork, integer *info );

void FFUNC(sgtsv,SGTSV)( integer* n, integer* nrhs, float* dl, float* d,
                   float* du, float* b, integer* ldb, integer *info );

void FFUNC(dgtsv,DGTSV)( integer* n, integer* nrhs, double* dl, double* d,
                   double* du, double* b, integer* ldb, integer *info );

void FFUNC(cgtsv,CGTSV)( integer* n, integer* nrhs, scomplex* dl, scomplex* d,
                    scomplex* du, scomplex* b, integer* ldb, integer *info );

void FFUNC(zgtsv,ZGTSV)( integer* n, integer* nrhs, dcomplex* dl, dcomplex* d,
                    dcomplex* du, dcomplex* b, integer* ldb, integer *info );

void FFUNC(sgtsvx,SGTSVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    const float* dl, const float* d, const float* du,
                    float* dlf, float* df, float* duf, float* du2,
                    integer* ipiv, const float* b, integer* ldb,
                    float* x, integer* ldx, float* rcond, float* ferr,
                    float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dgtsvx,DGTSVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    const double* dl, const double* d, const double* du,
                    double* dlf, double* df, double* duf, double* du2,
                    integer* ipiv, const double* b, integer* ldb,
                    double* x, integer* ldx, double* rcond, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(cgtsvx,CGTSVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    const scomplex* dl, const scomplex* d, const scomplex* du,
                    scomplex* dlf, scomplex* df, scomplex* duf, scomplex* du2,
                    integer* ipiv, const scomplex* b, integer* ldb,
                    scomplex* x, integer* ldx, float* rcond, float* ferr,
                    float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zgtsvx,ZGTSVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    const dcomplex* dl, const dcomplex* d, const dcomplex* du,
                    dcomplex* dlf, dcomplex* df, dcomplex* duf, dcomplex* du2,
                    integer* ipiv, const dcomplex* b, integer* ldb,
                    dcomplex* x, integer* ldx, double* rcond, double* ferr,
                    double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(sposv,SPOSV)( char* uplo, integer* n, integer* nrhs, float* a,
                   integer* lda, float* b, integer* ldb,
                   integer *info );

void FFUNC(dposv,DPOSV)( char* uplo, integer* n, integer* nrhs, double* a,
                   integer* lda, double* b, integer* ldb,
                   integer *info );

void FFUNC(cposv,CPOSV)( char* uplo, integer* n, integer* nrhs, scomplex* a,
                   integer* lda, scomplex* b, integer* ldb,
                   integer *info );

void FFUNC(zposv,ZPOSV)( char* uplo, integer* n, integer* nrhs, dcomplex* a,
                   integer* lda, dcomplex* b, integer* ldb,
                   integer *info );

void FFUNC(dsposv,DSPOSV)( char* uplo, integer* n, integer* nrhs, double* a,
                    integer* lda, double* b, integer* ldb, double* x,
                    integer* ldx, double* work, float* swork,
                    integer* iter, integer *info );

void FFUNC(zcposv,ZCPOSV)( char* uplo, integer* n, integer* nrhs, dcomplex* a,
                    integer* lda, dcomplex* b, integer* ldb, dcomplex* x,
                    integer* ldx, dcomplex* work, scomplex* swork, double* rwork,
                    integer* iter, integer *info );

void FFUNC(sposvx,SPOSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    float* a, integer* lda, float* af, integer* ldaf,
                    char* equed, float* s, float* b, integer* ldb,
                    float* x, integer* ldx, float* rcond, float* ferr,
                    float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dposvx,DPOSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    double* a, integer* lda, double* af, integer* ldaf,
                    char* equed, double* s, double* b, integer* ldb,
                    double* x, integer* ldx, double* rcond, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(cposvx,CPOSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    scomplex* a, integer* lda, scomplex* af, integer* ldaf,
                    char* equed, float* s, scomplex* b, integer* ldb,
                    scomplex* x, integer* ldx, float* rcond, float* ferr,
                    float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zposvx,ZPOSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    dcomplex* a, integer* lda, dcomplex* af, integer* ldaf,
                    char* equed, double* s, dcomplex* b, integer* ldb,
                    dcomplex* x, integer* ldx, double* rcond, double* ferr,
                    double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(sposvxx,SPOSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     float* a, integer* lda, float* af, integer* ldaf,
                     char* equed, float* s, float* b, integer* ldb,
                     float* x, integer* ldx, float* rcond, float* rpvgrw,
                     float* berr, integer* n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer* nparams, float* params, float* work,
                     integer* iwork, integer *info );

void FFUNC(dposvxx,DPOSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     double* a, integer* lda, double* af, integer* ldaf,
                     char* equed, double* s, double* b, integer* ldb,
                     double* x, integer* ldx, double* rcond, double* rpvgrw,
                     double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, double* work,
                     integer* iwork, integer *info );

void FFUNC(cposvxx,CPOSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    scomplex* a, integer* lda, scomplex* af, integer* ldaf,
                     char* equed, float* s, scomplex* b, integer* ldb,
                     scomplex* x, integer* ldx, float* rcond, float* rpvgrw,
                     float* berr, integer* n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer* nparams, float* params, scomplex* work,
                     float* rwork, integer *info );

void FFUNC(zposvxx,ZPOSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     dcomplex* a, integer* lda, dcomplex* af, integer* ldaf,
                     char* equed, double* s, dcomplex* b, integer* ldb,
                     dcomplex* x, integer* ldx, double* rcond, double* rpvgrw,
                     double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, dcomplex* work,
                     double* rwork, integer *info );

void FFUNC(sppsv,SPPSV)( char* uplo, integer* n, integer* nrhs, float* ap,
                   float* b, integer* ldb, integer *info );

void FFUNC(dppsv,DPPSV)( char* uplo, integer* n, integer* nrhs, double* ap,
                   double* b, integer* ldb, integer *info );

void FFUNC(cppsv,CPPSV)( char* uplo, integer* n, integer* nrhs, scomplex* ap,
                   scomplex* b, integer* ldb, integer *info );

void FFUNC(zppsv,ZPPSV)( char* uplo, integer* n, integer* nrhs, dcomplex* ap,
                   dcomplex* b, integer* ldb, integer *info );

void FFUNC(sppsvx,SPPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    float* ap, float* afp, char* equed, float* s, float* b,
                    integer* ldb, float* x, integer* ldx, float* rcond,
                    float* ferr, float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dppsvx,DPPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    double* ap, double* afp, char* equed, double* s, double* b,
                    integer* ldb, double* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(cppsvx,CPPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    scomplex* ap, scomplex* afp, char* equed, float* s, scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, float* rcond,
                    float* ferr, float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zppsvx,ZPPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    dcomplex* ap, dcomplex* afp, char* equed, double* s, dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(spbsv,SPBSV)( char* uplo, integer* n, integer* kd, integer* nrhs,
                   float* ab, integer* ldab, float* b, integer* ldb,
                   integer *info );

void FFUNC(dpbsv,DPBSV)( char* uplo, integer* n, integer* kd, integer* nrhs,
                   double* ab, integer* ldab, double* b, integer* ldb,
                   integer *info );

void FFUNC(cpbsv,CPBSV)( char* uplo, integer* n, integer* kd, integer* nrhs,
                   scomplex* ab, integer* ldab, scomplex* b, integer* ldb,
                   integer *info );

void FFUNC(zpbsv,ZPBSV)( char* uplo, integer* n, integer* kd, integer* nrhs,
                   dcomplex* ab, integer* ldab, dcomplex* b, integer* ldb,
                   integer *info );

void FFUNC(spbsvx,SPBSVX)( char* fact, char* uplo, integer* n, integer* kd,
                    integer* nrhs, float* ab, integer* ldab, float* afb,
                    integer* ldafb, char* equed, float* s, float* b,
                    integer* ldb, float* x, integer* ldx, float* rcond,
                    float* ferr, float* berr, float* work, integer* iwork,
                    integer *info );

void FFUNC(dpbsvx,DPBSVX)( char* fact, char* uplo, integer* n, integer* kd,
                    integer* nrhs, double* ab, integer* ldab, double* afb,
                    integer* ldafb, char* equed, double* s, double* b,
                    integer* ldb, double* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(cpbsvx,CPBSVX)( char* fact, char* uplo, integer* n, integer* kd,
                    integer* nrhs, scomplex* ab, integer* ldab, scomplex* afb,
                    integer* ldafb, char* equed, float* s, scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, float* rcond,
                    float* ferr, float* berr, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zpbsvx,ZPBSVX)( char* fact, char* uplo, integer* n, integer* kd,
                    integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* afb,
                    integer* ldafb, char* equed, double* s, dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(sptsv,SPTSV)( integer* n, integer* nrhs, float* d, float* e,
                    float* b, integer* ldb, integer *info );

void FFUNC(dptsv,DPTSV)( integer* n, integer* nrhs, double* d, double* e,
                   double* b, integer* ldb, integer *info );

void FFUNC(cptsv,CPTSV)( integer* n, integer* nrhs, float* d, scomplex* e,
                   scomplex* b, integer* ldb, integer *info );

void FFUNC(zptsv,ZPTSV)( integer* n, integer* nrhs, double* d, dcomplex* e,
                   dcomplex* b, integer* ldb, integer *info );

void FFUNC(sptsvx,SPTSVX)( char* fact, integer* n, integer* nrhs,
                    const float* d, const float* e, float* df, float* ef,
                    const float* b, integer* ldb, float* x,
                    integer* ldx, float* rcond, float* ferr, float* berr,
                    float* work, integer *info );

void FFUNC(dptsvx,DPTSVX)( char* fact, integer* n, integer* nrhs,
                    const double* d, const double* e, double* df, double* ef,
                    const double* b, integer* ldb, double* x,
                    integer* ldx, double* rcond, double* ferr, double* berr,
                    double* work, integer *info );

void FFUNC(cptsvx,CPTSVX)( char* fact, integer* n, integer* nrhs,
                    const float* d, const scomplex* e, float* df, scomplex* ef,
                    const scomplex* b, integer* ldb, scomplex* x,
                    integer* ldx, float* rcond, float* ferr, float* berr,
                    scomplex* work, float* rwork, integer *info );

void FFUNC(zptsvx,ZPTSVX)( char* fact, integer* n, integer* nrhs,
                    const double* d, const dcomplex* e, double* df, dcomplex* ef,
                    const dcomplex* b, integer* ldb, dcomplex* x,
                    integer* ldx, double* rcond, double* ferr, double* berr,
                    dcomplex* work, double* rwork, integer *info );

void FFUNC(ssysv,SSYSV)( char* uplo, integer* n, integer* nrhs, float* a,
                   integer* lda, integer* ipiv, float* b,
                   integer* ldb, float* work, integer* lwork,
                   integer *info );

void FFUNC(dsysv,DSYSV)( char* uplo, integer* n, integer* nrhs, double* a,
                   integer* lda, integer* ipiv, double* b,
                   integer* ldb, double* work, integer* lwork,
                   integer *info );

void FFUNC(csysv,CSYSV)( char* uplo, integer* n, integer* nrhs, scomplex* a,
                   integer* lda, integer* ipiv, scomplex* b,
                   integer* ldb, scomplex* work, integer* lwork,
                   integer *info );

void FFUNC(zsysv,ZSYSV)( char* uplo, integer* n, integer* nrhs, dcomplex* a,
                   integer* lda, integer* ipiv, dcomplex* b,
                   integer* ldb, dcomplex* work, integer* lwork,
                   integer *info );

void FFUNC(chesv,CHESV)( char* uplo, integer* n, integer* nrhs, scomplex* a,
                   integer* lda, integer* ipiv, scomplex* b,
                   integer* ldb, scomplex* work, integer* lwork,
                   integer *info );

void FFUNC(zhesv,ZHESV)( char* uplo, integer* n, integer* nrhs, dcomplex* a,
                   integer* lda, integer* ipiv, dcomplex* b,
                   integer* ldb, dcomplex* work, integer* lwork,
                   integer *info );

void FFUNC(ssysvx,SSYSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const float* a, integer* lda, float* af,
                    integer* ldaf, integer* ipiv, const float* b,
                    integer* ldb, float* x, integer* ldx, float* rcond,
                    float* ferr, float* berr, float* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(dsysvx,DSYSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, double* af,
                    integer* ldaf, integer* ipiv, const double* b,
                    integer* ldb, double* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(csysvx,CSYSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, scomplex* af,
                    integer* ldaf, integer* ipiv, const scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, float* rcond,
                    float* ferr, float* berr, scomplex* work, integer* lwork,
                    float* rwork, integer *info );

void FFUNC(zsysvx,ZSYSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, dcomplex* af,
                    integer* ldaf, integer* ipiv, const dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, dcomplex* work, integer* lwork,
                    double* rwork, integer *info );

void FFUNC(chesvx,CHESVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const scomplex* a, integer* lda, scomplex* af,
                    integer* ldaf, integer* ipiv, const scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, float* rcond,
                    float* ferr, float* berr, scomplex* work, integer* lwork,
                    float* rwork, integer *info );

void FFUNC(zhesvx,ZHESVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const dcomplex* a, integer* lda, dcomplex* af,
                    integer* ldaf, integer* ipiv, const dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, dcomplex* work, integer* lwork,
                    double* rwork, integer *info );

void FFUNC(ssysvxx,SSYSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    float* a, integer* lda, float* af, integer* ldaf,
                     integer* ipiv, char* equed, float* s, float* b,
                     integer* ldb, float* x, integer* ldx, float* rcond,
                     float* rpvgrw, float* berr, integer* n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer* nparams, float* params, float* work,
                     integer* iwork, integer *info );

void FFUNC(dsysvxx,DSYSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     double* a, integer* lda, double* af, integer* ldaf,
                     integer* ipiv, char* equed, double* s, double* b,
                     integer* ldb, double* x, integer* ldx, double* rcond,
                     double* rpvgrw, double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, double* work,
                     integer* iwork, integer *info );

void FFUNC(csysvxx,CSYSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     scomplex* a, integer* lda, scomplex* af, integer* ldaf,
                     integer* ipiv, char* equed, float* s, scomplex* b,
                     integer* ldb, scomplex* x, integer* ldx, float* rcond,
                     float* rpvgrw, float* berr, integer* n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer* nparams, float* params, scomplex* work,
                     float* rwork, integer *info );

void FFUNC(zsysvxx,ZSYSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     dcomplex* a, integer* lda, dcomplex* af, integer* ldaf,
                     integer* ipiv, char* equed, double* s, dcomplex* b,
                     integer* ldb, dcomplex* x, integer* ldx, double* rcond,
                     double* rpvgrw, double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, dcomplex* work,
                     double* rwork, integer *info );

void FFUNC(chesvxx,CHESVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     scomplex* a, integer* lda, scomplex* af, integer* ldaf,
                     integer* ipiv, char* equed, float* s, scomplex* b,
                     integer* ldb, scomplex* x, integer* ldx, float* rcond,
                     float* rpvgrw, float* berr, integer* n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer* nparams, float* params, scomplex* work,
                     float* rwork, integer *info );

void FFUNC(zhesvxx,ZHESVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     dcomplex* a, integer* lda, dcomplex* af, integer* ldaf,
                     integer* ipiv, char* equed, double* s, dcomplex* b,
                     integer* ldb, dcomplex* x, integer* ldx, double* rcond,
                     double* rpvgrw, double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, dcomplex* work,
                     double* rwork, integer *info );

void FFUNC(sspsv,SSPSV)( char* uplo, integer* n, integer* nrhs, float* ap,
                   integer* ipiv, float* b, integer* ldb,
                   integer *info );

void FFUNC(dspsv,DSPSV)( char* uplo, integer* n, integer* nrhs, double* ap,
                   integer* ipiv, double* b, integer* ldb,
                   integer *info );

void FFUNC(cspsv,CSPSV)( char* uplo, integer* n, integer* nrhs, scomplex* ap,
                   integer* ipiv, scomplex* b, integer* ldb,
                   integer *info );

void FFUNC(zspsv,ZSPSV)( char* uplo, integer* n, integer* nrhs, dcomplex* ap,
                   integer* ipiv, dcomplex* b, integer* ldb,
                   integer *info );

void FFUNC(chpsv,CHPSV)( char* uplo, integer* n, integer* nrhs, scomplex* ap,
                   integer* ipiv, scomplex* b, integer* ldb,
                   integer *info );

void FFUNC(zhpsv,ZHPSV)( char* uplo, integer* n, integer* nrhs, dcomplex* ap,
                   integer* ipiv, dcomplex* b, integer* ldb,
                   integer *info );

void FFUNC(sspsvx,SSPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const float* ap, float* afp, integer* ipiv,
                    const float* b, integer* ldb, float* x,
                    integer* ldx, float* rcond, float* ferr, float* berr,
                    float* work, integer* iwork, integer *info );

void FFUNC(dspsvx,DSPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const double* ap, double* afp, integer* ipiv,
                    const double* b, integer* ldb, double* x,
                    integer* ldx, double* rcond, double* ferr, double* berr,
                    double* work, integer* iwork, integer *info );

void FFUNC(cspsvx,CSPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const scomplex* ap, scomplex* afp, integer* ipiv,
                    const scomplex* b, integer* ldb, scomplex* x,
                    integer* ldx, float* rcond, float* ferr, float* berr,
                    scomplex* work, float* rwork, integer *info );

void FFUNC(zspsvx,ZSPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const dcomplex* ap, dcomplex* afp, integer* ipiv,
                    const dcomplex* b, integer* ldb, dcomplex* x,
                    integer* ldx, double* rcond, double* ferr, double* berr,
                    dcomplex* work, double* rwork, integer *info );

void FFUNC(chpsvx,CHPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const scomplex* ap, scomplex* afp, integer* ipiv,
                    const scomplex* b, integer* ldb, scomplex* x,
                    integer* ldx, float* rcond, float* ferr, float* berr,
                    scomplex* work, float* rwork, integer *info );

void FFUNC(zhpsvx,ZHPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const dcomplex* ap, dcomplex* afp, integer* ipiv,
                    const dcomplex* b, integer* ldb, dcomplex* x,
                    integer* ldx, double* rcond, double* ferr, double* berr,
                    dcomplex* work, double* rwork, integer *info );

void FFUNC(sgeqrf,SGEQRF)( integer* m, integer* n, float* a, integer* lda,
                    float* tau, float* work, integer* lwork,
                    integer *info );

void FFUNC(dgeqrf,DGEQRF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(cgeqrf,CGEQRF)( integer* m, integer* n, scomplex* a, integer* lda,
                    scomplex* tau, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zgeqrf,ZGEQRF)( integer* m, integer* n, dcomplex* a, integer* lda,
                    dcomplex* tau, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sgeqpf,SGEQPF)( integer* m, integer* n, float* a, integer* lda,
                    integer* jpvt, float* tau, float* work,
                    integer *info );

void FFUNC(dgeqpf,DGEQPF)( integer* m, integer* n, double* a, integer* lda,
                    integer* jpvt, double* tau, double* work,
                    integer *info );

void FFUNC(cgeqpf,CGEQPF)( integer* m, integer* n, scomplex* a, integer* lda,
                    integer* jpvt, scomplex* tau, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(zgeqpf,ZGEQPF)( integer* m, integer* n, dcomplex* a, integer* lda,
                    integer* jpvt, dcomplex* tau, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(sgeqp3,SGEQP3)( integer* m, integer* n, float* a, integer* lda,
                    integer* jpvt, float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dgeqp3,DGEQP3)( integer* m, integer* n, double* a, integer* lda,
                    integer* jpvt, double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(cgeqp3,CGEQP3)( integer* m, integer* n, scomplex* a, integer* lda,
                    integer* jpvt, scomplex* tau, scomplex* work,
                    integer* lwork, float* rwork, integer *info );

void FFUNC(zgeqp3,ZGEQP3)( integer* m, integer* n, dcomplex* a, integer* lda,
                    integer* jpvt, dcomplex* tau, dcomplex* work,
                    integer* lwork, double* rwork, integer *info );

void FFUNC(sorgqr,SORGQR)( integer* m, integer* n, integer* k, float* a,
                    integer* lda, const float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dorgqr,DORGQR)( integer* m, integer* n, integer* k, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(cungqr,CUNGQR)( integer* m, integer* n, integer* k, scomplex* a,
                    integer* lda, const scomplex* tau, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zungqr,ZUNGQR)( integer* m, integer* n, integer* k, dcomplex* a,
                    integer* lda, const dcomplex* tau, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sormqr,SORMQR)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormqr,DORMQR)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(cunmqr,CUNMQR)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const scomplex* a, integer* lda,
                    const scomplex* tau, scomplex* c, integer* ldc, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zunmqr,ZUNMQR)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const dcomplex* a, integer* lda,
                    const dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sgelqf,SGELQF)( integer* m, integer* n, float* a, integer* lda,
                    float* tau, float* work, integer* lwork,
                    integer *info );

void FFUNC(dgelqf,DGELQF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(cgelqf,CGELQF)( integer* m, integer* n, scomplex* a, integer* lda,
                    scomplex* tau, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zgelqf,ZGELQF)( integer* m, integer* n, dcomplex* a, integer* lda,
                    dcomplex* tau, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sorglq,SORGLQ)( integer* m, integer* n, integer* k, float* a,
                    integer* lda, const float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dorglq,DORGLQ)( integer* m, integer* n, integer* k, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(cunglq,CUNGLQ)( integer* m, integer* n, integer* k, scomplex* a,
                    integer* lda, const scomplex* tau, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zunglq,ZUNGLQ)( integer* m, integer* n, integer* k, dcomplex* a,
                    integer* lda, const dcomplex* tau, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sormlq,SORMLQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormlq,DORMLQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(cunmlq,CUNMLQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const scomplex* a, integer* lda,
                    const scomplex* tau, scomplex* c, integer* ldc, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zunmlq,ZUNMLQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const dcomplex* a, integer* lda,
                    const dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sgeqlf,SGEQLF)( integer* m, integer* n, float* a, integer* lda,
                    float* tau, float* work, integer* lwork,
                    integer *info );

void FFUNC(dgeqlf,DGEQLF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(cgeqlf,CGEQLF)( integer* m, integer* n, scomplex* a, integer* lda,
                    scomplex* tau, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zgeqlf,ZGEQLF)( integer* m, integer* n, dcomplex* a, integer* lda,
                    dcomplex* tau, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sorgql,SORGQL)( integer* m, integer* n, integer* k, float* a,
                    integer* lda, const float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dorgql,DORGQL)( integer* m, integer* n, integer* k, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(cungql,CUNGQL)( integer* m, integer* n, integer* k, scomplex* a,
                    integer* lda, const scomplex* tau, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zungql,ZUNGQL)( integer* m, integer* n, integer* k, dcomplex* a,
                    integer* lda, const dcomplex* tau, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sormql,SORMQL)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormql,DORMQL)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(cunmql,CUNMQL)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const scomplex* a, integer* lda,
                    const scomplex* tau, scomplex* c, integer* ldc, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zunmql,ZUNMQL)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const dcomplex* a, integer* lda,
                    const dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sgerqf,SGERQF)( integer* m, integer* n, float* a, integer* lda,
                    float* tau, float* work, integer* lwork,
                    integer *info );

void FFUNC(dgerqf,DGERQF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(cgerqf,CGERQF)( integer* m, integer* n, scomplex* a, integer* lda,
                    scomplex* tau, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zgerqf,ZGERQF)( integer* m, integer* n, dcomplex* a, integer* lda,
                    dcomplex* tau, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sorgrq,SORGRQ)( integer* m, integer* n, integer* k, float* a,
                    integer* lda, const float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dorgrq,DORGRQ)( integer* m, integer* n, integer* k, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(cungrq,CUNGRQ)( integer* m, integer* n, integer* k, scomplex* a,
                    integer* lda, const scomplex* tau, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zungrq,ZUNGRQ)( integer* m, integer* n, integer* k, dcomplex* a,
                    integer* lda, const dcomplex* tau, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sormrq,SORMRQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormrq,DORMRQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(cunmrq,CUNMRQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const scomplex* a, integer* lda,
                    const scomplex* tau, scomplex* c, integer* ldc, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zunmrq,ZUNMRQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const dcomplex* a, integer* lda,
                    const dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(stzrzf,STZRZF)( integer* m, integer* n, float* a, integer* lda,
                    float* tau, float* work, integer* lwork,
                    integer *info );

void FFUNC(dtzrzf,DTZRZF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(ctzrzf,CTZRZF)( integer* m, integer* n, scomplex* a, integer* lda,
                    scomplex* tau, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(ztzrzf,ZTZRZF)( integer* m, integer* n, dcomplex* a, integer* lda,
                    dcomplex* tau, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sormrz,SORMRZ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, integer* l, const float* a,
                    integer* lda, const float* tau, float* c,
                    integer* ldc, float* work, integer* lwork,
                    integer *info );

void FFUNC(dormrz,DORMRZ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, integer* l, const double* a,
                    integer* lda, const double* tau, double* c,
                    integer* ldc, double* work, integer* lwork,
                    integer *info );

void FFUNC(cunmrz,CUNMRZ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, integer* l, const scomplex* a,
                    integer* lda, const scomplex* tau, scomplex* c,
                    integer* ldc, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zunmrz,ZUNMRZ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, integer* l, const dcomplex* a,
                    integer* lda, const dcomplex* tau, dcomplex* c,
                    integer* ldc, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sggqrf,SGGQRF)( integer* n, integer* m, integer* p, float* a,
                    integer* lda, float* taua, float* b, integer* ldb,
                    float* taub, float* work, integer* lwork,
                    integer *info );

void FFUNC(dggqrf,DGGQRF)( integer* n, integer* m, integer* p, double* a,
                    integer* lda, double* taua, double* b, integer* ldb,
                    double* taub, double* work, integer* lwork,
                    integer *info );

void FFUNC(cggqrf,CGGQRF)( integer* n, integer* m, integer* p, scomplex* a,
                    integer* lda, scomplex* taua, scomplex* b, integer* ldb,
                    scomplex* taub, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zggqrf,ZGGQRF)( integer* n, integer* m, integer* p, dcomplex* a,
                    integer* lda, dcomplex* taua, dcomplex* b, integer* ldb,
                    dcomplex* taub, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sggrqf,SGGRQF)( integer* m, integer* p, integer* n, float* a,
                    integer* lda, float* taua, float* b, integer* ldb,
                    float* taub, float* work, integer* lwork,
                    integer *info );

void FFUNC(dggrqf,DGGRQF)( integer* m, integer* p, integer* n, double* a,
                    integer* lda, double* taua, double* b, integer* ldb,
                    double* taub, double* work, integer* lwork,
                    integer *info );

void FFUNC(cggrqf,CGGRQF)( integer* m, integer* p, integer* n, scomplex* a,
                    integer* lda, scomplex* taua, scomplex* b, integer* ldb,
                    scomplex* taub, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zggrqf,ZGGRQF)( integer* m, integer* p, integer* n, dcomplex* a,
                    integer* lda, dcomplex* taua, dcomplex* b, integer* ldb,
                    dcomplex* taub, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sgebrd,SGEBRD)( integer* m, integer* n, float* a, integer* lda,
                    float* d, float* e, float* tauq, float* taup,
                    float* work, integer* lwork, integer *info );

void FFUNC(dgebrd,DGEBRD)( integer* m, integer* n, double* a, integer* lda,
                    double* d, double* e, double* tauq, double* taup,
                    double* work, integer* lwork, integer *info );

void FFUNC(cgebrd,CGEBRD)( integer* m, integer* n, scomplex* a, integer* lda,
                    float* d, float* e, scomplex* tauq, scomplex* taup,
                    scomplex* work, integer* lwork, integer *info );

void FFUNC(zgebrd,ZGEBRD)( integer* m, integer* n, dcomplex* a, integer* lda,
                    double* d, double* e, dcomplex* tauq, dcomplex* taup,
                    dcomplex* work, integer* lwork, integer *info );

void FFUNC(sgbbrd,SGBBRD)( char* vect, integer* m, integer* n, integer* ncc,
                    integer* kl, integer* ku, float* ab,
                    integer* ldab, float* d, float* e, float* q,
                    integer* ldq, float* pt, integer* ldpt, float* c,
                    integer* ldc, float* work, integer *info );

void FFUNC(dgbbrd,DGBBRD)( char* vect, integer* m, integer* n, integer* ncc,
                    integer* kl, integer* ku, double* ab,
                    integer* ldab, double* d, double* e, double* q,
                    integer* ldq, double* pt, integer* ldpt, double* c,
                    integer* ldc, double* work, integer *info );

void FFUNC(cgbbrd,CGBBRD)( char* vect, integer* m, integer* n, integer* ncc,
                    integer* kl, integer* ku, scomplex* ab,
                    integer* ldab, float* d, float* e, scomplex* q,
                    integer* ldq, scomplex* pt, integer* ldpt, scomplex* c,
                    integer* ldc, scomplex* work, float* rwork, integer *info );

void FFUNC(zgbbrd,ZGBBRD)( char* vect, integer* m, integer* n, integer* ncc,
                    integer* kl, integer* ku, dcomplex* ab,
                    integer* ldab, double* d, double* e, dcomplex* q,
                    integer* ldq, dcomplex* pt, integer* ldpt, dcomplex* c,
                    integer* ldc, dcomplex* work, double* rwork, integer *info );

void FFUNC(sorgbr,SORGBR)( char* vect, integer* m, integer* n, integer* k,
                    float* a, integer* lda, const float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dorgbr,DORGBR)( char* vect, integer* m, integer* n, integer* k,
                    double* a, integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(cungbr,CUNGBR)( char* vect, integer* m, integer* n, integer* k,
                    scomplex* a, integer* lda, const scomplex* tau, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zungbr,ZUNGBR)( char* vect, integer* m, integer* n, integer* k,
                    dcomplex* a, integer* lda, const dcomplex* tau, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sormbr,SORMBR)( char* vect, char* side, char* trans, integer* m,
                    integer* n, integer* k, const float* a,
                    integer* lda, const float* tau, float* c,
                    integer* ldc, float* work, integer* lwork,
                    integer *info );

void FFUNC(dormbr,DORMBR)( char* vect, char* side, char* trans, integer* m,
                    integer* n, integer* k, const double* a,
                    integer* lda, const double* tau, double* c,
                    integer* ldc, double* work, integer* lwork,
                    integer *info );

void FFUNC(cunmbr,CUNMBR)( char* vect, char* side, char* trans, integer* m,
                    integer* n, integer* k, const scomplex* a,
                    integer* lda, const scomplex* tau, scomplex* c,
                    integer* ldc, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zunmbr,ZUNMBR)( char* vect, char* side, char* trans, integer* m,
                    integer* n, integer* k, const dcomplex* a,
                    integer* lda, const dcomplex* tau, dcomplex* c,
                    integer* ldc, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sbdsqr,SBDSQR)( char* uplo, integer* n, integer* ncvt,
                    integer* nru, integer* ncc, float* d, float* e,
                    float* vt, integer* ldvt, float* u, integer* ldu,
                    float* c, integer* ldc, float* work,
                    integer *info );

void FFUNC(dbdsqr,DBDSQR)( char* uplo, integer* n, integer* ncvt,
                    integer* nru, integer* ncc, double* d, double* e,
                    double* vt, integer* ldvt, double* u, integer* ldu,
                    double* c, integer* ldc, double* work,
                    integer *info );

void FFUNC(cbdsqr,CBDSQR)( char* uplo, integer* n, integer* ncvt,
                    integer* nru, integer* ncc, float* d, float* e,
                    scomplex* vt, integer* ldvt, scomplex* u, integer* ldu,
                    scomplex* c, integer* ldc, float* rwork,
                    integer *info );

void FFUNC(zbdsqr,ZBDSQR)( char* uplo, integer* n, integer* ncvt,
                    integer* nru, integer* ncc, double* d, double* e,
                    dcomplex* vt, integer* ldvt, dcomplex* u, integer* ldu,
                    dcomplex* c, integer* ldc, double* rwork,
                    integer *info );

void FFUNC(sbdsdc,SBDSDC)( char* uplo, char* compq, integer* n, float* d,
                    float* e, float* u, integer* ldu, float* vt,
                    integer* ldvt, float* q, integer* iq, float* work,
                    integer* iwork, integer *info );

void FFUNC(dbdsdc,DBDSDC)( char* uplo, char* compq, integer* n, double* d,
                    double* e, double* u, integer* ldu, double* vt,
                    integer* ldvt, double* q, integer* iq, double* work,
                    integer* iwork, integer *info );

void FFUNC(ssytrd,SSYTRD)( char* uplo, integer* n, float* a, integer* lda,
                    float* d, float* e, float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dsytrd,DSYTRD)( char* uplo, integer* n, double* a, integer* lda,
                    double* d, double* e, double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(chetrd,CHETRD)( char* uplo, integer* n, scomplex* a, integer* lda,
                    float* d, float* e, scomplex* tau, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zhetrd,ZHETRD)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    double* d, double* e, dcomplex* tau, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sorgtr,SORGTR)( char* uplo, integer* n, float* a, integer* lda,
                    const float* tau, float* work, integer* lwork,
                    integer *info );

void FFUNC(dorgtr,DORGTR)( char* uplo, integer* n, double* a, integer* lda,
                    const double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(cungtr,CUNGTR)( char* uplo, integer* n, scomplex* a, integer* lda,
                    const scomplex* tau, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zungtr,ZUNGTR)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    const dcomplex* tau, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sormtr,SORMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormtr,DORMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(cunmtr,CUNMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const scomplex* a, integer* lda,
                    const scomplex* tau, scomplex* c, integer* ldc, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zunmtr,ZUNMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const dcomplex* a, integer* lda,
                    const dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(ssptrd,SSPTRD)( char* uplo, integer* n, float* ap, float* d, float* e,
                    float* tau, integer *info );

void FFUNC(dsptrd,DSPTRD)( char* uplo, integer* n, double* ap, double* d, double* e,
                    double* tau, integer *info );

void FFUNC(chptrd,CHPTRD)( char* uplo, integer* n, scomplex* ap, float* d, float* e,
                    scomplex* tau, integer *info );

void FFUNC(zhptrd,ZHPTRD)( char* uplo, integer* n, dcomplex* ap, double* d, double* e,
                    dcomplex* tau, integer *info );

void FFUNC(sopgtr,SOPGTR)( char* uplo, integer* n, const float* ap,
                    const float* tau, float* q, integer* ldq, float* work,
                    integer *info );

void FFUNC(dopgtr,DOPGTR)( char* uplo, integer* n, const double* ap,
                    const double* tau, double* q, integer* ldq, double* work,
                    integer *info );

void FFUNC(cupgtr,CUPGTR)( char* uplo, integer* n, const scomplex* ap,
                    const scomplex* tau, scomplex* q, integer* ldq, scomplex* work,
                    integer *info );

void FFUNC(zupgtr,ZUPGTR)( char* uplo, integer* n, const dcomplex* ap,
                    const dcomplex* tau, dcomplex* q, integer* ldq, dcomplex* work,
                    integer *info );

void FFUNC(sopmtr,SOPMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const float* ap, const float* tau,
                    float* c, integer* ldc, float* work,
                    integer *info );

void FFUNC(dopmtr,DOPMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const double* ap, const double* tau,
                    double* c, integer* ldc, double* work,
                    integer *info );

void FFUNC(cupmtr,CUPMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const scomplex* ap, const scomplex* tau,
                    scomplex* c, integer* ldc, scomplex* work,
                    integer *info );

void FFUNC(zupmtr,ZUPMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const dcomplex* ap, const dcomplex* tau,
                    dcomplex* c, integer* ldc, dcomplex* work,
                    integer *info );

void FFUNC(ssbtrd,SSBTRD)( char* vect, char* uplo, integer* n, integer* kd,
                    float* ab, integer* ldab, float* d, float* e,
                    float* q, integer* ldq, float* work,
                    integer *info );

void FFUNC(dsbtrd,DSBTRD)( char* vect, char* uplo, integer* n, integer* kd,
                    double* ab, integer* ldab, double* d, double* e,
                    double* q, integer* ldq, double* work,
                    integer *info );

void FFUNC(chbtrd,CHBTRD)( char* vect, char* uplo, integer* n, integer* kd,
                    scomplex* ab, integer* ldab, float* d, float* e,
                    scomplex* q, integer* ldq, scomplex* work,
                    integer *info );

void FFUNC(zhbtrd,ZHBTRD)( char* vect, char* uplo, integer* n, integer* kd,
                    dcomplex* ab, integer* ldab, double* d, double* e,
                    dcomplex* q, integer* ldq, dcomplex* work,
                    integer *info );

void FFUNC(ssterf,SSTERF)( integer* n, float* d, float* e, integer *info );

void FFUNC(dsterf,DSTERF)( integer* n, double* d, double* e, integer *info );

void FFUNC(ssteqr,SSTEQR)( char* compz, integer* n, float* d, float* e, float* z,
                    integer* ldz, float* work, integer *info );

void FFUNC(dsteqr,DSTEQR)( char* compz, integer* n, double* d, double* e, double* z,
                    integer* ldz, double* work, integer *info );

void FFUNC(csteqr,CSTEQR)( char* compz, integer* n, float* d, float* e, scomplex* z,
                    integer* ldz, float* work, integer *info );

void FFUNC(zsteqr,ZSTEQR)( char* compz, integer* n, double* d, double* e, dcomplex* z,
                    integer* ldz, double* work, integer *info );

void FFUNC(sstemr,SSTEMR)( char* jobz, char* range, integer* n, float* d,
                    float* e, float* vl, float* vu, integer* il,
                    integer* iu, integer* m, float* w, float* z,
                    integer* ldz, integer* nzc, integer* isuppz,
                    logical* tryrac, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dstemr,DSTEMR)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, integer* m, double* w, double* z,
                    integer* ldz, integer* nzc, integer* isuppz,
                    logical* tryrac, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(cstemr,CSTEMR)( char* jobz, char* range, integer* n, float* d,
                    float* e, float* vl, float* vu, integer* il,
                    integer* iu, integer* m, float* w, scomplex* z,
                    integer* ldz, integer* nzc, integer* isuppz,
                    logical* tryrac, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(zstemr,ZSTEMR)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, integer* m, double* w, dcomplex* z,
                    integer* ldz, integer* nzc, integer* isuppz,
                    logical* tryrac, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(sstedc,SSTEDC)( char* compz, integer* n, float* d, float* e, float* z,
                    integer* ldz, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dstedc,DSTEDC)( char* compz, integer* n, double* d, double* e, double* z,
                    integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(cstedc,CSTEDC)( char* compz, integer* n, float* d, float* e, scomplex* z,
                    integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(zstedc,ZSTEDC)( char* compz, integer* n, double* d, double* e, dcomplex* z,
                    integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(sstegr,SSTEGR)( char* jobz, char* range, integer* n, float* d,
                    float* e, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    float* z, integer* ldz, integer* isuppz,
                    float* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(dstegr,DSTEGR)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, integer* isuppz,
                    double* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(cstegr,CSTEGR)( char* jobz, char* range, integer* n, float* d,
                    float* e, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    scomplex* z, integer* ldz, integer* isuppz,
                    float* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(zstegr,ZSTEGR)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    dcomplex* z, integer* ldz, integer* isuppz,
                    double* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(spteqr,SPTEQR)( char* compz, integer* n, float* d, float* e, float* z,
                    integer* ldz, float* work, integer *info );

void FFUNC(dpteqr,DPTEQR)( char* compz, integer* n, double* d, double* e, double* z,
                    integer* ldz, double* work, integer *info );

void FFUNC(cpteqr,CPTEQR)( char* compz, integer* n, float* d, float* e, scomplex* z,
                    integer* ldz, float* work, integer *info );

void FFUNC(zpteqr,ZPTEQR)( char* compz, integer* n, double* d, double* e, dcomplex* z,
                    integer* ldz, double* work, integer *info );

void FFUNC(sstebz,SSTEBZ)( char* range, char* order, integer* n, float* vl,
                    float* vu, integer* il, integer* iu, float* abstol,
                    const float* d, const float* e, integer* m,
                    integer* nsplit, float* w, integer* iblock,
                    integer* isplit, float* work, integer* iwork,
                    integer *info );

void FFUNC(dstebz,DSTEBZ)( char* range, char* order, integer* n, double* vl,
                    double* vu, integer* il, integer* iu, double* abstol,
                    const double* d, const double* e, integer* m,
                    integer* nsplit, double* w, integer* iblock,
                    integer* isplit, double* work, integer* iwork,
                    integer *info );

void FFUNC(sstein,SSTEIN)( integer* n, const float* d, const float* e,
                    integer* m, const float* w, const integer* iblock,
                    const integer* isplit, float* z, integer* ldz,
                    float* work, integer* iwork, integer* ifailv,
                    integer *info );

void FFUNC(dstein,DSTEIN)( integer* n, const double* d, const double* e,
                    integer* m, const double* w, const integer* iblock,
                    const integer* isplit, double* z, integer* ldz,
                    double* work, integer* iwork, integer* ifailv,
                    integer *info );

void FFUNC(cstein,CSTEIN)( integer* n, const float* d, const float* e,
                    integer* m, const float* w, const integer* iblock,
                    const integer* isplit, scomplex* z, integer* ldz,
                    float* work, integer* iwork, integer* ifailv,
                    integer *info );

void FFUNC(zstein,ZSTEIN)( integer* n, const double* d, const double* e,
                    integer* m, const double* w, const integer* iblock,
                    const integer* isplit, dcomplex* z, integer* ldz,
                    double* work, integer* iwork, integer* ifailv,
                    integer *info );

void FFUNC(sdisna,SDISNA)( char* job, integer* m, integer* n, const float* d,
                    float* sep, integer *info );

void FFUNC(ddisna,DDISNA)( char* job, integer* m, integer* n, const double* d,
                    double* sep, integer *info );

void FFUNC(ssygst,SSYGST)( integer* itype, char* uplo, integer* n, float* a,
                    integer* lda, const float* b, integer* ldb,
                    integer *info );

void FFUNC(dsygst,DSYGST)( integer* itype, char* uplo, integer* n, double* a,
                    integer* lda, const double* b, integer* ldb,
                    integer *info );

void FFUNC(chegst,CHEGST)( integer* itype, char* uplo, integer* n, scomplex* a,
                    integer* lda, const scomplex* b, integer* ldb,
                    integer *info );

void FFUNC(zhegst,ZHEGST)( integer* itype, char* uplo, integer* n, dcomplex* a,
                    integer* lda, const dcomplex* b, integer* ldb,
                    integer *info );

void FFUNC(sspgst,SSPGST)( integer* itype, char* uplo, integer* n, float* ap,
                    const float* bp, integer *info );

void FFUNC(dspgst,DSPGST)( integer* itype, char* uplo, integer* n, double* ap,
                    const double* bp, integer *info );

void FFUNC(chpgst,CHPGST)( integer* itype, char* uplo, integer* n, scomplex* ap,
                    const scomplex* bp, integer *info );

void FFUNC(zhpgst,ZHPGST)( integer* itype, char* uplo, integer* n, dcomplex* ap,
                    const dcomplex* bp, integer *info );

void FFUNC(ssbgst,SSBGST)( char* vect, char* uplo, integer* n, integer* ka,
                    integer* kb, float* ab, integer* ldab,
                    const float* bb, integer* ldbb, float* x,
                    integer* ldx, float* work, integer *info );

void FFUNC(dsbgst,DSBGST)( char* vect, char* uplo, integer* n, integer* ka,
                    integer* kb, double* ab, integer* ldab,
                    const double* bb, integer* ldbb, double* x,
                    integer* ldx, double* work, integer *info );

void FFUNC(chbgst,CHBGST)( char* vect, char* uplo, integer* n, integer* ka,
                    integer* kb, scomplex* ab, integer* ldab,
                    const scomplex* bb, integer* ldbb, scomplex* x,
                    integer* ldx, scomplex* work, float* rwork, integer *info );

void FFUNC(zhbgst,ZHBGST)( char* vect, char* uplo, integer* n, integer* ka,
                    integer* kb, dcomplex* ab, integer* ldab,
                    const dcomplex* bb, integer* ldbb, dcomplex* x,
                    integer* ldx, dcomplex* work, double* rwork, integer *info );

void FFUNC(spbstf,SPBSTF)( char* uplo, integer* n, integer* kb, float* bb,
                    integer* ldbb, integer *info );

void FFUNC(dpbstf,DPBSTF)( char* uplo, integer* n, integer* kb, double* bb,
                    integer* ldbb, integer *info );

void FFUNC(cpbstf,CPBSTF)( char* uplo, integer* n, integer* kb, scomplex* bb,
                    integer* ldbb, integer *info );

void FFUNC(zpbstf,ZPBSTF)( char* uplo, integer* n, integer* kb, dcomplex* bb,
                    integer* ldbb, integer *info );

void FFUNC(sgehrd,SGEHRD)( integer* n, integer* ilo, integer* ihi, float* a,
                    integer* lda, float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dgehrd,DGEHRD)( integer* n, integer* ilo, integer* ihi, double* a,
                    integer* lda, double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(cgehrd,CGEHRD)( integer* n, integer* ilo, integer* ihi, scomplex* a,
                    integer* lda, scomplex* tau, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zgehrd,ZGEHRD)( integer* n, integer* ilo, integer* ihi, dcomplex* a,
                    integer* lda, dcomplex* tau, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sorghr,SORGHR)( integer* n, integer* ilo, integer* ihi, float* a,
                    integer* lda, const float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dorghr,DORGHR)( integer* n, integer* ilo, integer* ihi, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(cunghr,CUNGHR)( integer* n, integer* ilo, integer* ihi, scomplex* a,
                    integer* lda, const scomplex* tau, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zunghr,ZUNGHR)( integer* n, integer* ilo, integer* ihi, dcomplex* a,
                    integer* lda, const dcomplex* tau, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sormhr,SORMHR)( char* side, char* trans, integer* m, integer* n,
                    integer* ilo, integer* ihi, const float* a,
                    integer* lda, const float* tau, float* c,
                    integer* ldc, float* work, integer* lwork,
                    integer *info );

void FFUNC(dormhr,DORMHR)( char* side, char* trans, integer* m, integer* n,
                    integer* ilo, integer* ihi, const double* a,
                    integer* lda, const double* tau, double* c,
                    integer* ldc, double* work, integer* lwork,
                    integer *info );

void FFUNC(cunmhr,CUNMHR)( char* side, char* trans, integer* m, integer* n,
                    integer* ilo, integer* ihi, const scomplex* a,
                    integer* lda, const scomplex* tau, scomplex* c,
                    integer* ldc, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zunmhr,ZUNMHR)( char* side, char* trans, integer* m, integer* n,
                    integer* ilo, integer* ihi, const dcomplex* a,
                    integer* lda, const dcomplex* tau, dcomplex* c,
                    integer* ldc, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sgebal,SGEBAL)( char* job, integer* n, float* a, integer* lda,
                    integer* ilo, integer* ihi, float* scale,
                    integer *info );

void FFUNC(dgebal,DGEBAL)( char* job, integer* n, double* a, integer* lda,
                    integer* ilo, integer* ihi, double* scale,
                    integer *info );

void FFUNC(cgebal,CGEBAL)( char* job, integer* n, scomplex* a, integer* lda,
                    integer* ilo, integer* ihi, float* scale,
                    integer *info );

void FFUNC(zgebal,ZGEBAL)( char* job, integer* n, dcomplex* a, integer* lda,
                    integer* ilo, integer* ihi, double* scale,
                    integer *info );

void FFUNC(sgebak,SGEBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const float* scale, integer* m,
                    float* v, integer* ldv, integer *info );

void FFUNC(dgebak,DGEBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const double* scale, integer* m,
                    double* v, integer* ldv, integer *info );

void FFUNC(cgebak,CGEBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const float* scale, integer* m,
                    scomplex* v, integer* ldv, integer *info );

void FFUNC(zgebak,ZGEBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const double* scale, integer* m,
                    dcomplex* v, integer* ldv, integer *info );

void FFUNC(shseqr,SHSEQR)( char* job, char* compz, integer* n, integer* ilo,
                    integer* ihi, float* h, integer* ldh, float* wr,
                    float* wi, float* z, integer* ldz, float* work,
                    integer* lwork, integer *info );

void FFUNC(dhseqr,DHSEQR)( char* job, char* compz, integer* n, integer* ilo,
                    integer* ihi, double* h, integer* ldh, double* wr,
                    double* wi, double* z, integer* ldz, double* work,
                    integer* lwork, integer *info );

void FFUNC(chseqr,CHSEQR)( char* job, char* compz, integer* n, integer* ilo,
                    integer* ihi, scomplex* h, integer* ldh, scomplex* w,
                    scomplex* z, integer* ldz, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zhseqr,ZHSEQR)( char* job, char* compz, integer* n, integer* ilo,
                    integer* ihi, dcomplex* h, integer* ldh, dcomplex* w,
                    dcomplex* z, integer* ldz, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(shsein,SHSEIN)( char* job, char* eigsrc, char* initv,
                    logical* select, integer* n, const float* h,
                    integer* ldh, float* wr, const float* wi, float* vl,
                    integer* ldvl, float* vr, integer* ldvr,
                    integer* mm, integer* m, float* work,
                    integer* ifaill, integer* ifailr, integer *info );

void FFUNC(dhsein,DHSEIN)( char* job, char* eigsrc, char* initv,
                    logical* select, integer* n, const double* h,
                    integer* ldh, double* wr, const double* wi, double* vl,
                    integer* ldvl, double* vr, integer* ldvr,
                    integer* mm, integer* m, double* work,
                    integer* ifaill, integer* ifailr, integer *info );

void FFUNC(chsein,CHSEIN)( char* job, char* eigsrc, char* initv,
                    logical* select, integer* n, const scomplex* h,
                    integer* ldh, scomplex* w, scomplex* vl,
                    integer* ldvl, scomplex* vr, integer* ldvr,
                    integer* mm, integer* m, scomplex* work, float* rwork,
                    integer* ifaill, integer* ifailr, integer *info );

void FFUNC(zhsein,ZHSEIN)( char* job, char* eigsrc, char* initv,
                    logical* select, integer* n, const dcomplex* h,
                    integer* ldh, dcomplex* w, dcomplex* vl,
                    integer* ldvl, dcomplex* vr, integer* ldvr,
                    integer* mm, integer* m, dcomplex* work, double* rwork,
                    integer* ifaill, integer* ifailr, integer *info );

void FFUNC(strevc,STREVC)( char* side, char* howmny, const logical* select,
                    integer* n, const float* t, integer* ldt, float* vl,
                    integer* ldvl, float* vr, integer* ldvr,
                    integer* mm, integer* m, float* work,
                    integer *info );

void FFUNC(dtrevc,DTREVC)( char* side, char* howmny, const logical* select,
                    integer* n, const double* t, integer* ldt, double* vl,
                    integer* ldvl, double* vr, integer* ldvr,
                    integer* mm, integer* m, double* work,
                    integer *info );

void FFUNC(ctrevc,CTREVC)( char* side, char* howmny, const logical* select,
                    integer* n, const scomplex* t, integer* ldt, scomplex* vl,
                    integer* ldvl, scomplex* vr, integer* ldvr,
                    integer* mm, integer* m, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(ztrevc,ZTREVC)( char* side, char* howmny, const logical* select,
                    integer* n, const dcomplex* t, integer* ldt, dcomplex* vl,
                    integer* ldvl, dcomplex* vr, integer* ldvr,
                    integer* mm, integer* m, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(strsna,STRSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const float* t, integer* ldt,
                    const float* vl, integer* ldvl, const float* vr,
                    integer* ldvr, float* s, float* sep, integer* mm,
                    integer* m, float* work, integer* ldwork,
                    integer* iwork, integer *info );

void FFUNC(dtrsna,DTRSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const double* t, integer* ldt,
                    const double* vl, integer* ldvl, const double* vr,
                    integer* ldvr, double* s, double* sep, integer* mm,
                    integer* m, double* work, integer* ldwork,
                    integer* iwork, integer *info );

void FFUNC(ctrsna,CTRSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const scomplex* t, integer* ldt,
                    const scomplex* vl, integer* ldvl, const scomplex* vr,
                    integer* ldvr, float* s, float* sep, integer* mm,
                    integer* m, scomplex* work, integer* ldwork,
                    float* rwork, integer *info );

void FFUNC(ztrsna,ZTRSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const dcomplex* t, integer* ldt,
                    const dcomplex* vl, integer* ldvl, const dcomplex* vr,
                    integer* ldvr, double* s, double* sep, integer* mm,
                    integer* m, dcomplex* work, integer* ldwork,
                    double* rwork, integer *info );

void FFUNC(strexc,STREXC)( char* compq, integer* n, float* t, integer* ldt,
                    float* q, integer* ldq, integer* ifst,
                    integer* ilst, float* work, integer *info );

void FFUNC(dtrexc,DTREXC)( char* compq, integer* n, double* t, integer* ldt,
                    double* q, integer* ldq, integer* ifst,
                    integer* ilst, double* work, integer *info );

void FFUNC(ctrexc,CTREXC)( char* compq, integer* n, scomplex* t, integer* ldt,
                    scomplex* q, integer* ldq, integer* ifst,
                    integer* ilst, integer *info );

void FFUNC(ztrexc,ZTREXC)( char* compq, integer* n, dcomplex* t, integer* ldt,
                    dcomplex* q, integer* ldq, integer* ifst,
                    integer* ilst, integer *info );

void FFUNC(strsen,STRSEN)( char* job, char* compq, const logical* select,
                    integer* n, float* t, integer* ldt, float* q,
                    integer* ldq, float* wr, float* wi, integer* m,
                    float* s, float* sep, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dtrsen,DTRSEN)( char* job, char* compq, const logical* select,
                    integer* n, double* t, integer* ldt, double* q,
                    integer* ldq, double* wr, double* wi, integer* m,
                    double* s, double* sep, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(ctrsen,CTRSEN)( char* job, char* compq, const logical* select,
                    integer* n, scomplex* t, integer* ldt, scomplex* q,
                    integer* ldq, scomplex* w, integer* m,
                    float* s, float* sep, scomplex* work, integer* lwork, integer *info );

void FFUNC(ztrsen,ZTRSEN)( char* job, char* compq, const logical* select,
                    integer* n, dcomplex* t, integer* ldt, dcomplex* q,
                    integer* ldq, dcomplex* w, integer* m,
                    double* s, double* sep, dcomplex* work, integer* lwork, integer *info );

void FFUNC(strsyl,STRSYL)( char* trana, char* tranb, integer* isgn, integer* m,
                    integer* n, const float* a, integer* lda,
                    const float* b, integer* ldb, float* c,
                    integer* ldc, float* scale, integer *info );

void FFUNC(dtrsyl,DTRSYL)( char* trana, char* tranb, integer* isgn, integer* m,
                    integer* n, const double* a, integer* lda,
                    const double* b, integer* ldb, double* c,
                    integer* ldc, double* scale, integer *info );

void FFUNC(ctrsyl,CTRSYL)( char* trana, char* tranb, integer* isgn, integer* m,
                    integer* n, const scomplex* a, integer* lda,
                    const scomplex* b, integer* ldb, scomplex* c,
                    integer* ldc, float* scale, integer *info );

void FFUNC(ztrsyl,ZTRSYL)( char* trana, char* tranb, integer* isgn, integer* m,
                    integer* n, const dcomplex* a, integer* lda,
                    const dcomplex* b, integer* ldb, dcomplex* c,
                    integer* ldc, double* scale, integer *info );

void FFUNC(sgghrd,SGGHRD)( char* compq, char* compz, integer* n, integer* ilo,
                    integer* ihi, float* a, integer* lda, float* b,
                    integer* ldb, float* q, integer* ldq, float* z,
                    integer* ldz, integer *info );

void FFUNC(dgghrd,DGGHRD)( char* compq, char* compz, integer* n, integer* ilo,
                    integer* ihi, double* a, integer* lda, double* b,
                    integer* ldb, double* q, integer* ldq, double* z,
                    integer* ldz, integer *info );

void FFUNC(cgghrd,CGGHRD)( char* compq, char* compz, integer* n, integer* ilo,
                    integer* ihi, scomplex* a, integer* lda, scomplex* b,
                    integer* ldb, scomplex* q, integer* ldq, scomplex* z,
                    integer* ldz, integer *info );

void FFUNC(zgghrd,ZGGHRD)( char* compq, char* compz, integer* n, integer* ilo,
                    integer* ihi, dcomplex* a, integer* lda, dcomplex* b,
                    integer* ldb, dcomplex* q, integer* ldq, dcomplex* z,
                    integer* ldz, integer *info );

void FFUNC(sggbal,SGGBAL)( char* job, integer* n, float* a, integer* lda,
                    float* b, integer* ldb, integer* ilo,
                    integer* ihi, float* lscale, float* rscale,
                    float* work, integer *info );

void FFUNC(dggbal,DGGBAL)( char* job, integer* n, double* a, integer* lda,
                    double* b, integer* ldb, integer* ilo,
                    integer* ihi, double* lscale, double* rscale,
                    double* work, integer *info );

void FFUNC(cggbal,CGGBAL)( char* job, integer* n, scomplex* a, integer* lda,
                    scomplex* b, integer* ldb, integer* ilo,
                    integer* ihi, float* lscale, float* rscale,
                    float* work, integer *info );

void FFUNC(zggbal,ZGGBAL)( char* job, integer* n, dcomplex* a, integer* lda,
                    dcomplex* b, integer* ldb, integer* ilo,
                    integer* ihi, double* lscale, double* rscale,
                    double* work, integer *info );

void FFUNC(sggbak,SGGBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const float* lscale, const float* rscale,
                    integer* m, float* v, integer* ldv,
                    integer *info );

void FFUNC(dggbak,DGGBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const double* lscale, const double* rscale,
                    integer* m, double* v, integer* ldv,
                    integer *info );

void FFUNC(cggbak,CGGBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const float* lscale, const float* rscale,
                    integer* m, scomplex* v, integer* ldv,
                    integer *info );

void FFUNC(zggbak,ZGGBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const double* lscale, const double* rscale,
                    integer* m, dcomplex* v, integer* ldv,
                    integer *info );

void FFUNC(shgeqz,SHGEQZ)( char* job, char* compq, char* compz, integer* n,
                    integer* ilo, integer* ihi, float* h,
                    integer* ldh, float* t, integer* ldt, float* alphar,
                    float* alphai, float* beta, float* q, integer* ldq,
                    float* z, integer* ldz, float* work, integer* lwork,
                    integer *info );

void FFUNC(dhgeqz,DHGEQZ)( char* job, char* compq, char* compz, integer* n,
                    integer* ilo, integer* ihi, double* h,
                    integer* ldh, double* t, integer* ldt, double* alphar,
                    double* alphai, double* beta, double* q, integer* ldq,
                    double* z, integer* ldz, double* work, integer* lwork,
                    integer *info );

void FFUNC(chgeqz,CHGEQZ)( char* job, char* compq, char* compz, integer* n,
                    integer* ilo, integer* ihi, scomplex* h,
                    integer* ldh, scomplex* t, integer* ldt, scomplex* alpha,
                    scomplex* beta, scomplex* q, integer* ldq,
                    scomplex* z, integer* ldz, scomplex* work, integer* lwork,
                    float* rwork, integer *info );

void FFUNC(zhgeqz,ZHGEQZ)( char* job, char* compq, char* compz, integer* n,
                    integer* ilo, integer* ihi, dcomplex* h,
                    integer* ldh, dcomplex* t, integer* ldt, dcomplex* alpha,
                    dcomplex* beta, dcomplex* q, integer* ldq,
                    dcomplex* z, integer* ldz, dcomplex* work, integer* lwork,
                    double* rwork, integer *info );

void FFUNC(stgevc,STGEVC)( char* side, char* howmny, const logical* select,
                    integer* n, const float* s, integer* lds,
                    const float* p, integer* ldp, float* vl,
                    integer* ldvl, float* vr, integer* ldvr,
                    integer* mm, integer* m, float* work,
                    integer *info );

void FFUNC(dtgevc,DTGEVC)( char* side, char* howmny, const logical* select,
                    integer* n, const double* s, integer* lds,
                    const double* p, integer* ldp, double* vl,
                    integer* ldvl, double* vr, integer* ldvr,
                    integer* mm, integer* m, double* work,
                    integer *info );

void FFUNC(ctgevc,CTGEVC)( char* side, char* howmny, const logical* select,
                    integer* n, const scomplex* s, integer* lds,
                    const scomplex* p, integer* ldp, scomplex* vl,
                    integer* ldvl, scomplex* vr, integer* ldvr,
                    integer* mm, integer* m, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(ztgevc,ZTGEVC)( char* side, char* howmny, const logical* select,
                    integer* n, const dcomplex* s, integer* lds,
                    const dcomplex* p, integer* ldp, dcomplex* vl,
                    integer* ldvl, dcomplex* vr, integer* ldvr,
                    integer* mm, integer* m, dcomplex* work, double* rwork,
                    integer *info );

void FFUNC(stgexc,STGEXC)( logical* wantq, logical* wantz, integer* n,
                    float* a, integer* lda, float* b, integer* ldb,
                    float* q, integer* ldq, float* z, integer* ldz,
                    integer* ifst, integer* ilst, float* work,
                    integer* lwork, integer *info );

void FFUNC(dtgexc,DTGEXC)( logical* wantq, logical* wantz, integer* n,
                    double* a, integer* lda, double* b, integer* ldb,
                    double* q, integer* ldq, double* z, integer* ldz,
                    integer* ifst, integer* ilst, double* work,
                    integer* lwork, integer *info );

void FFUNC(ctgexc,CTGEXC)( logical* wantq, logical* wantz, integer* n,
                    scomplex* a, integer* lda, scomplex* b, integer* ldb,
                    scomplex* q, integer* ldq, scomplex* z, integer* ldz,
                    integer* ifst, integer* ilst, integer *info );

void FFUNC(ztgexc,ZTGEXC)( logical* wantq, logical* wantz, integer* n,
                    dcomplex* a, integer* lda, dcomplex* b, integer* ldb,
                    dcomplex* q, integer* ldq, dcomplex* z, integer* ldz,
                    integer* ifst, integer* ilst, integer *info );

void FFUNC(stgsen,STGSEN)( integer* ijob, logical* wantq,
                    logical* wantz, const logical* select,
                    integer* n, float* a, integer* lda, float* b,
                    integer* ldb, float* alphar, float* alphai,
                    float* beta, float* q, integer* ldq, float* z,
                    integer* ldz, integer* m, float* pl, float* pr,
                    float* dif, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dtgsen,DTGSEN)( integer* ijob, logical* wantq,
                    logical* wantz, const logical* select,
                    integer* n, double* a, integer* lda, double* b,
                    integer* ldb, double* alphar, double* alphai,
                    double* beta, double* q, integer* ldq, double* z,
                    integer* ldz, integer* m, double* pl, double* pr,
                    double* dif, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(ctgsen,CTGSEN)( integer* ijob, logical* wantq,
                    logical* wantz, const logical* select,
                    integer* n, scomplex* a, integer* lda, scomplex* b,
                    integer* ldb, scomplex* alpha,
                    scomplex* beta, scomplex* q, integer* ldq, scomplex* z,
                    integer* ldz, integer* m, float* pl, float* pr,
                    float* dif, scomplex* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(ztgsen,ZTGSEN)( integer* ijob, logical* wantq,
                    logical* wantz, const logical* select,
                    integer* n, dcomplex* a, integer* lda, dcomplex* b,
                    integer* ldb, dcomplex* alpha,
                    dcomplex* beta, dcomplex* q, integer* ldq, dcomplex* z,
                    integer* ldz, integer* m, double* pl, double* pr,
                    double* dif, dcomplex* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(stgsyl,STGSYL)( char* trans, integer* ijob, integer* m, integer* n,
                    const float* a, integer* lda, const float* b,
                    integer* ldb, float* c, integer* ldc,
                    const float* d, integer* ldd, const float* e,
                    integer* lde, float* f, integer* ldf, float* scale,
                    float* dif, float* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(dtgsyl,DTGSYL)( char* trans, integer* ijob, integer* m, integer* n,
                    const double* a, integer* lda, const double* b,
                    integer* ldb, double* c, integer* ldc,
                    const double* d, integer* ldd, const double* e,
                    integer* lde, double* f, integer* ldf, double* scale,
                    double* dif, double* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(ctgsyl,CTGSYL)( char* trans, integer* ijob, integer* m, integer* n,
                    const scomplex* a, integer* lda, const scomplex* b,
                    integer* ldb, scomplex* c, integer* ldc,
                    const scomplex* d, integer* ldd, const scomplex* e,
                    integer* lde, scomplex* f, integer* ldf, float* scale,
                    float* dif, scomplex* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(ztgsyl,ZTGSYL)( char* trans, integer* ijob, integer* m, integer* n,
                    const dcomplex* a, integer* lda, const dcomplex* b,
                    integer* ldb, dcomplex* c, integer* ldc,
                    const dcomplex* d, integer* ldd, const dcomplex* e,
                    integer* lde, dcomplex* f, integer* ldf, double* scale,
                    double* dif, dcomplex* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(stgsna,STGSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const float* a, integer* lda,
                    const float* b, integer* ldb, const float* vl,
                    integer* ldvl, const float* vr, integer* ldvr,
                    float* s, float* dif, integer* mm, integer* m,
                    float* work, integer* lwork, integer* iwork,
                    integer *info );

void FFUNC(dtgsna,DTGSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const double* a, integer* lda,
                    const double* b, integer* ldb, const double* vl,
                    integer* ldvl, const double* vr, integer* ldvr,
                    double* s, double* dif, integer* mm, integer* m,
                    double* work, integer* lwork, integer* iwork,
                    integer *info );

void FFUNC(ctgsna,CTGSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const scomplex* a, integer* lda,
                    const scomplex* b, integer* ldb, const scomplex* vl,
                    integer* ldvl, const scomplex* vr, integer* ldvr,
                    float* s, float* dif, integer* mm, integer* m,
                    scomplex* work, integer* lwork, integer* iwork,
                    integer *info );

void FFUNC(ztgsna,ZTGSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const dcomplex* a, integer* lda,
                    const dcomplex* b, integer* ldb, const dcomplex* vl,
                    integer* ldvl, const dcomplex* vr, integer* ldvr,
                    double* s, double* dif, integer* mm, integer* m,
                    dcomplex* work, integer* lwork, integer* iwork,
                    integer *info );

void FFUNC(sggsvp,SGGSVP)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, float* a, integer* lda,
                    float* b, integer* ldb, float* tola, float* tolb,
                    integer* k, integer* l, float* u, integer* ldu,
                    float* v, integer* ldv, float* q, integer* ldq,
                    integer* iwork, float* tau, float* work,
                    integer *info );

void FFUNC(dggsvp,DGGSVP)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, double* a, integer* lda,
                    double* b, integer* ldb, double* tola, double* tolb,
                    integer* k, integer* l, double* u, integer* ldu,
                    double* v, integer* ldv, double* q, integer* ldq,
                    integer* iwork, double* tau, double* work,
                    integer *info );

void FFUNC(cggsvp,CGGSVP)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, scomplex* a, integer* lda,
                    scomplex* b, integer* ldb, float* tola, float* tolb,
                    integer* k, integer* l, scomplex* u, integer* ldu,
                    scomplex* v, integer* ldv, scomplex* q, integer* ldq,
                    integer* iwork, float* rwork, scomplex* tau, scomplex* work,
                    integer *info );

void FFUNC(zggsvp,ZGGSVP)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, dcomplex* a, integer* lda,
                    dcomplex* b, integer* ldb, double* tola, double* tolb,
                    integer* k, integer* l, dcomplex* u, integer* ldu,
                    dcomplex* v, integer* ldv, dcomplex* q, integer* ldq,
                    integer* iwork, double* rwork, dcomplex* tau, dcomplex* work,
                    integer *info );

void FFUNC(stgsja,STGSJA)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, integer* k, integer* l,
                    float* a, integer* lda, float* b, integer* ldb,
                    float* tola, float* tolb, float* alpha, float* beta,
                    float* u, integer* ldu, float* v, integer* ldv,
                    float* q, integer* ldq, float* work,
                    integer* ncycle, integer *info );

void FFUNC(dtgsja,DTGSJA)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, integer* k, integer* l,
                    double* a, integer* lda, double* b, integer* ldb,
                    double* tola, double* tolb, double* alpha, double* beta,
                    double* u, integer* ldu, double* v, integer* ldv,
                    double* q, integer* ldq, double* work,
                    integer* ncycle, integer *info );

void FFUNC(ctgsja,CTGSJA)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, integer* k, integer* l,
                    scomplex* a, integer* lda, scomplex* b, integer* ldb,
                    float* tola, float* tolb, float* alpha, float* beta,
                    scomplex* u, integer* ldu, scomplex* v, integer* ldv,
                    scomplex* q, integer* ldq, scomplex* work,
                    integer* ncycle, integer *info );

void FFUNC(ztgsja,ZTGSJA)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, integer* k, integer* l,
                    dcomplex* a, integer* lda, dcomplex* b, integer* ldb,
                    double* tola, double* tolb, double* alpha, double* beta,
                    dcomplex* u, integer* ldu, dcomplex* v, integer* ldv,
                    dcomplex* q, integer* ldq, dcomplex* work,
                    integer* ncycle, integer *info );

void FFUNC(sgels,SGELS)( char* trans, integer* m, integer* n, integer* nrhs,
                   float* a, integer* lda, float* b, integer* ldb,
                   float* work, integer* lwork, integer *info );

void FFUNC(dgels,DGELS)( char* trans, integer* m, integer* n, integer* nrhs,
                   double* a, integer* lda, double* b, integer* ldb,
                   double* work, integer* lwork, integer *info );

void FFUNC(cgels,CGELS)( char* trans, integer* m, integer* n, integer* nrhs,
                    scomplex* a, integer* lda, scomplex* b, integer* ldb,
                    scomplex* work, integer* lwork, integer *info );

void FFUNC(zgels,ZGELS)( char* trans, integer* m, integer* n, integer* nrhs,
                   dcomplex* a, integer* lda, dcomplex* b, integer* ldb,
                   dcomplex* work, integer* lwork, integer *info );

void FFUNC(sgelsy,SGELSY)( integer* m, integer* n, integer* nrhs, float* a,
                    integer* lda, float* b, integer* ldb,
                    integer* jpvt, float* rcond, integer* rank,
                    float* work, integer* lwork, integer *info );

void FFUNC(dgelsy,DGELSY)( integer* m, integer* n, integer* nrhs, double* a,
                    integer* lda, double* b, integer* ldb,
                    integer* jpvt, double* rcond, integer* rank,
                    double* work, integer* lwork, integer *info );

void FFUNC(cgelsy,CGELSY)( integer* m, integer* n, integer* nrhs, scomplex* a,
                    integer* lda, scomplex* b, integer* ldb,
                    integer* jpvt, float* rcond, integer* rank,
                    scomplex* work, integer* lwork, float* rwork, integer *info );

void FFUNC(zgelsy,ZGELSY)( integer* m, integer* n, integer* nrhs, dcomplex* a,
                    integer* lda, dcomplex* b, integer* ldb,
                    integer* jpvt, double* rcond, integer* rank,
                    dcomplex* work, integer* lwork, double* rwork, integer *info );

void FFUNC(sgelss,SGELSS)( integer* m, integer* n, integer* nrhs, float* a,
                    integer* lda, float* b, integer* ldb, float* s,
                    float* rcond, integer* rank, float* work,
                    integer* lwork, integer *info );

void FFUNC(dgelss,DGELSS)( integer* m, integer* n, integer* nrhs, double* a,
                    integer* lda, double* b, integer* ldb, double* s,
                    double* rcond, integer* rank, double* work,
                    integer* lwork, integer *info );

void FFUNC(cgelss,CGELSS)( integer* m, integer* n, integer* nrhs, scomplex* a,
                    integer* lda, scomplex* b, integer* ldb, float* s,
                    float* rcond, integer* rank, scomplex* work,
                    integer* lwork, float* rwork, integer *info );

void FFUNC(zgelss,ZGELSS)( integer* m, integer* n, integer* nrhs, dcomplex* a,
                    integer* lda, dcomplex* b, integer* ldb, double* s,
                    double* rcond, integer* rank, dcomplex* work,
                    integer* lwork, double* rwork, integer *info );

void FFUNC(sgglse,SGGLSE)( integer* m, integer* n, integer* p, float* a,
                    integer* lda, float* b, integer* ldb, float* c,
                    float* d, float* x, float* work, integer* lwork,
                    integer *info );

void FFUNC(dgglse,DGGLSE)( integer* m, integer* n, integer* p, double* a,
                    integer* lda, double* b, integer* ldb, double* c,
                    double* d, double* x, double* work, integer* lwork,
                    integer *info );

void FFUNC(cgglse,CGGLSE)( integer* m, integer* n, integer* p, scomplex* a,
                    integer* lda, scomplex* b, integer* ldb, scomplex* c,
                    scomplex* d, scomplex* x, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zgglse,ZGGLSE)( integer* m, integer* n, integer* p, dcomplex* a,
                    integer* lda, dcomplex* b, integer* ldb, dcomplex* c,
                    dcomplex* d, dcomplex* x, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(sggglm,SGGGLM)( integer* n, integer* m, integer* p, float* a,
                    integer* lda, float* b, integer* ldb, float* d,
                    float* x, float* y, float* work, integer* lwork,
                    integer *info );

void FFUNC(dggglm,DGGGLM)( integer* n, integer* m, integer* p, double* a,
                    integer* lda, double* b, integer* ldb, double* d,
                    double* x, double* y, double* work, integer* lwork,
                    integer *info );

void FFUNC(cggglm,CGGGLM)( integer* n, integer* m, integer* p, scomplex* a,
                    integer* lda, scomplex* b, integer* ldb, scomplex* d,
                    scomplex* x, scomplex* y, scomplex* work, integer* lwork,
                    integer *info );

void FFUNC(zggglm,ZGGGLM)( integer* n, integer* m, integer* p, dcomplex* a,
                    integer* lda, dcomplex* b, integer* ldb, dcomplex* d,
                    dcomplex* x, dcomplex* y, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(ssyev,SSYEV)( char* jobz, char* uplo, integer* n, float* a,
                   integer* lda, float* w, float* work, integer* lwork,
                   integer *info );

void FFUNC(dsyev,DSYEV)( char* jobz, char* uplo, integer* n, double* a,
                   integer* lda, double* w, double* work, integer* lwork,
                   integer *info );

void FFUNC(cheev,CHEEV)( char* jobz, char* uplo, integer* n, scomplex* a,
                   integer* lda, float* w, scomplex* work, integer* lwork,
                   float* rwork, integer *info );

void FFUNC(zheev,ZHEEV)( char* jobz, char* uplo, integer* n, dcomplex* a,
                   integer* lda, double* w, dcomplex* work, integer* lwork,
                   double* rwork, integer *info );

void FFUNC(ssyevd,SSYEVD)( char* jobz, char* uplo, integer* n, float* a,
                    integer* lda, float* w, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dsyevd,DSYEVD)( char* jobz, char* uplo, integer* n, double* a,
                    integer* lda, double* w, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(cheevd,CHEEVD)( char* jobz, char* uplo, integer* n, scomplex* a,
                    integer* lda, float* w, scomplex* work, integer* lwork,
                    float* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(zheevd,ZHEEVD)( char* jobz, char* uplo, integer* n, dcomplex* a,
                    integer* lda, double* w, dcomplex* work, integer* lwork,
                    double* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(ssyevx,SSYEVX)( char* jobz, char* range, char* uplo, integer* n,
                    float* a, integer* lda, float* vl, float* vu,
                    integer* il, integer* iu, float* abstol,
                    integer* m, float* w, float* z, integer* ldz,
                    float* work, integer* lwork, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dsyevx,DSYEVX)( char* jobz, char* range, char* uplo, integer* n,
                    double* a, integer* lda, double* vl, double* vu,
                    integer* il, integer* iu, double* abstol,
                    integer* m, double* w, double* z, integer* ldz,
                    double* work, integer* lwork, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(cheevx,CHEEVX)( char* jobz, char* range, char* uplo, integer* n,
                    scomplex* a, integer* lda, float* vl, float* vu,
                    integer* il, integer* iu, float* abstol,
                    integer* m, float* w, scomplex* z, integer* ldz,
                    scomplex* work, integer* lwork, float* rwork, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(zheevx,ZHEEVX)( char* jobz, char* range, char* uplo, integer* n,
                    dcomplex* a, integer* lda, double* vl, double* vu,
                    integer* il, integer* iu, double* abstol,
                    integer* m, double* w, dcomplex* z, integer* ldz,
                    dcomplex* work, integer* lwork, double* rwork, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(ssyevr,SSYEVR)( char* jobz, char* range, char* uplo, integer* n,
                    float* a, integer* lda, float* vl, float* vu,
                    integer* il, integer* iu, float* abstol,
                    integer* m, float* w, float* z, integer* ldz,
                    integer* isuppz, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dsyevr,DSYEVR)( char* jobz, char* range, char* uplo, integer* n,
                    double* a, integer* lda, double* vl, double* vu,
                    integer* il, integer* iu, double* abstol,
                    integer* m, double* w, double* z, integer* ldz,
                    integer* isuppz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(cheevr,CHEEVR)( char* jobz, char* range, char* uplo, integer* n,
                    scomplex* a, integer* lda, float* vl, float* vu,
                    integer* il, integer* iu, float* abstol,
                    integer* m, float* w, scomplex* z, integer* ldz,
                    integer* isuppz, scomplex* work, integer* lwork,
                    float* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(zheevr,ZHEEVR)( char* jobz, char* range, char* uplo, integer* n,
                    dcomplex* a, integer* lda, double* vl, double* vu,
                    integer* il, integer* iu, double* abstol,
                    integer* m, double* w, dcomplex* z, integer* ldz,
                    integer* isuppz, dcomplex* work, integer* lwork,
                    double* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(sspev,SSPEV)( char* jobz, char* uplo, integer* n, float* ap, float* w,
                    float* z, integer* ldz, float* work, integer *info );

void FFUNC(dspev,DSPEV)( char* jobz, char* uplo, integer* n, double* ap, double* w,
                   double* z, integer* ldz, double* work, integer *info );

void FFUNC(chpev,CHPEV)( char* jobz, char* uplo, integer* n, scomplex* ap, float* w,
                   scomplex* z, integer* ldz, scomplex* work, float* rwork, integer *info );

void FFUNC(zhpev,ZHPEV)( char* jobz, char* uplo, integer* n, dcomplex* ap, double* w,
                   dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer *info );

void FFUNC(sspevd,SSPEVD)( char* jobz, char* uplo, integer* n, float* ap,
                    float* w, float* z, integer* ldz, float* work,
                    integer* lwork, integer* iwork, integer* liwork,
                    integer *info );

void FFUNC(dspevd,DSPEVD)( char* jobz, char* uplo, integer* n, double* ap,
                    double* w, double* z, integer* ldz, double* work,
                    integer* lwork, integer* iwork, integer* liwork,
                    integer *info );

void FFUNC(chpevd,CHPEVD)( char* jobz, char* uplo, integer* n, scomplex* ap,
                    float* w, scomplex* z, integer* ldz, scomplex* work,
                    integer* lwork, float* rwork, integer* lrwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(zhpevd,ZHPEVD)( char* jobz, char* uplo, integer* n, dcomplex* ap,
                    double* w, dcomplex* z, integer* ldz, dcomplex* work,
                    integer* lwork, double* rwork, integer* lrwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(sspevx,SSPEVX)( char* jobz, char* range, char* uplo, integer* n,
                    float* ap, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    float* z, integer* ldz, float* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dspevx,DSPEVX)( char* jobz, char* range, char* uplo, integer* n,
                    double* ap, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(chpevx,CHPEVX)( char* jobz, char* range, char* uplo, integer* n,
                    scomplex* ap, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    scomplex* z, integer* ldz, scomplex* work, float* rwork,
                    integer* iwork, integer* ifail, integer *info );

void FFUNC(zhpevx,ZHPEVX)( char* jobz, char* range, char* uplo, integer* n,
                    dcomplex* ap, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    dcomplex* z, integer* ldz, dcomplex* work, double* rwork,
                    integer* iwork, integer* ifail, integer *info );

void FFUNC(ssbev,SSBEV)( char* jobz, char* uplo, integer* n, integer* kd,
                   float* ab, integer* ldab, float* w, float* z,
                   integer* ldz, float* work, integer *info );

void FFUNC(dsbev,DSBEV)( char* jobz, char* uplo, integer* n, integer* kd,
                   double* ab, integer* ldab, double* w, double* z,
                   integer* ldz, double* work, integer *info );

void FFUNC(chbev,CHBEV)( char* jobz, char* uplo, integer* n, integer* kd,
                    scomplex* ab, integer* ldab, float* w, scomplex* z,
                   integer* ldz, scomplex* work, float* rwork, integer *info );

void FFUNC(zhbev,ZHBEV)( char* jobz, char* uplo, integer* n, integer* kd,
                    dcomplex* ab, integer* ldab, double* w, dcomplex* z,
                   integer* ldz, dcomplex* work, double* rwork, integer *info );

void FFUNC(ssbevd,SSBEVD)( char* jobz, char* uplo, integer* n, integer* kd,
                    float* ab, integer* ldab, float* w, float* z,
                    integer* ldz, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dsbevd,DSBEVD)( char* jobz, char* uplo, integer* n, integer* kd,
                    double* ab, integer* ldab, double* w, double* z,
                    integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(chbevd,CHBEVD)( char* jobz, char* uplo, integer* n, integer* kd,
                    scomplex* ab, integer* ldab, float* w, scomplex* z,
                    integer* ldz, scomplex* work, integer* lwork,
                    float* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(zhbevd,ZHBEVD)( char* jobz, char* uplo, integer* n, integer* kd,
                    dcomplex* ab, integer* ldab, double* w, dcomplex* z,
                    integer* ldz, dcomplex* work, integer* lwork,
                    double* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(ssbevx,SSBEVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* kd, float* ab, integer* ldab, float* q,
                    integer* ldq, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    float* z, integer* ldz, float* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dsbevx,DSBEVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* kd, double* ab, integer* ldab, double* q,
                    integer* ldq, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(chbevx,CHBEVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* kd, scomplex* ab, integer* ldab, scomplex* q,
                    integer* ldq, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    scomplex* z, integer* ldz, scomplex* work, float* rwork,
                    integer* iwork, integer* ifail, integer *info );

void FFUNC(zhbevx,ZHBEVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* kd, dcomplex* ab, integer* ldab, dcomplex* q,
                    integer* ldq, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    dcomplex* z, integer* ldz, dcomplex* work, double* rwork,
                    integer* iwork, integer* ifail, integer *info );

void FFUNC(sstev,SSTEV)( char* jobz, integer* n, float* d, float* e, float* z,
                   integer* ldz, float* work, integer *info );

void FFUNC(dstev,DSTEV)( char* jobz, integer* n, double* d, double* e, double* z,
                   integer* ldz, double* work, integer *info );

void FFUNC(sstevd,SSTEVD)( char* jobz, integer* n, float* d, float* e, float* z,
                    integer* ldz, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dstevd,DSTEVD)( char* jobz, integer* n, double* d, double* e, double* z,
                    integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(sstevx,SSTEVX)( char* jobz, char* range, integer* n, float* d,
                    float* e, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    float* z, integer* ldz, float* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dstevx,DSTEVX)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(sstevr,SSTEVR)( char* jobz, char* range, integer* n, float* d,
                    float* e, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    float* z, integer* ldz, integer* isuppz,
                    float* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(dstevr,DSTEVR)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, integer* isuppz,
                    double* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(sgeev,SGEEV)( char* jobvl, char* jobvr, integer* n, float* a,
                   integer* lda, float* wr, float* wi, float* vl,
                   integer* ldvl, float* vr, integer* ldvr, float* work,
                   integer* lwork, integer *info );

void FFUNC(dgeev,DGEEV)( char* jobvl, char* jobvr, integer* n, double* a,
                   integer* lda, double* wr, double* wi, double* vl,
                   integer* ldvl, double* vr, integer* ldvr, double* work,
                   integer* lwork, integer *info );

void FFUNC(cgeev,CGEEV)( char* jobvl, char* jobvr, integer* n, scomplex* a,
                   integer* lda, scomplex* w, scomplex* vl,
                   integer* ldvl, scomplex* vr, integer* ldvr, scomplex* work,
                   integer* lwork, float* rwork, integer *info );

void FFUNC(zgeev,ZGEEV)( char* jobvl, char* jobvr, integer* n, dcomplex* a,
                   integer* lda, dcomplex* w, dcomplex* vl,
                   integer* ldvl, dcomplex* vr, integer* ldvr, dcomplex* work,
                   integer* lwork, double* rwork, integer *info );

void FFUNC(sgeevx,SGEEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, float* a, integer* lda, float* wr,
                    float* wi, float* vl, integer* ldvl, float* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    float* scale, float* abnrm, float* rconde,
                    float* rcondv, float* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(dgeevx,DGEEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, double* a, integer* lda, double* wr,
                    double* wi, double* vl, integer* ldvl, double* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    double* scale, double* abnrm, double* rconde,
                    double* rcondv, double* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(cgeevx,CGEEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, scomplex* a, integer* lda, scomplex* w,
                    scomplex* vl, integer* ldvl, scomplex* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    float* scale, float* abnrm, float* rconde,
                    float* rcondv, scomplex* work, integer* lwork,
                    float* rwork, integer *info );

void FFUNC(zgeevx,ZGEEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, dcomplex* a, integer* lda, dcomplex* w,
                    dcomplex* vl, integer* ldvl, dcomplex* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    double* scale, double* abnrm, double* rconde,
                    double* rcondv, dcomplex* work, integer* lwork,
                    double* rwork, integer *info );

void FFUNC(sgesvd,SGESVD)( char* jobu, char* jobvt, integer* m, integer* n,
                    float* a, integer* lda, float* s, float* u,
                    integer* ldu, float* vt, integer* ldvt, float* work,
                    integer* lwork, integer *info );

void FFUNC(dgesvd,DGESVD)( char* jobu, char* jobvt, integer* m, integer* n,
                    double* a, integer* lda, double* s, double* u,
                    integer* ldu, double* vt, integer* ldvt, double* work,
                    integer* lwork, integer *info );

void FFUNC(cgesvd,CGESVD)( char* jobu, char* jobvt, integer* m, integer* n,
                    scomplex* a, integer* lda, float* s, scomplex* u,
                    integer* ldu, scomplex* vt, integer* ldvt, scomplex* work,
                    integer* lwork, float* rwork, integer *info );

void FFUNC(zgesvd,ZGESVD)( char* jobu, char* jobvt, integer* m, integer* n,
                    dcomplex* a, integer* lda, double* s, dcomplex* u,
                    integer* ldu, dcomplex* vt, integer* ldvt, dcomplex* work,
                    integer* lwork, double* rwork, integer *info );

void FFUNC(sgesdd,SGESDD)( char* jobz, integer* m, integer* n, float* a,
                    integer* lda, float* s, float* u, integer* ldu,
                    float* vt, integer* ldvt, float* work,
                    integer* lwork, integer* iwork, integer *info );

void FFUNC(dgesdd,DGESDD)( char* jobz, integer* m, integer* n, double* a,
                    integer* lda, double* s, double* u, integer* ldu,
                    double* vt, integer* ldvt, double* work,
                    integer* lwork, integer* iwork, integer *info );

void FFUNC(cgesdd,CGESDD)( char* jobz, integer* m, integer* n, scomplex* a,
                    integer* lda, float* s, scomplex* u, integer* ldu,
                    scomplex* vt, integer* ldvt, scomplex* work,
                    integer* lwork, float* rwork, integer* iwork, integer *info );

void FFUNC(zgesdd,ZGESDD)( char* jobz, integer* m, integer* n, dcomplex* a,
                    integer* lda, double* s, dcomplex* u, integer* ldu,
                    dcomplex* vt, integer* ldvt, dcomplex* work,
                    integer* lwork, double* rwork, integer* iwork, integer *info );

void FFUNC(sgejsv,SGEJSV)( char* joba, char* jobu, char* jobv, char* jobr, char* jobt,
                    char* jobp, integer* m, integer* n, float* a,
                    integer* lda, float* sva, float* u, integer* ldu,
                    float* v, integer* ldv, float* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(dgejsv,DGEJSV)( char* joba, char* jobu, char* jobv, char* jobr, char* jobt,
                    char* jobp, integer* m, integer* n, double* a,
                    integer* lda, double* sva, double* u, integer* ldu,
                    double* v, integer* ldv, double* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(sgesvj,SGESVJ)( char* joba, char* jobu, char* jobv, integer* m,
                    integer* n, float* a, integer* lda, float* sva,
                    integer* mv, float* v, integer* ldv, float* work,
                    integer* lwork, integer *info );

void FFUNC(dgesvj,DGESVJ)( char* joba, char* jobu, char* jobv, integer* m,
                    integer* n, double* a, integer* lda, double* sva,
                    integer* mv, double* v, integer* ldv, double* work,
                    integer* lwork, integer *info );

void FFUNC(sggsvd,SGGSVD)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* n, integer* p, integer* k, integer* l,
                    float* a, integer* lda, float* b, integer* ldb,
                    float* alpha, float* beta, float* u, integer* ldu,
                    float* v, integer* ldv, float* q, integer* ldq,
                    float* work, integer* iwork, integer *info );

void FFUNC(dggsvd,DGGSVD)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* n, integer* p, integer* k, integer* l,
                    double* a, integer* lda, double* b, integer* ldb,
                    double* alpha, double* beta, double* u, integer* ldu,
                    double* v, integer* ldv, double* q, integer* ldq,
                    double* work, integer* iwork, integer *info );

void FFUNC(cggsvd,CGGSVD)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* n, integer* p, integer* k, integer* l,
                    scomplex* a, integer* lda, scomplex* b, integer* ldb,
                    float* alpha, float* beta, scomplex* u, integer* ldu,
                    scomplex* v, integer* ldv, scomplex* q, integer* ldq,
                    scomplex* work, float* rwork, integer* iwork, integer *info );

void FFUNC(zggsvd,ZGGSVD)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* n, integer* p, integer* k, integer* l,
                    dcomplex* a, integer* lda, dcomplex* b, integer* ldb,
                    double* alpha, double* beta, dcomplex* u, integer* ldu,
                    dcomplex* v, integer* ldv, dcomplex* q, integer* ldq,
                    dcomplex* work, double* rwork, integer* iwork, integer *info );

void FFUNC(ssygv,SSYGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   float* a, integer* lda, float* b, integer* ldb,
                   float* w, float* work, integer* lwork,
                   integer *info );

void FFUNC(dsygv,DSYGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   double* a, integer* lda, double* b, integer* ldb,
                   double* w, double* work, integer* lwork,
                   integer *info );

void FFUNC(chegv,CHEGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   scomplex* a, integer* lda, scomplex* b, integer* ldb,
                   float* w, scomplex* work, integer* lwork,
                   float* rwork, integer *info );

void FFUNC(zhegv,ZHEGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   dcomplex* a, integer* lda, dcomplex* b, integer* ldb,
                   double* w, dcomplex* work, integer* lwork,
                   double* rwork, integer *info );

void FFUNC(ssygvd,SSYGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    float* a, integer* lda, float* b, integer* ldb,
                    float* w, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dsygvd,DSYGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    double* a, integer* lda, double* b, integer* ldb,
                    double* w, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(chegvd,CHEGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    scomplex* a, integer* lda, scomplex* b, integer* ldb,
                    float* w, scomplex* work, integer* lwork,
                    float* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(zhegvd,ZHEGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    dcomplex* a, integer* lda, dcomplex* b, integer* ldb,
                    double* w, dcomplex* work, integer* lwork,
                    double* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(ssygvx,SSYGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, float* a, integer* lda, float* b,
                    integer* ldb, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    float* z, integer* ldz, float* work, integer* lwork,
                    integer* iwork, integer* ifail, integer *info );

void FFUNC(dsygvx,DSYGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, double* a, integer* lda, double* b,
                    integer* ldb, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* ifail, integer *info );

void FFUNC(chegvx,CHEGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, scomplex* a, integer* lda, scomplex* b,
                    integer* ldb, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    scomplex* z, integer* ldz, scomplex* work, integer* lwork,
                    float* rwork, integer* iwork, integer* ifail, integer *info );

void FFUNC(zhegvx,ZHEGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, dcomplex* a, integer* lda, dcomplex* b,
                    integer* ldb, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    dcomplex* z, integer* ldz, dcomplex* work, integer* lwork,
                    double* rwork, integer* iwork, integer* ifail, integer *info );

void FFUNC(sspgv,SSPGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   float* ap, float* bp, float* w, float* z,
                   integer* ldz, float* work, integer *info );

void FFUNC(dspgv,DSPGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   double* ap, double* bp, double* w, double* z,
                   integer* ldz, double* work, integer *info );

void FFUNC(chpgv,CHPGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   scomplex* ap, scomplex* bp, float* w, scomplex* z,
                   integer* ldz, scomplex* work, float* rwork, integer *info );

void FFUNC(zhpgv,ZHPGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   dcomplex* ap, dcomplex* bp, double* w, dcomplex* z,
                   integer* ldz, dcomplex* work, double* rwork, integer *info );

void FFUNC(sspgvd,SSPGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    float* ap, float* bp, float* w, float* z,
                    integer* ldz, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dspgvd,DSPGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    double* ap, double* bp, double* w, double* z,
                    integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(chpgvd,CHPGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    scomplex* ap, scomplex* bp, float* w, scomplex* z,
                    integer* ldz, scomplex* work, integer* lwork,
                    float* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(zhpgvd,ZHPGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    dcomplex* ap, dcomplex* bp, double* w, dcomplex* z,
                    integer* ldz, dcomplex* work, integer* lwork,
                    double* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(sspgvx,SSPGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, float* ap, float* bp, float* vl,
                    float* vu, integer* il, integer* iu, float* abstol,
                    integer* m, float* w, float* z, integer* ldz,
                    float* work, integer* iwork, integer* ifail,
                    integer *info );

void FFUNC(dspgvx,DSPGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, double* ap, double* bp, double* vl,
                    double* vu, integer* il, integer* iu, double* abstol,
                    integer* m, double* w, double* z, integer* ldz,
                    double* work, integer* iwork, integer* ifail,
                    integer *info );

void FFUNC(chpgvx,CHPGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, scomplex* ap, scomplex* bp, float* vl,
                    float* vu, integer* il, integer* iu, float* abstol,
                    integer* m, float* w, scomplex* z, integer* ldz,
                    scomplex* work, float* rwork, integer* iwork, integer* ifail,
                    integer *info );

void FFUNC(zhpgvx,ZHPGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, dcomplex* ap, dcomplex* bp, double* vl,
                    double* vu, integer* il, integer* iu, double* abstol,
                    integer* m, double* w, dcomplex* z, integer* ldz,
                    dcomplex* work, double* rwork, integer* iwork, integer* ifail,
                    integer *info );

void FFUNC(ssbgv,SSBGV)( char* jobz, char* uplo, integer* n, integer* ka,
                   integer* kb, float* ab, integer* ldab, float* bb,
                   integer* ldbb, float* w, float* z, integer* ldz,
                   float* work, integer *info );

void FFUNC(dsbgv,DSBGV)( char* jobz, char* uplo, integer* n, integer* ka,
                   integer* kb, double* ab, integer* ldab, double* bb,
                   integer* ldbb, double* w, double* z, integer* ldz,
                   double* work, integer *info );

void FFUNC(chbgv,CHBGV)( char* jobz, char* uplo, integer* n, integer* ka,
                   integer* kb, scomplex* ab, integer* ldab, scomplex* bb,
                   integer* ldbb, float* w, scomplex* z, integer* ldz,
                   scomplex* work, float* rwork, integer *info );

void FFUNC(zhbgv,ZHBGV)( char* jobz, char* uplo, integer* n, integer* ka,
                   integer* kb, dcomplex* ab, integer* ldab, dcomplex* bb,
                   integer* ldbb, double* w, dcomplex* z, integer* ldz,
                   dcomplex* work, double* rwork, integer *info );

void FFUNC(ssbgvd,SSBGVD)( char* jobz, char* uplo, integer* n, integer* ka,
                    integer* kb, float* ab, integer* ldab, float* bb,
                    integer* ldbb, float* w, float* z, integer* ldz,
                    float* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(dsbgvd,DSBGVD)( char* jobz, char* uplo, integer* n, integer* ka,
                    integer* kb, double* ab, integer* ldab, double* bb,
                    integer* ldbb, double* w, double* z, integer* ldz,
                    double* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(chbgvd,CHBGVD)( char* jobz, char* uplo, integer* n, integer* ka,
                    integer* kb, scomplex* ab, integer* ldab, scomplex* bb,
                    integer* ldbb, float* w, scomplex* z, integer* ldz,
                    scomplex* work, integer* lwork, float* rwork, integer* lrwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(zhbgvd,ZHBGVD)( char* jobz, char* uplo, integer* n, integer* ka,
                    integer* kb, dcomplex* ab, integer* ldab, dcomplex* bb,
                    integer* ldbb, double* w, dcomplex* z, integer* ldz,
                    dcomplex* work, integer* lwork, double* rwork, integer* lrwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(ssbgvx,SSBGVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* ka, integer* kb, float* ab,
                    integer* ldab, float* bb, integer* ldbb, float* q,
                    integer* ldq, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    float* z, integer* ldz, float* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dsbgvx,DSBGVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* ka, integer* kb, double* ab,
                    integer* ldab, double* bb, integer* ldbb, double* q,
                    integer* ldq, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(chbgvx,CHBGVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* ka, integer* kb, scomplex* ab,
                    integer* ldab, scomplex* bb, integer* ldbb, scomplex* q,
                    integer* ldq, float* vl, float* vu, integer* il,
                    integer* iu, float* abstol, integer* m, float* w,
                    scomplex* z, integer* ldz, scomplex* work, float* rwork,
                    integer* iwork, integer* ifail, integer *info );

void FFUNC(zhbgvx,ZHBGVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* ka, integer* kb, dcomplex* ab,
                    integer* ldab, dcomplex* bb, integer* ldbb, dcomplex* q,
                    integer* ldq, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    dcomplex* z, integer* ldz, dcomplex* work, double* rwork,
                    integer* iwork, integer* ifail, integer *info );

void FFUNC(sggev,SGGEV)( char* jobvl, char* jobvr, integer* n, float* a,
                   integer* lda, float* b, integer* ldb, float* alphar,
                   float* alphai, float* beta, float* vl, integer* ldvl,
                   float* vr, integer* ldvr, float* work,
                   integer* lwork, integer *info );

void FFUNC(dggev,DGGEV)( char* jobvl, char* jobvr, integer* n, double* a,
                   integer* lda, double* b, integer* ldb, double* alphar,
                   double* alphai, double* beta, double* vl, integer* ldvl,
                   double* vr, integer* ldvr, double* work,
                   integer* lwork, integer *info );

void FFUNC(cggev,CGGEV)( char* jobvl, char* jobvr, integer* n, scomplex* a,
                   integer* lda, scomplex* b, integer* ldb, scomplex* alpha,
                   scomplex* beta, scomplex* vl, integer* ldvl,
                   scomplex* vr, integer* ldvr, scomplex* work,
                   integer* lwork, float* rwork, integer *info );

void FFUNC(zggev,ZGGEV)( char* jobvl, char* jobvr, integer* n, dcomplex* a,
                   integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha,
                   dcomplex* beta, dcomplex* vl, integer* ldvl,
                   dcomplex* vr, integer* ldvr, dcomplex* work,
                   integer* lwork, double* rwork, integer *info );

void FFUNC(sggevx,SGGEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, float* a, integer* lda, float* b,
                    integer* ldb, float* alphar, float* alphai,
                    float* beta, float* vl, integer* ldvl, float* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    float* lscale, float* rscale, float* abnrm,
                    float* bbnrm, float* rconde, float* rcondv, float* work,
                    integer* lwork, integer* iwork, logical* bwork,
                    integer *info );

void FFUNC(dggevx,DGGEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, double* a, integer* lda, double* b,
                    integer* ldb, double* alphar, double* alphai,
                    double* beta, double* vl, integer* ldvl, double* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    double* lscale, double* rscale, double* abnrm,
                    double* bbnrm, double* rconde, double* rcondv, double* work,
                    integer* lwork, integer* iwork, logical* bwork,
                    integer *info );

void FFUNC(cggevx,CGGEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, scomplex* a, integer* lda, scomplex* b,
                    integer* ldb, scomplex* alpha,
                    scomplex* beta, scomplex* vl, integer* ldvl, scomplex* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    float* lscale, float* rscale, float* abnrm,
                    float* bbnrm, float* rconde, float* rcondv, scomplex* work,
                    integer* lwork, float* rwork, integer* iwork, logical* bwork,
                    integer *info );

void FFUNC(zggevx,ZGGEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, dcomplex* a, integer* lda, dcomplex* b,
                    integer* ldb, dcomplex* alpha,
                    dcomplex* beta, dcomplex* vl, integer* ldvl, dcomplex* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    double* lscale, double* rscale, double* abnrm,
                    double* bbnrm, double* rconde, double* rcondv, dcomplex* work,
                    integer* lwork, double* rwork, integer* iwork, logical* bwork,
                    integer *info );

void FFUNC(ssfrk,SSFRK)( char* transr, char* uplo, char* trans, integer* n,
                   integer* k, float* alpha, const float* a,
                   integer* lda, float* beta, float* c );

void FFUNC(dsfrk,DSFRK)( char* transr, char* uplo, char* trans, integer* n,
                   integer* k, double* alpha, const double* a,
                   integer* lda, double* beta, double* c );

void FFUNC(chfrk,CHFRK)( char* transr, char* uplo, char* trans, integer* n,
                   integer* k, float* alpha, const scomplex* a,
                   integer* lda, float* beta, scomplex* c );

void FFUNC(zhfrk,ZHFRK)( char* transr, char* uplo, char* trans, integer* n,
                   integer* k, double* alpha, const dcomplex* a,
                   integer* lda, double* beta, dcomplex* c );

void FFUNC(stfsm,STFSM)( char* transr, char* side, char* uplo, char* trans,
                   char* diag, integer* m, integer* n, float* alpha,
                   const float* a, float* b, integer* ldb );

void FFUNC(dtfsm,DTFSM)( char* transr, char* side, char* uplo, char* trans,
                   char* diag, integer* m, integer* n, double* alpha,
                   const double* a, double* b, integer* ldb );

void FFUNC(ctfsm,CTFSM)( char* transr, char* side, char* uplo, char* trans,
                   char* diag, integer* m, integer* n, scomplex* alpha,
                   const scomplex* a, scomplex* b, integer* ldb );

void FFUNC(ztfsm,ZTFSM)( char* transr, char* side, char* uplo, char* trans,
                   char* diag, integer* m, integer* n, dcomplex* alpha,
                   const dcomplex* a, dcomplex* b, integer* ldb );

void FFUNC(stfttp,STFTTP)( char* transr, char* uplo, integer* n, const float* arf,
                    float* ap, integer *info );

void FFUNC(dtfttp,DTFTTP)( char* transr, char* uplo, integer* n, const double* arf,
                    double* ap, integer *info );

void FFUNC(ctfttp,CTFTTP)( char* transr, char* uplo, integer* n, const scomplex* arf,
                    scomplex* ap, integer *info );

void FFUNC(ztfttp,ZTFTTP)( char* transr, char* uplo, integer* n, const dcomplex* arf,
                    dcomplex* ap, integer *info );

void FFUNC(stfttr,STFTTR)( char* transr, char* uplo, integer* n, const float* arf,
                    float* a, integer* lda, integer *info );

void FFUNC(dtfttr,DTFTTR)( char* transr, char* uplo, integer* n, const double* arf,
                    double* a, integer* lda, integer *info );

void FFUNC(ctfttr,CTFTTR)( char* transr, char* uplo, integer* n, const scomplex* arf,
                    scomplex* a, integer* lda, integer *info );

void FFUNC(ztfttr,ZTFTTR)( char* transr, char* uplo, integer* n, const dcomplex* arf,
                    dcomplex* a, integer* lda, integer *info );

void FFUNC(stpttf,STPTTF)( char* transr, char* uplo, integer* n, const float* ap,
                    float* arf, integer *info );

void FFUNC(dtpttf,DTPTTF)( char* transr, char* uplo, integer* n, const double* ap,
                    double* arf, integer *info );

void FFUNC(ctpttf,CTPTTF)( char* transr, char* uplo, integer* n, const scomplex* ap,
                    scomplex* arf, integer *info );

void FFUNC(ztpttf,ZTPTTF)( char* transr, char* uplo, integer* n, const dcomplex* ap,
                    dcomplex* arf, integer *info );

void FFUNC(stpttr,STPTTR)( char* uplo, integer* n, const float* ap, float* a,
                    integer* lda, integer *info );

void FFUNC(dtpttr,DTPTTR)( char* uplo, integer* n, const double* ap, double* a,
                    integer* lda, integer *info );

void FFUNC(ctpttr,CTPTTR)( char* uplo, integer* n, const scomplex* ap, scomplex* a,
                    integer* lda, integer *info );

void FFUNC(ztpttr,ZTPTTR)( char* uplo, integer* n, const dcomplex* ap, dcomplex* a,
                    integer* lda, integer *info );

void FFUNC(strttf,STRTTF)( char* transr, char* uplo, integer* n, const float* a,
                    integer* lda, float* arf, integer *info );

void FFUNC(dtrttf,DTRTTF)( char* transr, char* uplo, integer* n, const double* a,
                    integer* lda, double* arf, integer *info );

void FFUNC(ctrttf,CTRTTF)( char* transr, char* uplo, integer* n, const scomplex* a,
                    integer* lda, scomplex* arf, integer *info );

void FFUNC(ztrttf,ZTRTTF)( char* transr, char* uplo, integer* n, const dcomplex* a,
                    integer* lda, dcomplex* arf, integer *info );

void FFUNC(strttp,STRTTP)( char* uplo, integer* n, const float* a, integer* lda,
                    float* ap, integer *info );

void FFUNC(dtrttp,DTRTTP)( char* uplo, integer* n, const double* a, integer* lda,
                    double* ap, integer *info );

void FFUNC(ctrttp,CTRTTP)( char* uplo, integer* n, const scomplex* a, integer* lda,
                    scomplex* ap, integer *info );

void FFUNC(ztrttp,ZTRTTP)( char* uplo, integer* n, const dcomplex* a, integer* lda,
                    dcomplex* ap, integer *info );

void FFUNC(sgeqrfp,SGEQRFP)( integer* m, integer* n, float* a, integer* lda,
                     float* tau, float* work, integer* lwork,
                     integer *info );

void FFUNC(dgeqrfp,DGEQRFP)( integer* m, integer* n, double* a, integer* lda,
                     double* tau, double* work, integer* lwork,
                     integer *info );

void FFUNC(cgeqrfp,CGEQRFP)( integer* m, integer* n, scomplex* a, integer* lda,
                     scomplex* tau, scomplex* work, integer* lwork,
                     integer *info );

void FFUNC(zgeqrfp,ZGEQRFP)( integer* m, integer* n, dcomplex* a, integer* lda,
                     dcomplex* tau, dcomplex* work, integer* lwork,
                     integer *info );

/*
 * C wrappers
 */
static inline integer c_sgetrf( integer m, integer n, float* a, integer lda,
                    integer* ipiv )
{
    integer info;
    FFUNC(sgetrf,SGETRF)(&m, &n, a, &lda, ipiv, &info);
    return info;
}

static inline integer c_dgetrf( integer m, integer n, double* a, integer lda,
                    integer* ipiv )
{
    integer info;
    FFUNC(dgetrf,DGETRF)(&m, &n, a, &lda, ipiv, &info);
    return info;
}

static inline integer c_cgetrf( integer m, integer n, scomplex* a, integer lda,
                    integer* ipiv )
{
    integer info;
    FFUNC(cgetrf,CGETRF)(&m, &n, a, &lda, ipiv, &info);
    return info;
}

static inline integer c_zgetrf( integer m, integer n, dcomplex* a, integer lda,
                    integer* ipiv )
{
    integer info;
    FFUNC(zgetrf,ZGETRF)(&m, &n, a, &lda, ipiv, &info);
    return info;
}

static inline integer c_sgbtrf( integer m, integer n, integer kl,
                    integer ku, float* ab, integer ldab,
                    integer* ipiv )
{
    integer info;
    FFUNC(sgbtrf,SGBTRF)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

static inline integer c_dgbtrf( integer m, integer n, integer kl,
                    integer ku, double* ab, integer ldab,
                    integer* ipiv )
{
    integer info;
    FFUNC(dgbtrf,DGBTRF)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

static inline integer c_cgbtrf( integer m, integer n, integer kl,
                    integer ku, scomplex* ab, integer ldab,
                    integer* ipiv )
{
    integer info;
    FFUNC(cgbtrf,CGBTRF)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

static inline integer c_zgbtrf( integer m, integer n, integer kl,
                    integer ku, dcomplex* ab, integer ldab,
                    integer* ipiv )
{
    integer info;
    FFUNC(zgbtrf,ZGBTRF)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

static inline integer c_sgttrf( integer n, float* dl, float* d, float* du,
                    float* du2, integer* ipiv )
{
    integer info;
    FFUNC(sgttrf,SGTTRF)(&n, dl, d, du, du2, ipiv, &info);
    return info;
}

static inline integer c_dgttrf( integer n, double* dl, double* d, double* du,
                    double* du2, integer* ipiv )
{
    integer info;
    FFUNC(dgttrf,DGTTRF)(&n, dl, d, du, du2, ipiv, &info);
    return info;
}

static inline integer c_cgttrf( integer n, scomplex* dl, scomplex* d, scomplex* du,
                    scomplex* du2, integer* ipiv )
{
    integer info;
    FFUNC(cgttrf,CGTTRF)(&n, dl, d, du, du2, ipiv, &info);
    return info;
}

static inline integer c_zgttrf( integer n, dcomplex* dl, dcomplex* d, dcomplex* du,
                    dcomplex* du2, integer* ipiv )
{
    integer info;
    FFUNC(zgttrf,ZGTTRF)(&n, dl, d, du, du2, ipiv, &info);
    return info;
}

static inline integer c_spotrf( char uplo, integer n, float* a, integer lda )
{
    integer info;
    FFUNC(spotrf,SPOTRF)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_dpotrf( char uplo, integer n, double* a, integer lda )
{
    integer info;
    FFUNC(dpotrf,DPOTRF)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_cpotrf( char uplo, integer n, scomplex* a, integer lda )
{
    integer info;
    FFUNC(cpotrf,CPOTRF)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_zpotrf( char uplo, integer n, dcomplex* a, integer lda )
{
    integer info;
    FFUNC(zpotrf,ZPOTRF)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_spstrf( char uplo, integer n, float* a, integer lda,
                    integer* piv, integer* rank, float tol)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    FFUNC(spstrf,SPSTRF)(&uplo, &n, a, &lda, piv, rank, &tol, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dpstrf( char uplo, integer n, double* a, integer lda,
                    integer* piv, integer* rank, double tol)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    FFUNC(dpstrf,DPSTRF)(&uplo, &n, a, &lda, piv, rank, &tol, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cpstrf( char uplo, integer n, scomplex* a, integer lda,
                    integer* piv, integer* rank, float tol)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    FFUNC(cpstrf,CPSTRF)(&uplo, &n, a, &lda, piv, rank, &tol, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zpstrf( char uplo, integer n, dcomplex* a, integer lda,
                    integer* piv, integer* rank, double tol)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    FFUNC(zpstrf,ZPSTRF)(&uplo, &n, a, &lda, piv, rank, &tol, work, &info);
    FREE(work);
    return info;
}

static inline integer c_spftrf( char transr, char uplo, integer n, float* a )
{
    integer info;
    FFUNC(spftrf,SPFTRF)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_dpftrf( char transr, char uplo, integer n, double* a )
{
    integer info;
    FFUNC(dpftrf,DPFTRF)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_cpftrf( char transr, char uplo, integer n, scomplex* a )
{
    integer info;
    FFUNC(cpftrf,CPFTRF)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_zpftrf( char transr, char uplo, integer n, dcomplex* a )
{
    integer info;
    FFUNC(zpftrf,ZPFTRF)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_spptrf( char uplo, integer n, float* ap )
{
    integer info;
    FFUNC(spptrf,SPPTRF)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_dpptrf( char uplo, integer n, double* ap )
{
    integer info;
    FFUNC(dpptrf,DPPTRF)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_cpptrf( char uplo, integer n, scomplex* ap )
{
    integer info;
    FFUNC(cpptrf,CPPTRF)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_zpptrf( char uplo, integer n, dcomplex* ap )
{
    integer info;
    FFUNC(zpptrf,ZPPTRF)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_spbtrf( char uplo, integer n, integer kd, float* ab,
                    integer ldab )
{
    integer info;
    FFUNC(spbtrf,SPBTRF)(&uplo, &n, &kd, ab, &ldab, &info);
    return info;
}

static inline integer c_dpbtrf( char uplo, integer n, integer kd, double* ab,
                    integer ldab )
{
    integer info;
    FFUNC(dpbtrf,DPBTRF)(&uplo, &n, &kd, ab, &ldab, &info);
    return info;
}

static inline integer c_cpbtrf( char uplo, integer n, integer kd, scomplex* ab,
                    integer ldab )
{
    integer info;
    FFUNC(cpbtrf,CPBTRF)(&uplo, &n, &kd, ab, &ldab, &info);
    return info;
}

static inline integer c_zpbtrf( char uplo, integer n, integer kd, dcomplex* ab,
                    integer ldab )
{
    integer info;
    FFUNC(zpbtrf,ZPBTRF)(&uplo, &n, &kd, ab, &ldab, &info);
    return info;
}

static inline integer c_spttrf( integer n, float* d, float* e )
{
    integer info;
    FFUNC(spttrf,SPTTRF)(&n, d, e, &info);
    return info;
}

static inline integer c_dpttrf( integer n, double* d, double* e )
{
    integer info;
    FFUNC(dpttrf,DPTTRF)(&n, d, e, &info);
    return info;
}

static inline integer c_cpttrf( integer n, float* d, scomplex* e )
{
    integer info;
    FFUNC(cpttrf,CPTTRF)(&n, d, e, &info);
    return info;
}

static inline integer c_zpttrf( integer n, double* d, dcomplex* e )
{
    integer info;
    FFUNC(zpttrf,ZPTTRF)(&n, d, e, &info);
    return info;
}

static inline integer c_ssytrf( char uplo, integer n, float* a, integer lda,
                    integer* ipiv)
{
    integer info, lwork = -1;
    float query, *work;
    FFUNC(ssytrf,SSYTRF)(&uplo, &n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(ssytrf,SSYTRF)(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dsytrf( char uplo, integer n, double* a, integer lda,
                    integer* ipiv)
{
    integer info, lwork = -1;
    double query, *work;
    FFUNC(dsytrf,DSYTRF)(&uplo, &n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dsytrf,DSYTRF)(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_csytrf( char uplo, integer n, scomplex* a, integer lda,
                    integer* ipiv)
{
    integer info, lwork = -1;
    scomplex query, *work;
    FFUNC(csytrf,CSYTRF)(&uplo, &n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(csytrf,CSYTRF)(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zsytrf( char uplo, integer n, dcomplex* a, integer lda,
                    integer* ipiv)
{
    integer info, lwork = -1;
    dcomplex query, *work;
    FFUNC(zsytrf,ZSYTRF)(&uplo, &n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zsytrf,ZSYTRF)(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_chetrf( char uplo, integer n, scomplex* a, integer lda,
                    integer* ipiv)
{
    integer info, lwork = -1;
    scomplex query, *work;
    FFUNC(chetrf,CHETRF)(&uplo, &n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(chetrf,CHETRF)(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zhetrf( char uplo, integer n, dcomplex* a, integer lda,
                    integer* ipiv)
{
    integer info, lwork = -1;
    dcomplex query, *work;
    FFUNC(zhetrf,ZHETRF)(&uplo, &n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zhetrf,ZHETRF)(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_ssptrf( char uplo, integer n, float* ap, integer* ipiv )
{
    integer info;
    FFUNC(ssptrf,SSPTRF)(&uplo, &n, ap, ipiv, &info);
    return info;
}

static inline integer c_dsptrf( char uplo, integer n, double* ap, integer* ipiv )
{
    integer info;
    FFUNC(dsptrf,DSPTRF)(&uplo, &n, ap, ipiv, &info);
    return info;
}

static inline integer c_csptrf( char uplo, integer n, scomplex* ap, integer* ipiv )
{
    integer info;
    FFUNC(csptrf,CSPTRF)(&uplo, &n, ap, ipiv, &info);
    return info;
}

static inline integer c_zsptrf( char uplo, integer n, dcomplex* ap, integer* ipiv )
{
    integer info;
    FFUNC(zsptrf,ZSPTRF)(&uplo, &n, ap, ipiv, &info);
    return info;
}

static inline integer c_chptrf( char uplo, integer n, scomplex* ap, integer* ipiv )
{
    integer info;
    FFUNC(chptrf,CHPTRF)(&uplo, &n, ap, ipiv, &info);
    return info;
}

static inline integer c_zhptrf( char uplo, integer n, dcomplex* ap, integer* ipiv )
{
    integer info;
    FFUNC(zhptrf,ZHPTRF)(&uplo, &n, ap, ipiv, &info);
    return info;
}

static inline integer c_sgetrs( char trans, integer n, integer nrhs,
                    const float* a, integer lda, const integer* ipiv,
                    float* b, integer ldb )
{
    integer info;
    FFUNC(sgetrs,SGETRS)(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dgetrs( char trans, integer n, integer nrhs,
                    const double* a, integer lda, const integer* ipiv,
                    double* b, integer ldb )
{
    integer info;
    FFUNC(dgetrs,DGETRS)(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_cgetrs( char trans, integer n, integer nrhs,
                    const scomplex* a, integer lda, const integer* ipiv,
                    scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cgetrs,CGETRS)(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zgetrs( char trans, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const integer* ipiv,
                    dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zgetrs,ZGETRS)(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_sgbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const float* ab, integer ldab,
                    const integer* ipiv, float* b, integer ldb )
{
    integer info;
    FFUNC(sgbtrs,SGBTRS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dgbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const double* ab, integer ldab,
                    const integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dgbtrs,DGBTRS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_cgbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const scomplex* ab, integer ldab,
                    const integer* ipiv, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cgbtrs,CGBTRS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zgbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const dcomplex* ab, integer ldab,
                    const integer* ipiv, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zgbtrs,ZGBTRS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_sgttrs( char trans, integer n, integer nrhs,
                    const float* dl, const float* d, const float* du,
                    const float* du2, const integer* ipiv, float* b,
                    integer ldb )
{
    integer info;
    FFUNC(sgttrs,SGTTRS)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dgttrs( char trans, integer n, integer nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* du2, const integer* ipiv, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dgttrs,DGTTRS)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_cgttrs( char trans, integer n, integer nrhs,
                    const scomplex* dl, const scomplex* d, const scomplex* du,
                    const scomplex* du2, const integer* ipiv, scomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(cgttrs,CGTTRS)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zgttrs( char trans, integer n, integer nrhs,
                    const dcomplex* dl, const dcomplex* d, const dcomplex* du,
                    const dcomplex* du2, const integer* ipiv, dcomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(zgttrs,ZGTTRS)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_spotrs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, float* b,
                    integer ldb )
{
    integer info;
    FFUNC(spotrs,SPOTRS)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_dpotrs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dpotrs,DPOTRS)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_cpotrs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, scomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(cpotrs,CPOTRS)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_zpotrs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, dcomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(zpotrs,ZPOTRS)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_spftrs( char transr, char uplo, integer n, integer nrhs,
                    const float* a, float* b, integer ldb )
{
    integer info;
    FFUNC(spftrs,SPFTRS)(&transr, &uplo, &n, &nrhs, a, b, &ldb, &info);
    return info;
}

static inline integer c_dpftrs( char transr, char uplo, integer n, integer nrhs,
                    const double* a, double* b, integer ldb )
{
    integer info;
    FFUNC(dpftrs,DPFTRS)(&transr, &uplo, &n, &nrhs, a, b, &ldb, &info);
    return info;
}

static inline integer c_cpftrs( char transr, char uplo, integer n, integer nrhs,
                    const scomplex* a, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cpftrs,CPFTRS)(&transr, &uplo, &n, &nrhs, a, b, &ldb, &info);
    return info;
}

static inline integer c_zpftrs( char transr, char uplo, integer n, integer nrhs,
                    const dcomplex* a, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zpftrs,ZPFTRS)(&transr, &uplo, &n, &nrhs, a, b, &ldb, &info);
    return info;
}

static inline integer c_spptrs( char uplo, integer n, integer nrhs,
                    const float* ap, float* b, integer ldb )
{
    integer info;
    FFUNC(spptrs,SPPTRS)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_dpptrs( char uplo, integer n, integer nrhs,
                    const double* ap, double* b, integer ldb )
{
    integer info;
    FFUNC(dpptrs,DPPTRS)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_cpptrs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cpptrs,CPPTRS)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_zpptrs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zpptrs,ZPPTRS)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_spbtrs( char uplo, integer n, integer kd, integer nrhs,
                    const float* ab, integer ldab, float* b,
                    integer ldb )
{
    integer info;
    FFUNC(spbtrs,SPBTRS)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_dpbtrs( char uplo, integer n, integer kd, integer nrhs,
                    const double* ab, integer ldab, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dpbtrs,DPBTRS)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_cpbtrs( char uplo, integer n, integer kd, integer nrhs,
                    const scomplex* ab, integer ldab, scomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(cpbtrs,CPBTRS)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_zpbtrs( char uplo, integer n, integer kd, integer nrhs,
                    const dcomplex* ab, integer ldab, dcomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(zpbtrs,ZPBTRS)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_spttrs( integer n, integer nrhs, const float* d,
                    const float* e, float* b, integer ldb )
{
    integer info;
    FFUNC(spttrs,SPTTRS)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

static inline integer c_dpttrs( integer n, integer nrhs, const double* d,
                    const double* e, double* b, integer ldb )
{
    integer info;
    FFUNC(dpttrs,DPTTRS)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

static inline integer c_cpttrs( integer n, integer nrhs, const float* d,
                    const scomplex* e, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cpttrs,CPTTRS)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

static inline integer c_zpttrs( integer n, integer nrhs, const double* d,
                    const dcomplex* e, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zpttrs,ZPTTRS)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

static inline integer c_ssytrs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, const integer* ipiv,
                    float* b, integer ldb )
{
    integer info;
    FFUNC(ssytrs,SSYTRS)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dsytrs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const integer* ipiv,
                    double* b, integer ldb )
{
    integer info;
    FFUNC(dsytrs,DSYTRS)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_csytrs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const integer* ipiv,
                    scomplex* b, integer ldb )
{
    integer info;
    FFUNC(csytrs,CSYTRS)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zsytrs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const integer* ipiv,
                    dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zsytrs,ZSYTRS)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_chetrs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const integer* ipiv,
                    scomplex* b, integer ldb )
{
    integer info;
    FFUNC(chetrs,CHETRS)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zhetrs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const integer* ipiv,
                    dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zhetrs,ZHETRS)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_ssptrs( char uplo, integer n, integer nrhs,
                    const float* ap, const integer* ipiv, float* b,
                    integer ldb )
{
    integer info;
    FFUNC(ssptrs,SSPTRS)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dsptrs( char uplo, integer n, integer nrhs,
                    const double* ap, const integer* ipiv, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dsptrs,DSPTRS)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_csptrs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const integer* ipiv, scomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(csptrs,CSPTRS)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zsptrs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const integer* ipiv, dcomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(zsptrs,ZSPTRS)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_chptrs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const integer* ipiv, scomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(chptrs,CHPTRS)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zhptrs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const integer* ipiv, dcomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(zhptrs,ZHPTRS)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_strtrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const float* a, integer lda,
                    float* b, integer ldb )
{
    integer info;
    FFUNC(strtrs,STRTRS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_dtrtrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* a, integer lda,
                    double* b, integer ldb )
{
    integer info;
    FFUNC(dtrtrs,DTRTRS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_ctrtrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const scomplex* a, integer lda,
                    scomplex* b, integer ldb )
{
    integer info;
    FFUNC(ctrtrs,CTRTRS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_ztrtrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const dcomplex* a, integer lda,
                    dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(ztrtrs,ZTRTRS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_stptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const float* ap, float* b,
                    integer ldb )
{
    integer info;
    FFUNC(stptrs,STPTRS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_dtptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* ap, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dtptrs,DTPTRS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_ctptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const scomplex* ap, scomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(ctptrs,CTPTRS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_ztptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const dcomplex* ap, dcomplex* b,
                    integer ldb )
{
    integer info;
    FFUNC(ztptrs,ZTPTRS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_stbtrs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const float* ab,
                    integer ldab, float* b, integer ldb )
{
    integer info;
    FFUNC(stbtrs,STBTRS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_dtbtrs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const double* ab,
                    integer ldab, double* b, integer ldb )
{
    integer info;
    FFUNC(dtbtrs,DTBTRS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_ctbtrs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const scomplex* ab,
                    integer ldab, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(ctbtrs,CTBTRS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_ztbtrs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const dcomplex* ab,
                    integer ldab, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(ztbtrs,ZTBTRS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_sgecon( char norm, integer n, const float* a, integer lda,
                    float anorm, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgecon,SGECON)(&norm, &n, a, &lda, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgecon( char norm, integer n, const double* a, integer lda,
                    double anorm, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgecon,DGECON)(&norm, &n, a, &lda, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgecon( char norm, integer n, const scomplex* a, integer lda,
                    float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgecon,CGECON)(&norm, &n, a, &lda, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgecon( char norm, integer n, const dcomplex* a, integer lda,
                    double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgecon,ZGECON)(&norm, &n, a, &lda, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgbcon( char norm, integer n, integer kl, integer ku,
                    const float* ab, integer ldab, const integer* ipiv,
                    float anorm, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgbcon,SGBCON)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgbcon( char norm, integer n, integer kl, integer ku,
                    const double* ab, integer ldab, const integer* ipiv,
                    double anorm, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgbcon,DGBCON)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgbcon( char norm, integer n, integer kl, integer ku,
                    const scomplex* ab, integer ldab, const integer* ipiv,
                    float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cgbcon,CGBCON)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgbcon( char norm, integer n, integer kl, integer ku,
                    const dcomplex* ab, integer ldab, const integer* ipiv,
                    double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zgbcon,ZGBCON)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgtcon( char norm, integer n, const float* dl,
                    const float* d, const float* du, const float* du2,
                    const integer* ipiv, float anorm, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgtcon,SGTCON)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgtcon( char norm, integer n, const double* dl,
                    const double* d, const double* du, const double* du2,
                    const integer* ipiv, double anorm, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgtcon,DGTCON)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgtcon( char norm, integer n, const scomplex* dl,
                    const scomplex* d, const scomplex* du, const scomplex* du2,
                    const integer* ipiv, float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    FFUNC(cgtcon,CGTCON)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zgtcon( char norm, integer n, const dcomplex* dl,
                    const dcomplex* d, const dcomplex* du, const dcomplex* du2,
                    const integer* ipiv, double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    FFUNC(zgtcon,ZGTCON)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_spocon( char uplo, integer n, const float* a, integer lda,
                    float anorm, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(spocon,SPOCON)(&uplo, &n, a, &lda, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dpocon( char uplo, integer n, const double* a, integer lda,
                    double anorm, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dpocon,DPOCON)(&uplo, &n, a, &lda, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cpocon( char uplo, integer n, const scomplex* a, integer lda,
                    float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cpocon,CPOCON)(&uplo, &n, a, &lda, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zpocon( char uplo, integer n, const dcomplex* a, integer lda,
                    double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zpocon,ZPOCON)(&uplo, &n, a, &lda, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sppcon( char uplo, integer n, const float* ap, float anorm, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sppcon,SPPCON)(&uplo, &n, ap, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dppcon( char uplo, integer n, const double* ap, double anorm, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dppcon,DPPCON)(&uplo, &n, ap, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cppcon( char uplo, integer n, const scomplex* ap, float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cppcon,CPPCON)(&uplo, &n, ap, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zppcon( char uplo, integer n, const dcomplex* ap, double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zppcon,ZPPCON)(&uplo, &n, ap, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_spbcon( char uplo, integer n, integer kd, const float* ab,
                    integer ldab, float anorm, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(spbcon,SPBCON)(&uplo, &n, &kd, ab, &ldab, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dpbcon( char uplo, integer n, integer kd, const double* ab,
                    integer ldab, double anorm, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dpbcon,DPBCON)(&uplo, &n, &kd, ab, &ldab, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cpbcon( char uplo, integer n, integer kd, const scomplex* ab,
                    integer ldab, float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cpbcon,CPBCON)(&uplo, &n, &kd, ab, &ldab, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zpbcon( char uplo, integer n, integer kd, const dcomplex* ab,
                    integer ldab, double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zpbcon,ZPBCON)(&uplo, &n, &kd, ab, &ldab, &anorm, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sptcon( integer n, const float* d, const float* e,
                    float anorm, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,n);
    FFUNC(sptcon,SPTCON)(&n, d, e, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dptcon( integer n, const double* d, const double* e,
                    double anorm, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,n);
    FFUNC(dptcon,DPTCON)(&n, d, e, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cptcon( integer n, const float* d, const scomplex* e,
                    float anorm, float* rcond)
{
    integer info;
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cptcon,CPTCON)(&n, d, e, &anorm, rcond, rwork, &info);
    FREE(rwork);
    return info;
}

static inline integer c_zptcon( integer n, const double* d, const dcomplex* e,
                    double anorm, double* rcond)
{
    integer info;
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zptcon,ZPTCON)(&n, d, e, &anorm, rcond, rwork, &info);
    FREE(rwork);
    return info;
}

static inline integer c_ssycon( char uplo, integer n, const float* a, integer lda,
                    const integer* ipiv, float anorm, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(ssycon,SSYCON)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsycon( char uplo, integer n, const double* a, integer lda,
                    const integer* ipiv, double anorm, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dsycon,DSYCON)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_csycon( char uplo, integer n, const scomplex* a, integer lda,
                    const integer* ipiv, float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    FFUNC(csycon,CSYCON)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zsycon( char uplo, integer n, const dcomplex* a, integer lda,
                    const integer* ipiv, double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    FFUNC(zsycon,ZSYCON)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_checon( char uplo, integer n, const scomplex* a, integer lda,
                    const integer* ipiv, float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    FFUNC(checon,CHECON)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zhecon( char uplo, integer n, const dcomplex* a, integer lda,
                    const integer* ipiv, double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    FFUNC(zhecon,ZHECON)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_sspcon( char uplo, integer n, const float* ap,
                    const integer* ipiv, float anorm, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sspcon,SSPCON)(&uplo, &n, ap, ipiv, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dspcon( char uplo, integer n, const double* ap,
                    const integer* ipiv, double anorm, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dspcon,DSPCON)(&uplo, &n, ap, ipiv, &anorm, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cspcon( char uplo, integer n, const scomplex* ap,
                    const integer* ipiv, float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    FFUNC(cspcon,CSPCON)(&uplo, &n, ap, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zspcon( char uplo, integer n, const dcomplex* ap,
                    const integer* ipiv, double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    FFUNC(zspcon,ZSPCON)(&uplo, &n, ap, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_chpcon( char uplo, integer n, const scomplex* ap,
                    const integer* ipiv, float anorm, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    FFUNC(chpcon,CHPCON)(&uplo, &n, ap, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zhpcon( char uplo, integer n, const dcomplex* ap,
                    const integer* ipiv, double anorm, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    FFUNC(zhpcon,ZHPCON)(&uplo, &n, ap, ipiv, &anorm, rcond, work, &info);
    FREE(work);
    return info;
}

static inline integer c_strcon( char norm, char uplo, char diag, integer n,
                    const float* a, integer lda, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(strcon,STRCON)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtrcon( char norm, char uplo, char diag, integer n,
                    const double* a, integer lda, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dtrcon,DTRCON)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctrcon( char norm, char uplo, char diag, integer n,
                    const scomplex* a, integer lda, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(ctrcon,CTRCON)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ztrcon( char norm, char uplo, char diag, integer n,
                    const dcomplex* a, integer lda, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(ztrcon,ZTRCON)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_stpcon( char norm, char uplo, char diag, integer n,
                    const float* ap, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(stpcon,STPCON)(&norm, &uplo, &diag, &n, ap, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtpcon( char norm, char uplo, char diag, integer n,
                    const double* ap, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dtpcon,DTPCON)(&norm, &uplo, &diag, &n, ap, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctpcon( char norm, char uplo, char diag, integer n,
                    const scomplex* ap, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(ctpcon,CTPCON)(&norm, &uplo, &diag, &n, ap, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ztpcon( char norm, char uplo, char diag, integer n,
                    const dcomplex* ap, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(ztpcon,ZTPCON)(&norm, &uplo, &diag, &n, ap, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_stbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const float* ab, integer ldab, float* rcond)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(stbcon,STBCON)(&norm, &uplo, &diag, &n, &kd, ab, &ldab, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const double* ab, integer ldab, double* rcond)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dtbcon,DTBCON)(&norm, &uplo, &diag, &n, &kd, ab, &ldab, rcond, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const scomplex* ab, integer ldab, float* rcond)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(ctbcon,CTBCON)(&norm, &uplo, &diag, &n, &kd, ab, &ldab, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ztbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const dcomplex* ab, integer ldab, double* rcond)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(ztbcon,ZTBCON)(&norm, &uplo, &diag, &n, &kd, ab, &ldab, rcond, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgerfs( char trans, integer n, integer nrhs,
                    const float* a, integer lda, const float* af,
                    integer ldaf, const integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgerfs,SGERFS)(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgerfs( char trans, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgerfs,DGERFS)(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgerfs( char trans, integer n, integer nrhs,
                    const scomplex* a, integer lda, const scomplex* af,
                    integer ldaf, const integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cgerfs,CGERFS)(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgerfs( char trans, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const dcomplex* af,
                    integer ldaf, const integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zgerfs,ZGERFS)(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgerfsx( char trans, char equed, integer n, integer nrhs,
                     const float* a, integer lda, const float* af,
                     integer ldaf, const integer* ipiv, const float* r,
                     const float* c, const float* b, integer ldb,
                     float* x, integer ldx, float* rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgerfsx,SGERFSX)(&trans, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgerfsx( char trans, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const integer* ipiv, const double* r,
                     const double* c, const double* b, integer ldb,
                     double* x, integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgerfsx,DGERFSX)(&trans, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgerfsx( char trans, char equed, integer n, integer nrhs,
                     const scomplex* a, integer lda, const scomplex* af,
                     integer ldaf, const integer* ipiv, const float* r,
                     const float* c, const scomplex* b, integer ldb,
                     scomplex* x, integer ldx, float* rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,4*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cgerfsx,CGERFSX)(&trans, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgerfsx( char trans, char equed, integer n, integer nrhs,
                     const dcomplex* a, integer lda, const dcomplex* af,
                     integer ldaf, const integer* ipiv, const double* r,
                     const double* c, const dcomplex* b, integer ldb,
                     dcomplex* x, integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,4*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zgerfsx,ZGERFSX)(&trans, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgbrfs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const float* ab, integer ldab,
                    const float* afb, integer ldafb,
                    const integer* ipiv, const float* b, integer ldb,
                    float* x, integer ldx, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgbrfs,SGBRFS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb, x, &ldx,
                            ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgbrfs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const double* ab, integer ldab,
                    const double* afb, integer ldafb,
                    const integer* ipiv, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgbrfs,DGBRFS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb, x, &ldx,
                            ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgbrfs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const scomplex* ab, integer ldab,
                    const scomplex* afb, integer ldafb,
                    const integer* ipiv, const scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cgbrfs,CGBRFS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb, x, &ldx,
                            ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgbrfs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const dcomplex* ab, integer ldab,
                    const dcomplex* afb, integer ldafb,
                    const integer* ipiv, const dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zgbrfs,ZGBRFS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb, x, &ldx,
                            ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgbrfsx( char trans, char equed, integer n, integer kl,
                     integer ku, integer nrhs, const float* ab,
                     integer ldab, const float* afb, integer ldafb,
                     const integer* ipiv, const float* r, const float* c,
                     const float* b, integer ldb, float* x,
                     integer ldx, float* rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgbrfsx,SGBRFSX)(&trans, &equed, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgbrfsx( char trans, char equed, integer n, integer kl,
                     integer ku, integer nrhs, const double* ab,
                     integer ldab, const double* afb, integer ldafb,
                     const integer* ipiv, const double* r, const double* c,
                     const double* b, integer ldb, double* x,
                     integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgbrfsx,DGBRFSX)(&trans, &equed, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgbrfsx( char trans, char equed, integer n, integer kl,
                     integer ku, integer nrhs, const scomplex* ab,
                     integer ldab, const scomplex* afb, integer ldafb,
                     const integer* ipiv, const float* r, const float* c,
                     const scomplex* b, integer ldb, scomplex* x,
                     integer ldx, float* rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgbrfsx,CGBRFSX)(&trans, &equed, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgbrfsx( char trans, char equed, integer n, integer kl,
                     integer ku, integer nrhs, const dcomplex* ab,
                     integer ldab, const dcomplex* afb, integer ldafb,
                     const integer* ipiv, const double* r, const double* c,
                     const dcomplex* b, integer ldb, dcomplex* x,
                     integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgbrfsx,ZGBRFSX)(&trans, &equed, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgtrfs( char trans, integer n, integer nrhs,
                    const float* dl, const float* d, const float* du,
                    const float* dlf, const float* df, const float* duf,
                    const float* du2, const integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgtrfs,SGTRFS)(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgtrfs( char trans, integer n, integer nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* dlf, const double* df, const double* duf,
                    const double* du2, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgtrfs,DGTRFS)(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgtrfs( char trans, integer n, integer nrhs,
                    const scomplex* dl, const scomplex* d, const scomplex* du,
                    const scomplex* dlf, const scomplex* df, const scomplex* duf,
                    const scomplex* du2, const integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cgtrfs,CGTRFS)(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr,
                            work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgtrfs( char trans, integer n, integer nrhs,
                    const dcomplex* dl, const dcomplex* d, const dcomplex* du,
                    const dcomplex* dlf, const dcomplex* df, const dcomplex* duf,
                    const dcomplex* du2, const integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zgtrfs,ZGTRFS)(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr,
                            work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sporfs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, const float* af,
                    integer ldaf, const float* b, integer ldb,
                    float* x, integer ldx, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sporfs,SPORFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dporfs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dporfs,DPORFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cporfs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const scomplex* af,
                    integer ldaf, const scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cporfs,CPORFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zporfs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const dcomplex* af,
                    integer ldaf, const dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zporfs,ZPORFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sporfsx( char uplo, char equed, integer n, integer nrhs,
                     const float* a, integer lda, const float* af,
                     integer ldaf, const float* s, const float* b,
                     integer ldb, float* x, integer ldx, float* rcond,
                     float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sporfsx,SPORFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, s, b, &ldb, x, &ldx, rcond, berr, &n_err_bnds,
                              err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dporfsx( char uplo, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const double* s, const double* b,
                     integer ldb, double* x, integer ldx, double* rcond,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dporfsx,DPORFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, s, b, &ldb, x, &ldx, rcond, berr, &n_err_bnds,
                              err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cporfsx( char uplo, char equed, integer n, integer nrhs,
                     const scomplex* a, integer lda, const scomplex* af,
                     integer ldaf, const float* s, const scomplex* b,
                     integer ldb, scomplex* x, integer ldx, float* rcond,
                     float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cporfsx,CPORFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, s, b, &ldb, x, &ldx, rcond, berr, &n_err_bnds,
                              err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zporfsx( char uplo, char equed, integer n, integer nrhs,
                     const dcomplex* a, integer lda, const dcomplex* af,
                     integer ldaf, const double* s, const dcomplex* b,
                     integer ldb, dcomplex* x, integer ldx, double* rcond,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zporfsx,ZPORFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, s, b, &ldb, x, &ldx, rcond, berr, &n_err_bnds,
                              err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_spprfs( char uplo, integer n, integer nrhs,
                    const float* ap, const float* afp, const float* b,
                    integer ldb, float* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(spprfs,SPPRFS)(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dpprfs( char uplo, integer n, integer nrhs,
                    const double* ap, const double* afp, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dpprfs,DPPRFS)(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cpprfs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const scomplex* afp, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cpprfs,CPPRFS)(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zpprfs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const dcomplex* afp, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zpprfs,ZPPRFS)(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_spbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const float* ab, integer ldab, const float* afb,
                    integer ldafb, const float* b, integer ldb,
                    float* x, integer ldx, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(spbrfs,SPBRFS)(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dpbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const double* ab, integer ldab, const double* afb,
                    integer ldafb, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dpbrfs,DPBRFS)(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cpbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const scomplex* ab, integer ldab, const scomplex* afb,
                    integer ldafb, const scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cpbrfs,CPBRFS)(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zpbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const dcomplex* ab, integer ldab, const dcomplex* afb,
                    integer ldafb, const dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zpbrfs,ZPBRFS)(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sptrfs( integer n, integer nrhs, const float* d,
                    const float* e, const float* df, const float* ef,
                    const float* b, integer ldb, float* x,
                    integer ldx, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    FFUNC(sptrfs,SPTRFS)(&n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dptrfs( integer n, integer nrhs, const double* d,
                    const double* e, const double* df, const double* ef,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    FFUNC(dptrfs,DPTRFS)(&n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cptrfs( integer n, integer nrhs, const float* d,
                    const scomplex* e, const float* df, const scomplex* ef,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cptrfs,CPTRFS)(&n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zptrfs( integer n, integer nrhs, const double* d,
                    const dcomplex* e, const double* df, const dcomplex* ef,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zptrfs,ZPTRFS)(&n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ssyrfs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, const float* af,
                    integer ldaf, const integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(ssyrfs,SSYRFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsyrfs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dsyrfs,DSYRFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_csyrfs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const scomplex* af,
                    integer ldaf, const integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(csyrfs,CSYRFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zsyrfs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const dcomplex* af,
                    integer ldaf, const integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zsyrfs,ZSYRFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_cherfs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const scomplex* af,
                    integer ldaf, const integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cherfs,CHERFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zherfs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const dcomplex* af,
                    integer ldaf, const integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zherfs,ZHERFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ssyrfsx( char uplo, char equed, integer n, integer nrhs,
                     const float* a, integer lda, const float* af,
                     integer ldaf, const integer* ipiv, const float* s,
                     const float* b, integer ldb, float* x,
                     integer ldx, float* rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(ssyrfsx,SSYRFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, s, b, &ldb, x, &ldx, rcond, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsyrfsx( char uplo, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const integer* ipiv, const double* s,
                     const double* b, integer ldb, double* x,
                     integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dsyrfsx,DSYRFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, s, b, &ldb, x, &ldx, rcond, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_csyrfsx( char uplo, char equed, integer n, integer nrhs,
                     const scomplex* a, integer lda, const scomplex* af,
                     integer ldaf, const integer* ipiv, const float* s,
                     const scomplex* b, integer ldb, scomplex* x,
                     integer ldx, float* rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(csyrfsx,CSYRFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, s, b, &ldb, x, &ldx, rcond, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zsyrfsx( char uplo, char equed, integer n, integer nrhs,
                     const dcomplex* a, integer lda, const dcomplex* af,
                     integer ldaf, const integer* ipiv, const double* s,
                     const dcomplex* b, integer ldb, dcomplex* x,
                     integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zsyrfsx,ZSYRFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, s, b, &ldb, x, &ldx, rcond, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_cherfsx( char uplo, char equed, integer n, integer nrhs,
                     const scomplex* a, integer lda, const scomplex* af,
                     integer ldaf, const integer* ipiv, const float* s,
                     const scomplex* b, integer ldb, scomplex* x,
                     integer ldx, float* rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cherfsx,CHERFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, s, b, &ldb, x, &ldx, rcond, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zherfsx( char uplo, char equed, integer n, integer nrhs,
                     const dcomplex* a, integer lda, const dcomplex* af,
                     integer ldaf, const integer* ipiv, const double* s,
                     const dcomplex* b, integer ldb, dcomplex* x,
                     integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zherfsx,ZHERFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, s, b, &ldb, x, &ldx, rcond, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ssprfs( char uplo, integer n, integer nrhs,
                    const float* ap, const float* afp, const integer* ipiv,
                    const float* b, integer ldb, float* x,
                    integer ldx, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(ssprfs,SSPRFS)(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsprfs( char uplo, integer n, integer nrhs,
                    const double* ap, const double* afp, const integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dsprfs,DSPRFS)(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_csprfs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const scomplex* afp, const integer* ipiv,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(csprfs,CSPRFS)(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zsprfs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const dcomplex* afp, const integer* ipiv,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zsprfs,ZSPRFS)(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_chprfs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const scomplex* afp, const integer* ipiv,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(chprfs,CHPRFS)(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhprfs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const dcomplex* afp, const integer* ipiv,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zhprfs,ZHPRFS)(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_strrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const float* a, integer lda,
                    const float* b, integer ldb, const float* x,
                    integer ldx, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(strrfs,STRRFS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtrrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* a, integer lda,
                    const double* b, integer ldb, const double* x,
                    integer ldx, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dtrrfs,DTRRFS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctrrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const scomplex* a, integer lda,
                    const scomplex* b, integer ldb, const scomplex* x,
                    integer ldx, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(ctrrfs,CTRRFS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ztrrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const dcomplex* a, integer lda,
                    const dcomplex* b, integer ldb, const dcomplex* x,
                    integer ldx, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(ztrrfs,ZTRRFS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_stprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const float* ap, const float* b,
                    integer ldb, const float* x, integer ldx,
                    float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(stprfs,STPRFS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* ap, const double* b,
                    integer ldb, const double* x, integer ldx,
                    double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dtprfs,DTPRFS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const scomplex* ap, const scomplex* b,
                    integer ldb, const scomplex* x, integer ldx,
                    float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(ctprfs,CTPRFS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ztprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const dcomplex* ap, const dcomplex* b,
                    integer ldb, const dcomplex* x, integer ldx,
                    double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(ztprfs,ZTPRFS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_stbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const float* ab,
                    integer ldab, const float* b, integer ldb,
                    const float* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(stbrfs,STBRFS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const double* ab,
                    integer ldab, const double* b, integer ldb,
                    const double* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dtbrfs,DTBRFS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const scomplex* ab,
                    integer ldab, const scomplex* b, integer ldb,
                    const scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(ctbrfs,CTBRFS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ztbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const dcomplex* ab,
                    integer ldab, const dcomplex* b, integer ldb,
                    const dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(ztbrfs,ZTBRFS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgetri( integer n, float* a, integer lda,
                    const integer* ipiv)
{
    integer info, lwork = -1;
    float query, *work;
    FFUNC(sgetri,SGETRI)(&n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgetri,SGETRI)(&n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgetri( integer n, double* a, integer lda,
                    const integer* ipiv)
{
    integer info, lwork = -1;
    double query, *work;
    FFUNC(dgetri,DGETRI)(&n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgetri,DGETRI)(&n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgetri( integer n, scomplex* a, integer lda,
                    const integer* ipiv)
{
    integer info, lwork = -1;
    scomplex query, *work;
    FFUNC(cgetri,CGETRI)(&n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgetri,CGETRI)(&n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgetri( integer n, dcomplex* a, integer lda,
                    const integer* ipiv)
{
    integer info, lwork = -1;
    dcomplex query, *work;
    FFUNC(zgetri,ZGETRI)(&n, a, &lda, ipiv, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgetri,ZGETRI)(&n, a, &lda, ipiv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_spotri( char uplo, integer n, float* a, integer lda )
{
    integer info;
    FFUNC(spotri,SPOTRI)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_dpotri( char uplo, integer n, double* a, integer lda )
{
    integer info;
    FFUNC(dpotri,DPOTRI)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_cpotri( char uplo, integer n, scomplex* a, integer lda )
{
    integer info;
    FFUNC(cpotri,CPOTRI)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_zpotri( char uplo, integer n, dcomplex* a, integer lda )
{
    integer info;
    FFUNC(zpotri,ZPOTRI)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_spftri( char transr, char uplo, integer n, float* a )
{
    integer info;
    FFUNC(spftri,SPFTRI)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_dpftri( char transr, char uplo, integer n, double* a )
{
    integer info;
    FFUNC(dpftri,DPFTRI)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_cpftri( char transr, char uplo, integer n, scomplex* a )
{
    integer info;
    FFUNC(cpftri,CPFTRI)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_zpftri( char transr, char uplo, integer n, dcomplex* a )
{
    integer info;
    FFUNC(zpftri,ZPFTRI)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_spptri( char uplo, integer n, float* ap )
{
    integer info;
    FFUNC(spptri,SPPTRI)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_dpptri( char uplo, integer n, double* ap )
{
    integer info;
    FFUNC(dpptri,DPPTRI)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_cpptri( char uplo, integer n, scomplex* ap )
{
    integer info;
    FFUNC(cpptri,CPPTRI)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_zpptri( char uplo, integer n, dcomplex* ap )
{
    integer info;
    FFUNC(zpptri,ZPPTRI)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_ssytri( char uplo, integer n, float* a, integer lda,
                    const integer* ipiv)
{
    integer info;
    float *work = SAFE_MALLOC(float,n);
    FFUNC(ssytri,SSYTRI)(&uplo, &n, a, &lda, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dsytri( char uplo, integer n, double* a, integer lda,
                    const integer* ipiv)
{
    integer info;
    double *work = SAFE_MALLOC(double,n);
    FFUNC(dsytri,DSYTRI)(&uplo, &n, a, &lda, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_csytri( char uplo, integer n, scomplex* a, integer lda,
                    const integer* ipiv)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    FFUNC(csytri,CSYTRI)(&uplo, &n, a, &lda, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zsytri( char uplo, integer n, dcomplex* a, integer lda,
                    const integer* ipiv)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    FFUNC(zsytri,ZSYTRI)(&uplo, &n, a, &lda, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_chetri( char uplo, integer n, scomplex* a, integer lda,
                    const integer* ipiv)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    FFUNC(chetri,CHETRI)(&uplo, &n, a, &lda, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zhetri( char uplo, integer n, dcomplex* a, integer lda,
                    const integer* ipiv)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    FFUNC(zhetri,ZHETRI)(&uplo, &n, a, &lda, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_ssptri( char uplo, integer n, float* ap,
                    const integer* ipiv)
{
    integer info;
    float *work = SAFE_MALLOC(float,n);
    FFUNC(ssptri,SSPTRI)(&uplo, &n, ap, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dsptri( char uplo, integer n, double* ap,
                    const integer* ipiv)
{
    integer info;
    double *work = SAFE_MALLOC(double,n);
    FFUNC(dsptri,DSPTRI)(&uplo, &n, ap, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_csptri( char uplo, integer n, scomplex* ap,
                    const integer* ipiv)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    FFUNC(csptri,CSPTRI)(&uplo, &n, ap, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zsptri( char uplo, integer n, dcomplex* ap,
                    const integer* ipiv)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    FFUNC(zsptri,ZSPTRI)(&uplo, &n, ap, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_chptri( char uplo, integer n, scomplex* ap,
                    const integer* ipiv)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    FFUNC(chptri,CHPTRI)(&uplo, &n, ap, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zhptri( char uplo, integer n, dcomplex* ap,
                    const integer* ipiv)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    FFUNC(zhptri,ZHPTRI)(&uplo, &n, ap, ipiv, work, &info);
    FREE(work);
    return info;
}

static inline integer c_strtri( char uplo, char diag, integer n, float* a,
                    integer lda )
{
    integer info;
    FFUNC(strtri,STRTRI)(&uplo, &diag, &n, a, &lda, &info);
    return info;
}

static inline integer c_dtrtri( char uplo, char diag, integer n, double* a,
                    integer lda )
{
    integer info;
    FFUNC(dtrtri,DTRTRI)(&uplo, &diag, &n, a, &lda, &info);
    return info;
}

static inline integer c_ctrtri( char uplo, char diag, integer n, scomplex* a,
                    integer lda )
{
    integer info;
    FFUNC(ctrtri,CTRTRI)(&uplo, &diag, &n, a, &lda, &info);
    return info;
}

static inline integer c_ztrtri( char uplo, char diag, integer n, dcomplex* a,
                    integer lda )
{
    integer info;
    FFUNC(ztrtri,ZTRTRI)(&uplo, &diag, &n, a, &lda, &info);
    return info;
}

static inline integer c_stftri( char transr, char uplo, char diag, integer n,
                    float* a )
{
    integer info;
    FFUNC(stftri,STFTRI)(&transr, &uplo, &diag, &n, a, &info);
    return info;
}

static inline integer c_dtftri( char transr, char uplo, char diag, integer n,
                    double* a )
{
    integer info;
    FFUNC(dtftri,DTFTRI)(&transr, &uplo, &diag, &n, a, &info);
    return info;
}

static inline integer c_ctftri( char transr, char uplo, char diag, integer n,
                    scomplex* a )
{
    integer info;
    FFUNC(ctftri,CTFTRI)(&transr, &uplo, &diag, &n, a, &info);
    return info;
}

static inline integer c_ztftri( char transr, char uplo, char diag, integer n,
                    dcomplex* a )
{
    integer info;
    FFUNC(ztftri,ZTFTRI)(&transr, &uplo, &diag, &n, a, &info);
    return info;
}

static inline integer c_stptri( char uplo, char diag, integer n, float* ap )
{
    integer info;
    FFUNC(stptri,STPTRI)(&uplo, &diag, &n, ap, &info);
    return info;
}

static inline integer c_dtptri( char uplo, char diag, integer n, double* ap )
{
    integer info;
    FFUNC(dtptri,DTPTRI)(&uplo, &diag, &n, ap, &info);
    return info;
}

static inline integer c_ctptri( char uplo, char diag, integer n, scomplex* ap )
{
    integer info;
    FFUNC(ctptri,CTPTRI)(&uplo, &diag, &n, ap, &info);
    return info;
}

static inline integer c_ztptri( char uplo, char diag, integer n, dcomplex* ap )
{
    integer info;
    FFUNC(ztptri,ZTPTRI)(&uplo, &diag, &n, ap, &info);
    return info;
}

static inline integer c_sgeequ( integer m, integer n, const float* a,
                    integer lda, float* r, float* c, float* rowcnd,
                    float* colcnd, float* amax )
{
    integer info;
    FFUNC(sgeequ,SGEEQU)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_dgeequ( integer m, integer n, const double* a,
                    integer lda, double* r, double* c, double* rowcnd,
                    double* colcnd, double* amax )
{
    integer info;
    FFUNC(dgeequ,DGEEQU)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_cgeequ( integer m, integer n, const scomplex* a,
                    integer lda, float* r, float* c, float* rowcnd,
                    float* colcnd, float* amax )
{
    integer info;
    FFUNC(cgeequ,CGEEQU)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_zgeequ( integer m, integer n, const dcomplex* a,
                    integer lda, double* r, double* c, double* rowcnd,
                    double* colcnd, double* amax )
{
    integer info;
    FFUNC(zgeequ,ZGEEQU)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_sgeequb( integer m, integer n, const float* a,
                     integer lda, float* r, float* c, float* rowcnd,
                     float* colcnd, float* amax )
{
    integer info;
    FFUNC(sgeequb,SGEEQUB)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_dgeequb( integer m, integer n, const double* a,
                     integer lda, double* r, double* c, double* rowcnd,
                     double* colcnd, double* amax )
{
    integer info;
    FFUNC(dgeequb,DGEEQUB)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_cgeequb( integer m, integer n, const scomplex* a,
                     integer lda, float* r, float* c, float* rowcnd,
                     float* colcnd, float* amax )
{
    integer info;
    FFUNC(cgeequb,CGEEQUB)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_zgeequb( integer m, integer n, const dcomplex* a,
                     integer lda, double* r, double* c, double* rowcnd,
                     double* colcnd, double* amax )
{
    integer info;
    FFUNC(zgeequb,ZGEEQUB)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_sgbequ( integer m, integer n, integer kl,
                    integer ku, const float* ab, integer ldab,
                    float* r, float* c, float* rowcnd, float* colcnd,
                    float* amax )
{
    integer info;
    FFUNC(sgbequ,SGBEQU)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_dgbequ( integer m, integer n, integer kl,
                    integer ku, const double* ab, integer ldab,
                    double* r, double* c, double* rowcnd, double* colcnd,
                    double* amax )
{
    integer info;
    FFUNC(dgbequ,DGBEQU)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_cgbequ( integer m, integer n, integer kl,
                    integer ku, const scomplex* ab, integer ldab,
                    float* r, float* c, float* rowcnd, float* colcnd,
                    float* amax )
{
    integer info;
    FFUNC(cgbequ,CGBEQU)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_zgbequ( integer m, integer n, integer kl,
                    integer ku, const dcomplex* ab, integer ldab,
                    double* r, double* c, double* rowcnd, double* colcnd,
                    double* amax )
{
    integer info;
    FFUNC(zgbequ,ZGBEQU)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_sgbequb( integer m, integer n, integer kl,
                     integer ku, const float* ab, integer ldab,
                     float* r, float* c, float* rowcnd, float* colcnd,
                     float* amax )
{
    integer info;
    FFUNC(sgbequb,SGBEQUB)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_dgbequb( integer m, integer n, integer kl,
                     integer ku, const double* ab, integer ldab,
                     double* r, double* c, double* rowcnd, double* colcnd,
                     double* amax )
{
    integer info;
    FFUNC(dgbequb,DGBEQUB)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_cgbequb( integer m, integer n, integer kl,
                     integer ku, const scomplex* ab, integer ldab,
                     float* r, float* c, float* rowcnd, float* colcnd,
                     float* amax )
{
    integer info;
    FFUNC(cgbequb,CGBEQUB)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_zgbequb( integer m, integer n, integer kl,
                     integer ku, const dcomplex* ab, integer ldab,
                     double* r, double* c, double* rowcnd, double* colcnd,
                     double* amax )
{
    integer info;
    FFUNC(zgbequb,ZGBEQUB)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_spoequ( integer n, const float* a, integer lda, float* s,
                    float* scond, float* amax )
{
    integer info;
    FFUNC(spoequ,SPOEQU)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_dpoequ( integer n, const double* a, integer lda, double* s,
                    double* scond, double* amax )
{
    integer info;
    FFUNC(dpoequ,DPOEQU)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_cpoequ( integer n, const scomplex* a, integer lda, float* s,
                    float* scond, float* amax )
{
    integer info;
    FFUNC(cpoequ,CPOEQU)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_zpoequ( integer n, const dcomplex* a, integer lda, double* s,
                    double* scond, double* amax )
{
    integer info;
    FFUNC(zpoequ,ZPOEQU)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_spoequb( integer n, const float* a, integer lda, float* s,
                     float* scond, float* amax )
{
    integer info;
    FFUNC(spoequb,SPOEQUB)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_dpoequb( integer n, const double* a, integer lda, double* s,
                     double* scond, double* amax )
{
    integer info;
    FFUNC(dpoequb,DPOEQUB)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_cpoequb( integer n, const scomplex* a, integer lda, float* s,
                    float* scond, float* amax )
{
    integer info;
    FFUNC(cpoequb,CPOEQUB)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_zpoequb( integer n, const dcomplex* a, integer lda, double* s,
                     double* scond, double* amax )
{
    integer info;
    FFUNC(zpoequb,ZPOEQUB)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_sppequ( char uplo, integer n, const float* ap, float* s,
                    float* scond, float* amax )
{
    integer info;
    FFUNC(sppequ,SPPEQU)(&uplo, &n, ap, s, scond, amax, &info);
    return info;
}

static inline integer c_dppequ( char uplo, integer n, const double* ap, double* s,
                    double* scond, double* amax )
{
    integer info;
    FFUNC(dppequ,DPPEQU)(&uplo, &n, ap, s, scond, amax, &info);
    return info;
}

static inline integer c_cppequ( char uplo, integer n, const scomplex* ap, float* s,
                    float* scond, float* amax )
{
    integer info;
    FFUNC(cppequ,CPPEQU)(&uplo, &n, ap, s, scond, amax, &info);
    return info;
}

static inline integer c_zppequ( char uplo, integer n, const dcomplex* ap, double* s,
                    double* scond, double* amax )
{
    integer info;
    FFUNC(zppequ,ZPPEQU)(&uplo, &n, ap, s, scond, amax, &info);
    return info;
}

static inline integer c_spbequ( char uplo, integer n, integer kd, const float* ab,
                    integer ldab, float* s, float* scond, float* amax )
{
    integer info;
    FFUNC(spbequ,SPBEQU)(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info);
    return info;
}

static inline integer c_dpbequ( char uplo, integer n, integer kd, const double* ab,
                    integer ldab, double* s, double* scond, double* amax )
{
    integer info;
    FFUNC(dpbequ,DPBEQU)(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info);
    return info;
}

static inline integer c_cpbequ( char uplo, integer n, integer kd, const scomplex* ab,
                    integer ldab, float* s, float* scond, float* amax )
{
    integer info;
    FFUNC(cpbequ,CPBEQU)(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info);
    return info;
}

static inline integer c_zpbequ( char uplo, integer n, integer kd, const dcomplex* ab,
                    integer ldab, double* s, double* scond, double* amax )
{
    integer info;
    FFUNC(zpbequ,ZPBEQU)(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info);
    return info;
}

static inline integer c_ssyequb( char uplo, integer n, const float* a,
                     integer lda, float* s, float* scond, float* amax)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    FFUNC(ssyequb,SSYEQUB)(&uplo, &n, a, &lda, s, scond, amax, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dsyequb( char uplo, integer n, const double* a,
                     integer lda, double* s, double* scond, double* amax)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    FFUNC(dsyequb,DSYEQUB)(&uplo, &n, a, &lda, s, scond, amax, work, &info);
    FREE(work);
    return info;
}

static inline integer c_csyequb( char uplo, integer n, const scomplex* a,
                     integer lda, float* s, float* scond, float* amax)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,3*n);
    FFUNC(csyequb,CSYEQUB)(&uplo, &n, a, &lda, s, scond, amax, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zsyequb( char uplo, integer n, const dcomplex* a,
                     integer lda, double* s, double* scond, double* amax)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,3*n);
    FFUNC(zsyequb,ZSYEQUB)(&uplo, &n, a, &lda, s, scond, amax, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cheequb( char uplo, integer n, const scomplex* a,
                     integer lda, float* s, float* scond, float* amax)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,3*n);
    FFUNC(cheequb,CHEEQUB)(&uplo, &n, a, &lda, s, scond, amax, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zheequb( char uplo, integer n, const dcomplex* a,
                     integer lda, double* s, double* scond, double* amax)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,3*n);
    FFUNC(zheequb,ZHEEQUB)(&uplo, &n, a, &lda, s, scond, amax, work, &info);
    FREE(work);
    return info;
}

static inline integer c_sgesv( integer n, integer nrhs, float* a, integer lda,
                   integer* ipiv, float* b, integer ldb )
{
    integer info;
    FFUNC(sgesv,SGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dgesv( integer n, integer nrhs, double* a, integer lda,
                   integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dgesv,DGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_cgesv( integer n, integer nrhs, scomplex* a, integer lda,
                   integer* ipiv, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cgesv,CGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zgesv( integer n, integer nrhs, dcomplex* a, integer lda,
                   integer* ipiv, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zgesv,ZGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dsgesv( integer n, integer nrhs, double* a, integer lda,
                    integer* ipiv, double* b, integer ldb, double* x,
                    integer ldx, integer* iter )
{
    integer info;
    double *work = SAFE_MALLOC(double,n*nrhs);
    float *swork = SAFE_MALLOC(float,n*(n+nrhs));
    FFUNC(dsgesv,DSGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, x, &ldx, work, swork, iter, &info);
    FREE(work);
    FREE(swork);
    return info;
}

static inline integer c_zcgesv( integer n, integer nrhs, dcomplex* a, integer lda,
                    integer* ipiv, dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, integer* iter )
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n*nrhs);
    scomplex *swork = SAFE_MALLOC(scomplex,n*(n+nrhs));
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zcgesv,ZCGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, x, &ldx, work, swork, rwork, iter, &info);
    FREE(work);
    FREE(swork);
    FREE(rwork);
    return info;
}

static inline integer c_sgesvx( char fact, char trans, integer n, integer nrhs,
                    float* a, integer lda, float* af, integer ldaf,
                    integer* ipiv, char* equed, float* r, float* c,
                    float* b, integer ldb, float* x, integer ldx,
                    float* rcond, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgesvx,SGESVX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldx, x, &ldx, rcond,
                           ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgesvx( char fact, char trans, integer n, integer nrhs,
                    double* a, integer lda, double* af, integer ldaf,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer ldb, double* x, integer ldx,
                    double* rcond, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgesvx,DGESVX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldx, x, &ldx, rcond,
                           ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgesvx( char fact, char trans, integer n, integer nrhs,
                    scomplex* a, integer lda, scomplex* af, integer ldaf,
                    integer* ipiv, char* equed, float* r, float* c,
                    scomplex* b, integer ldb, scomplex* x, integer ldx,
                    float* rcond, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgesvx,CGESVX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldx, x, &ldx, rcond,
                           ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgesvx( char fact, char trans, integer n, integer nrhs,
                    dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                    integer* ipiv, char* equed, double* r, double* c,
                    dcomplex* b, integer ldb, dcomplex* x, integer ldx,
                    double* rcond, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgesvx,ZGESVX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldx, x, &ldx, rcond,
                           ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgesvxx( char fact, char trans, integer n, integer nrhs,
                    float* a, integer lda, float* af, integer ldaf,
                     integer* ipiv, char* equed, float* r, float* c,
                     float* b, integer ldb, float* x, integer ldx,
                     float* rcond, float* rpvgrw, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgesvxx,SGESVXX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldb, x, &ldx, rcond, rpvgrw,
                              berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgesvxx( char fact, char trans, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer ldb, double* x, integer ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgesvxx,DGESVXX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldb, x, &ldx, rcond, rpvgrw,
                              berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgesvxx( char fact, char trans, integer n, integer nrhs,
                     scomplex* a, integer lda, scomplex* af, integer ldaf,
                     integer* ipiv, char* equed, float* r, float* c,
                     scomplex* b, integer ldb, scomplex* x, integer ldx,
                     float* rcond, float* rpvgrw, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,3*n);
    FFUNC(cgesvxx,CGESVXX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldb, x, &ldx, rcond, rpvgrw,
                              berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgesvxx( char fact, char trans, integer n, integer nrhs,
                    dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                     integer* ipiv, char* equed, double* r, double* c,
                     dcomplex* b, integer ldb, dcomplex* x, integer ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,3*n);
    FFUNC(zgesvxx,ZGESVXX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldb, x, &ldx, rcond, rpvgrw,
                              berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgbsv( integer n, integer kl, integer ku,
                   integer nrhs, float* ab, integer ldab,
                   integer* ipiv, float* b, integer ldb )
{
    integer info;
    FFUNC(sgbsv,SGBSV)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dgbsv( integer n, integer kl, integer ku,
                   integer nrhs, double* ab, integer ldab,
                   integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dgbsv,DGBSV)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_cgbsv( integer n, integer kl, integer ku,
                   integer nrhs, scomplex* ab, integer ldab,
                   integer* ipiv, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cgbsv,CGBSV)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zgbsv( integer n, integer kl, integer ku,
                   integer nrhs, dcomplex* ab, integer ldab,
                   integer* ipiv, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zgbsv,ZGBSV)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_sgbsvx( char fact, char trans, integer n, integer kl,
                    integer ku, integer nrhs, float* ab,
                    integer ldab, float* afb, integer ldafb,
                    integer* ipiv, char* equed, float* r, float* c,
                    float* b, integer ldb, float* x, integer ldx,
                    float* rcond, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgbsvx,SGBSVX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                           rcond, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgbsvx( char fact, char trans, integer n, integer kl,
                    integer ku, integer nrhs, double* ab,
                    integer ldab, double* afb, integer ldafb,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer ldb, double* x, integer ldx,
                    double* rcond, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgbsvx,DGBSVX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                           rcond, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgbsvx( char fact, char trans, integer n, integer kl,
                    integer ku, integer nrhs, scomplex* ab,
                    integer ldab, scomplex* afb, integer ldafb,
                    integer* ipiv, char* equed, float* r, float* c,
                    scomplex* b, integer ldb, scomplex* x, integer ldx,
                    float* rcond, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cgbsvx,CGBSVX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                           rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgbsvx( char fact, char trans, integer n, integer kl,
                    integer ku, integer nrhs, dcomplex* ab,
                    integer ldab, dcomplex* afb, integer ldafb,
                    integer* ipiv, char* equed, double* r, double* c,
                    dcomplex* b, integer ldb, dcomplex* x, integer ldx,
                    double* rcond, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zgbsvx,ZGBSVX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                           rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgbsvxx( char fact, char trans, integer n, integer kl,
                     integer ku, integer nrhs, float* ab,
                     integer ldab, float* afb, integer ldafb,
                     integer* ipiv, char* equed, float* r, float* c,
                     float* b, integer ldb, float* x, integer ldx,
                     float* rcond, float* rpvgrw, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgbsvxx,SGBSVXX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                             rcond, rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgbsvxx( char fact, char trans, integer n, integer kl,
                     integer ku, integer nrhs, double* ab,
                     integer ldab, double* afb, integer ldafb,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer ldb, double* x, integer ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgbsvxx,DGBSVXX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                             rcond, rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgbsvxx( char fact, char trans, integer n, integer kl,
                     integer ku, integer nrhs, scomplex* ab,
                     integer ldab, scomplex* afb, integer ldafb,
                     integer* ipiv, char* equed, float* r, float* c,
                     scomplex* b, integer ldb, scomplex* x, integer ldx,
                     float* rcond, float* rpvgrw, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgbsvxx,CGBSVXX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                             rcond, rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgbsvxx( char fact, char trans, integer n, integer kl,
                     integer ku, integer nrhs, dcomplex* ab,
                     integer ldab, dcomplex* afb, integer ldafb,
                     integer* ipiv, char* equed, double* r, double* c,
                     dcomplex* b, integer ldb, dcomplex* x, integer ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgbsvxx,ZGBSVXX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                             rcond, rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgtsv( integer n, integer nrhs, float* dl, float* d,
                   float* du, float* b, integer ldb )
{
    integer info;
    FFUNC(sgtsv,SGTSV)(&n, &nrhs, dl, d, du, b, &ldb, &info);
    return info;
}

static inline integer c_dgtsv( integer n, integer nrhs, double* dl, double* d,
                   double* du, double* b, integer ldb )
{
    integer info;
    FFUNC(dgtsv,DGTSV)(&n, &nrhs, dl, d, du, b, &ldb, &info);
    return info;
}

static inline integer c_cgtsv( integer n, integer nrhs, scomplex* dl, scomplex* d,
                   scomplex* du, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cgtsv,CGTSV)(&n, &nrhs, dl, d, du, b, &ldb, &info);
    return info;
}

static inline integer c_zgtsv( integer n, integer nrhs, dcomplex* dl, dcomplex* d,
                   dcomplex* du, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zgtsv,ZGTSV)(&n, &nrhs, dl, d, du, b, &ldb, &info);
    return info;
}

static inline integer c_sgtsvx( char fact, char trans, integer n, integer nrhs,
                    const float* dl, const float* d, const float* du,
                    float* dlf, float* df, float* duf, float* du2,
                    integer* ipiv, const float* b, integer ldb,
                    float* x, integer ldx, float* rcond, float* ferr,
                    float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sgtsvx,SGTSVX)(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgtsvx( char fact, char trans, integer n, integer nrhs,
                    const double* dl, const double* d, const double* du,
                    double* dlf, double* df, double* duf, double* du2,
                    integer* ipiv, const double* b, integer ldb,
                    double* x, integer ldx, double* rcond, double* ferr,
                    double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dgtsvx,DGTSVX)(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgtsvx( char fact, char trans, integer n, integer nrhs,
                    const scomplex* dl, const scomplex* d, const scomplex* du,
                    scomplex* dlf, scomplex* df, scomplex* duf, scomplex* du2,
                    integer* ipiv, const scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float* rcond, float* ferr,
                    float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cgtsvx,CGTSVX)(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgtsvx( char fact, char trans, integer n, integer nrhs,
                    const dcomplex* dl, const dcomplex* d, const dcomplex* du,
                    dcomplex* dlf, dcomplex* df, dcomplex* duf, dcomplex* du2,
                    integer* ipiv, const dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double* rcond, double* ferr,
                    double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zgtsvx,ZGTSVX)(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sposv( char uplo, integer n, integer nrhs, float* a,
                   integer lda, float* b, integer ldb )
{
    integer info;
    FFUNC(sposv,SPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_dposv( char uplo, integer n, integer nrhs, double* a,
                   integer lda, double* b, integer ldb )
{
    integer info;
    FFUNC(dposv,DPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_cposv( char uplo, integer n, integer nrhs, scomplex* a,
                   integer lda, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cposv,CPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_zposv( char uplo, integer n, integer nrhs, dcomplex* a,
                   integer lda, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zposv,ZPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_dsposv( char uplo, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* x,
                    integer ldx, integer* iter )
{
    integer info;
    double *work = SAFE_MALLOC(double,n*nrhs);
    float *swork = SAFE_MALLOC(float,n*(n+nrhs));
    FFUNC(dsposv,DSPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, work, swork, iter, &info);
    FREE(work);
    FREE(swork);
    return info;
}

static inline integer c_zcposv( char uplo, integer n, integer nrhs, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, integer* iter )
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n*nrhs);
    scomplex *swork = SAFE_MALLOC(scomplex,n*(n+nrhs));
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zcposv,ZCPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, work, swork, rwork, iter, &info);
    FREE(work);
    FREE(swork);
    FREE(rwork);
    return info;
}

static inline integer c_sposvx( char fact, char uplo, integer n, integer nrhs,
                    float* a, integer lda, float* af, integer ldaf,
                    char* equed, float* s, float* b, integer ldb,
                    float* x, integer ldx, float* rcond, float* ferr,
                    float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sposvx,SPOSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dposvx( char fact, char uplo, integer n, integer nrhs,
                    double* a, integer lda, double* af, integer ldaf,
                    char* equed, double* s, double* b, integer ldb,
                    double* x, integer ldx, double* rcond, double* ferr,
                    double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dposvx,DPOSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cposvx( char fact, char uplo, integer n, integer nrhs,
                    scomplex* a, integer lda, scomplex* af, integer ldaf,
                    char* equed, float* s, scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float* rcond, float* ferr,
                    float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cposvx,CPOSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zposvx( char fact, char uplo, integer n, integer nrhs,
                    dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                    char* equed, double* s, dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double* rcond, double* ferr,
                    double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zposvx,ZPOSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sposvxx( char fact, char uplo, integer n, integer nrhs,
                    float* a, integer lda, float* af, integer ldaf,
                     char* equed, float* s, float* b, integer ldb,
                     float* x, integer ldx, float* rcond, float* rpvgrw,
                     float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sposvxx,SPOSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, rpvgrw, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dposvxx( char fact, char uplo, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     char* equed, double* s, double* b, integer ldb,
                     double* x, integer ldx, double* rcond, double* rpvgrw,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dposvxx,DPOSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, rpvgrw, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cposvxx( char fact, char uplo, integer n, integer nrhs,
                     scomplex* a, integer lda, scomplex* af, integer ldaf,
                     char* equed, float* s, scomplex* b, integer ldb,
                     scomplex* x, integer ldx, float* rcond, float* rpvgrw,
                     float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cposvxx,CPOSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, rpvgrw, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zposvxx( char fact, char uplo, integer n, integer nrhs,
                     dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                     char* equed, double* s, dcomplex* b, integer ldb,
                     dcomplex* x, integer ldx, double* rcond, double* rpvgrw,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zposvxx,ZPOSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, rpvgrw, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sppsv( char uplo, integer n, integer nrhs, float* ap,
                    float* b, integer ldb )
{
    integer info;
    FFUNC(sppsv,SPPSV)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_dppsv( char uplo, integer n, integer nrhs, double* ap,
                   double* b, integer ldb )
{
    integer info;
    FFUNC(dppsv,DPPSV)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_cppsv( char uplo, integer n, integer nrhs, scomplex* ap,
                   scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cppsv,CPPSV)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_zppsv( char uplo, integer n, integer nrhs, dcomplex* ap,
                   dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zppsv,ZPPSV)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_sppsvx( char fact, char uplo, integer n, integer nrhs,
                    float* ap, float* afp, char* equed, float* s, float* b,
                    integer ldb, float* x, integer ldx, float* rcond,
                    float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sppsvx,SPPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dppsvx( char fact, char uplo, integer n, integer nrhs,
                    double* ap, double* afp, char* equed, double* s, double* b,
                    integer ldb, double* x, integer ldx, double* rcond,
                    double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dppsvx,DPPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cppsvx( char fact, char uplo, integer n, integer nrhs,
                    scomplex* ap, scomplex* afp, char* equed, float* s, scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* rcond,
                    float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cppsvx,CPPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zppsvx( char fact, char uplo, integer n, integer nrhs,
                    dcomplex* ap, dcomplex* afp, char* equed, double* s, dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* rcond,
                    double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zppsvx,ZPPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_spbsv( char uplo, integer n, integer kd, integer nrhs,
                   float* ab, integer ldab, float* b, integer ldb )
{
    integer info;
    FFUNC(spbsv,SPBSV)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_dpbsv( char uplo, integer n, integer kd, integer nrhs,
                   double* ab, integer ldab, double* b, integer ldb )
{
    integer info;
    FFUNC(dpbsv,DPBSV)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_cpbsv( char uplo, integer n, integer kd, integer nrhs,
                   scomplex* ab, integer ldab, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cpbsv,CPBSV)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_zpbsv( char uplo, integer n, integer kd, integer nrhs,
                   dcomplex* ab, integer ldab, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zpbsv,ZPBSV)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_spbsvx( char fact, char uplo, integer n, integer kd,
                    integer nrhs, float* ab, integer ldab, float* afb,
                    integer ldafb, char* equed, float* s, float* b,
                    integer ldb, float* x, integer ldx, float* rcond,
                    float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(spbsvx,SPBSVX)(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dpbsvx( char fact, char uplo, integer n, integer kd,
                    integer nrhs, double* ab, integer ldab, double* afb,
                    integer ldafb, char* equed, double* s, double* b,
                    integer ldb, double* x, integer ldx, double* rcond,
                    double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dpbsvx,DPBSVX)(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cpbsvx( char fact, char uplo, integer n, integer kd,
                    integer nrhs, scomplex* ab, integer ldab, scomplex* afb,
                    integer ldafb, char* equed, float* s, scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* rcond,
                    float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cpbsvx,CPBSVX)(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zpbsvx( char fact, char uplo, integer n, integer kd,
                    integer nrhs, dcomplex* ab, integer ldab, dcomplex* afb,
                    integer ldafb, char* equed, double* s, dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* rcond,
                    double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zpbsvx,ZPBSVX)(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sptsv( integer n, integer nrhs, float* d, float* e,
                   float* b, integer ldb )
{
    integer info;
    FFUNC(sptsv,SPTSV)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

static inline integer c_dptsv( integer n, integer nrhs, double* d, double* e,
                   double* b, integer ldb )
{
    integer info;
    FFUNC(dptsv,DPTSV)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

static inline integer c_cptsv( integer n, integer nrhs, float* d, scomplex* e,
                   scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cptsv,CPTSV)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

static inline integer c_zptsv( integer n, integer nrhs, double* d, dcomplex* e,
                   dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zptsv,ZPTSV)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

static inline integer c_sptsvx( char fact, integer n, integer nrhs,
                    const float* d, const float* e, float* df, float* ef,
                    const float* b, integer ldb, float* x,
                    integer ldx, float* rcond, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    FFUNC(sptsvx,SPTSVX)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dptsvx( char fact, integer n, integer nrhs,
                    const double* d, const double* e, double* df, double* ef,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* rcond, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    FFUNC(dptsvx,DPTSVX)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cptsvx( char fact, integer n, integer nrhs,
                    const float* d, const scomplex* e, float* df, scomplex* ef,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float* rcond, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cptsvx,CPTSVX)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zptsvx( char fact, integer n, integer nrhs,
                    const double* d, const dcomplex* e, double* df, dcomplex* ef,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double* rcond, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zptsvx,ZPTSVX)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ssysv( char uplo, integer n, integer nrhs, float* a,
                   integer lda, integer* ipiv, float* b,
                   integer ldb)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(ssysv,SSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(ssysv,SSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dsysv( char uplo, integer n, integer nrhs, double* a,
                   integer lda, integer* ipiv, double* b,
                   integer ldb)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dsysv,DSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dsysv,DSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_csysv( char uplo, integer n, integer nrhs, scomplex* a,
                   integer lda, integer* ipiv, scomplex* b,
                   integer ldb)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(csysv,CSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(csysv,CSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zsysv( char uplo, integer n, integer nrhs, dcomplex* a,
                   integer lda, integer* ipiv, dcomplex* b,
                   integer ldb)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zsysv,ZSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zsysv,ZSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_chesv( char uplo, integer n, integer nrhs, scomplex* a,
                   integer lda, integer* ipiv, scomplex* b,
                   integer ldb)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(chesv,CHESV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(chesv,CHESV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zhesv( char uplo, integer n, integer nrhs, dcomplex* a,
                   integer lda, integer* ipiv, dcomplex* b,
                   integer ldb)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zhesv,ZHESV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zhesv,ZHESV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_ssysvx( char fact, char uplo, integer n, integer nrhs,
                    const float* a, integer lda, float* af,
                    integer ldaf, integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float* rcond,
                    float* ferr, float* berr)
{
    integer info, lwork = -1;
    float *work, query;
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(ssysvx,SSYSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            &query, &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(ssysvx,SSYSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsysvx( char fact, char uplo, integer n, integer nrhs,
                    const double* a, integer lda, double* af,
                    integer ldaf, integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* rcond,
                    double* ferr, double* berr)
{
    integer info, lwork = -1;
    double *work, query;
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dsysvx,DSYSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            &query, &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dsysvx,DSYSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_csysvx( char fact, char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, scomplex* af,
                    integer ldaf, integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* rcond,
                    float* ferr, float* berr)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(csysvx,CSYSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(csysvx,CSYSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zsysvx( char fact, char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, dcomplex* af,
                    integer ldaf, integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* rcond,
                    double* ferr, double* berr)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zsysvx,ZSYSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zsysvx,ZSYSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_chesvx( char fact, char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, scomplex* af,
                    integer ldaf, integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* rcond,
                    float* ferr, float* berr)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(chesvx,CHESVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(chesvx,CHESVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhesvx( char fact, char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, dcomplex* af,
                    integer ldaf, integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* rcond,
                    double* ferr, double* berr)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zhesvx,ZHESVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zhesvx,ZHESVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ssysvxx( char fact, char uplo, integer n, integer nrhs,
                     float* a, integer lda, float* af, integer ldaf,
                     integer* ipiv, char* equed, float* s, float* b,
                     integer ldb, float* x, integer ldx, float* rcond,
                     float* rpvgrw, float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(ssysvxx,SSYSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, s, b, &ldb, x, &ldx, rcond,
                              rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsysvxx( char fact, char uplo, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     integer* ipiv, char* equed, double* s, double* b,
                     integer ldb, double* x, integer ldx, double* rcond,
                     double* rpvgrw, double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dsysvxx,DSYSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, s, b, &ldb, x, &ldx, rcond,
                              rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_csysvxx( char fact, char uplo, integer n, integer nrhs,
                     scomplex* a, integer lda, scomplex* af, integer ldaf,
                     integer* ipiv, char* equed, float* s, scomplex* b,
                     integer ldb, scomplex* x, integer ldx, float* rcond,
                     float* rpvgrw, float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,3*n);
    FFUNC(csysvxx,CSYSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, s, b, &ldb, x, &ldx, rcond,
                              rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zsysvxx( char fact, char uplo, integer n, integer nrhs,
                     dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                     integer* ipiv, char* equed, double* s, dcomplex* b,
                     integer ldb, dcomplex* x, integer ldx, double* rcond,
                     double* rpvgrw, double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,3*n);
    FFUNC(zsysvxx,ZSYSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, s, b, &ldb, x, &ldx, rcond,
                              rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_chesvxx( char fact, char uplo, integer n, integer nrhs,
                     scomplex* a, integer lda, scomplex* af, integer ldaf,
                     integer* ipiv, char* equed, float* s, scomplex* b,
                     integer ldb, scomplex* x, integer ldx, float* rcond,
                     float* rpvgrw, float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,3*n);
    FFUNC(chesvxx,CHESVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, s, b, &ldb, x, &ldx, rcond,
                              rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhesvxx( char fact, char uplo, integer n, integer nrhs,
                     dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                     integer* ipiv, char* equed, double* s, dcomplex* b,
                     integer ldb, dcomplex* x, integer ldx, double* rcond,
                     double* rpvgrw, double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,3*n);
    FFUNC(zhesvxx,ZHESVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, s, b, &ldb, x, &ldx, rcond,
                              rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sspsv( char uplo, integer n, integer nrhs, float* ap,
                   integer* ipiv, float* b, integer ldb )
{
    integer info;
    FFUNC(sspsv,SSPSV)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dspsv( char uplo, integer n, integer nrhs, double* ap,
                   integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dspsv,DSPSV)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_cspsv( char uplo, integer n, integer nrhs, scomplex* ap,
                   integer* ipiv, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(cspsv,CSPSV)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zspsv( char uplo, integer n, integer nrhs, dcomplex* ap,
                   integer* ipiv, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zspsv,ZSPSV)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_chpsv( char uplo, integer n, integer nrhs, scomplex* ap,
                   integer* ipiv, scomplex* b, integer ldb )
{
    integer info;
    FFUNC(chpsv,CHPSV)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_zhpsv( char uplo, integer n, integer nrhs, dcomplex* ap,
                   integer* ipiv, dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zhpsv,ZHPSV)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_sspsvx( char fact, char uplo, integer n, integer nrhs,
                    const float* ap, float* afp, integer* ipiv,
                    const float* b, integer ldb, float* x,
                    integer ldx, float* rcond, float* ferr, float* berr)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sspsvx,SSPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dspsvx( char fact, char uplo, integer n, integer nrhs,
                    const double* ap, double* afp, integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* rcond, double* ferr, double* berr)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dspsvx,DSPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cspsvx( char fact, char uplo, integer n, integer nrhs,
                    const scomplex* ap, scomplex* afp, integer* ipiv,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float* rcond, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(cspsvx,CSPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zspsvx( char fact, char uplo, integer n, integer nrhs,
                    const dcomplex* ap, dcomplex* afp, integer* ipiv,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double* rcond, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zspsvx,ZSPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_chpsvx( char fact, char uplo, integer n, integer nrhs,
                    const scomplex* ap, scomplex* afp, integer* ipiv,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float* rcond, float* ferr, float* berr)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(chpsvx,CHPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhpsvx( char fact, char uplo, integer n, integer nrhs,
                    const dcomplex* ap, dcomplex* afp, integer* ipiv,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double* rcond, double* ferr, double* berr)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zhpsvx,ZHPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgeqrf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgeqrf,SGEQRF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgeqrf,SGEQRF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgeqrf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgeqrf,DGEQRF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgeqrf,DGEQRF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgeqrf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cgeqrf,CGEQRF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgeqrf,CGEQRF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgeqrf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zgeqrf,ZGEQRF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgeqrf,ZGEQRF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sgeqpf( integer m, integer n, float* a, integer lda,
                    integer* jpvt, float* tau)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    FFUNC(sgeqpf,SGEQPF)(&m, &n, a, &lda, jpvt, tau, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dgeqpf( integer m, integer n, double* a, integer lda,
                    integer* jpvt, double* tau)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    FFUNC(dgeqpf,DGEQPF)(&m, &n, a, &lda, jpvt, tau, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cgeqpf( integer m, integer n, scomplex* a, integer lda,
                    integer* jpvt, scomplex* tau)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgeqpf,CGEQPF)(&m, &n, a, &lda, jpvt, tau, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgeqpf( integer m, integer n, dcomplex* a, integer lda,
                    integer* jpvt, dcomplex* tau)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgeqpf,ZGEQPF)(&m, &n, a, &lda, jpvt, tau, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgeqp3( integer m, integer n, float* a, integer lda,
                    integer* jpvt, float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgeqp3,SGEQP3)(&m, &n, a, &lda, jpvt, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgeqp3,SGEQP3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgeqp3( integer m, integer n, double* a, integer lda,
                    integer* jpvt, double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgeqp3,DGEQP3)(&m, &n, a, &lda, jpvt, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgeqp3,DGEQP3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgeqp3( integer m, integer n, scomplex* a, integer lda,
                    integer* jpvt, scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgeqp3,CGEQP3)(&m, &n, a, &lda, jpvt, tau, &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgeqp3,CGEQP3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgeqp3( integer m, integer n, dcomplex* a, integer lda,
                    integer* jpvt, dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgeqp3,ZGEQP3)(&m, &n, a, &lda, jpvt, tau, &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgeqp3,ZGEQP3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sorgqr( integer m, integer n, integer k, float* a,
                    integer lda, const float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sorgqr,SORGQR)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sorgqr,SORGQR)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dorgqr( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dorgqr,DORGQR)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dorgqr,DORGQR)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cungqr( integer m, integer n, integer k, scomplex* a,
                    integer lda, const scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cungqr,CUNGQR)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cungqr,CUNGQR)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zungqr( integer m, integer n, integer k, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zungqr,ZUNGQR)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zungqr,ZUNGQR)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sormqr( char side, char trans, integer m, integer n,
                    integer k, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sormqr,SORMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sormqr,SORMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dormqr( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dormqr,DORMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dormqr,DORMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunmqr( char side, char trans, integer m, integer n,
                    integer k, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunmqr,CUNMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunmqr,CUNMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunmqr( char side, char trans, integer m, integer n,
                    integer k, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunmqr,ZUNMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunmqr,ZUNMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sgelqf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgelqf,SGELQF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgelqf,SGELQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgelqf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgelqf,DGELQF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgelqf,DGELQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgelqf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cgelqf,CGELQF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgelqf,CGELQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgelqf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zgelqf,ZGELQF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgelqf,ZGELQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sorglq( integer m, integer n, integer k, float* a,
                    integer lda, const float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sorglq,SORGLQ)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sorglq,SORGLQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dorglq( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dorglq,DORGLQ)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dorglq,DORGLQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunglq( integer m, integer n, integer k, scomplex* a,
                    integer lda, const scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunglq,CUNGLQ)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunglq,CUNGLQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunglq( integer m, integer n, integer k, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunglq,ZUNGLQ)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunglq,ZUNGLQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sormlq( char side, char trans, integer m, integer n,
                    integer k, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sormlq,SORMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sormlq,SORMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dormlq( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dormlq,DORMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dormlq,DORMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunmlq( char side, char trans, integer m, integer n,
                    integer k, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunmlq,CUNMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunmlq,CUNMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunmlq( char side, char trans, integer m, integer n,
                    integer k, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunmlq,ZUNMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunmlq,ZUNMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sgeqlf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgeqlf,SGEQLF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgeqlf,SGEQLF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgeqlf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgeqlf,DGEQLF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgeqlf,DGEQLF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgeqlf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cgeqlf,CGEQLF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgeqlf,CGEQLF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgeqlf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zgeqlf,ZGEQLF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgeqlf,ZGEQLF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sorgql( integer m, integer n, integer k, float* a,
                    integer lda, const float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sorgql,SORGQL)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sorgql,SORGQL)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dorgql( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dorgql,DORGQL)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dorgql,DORGQL)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cungql( integer m, integer n, integer k, scomplex* a,
                    integer lda, const scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cungql,CUNGQL)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cungql,CUNGQL)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zungql( integer m, integer n, integer k, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zungql,ZUNGQL)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zungql,ZUNGQL)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sormql( char side, char trans, integer m, integer n,
                    integer k, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sormql,SORMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sormql,SORMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dormql( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dormql,DORMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dormql,DORMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunmql( char side, char trans, integer m, integer n,
                    integer k, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunmql,CUNMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunmql,CUNMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunmql( char side, char trans, integer m, integer n,
                    integer k, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunmql,ZUNMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunmql,ZUNMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sgerqf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgerqf,SGERQF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgerqf,SGERQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgerqf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgerqf,DGERQF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgerqf,DGERQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgerqf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cgerqf,CGERQF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgerqf,CGERQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgerqf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zgerqf,ZGERQF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgerqf,ZGERQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sorgrq( integer m, integer n, integer k, float* a,
                    integer lda, const float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sorgrq,SORGRQ)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sorgrq,SORGRQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dorgrq( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dorgrq,DORGRQ)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dorgrq,DORGRQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cungrq( integer m, integer n, integer k, scomplex* a,
                    integer lda, const scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cungrq,CUNGRQ)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cungrq,CUNGRQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zungrq( integer m, integer n, integer k, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zungrq,ZUNGRQ)(&m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zungrq,ZUNGRQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sormrq( char side, char trans, integer m, integer n,
                    integer k, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sormrq,SORMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sormrq,SORMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dormrq( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dormrq,DORMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dormrq,DORMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunmrq( char side, char trans, integer m, integer n,
                    integer k, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunmrq,CUNMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunmrq,CUNMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunmrq( char side, char trans, integer m, integer n,
                    integer k, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunmrq,ZUNMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunmrq,ZUNMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_stzrzf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(stzrzf,STZRZF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(stzrzf,STZRZF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dtzrzf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dtzrzf,DTZRZF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dtzrzf,DTZRZF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_ctzrzf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(ctzrzf,CTZRZF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(ctzrzf,CTZRZF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_ztzrzf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(ztzrzf,ZTZRZF)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(ztzrzf,ZTZRZF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sormrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const float* a,
                    integer lda, const float* tau, float* c,
                    integer ldc)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sormrz,SORMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sormrz,SORMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dormrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dormrz,DORMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dormrz,DORMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunmrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const scomplex* a,
                    integer lda, const scomplex* tau, scomplex* c,
                    integer ldc)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunmrz,CUNMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunmrz,CUNMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunmrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const dcomplex* a,
                    integer lda, const dcomplex* tau, dcomplex* c,
                    integer ldc)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunmrz,ZUNMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunmrz,ZUNMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sggqrf( integer n, integer m, integer p, float* a,
                    integer lda, float* taua, float* b, integer ldb,
                    float* taub)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sggqrf,SGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sggqrf,SGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dggqrf( integer n, integer m, integer p, double* a,
                    integer lda, double* taua, double* b, integer ldb,
                    double* taub)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dggqrf,DGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dggqrf,DGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cggqrf( integer n, integer m, integer p, scomplex* a,
                    integer lda, scomplex* taua, scomplex* b, integer ldb,
                    scomplex* taub)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cggqrf,CGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cggqrf,CGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zggqrf( integer n, integer m, integer p, dcomplex* a,
                    integer lda, dcomplex* taua, dcomplex* b, integer ldb,
                    dcomplex* taub)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zggqrf,ZGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zggqrf,ZGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sggrqf( integer m, integer p, integer n, float* a,
                    integer lda, float* taua, float* b, integer ldb,
                    float* taub)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sggrqf,SGGRQF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sggrqf,SGGRQF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dggrqf( integer m, integer p, integer n, double* a,
                    integer lda, double* taua, double* b, integer ldb,
                    double* taub)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dggrqf,DGGRQF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dggrqf,DGGRQF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cggrqf( integer m, integer p, integer n, scomplex* a,
                    integer lda, scomplex* taua, scomplex* b, integer ldb,
                    scomplex* taub)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cggrqf,CGGRQF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cggrqf,CGGRQF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zggrqf( integer m, integer p, integer n, dcomplex* a,
                    integer lda, dcomplex* taua, dcomplex* b, integer ldb,
                    dcomplex* taub)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zggrqf,ZGGRQF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zggrqf,ZGGRQF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sgebrd( integer m, integer n, float* a, integer lda,
                    float* d, float* e, float* tauq, float* taup)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgebrd,SGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgebrd,SGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgebrd( integer m, integer n, double* a, integer lda,
                    double* d, double* e, double* tauq, double* taup)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgebrd,DGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgebrd,DGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgebrd( integer m, integer n, scomplex* a, integer lda,
                    float* d, float* e, scomplex* tauq, scomplex* taup)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cgebrd,CGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgebrd,CGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgebrd( integer m, integer n, dcomplex* a, integer lda,
                    double* d, double* e, dcomplex* tauq, dcomplex* taup)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zgebrd,ZGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgebrd,ZGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sgbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, float* ab,
                    integer ldab, float* d, float* e, float* q,
                    integer ldq, float* pt, integer ldpt, float* c,
                    integer ldc)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*MAX(m,n));
    FFUNC(sgbbrd,SGBBRD)(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dgbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, double* ab,
                    integer ldab, double* d, double* e, double* q,
                    integer ldq, double* pt, integer ldpt, double* c,
                    integer ldc)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*MAX(m,n));
    FFUNC(dgbbrd,DGBBRD)(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cgbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, scomplex* ab,
                    integer ldab, float* d, float* e, scomplex* q,
                    integer ldq, scomplex* pt, integer ldpt, scomplex* c,
                    integer ldc)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,MAX(m,n));
    float *rwork = SAFE_MALLOC(float,MAX(m,n));
    FFUNC(cgbbrd,CGBBRD)(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, dcomplex* ab,
                    integer ldab, double* d, double* e, dcomplex* q,
                    integer ldq, dcomplex* pt, integer ldpt, dcomplex* c,
                    integer ldc)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,MAX(m,n));
    double *rwork = SAFE_MALLOC(double,MAX(m,n));
    FFUNC(zgbbrd,ZGBBRD)(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sorgbr( char vect, integer m, integer n, integer k,
                    float* a, integer lda, const float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sorgbr,SORGBR)(&vect, &m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sorgbr,SORGBR)(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dorgbr( char vect, integer m, integer n, integer k,
                    double* a, integer lda, const double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dorgbr,DORGBR)(&vect, &m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dorgbr,DORGBR)(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cungbr( char vect, integer m, integer n, integer k,
                    scomplex* a, integer lda, const scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cungbr,CUNGBR)(&vect, &m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cungbr,CUNGBR)(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zungbr( char vect, integer m, integer n, integer k,
                    dcomplex* a, integer lda, const dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zungbr,ZUNGBR)(&vect, &m, &n, &k, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zungbr,ZUNGBR)(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sormbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const float* a,
                    integer lda, const float* tau, float* c,
                    integer ldc)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sormbr,SORMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sormbr,SORMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dormbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dormbr,DORMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dormbr,DORMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunmbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const scomplex* a,
                    integer lda, const scomplex* tau, scomplex* c,
                    integer ldc)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunmbr,CUNMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunmbr,CUNMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunmbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const dcomplex* a,
                    integer lda, const dcomplex* tau, dcomplex* c,
                    integer ldc)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunmbr,ZUNMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunmbr,ZUNMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sbdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, float* d, float* e,
                    float* vt, integer ldvt, float* u, integer ldu,
                    float* c, integer ldc)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    FFUNC(sbdsqr,SBDSQR)(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dbdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, double* d, double* e,
                    double* vt, integer ldvt, double* u, integer ldu,
                    double* c, integer ldc)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    FFUNC(dbdsqr,DBDSQR)(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cbdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, float* d, float* e,
                    scomplex* vt, integer ldvt, scomplex* u, integer ldu,
                    scomplex* c, integer ldc)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    FFUNC(cbdsqr,CBDSQR)(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zbdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, double* d, double* e,
                    dcomplex* vt, integer ldvt, dcomplex* u, integer ldu,
                    dcomplex* c, integer ldc)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    FFUNC(zbdsqr,ZBDSQR)(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_sbdsdc( char uplo, char compq, integer n, float* d,
                    float* e, float* u, integer ldu, float* vt,
                    integer ldvt, float* q, integer* iq)
{
    integer info;
    float *work;
    if (compq == 'N')
    {
        work = SAFE_MALLOC(float,4*n);
    }
    else if (compq == 'P')
    {
        work = SAFE_MALLOC(float,6*n);
    }
    else
    {
        work = SAFE_MALLOC(float,3*n*n+4*n);
    }
    integer *iwork = SAFE_MALLOC(integer,8*n);
    FFUNC(sbdsdc,SBDSDC)(&uplo, &compq, &n, d, e, u, &ldu, vt, &ldvt, q, iq, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dbdsdc( char uplo, char compq, integer n, double* d,
                    double* e, double* u, integer ldu, double* vt,
                    integer ldvt, double* q, integer* iq)
{
    integer info;
    double *work;
    if (compq == 'N')
    {
        work = SAFE_MALLOC(double,4*n);
    }
    else if (compq == 'P')
    {
        work = SAFE_MALLOC(double,6*n);
    }
    else
    {
        work = SAFE_MALLOC(double,3*n*n+4*n);
    }
    integer *iwork = SAFE_MALLOC(integer,8*n);
    FFUNC(dbdsdc,DBDSDC)(&uplo, &compq, &n, d, e, u, &ldu, vt, &ldvt, q, iq, work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ssytrd( char uplo, integer n, float* a, integer lda,
                    float* d, float* e, float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(ssytrd,SSYTRD)(&uplo, &n, a, &lda, d, e, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(ssytrd,SSYTRD)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dsytrd( char uplo, integer n, double* a, integer lda,
                    double* d, double* e, double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dsytrd,DSYTRD)(&uplo, &n, a, &lda, d, e, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dsytrd,DSYTRD)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_chetrd( char uplo, integer n, scomplex* a, integer lda,
                    float* d, float* e, scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(chetrd,CHETRD)(&uplo, &n, a, &lda, d, e, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(chetrd,CHETRD)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zhetrd( char uplo, integer n, dcomplex* a, integer lda,
                    double* d, double* e, dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zhetrd,ZHETRD)(&uplo, &n, a, &lda, d, e, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zhetrd,ZHETRD)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sorgtr( char uplo, integer n, float* a, integer lda,
                    const float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sorgtr,SORGTR)(&uplo, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sorgtr,SORGTR)(&uplo, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dorgtr( char uplo, integer n, double* a, integer lda,
                    const double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dorgtr,DORGTR)(&uplo, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dorgtr,DORGTR)(&uplo, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cungtr( char uplo, integer n, scomplex* a, integer lda,
                    const scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cungtr,CUNGTR)(&uplo, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cungtr,CUNGTR)(&uplo, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zungtr( char uplo, integer n, dcomplex* a, integer lda,
                    const dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zungtr,ZUNGTR)(&uplo, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zungtr,ZUNGTR)(&uplo, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sormtr( char side, char uplo, char trans, integer m,
                    integer n, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sormtr,SORMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sormtr,SORMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dormtr( char side, char uplo, char trans, integer m,
                    integer n, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dormtr,DORMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dormtr,DORMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunmtr( char side, char uplo, char trans, integer m,
                    integer n, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunmtr,CUNMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunmtr,CUNMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunmtr( char side, char uplo, char trans, integer m,
                    integer n, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunmtr,ZUNMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunmtr,ZUNMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_ssptrd( char uplo, integer n, float* ap, float* d, float* e,
                    float* tau )
{
    integer info;
    FFUNC(ssptrd,SSPTRD)(&uplo, &n, ap, d, e, tau, &info);
    return info;
}

static inline integer c_dsptrd( char uplo, integer n, double* ap, double* d, double* e,
                    double* tau )
{
    integer info;
    FFUNC(dsptrd,DSPTRD)(&uplo, &n, ap, d, e, tau, &info);
    return info;
}

static inline integer c_chptrd( char uplo, integer n, scomplex* ap, float* d, float* e,
                    scomplex* tau )
{
    integer info;
    FFUNC(chptrd,CHPTRD)(&uplo, &n, ap, d, e, tau, &info);
    return info;
}

static inline integer c_zhptrd( char uplo, integer n, dcomplex* ap, double* d, double* e,
                    dcomplex* tau )
{
    integer info;
    FFUNC(zhptrd,ZHPTRD)(&uplo, &n, ap, d, e, tau, &info);
    return info;
}

static inline integer c_sopgtr( char uplo, integer n, const float* ap,
                    const float* tau, float* q, integer ldq)
{
    integer info;
    float *work = SAFE_MALLOC(float,n);
    FFUNC(sopgtr,SOPGTR)(&uplo, &n, ap, tau, q, &ldq, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dopgtr( char uplo, integer n, const double* ap,
                    const double* tau, double* q, integer ldq)
{
    integer info;
    double *work = SAFE_MALLOC(double,n);
    FFUNC(dopgtr,DOPGTR)(&uplo, &n, ap, tau, q, &ldq, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cupgtr( char uplo, integer n, const scomplex* ap,
                    const scomplex* tau, scomplex* q, integer ldq)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    FFUNC(cupgtr,CUPGTR)(&uplo, &n, ap, tau, q, &ldq, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zupgtr( char uplo, integer n, const dcomplex* ap,
                    const dcomplex* tau, dcomplex* q, integer ldq)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    FFUNC(zupgtr,ZUPGTR)(&uplo, &n, ap, tau, q, &ldq, work, &info);
    FREE(work);
    return info;
}

static inline integer c_sopmtr( char side, char uplo, char trans, integer m,
                    integer n, const float* ap, const float* tau,
                    float* c, integer ldc)
{
    integer info;
    float *work = SAFE_MALLOC(float,n);
    FFUNC(sopmtr,SOPMTR)(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dopmtr( char side, char uplo, char trans, integer m,
                    integer n, const double* ap, const double* tau,
                    double* c, integer ldc)
{
    integer info;
    double *work = SAFE_MALLOC(double,n);
    FFUNC(dopmtr,DOPMTR)(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cupmtr( char side, char uplo, char trans, integer m,
                    integer n, const scomplex* ap, const scomplex* tau,
                    scomplex* c, integer ldc)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    FFUNC(cupmtr,CUPMTR)(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zupmtr( char side, char uplo, char trans, integer m,
                    integer n, const dcomplex* ap, const dcomplex* tau,
                    dcomplex* c, integer ldc)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    FFUNC(zupmtr,ZUPMTR)(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info);
    FREE(work);
    return info;
}

static inline integer c_ssbtrd( char vect, char uplo, integer n, integer kd,
                    float* ab, integer ldab, float* d, float* e,
                    float* q, integer ldq)
{
    integer info;
    float *work = SAFE_MALLOC(float,n);
    FFUNC(ssbtrd,SSBTRD)(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dsbtrd( char vect, char uplo, integer n, integer kd,
                    double* ab, integer ldab, double* d, double* e,
                    double* q, integer ldq)
{
    integer info;
    double *work = SAFE_MALLOC(double,n);
    FFUNC(dsbtrd,DSBTRD)(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info);
    FREE(work);
    return info;
}

static inline integer c_chbtrd( char vect, char uplo, integer n, integer kd,
                    scomplex* ab, integer ldab, float* d, float* e,
                    scomplex* q, integer ldq)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    FFUNC(chbtrd,CHBTRD)(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zhbtrd( char vect, char uplo, integer n, integer kd,
                    dcomplex* ab, integer ldab, double* d, double* e,
                    dcomplex* q, integer ldq)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    FFUNC(zhbtrd,ZHBTRD)(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info);
    FREE(work);
    return info;
}

static inline integer c_ssterf( integer n, float* d, float* e )
{
    integer info;
    FFUNC(ssterf,SSTERF)(&n, d, e, &info);
    return info;
}

static inline integer c_dsterf( integer n, double* d, double* e )
{
    integer info;
    FFUNC(dsterf,DSTERF)(&n, d, e, &info);
    return info;
}

static inline integer c_ssteqr( char compz, integer n, float* d, float* e, float* z,
                    integer ldz)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    FFUNC(ssteqr,SSTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dsteqr( char compz, integer n, double* d, double* e, double* z,
                    integer ldz)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    FFUNC(dsteqr,DSTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_csteqr( char compz, integer n, float* d, float* e, scomplex* z,
                    integer ldz)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    FFUNC(csteqr,CSTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zsteqr( char compz, integer n, double* d, double* e, dcomplex* z,
                    integer ldz)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    FFUNC(zsteqr,ZSTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_sstemr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, integer* m, float* w, float* z,
                    integer ldz, integer nzc, integer* isuppz,
                    logical* tryrac)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(sstemr,SSTEMR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac,
                            &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(sstemr,SSTEMR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac,
                            work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dstemr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, integer* m, double* w, double* z,
                    integer ldz, integer nzc, integer* isuppz,
                    logical* tryrac)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dstemr,DSTEMR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac,
                            &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dstemr,DSTEMR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac,
                            work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cstemr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, integer* m, float* w, scomplex* z,
                    integer ldz, integer nzc, integer* isuppz,
                    logical* tryrac)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(cstemr,CSTEMR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac,
                            &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(cstemr,CSTEMR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac,
                            work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_zstemr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, integer* m, double* w, dcomplex* z,
                    integer ldz, integer nzc, integer* isuppz,
                    logical* tryrac)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(zstemr,ZSTEMR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac,
                            &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zstemr,ZSTEMR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac,
                            work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_sstedc( char compz, integer n, float* d, float* e, float* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(sstedc,SSTEDC)(&compz, &n, d, e, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(sstedc,SSTEDC)(&compz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dstedc( char compz, integer n, double* d, double* e, double* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dstedc,DSTEDC)(&compz, &n, d, e, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dstedc,DSTEDC)(&compz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cstedc( char compz, integer n, float* d, float* e, scomplex* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    scomplex *work, query;
    float *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(cstedc,CSTEDC)(&compz, &n, d, e, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)crealf(query);
    lwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(scomplex,lwork);
    rwork = SAFE_MALLOC(float,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(cstedc,CSTEDC)(&compz, &n, d, e, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zstedc( char compz, integer n, double* d, double* e, dcomplex* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    dcomplex *work, query;
    double *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(zstedc,ZSTEDC)(&compz, &n, d, e, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)creal(query);
    lwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(dcomplex,lwork);
    rwork = SAFE_MALLOC(double,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zstedc,ZSTEDC)(&compz, &n, d, e, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_sstegr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    float* z, integer ldz, integer* isuppz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(sstegr,SSTEGR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(sstegr,SSTEGR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dstegr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, integer* isuppz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dstegr,DSTEGR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dstegr,DSTEGR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cstegr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    scomplex* z, integer ldz, integer* isuppz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(cstegr,CSTEGR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(cstegr,CSTEGR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_zstegr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    dcomplex* z, integer ldz, integer* isuppz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(zstegr,ZSTEGR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zstegr,ZSTEGR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_spteqr( char compz, integer n, float* d, float* e, float* z,
                    integer ldz)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    FFUNC(spteqr,SPTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dpteqr( char compz, integer n, double* d, double* e, double* z,
                    integer ldz)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    FFUNC(dpteqr,DPTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cpteqr( char compz, integer n, float* d, float* e, scomplex* z,
                    integer ldz)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    FFUNC(cpteqr,CPTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zpteqr( char compz, integer n, double* d, double* e, dcomplex* z,
                    integer ldz)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    FFUNC(zpteqr,ZPTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_sstebz( char range, char order, integer n, float vl,
                    float vu, integer il, integer iu, float abstol,
                    const float* d, const float* e, integer* m,
                    integer* nsplit, float* w, integer* iblock,
                    integer* isplit)
{
    integer info;
    float *work = SAFE_MALLOC(float,4*n);
    integer *iwork = SAFE_MALLOC(integer,3*n);
    FFUNC(sstebz,SSTEBZ)(&range, &order, &n, &vl, &vu, &il, &iu, &abstol, d, e, m, nsplit, w, iblock, isplit,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dstebz( char range, char order, integer n, double vl,
                    double vu, integer il, integer iu, double abstol,
                    const double* d, const double* e, integer* m,
                    integer* nsplit, double* w, integer* iblock,
                    integer* isplit)
{
    integer info;
    double *work = SAFE_MALLOC(double,4*n);
    integer *iwork = SAFE_MALLOC(integer,3*n);
    FFUNC(dstebz,DSTEBZ)(&range, &order, &n, &vl, &vu, &il, &iu, &abstol, d, e, m, nsplit, w, iblock, isplit,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_sstein( integer n, const float* d, const float* e,
                    integer m, const float* w, const integer* iblock,
                    const integer* isplit, float* z, integer ldz,
                    integer* ifailv )
{
    integer info;
    float *work = SAFE_MALLOC(float,5*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sstein,SSTEIN)(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifailv, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dstein( integer n, const double* d, const double* e,
                    integer m, const double* w, const integer* iblock,
                    const integer* isplit, double* z, integer ldz,
                    integer* ifailv )
{
    integer info;
    double *work = SAFE_MALLOC(double,5*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dstein,DSTEIN)(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifailv, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cstein( integer n, const float* d, const float* e,
                    integer m, const float* w, const integer* iblock,
                    const integer* isplit, scomplex* z, integer ldz,
                    integer* ifailv )
{
    integer info;
    float *work = SAFE_MALLOC(float,5*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(cstein,CSTEIN)(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifailv, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_zstein( integer n, const double* d, const double* e,
                    integer m, const double* w, const integer* iblock,
                    const integer* isplit, dcomplex* z, integer ldz,
                    integer* ifailv )
{
    integer info;
    double *work = SAFE_MALLOC(double,5*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(zstein,ZSTEIN)(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifailv, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_sdisna( char job, integer m, integer n, const float* d,
                    float* sep )
{
    integer info;
    FFUNC(sdisna,SDISNA)(&job, &m, &n, d, sep, &info);
    return info;
}

static inline integer c_ddisna( char job, integer m, integer n, const double* d,
                    double* sep )
{
    integer info;
    FFUNC(ddisna,DDISNA)(&job, &m, &n, d, sep, &info);
    return info;
}

static inline integer c_ssygst( integer itype, char uplo, integer n, float* a,
                    integer lda, const float* b, integer ldb )
{
    integer info;
    FFUNC(ssygst,SSYGST)(&itype, &uplo, &n, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_dsygst( integer itype, char uplo, integer n, double* a,
                    integer lda, const double* b, integer ldb )
{
    integer info;
    FFUNC(dsygst,DSYGST)(&itype, &uplo, &n, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_chegst( integer itype, char uplo, integer n, scomplex* a,
                    integer lda, const scomplex* b, integer ldb )
{
    integer info;
    FFUNC(chegst,CHEGST)(&itype, &uplo, &n, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_zhegst( integer itype, char uplo, integer n, dcomplex* a,
                    integer lda, const dcomplex* b, integer ldb )
{
    integer info;
    FFUNC(zhegst,ZHEGST)(&itype, &uplo, &n, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_sspgst( integer itype, char uplo, integer n, float* ap,
                    const float* bp )
{
    integer info;
    FFUNC(sspgst,SSPGST)(&itype, &uplo, &n, ap, bp, &info);
    return info;
}

static inline integer c_dspgst( integer itype, char uplo, integer n, double* ap,
                    const double* bp )
{
    integer info;
    FFUNC(dspgst,DSPGST)(&itype, &uplo, &n, ap, bp, &info);
    return info;
}

static inline integer c_chpgst( integer itype, char uplo, integer n, scomplex* ap,
                    const scomplex* bp )
{
    integer info;
    FFUNC(chpgst,CHPGST)(&itype, &uplo, &n, ap, bp, &info);
    return info;
}

static inline integer c_zhpgst( integer itype, char uplo, integer n, dcomplex* ap,
                    const dcomplex* bp )
{
    integer info;
    FFUNC(zhpgst,ZHPGST)(&itype, &uplo, &n, ap, bp, &info);
    return info;
}

static inline integer c_ssbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, float* ab, integer ldab,
                    const float* bb, integer ldbb, float* x,
                    integer ldx)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    FFUNC(ssbgst,SSBGST)(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dsbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, double* ab, integer ldab,
                    const double* bb, integer ldbb, double* x,
                    integer ldx)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    FFUNC(dsbgst,DSBGST)(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, &info);
    FREE(work);
    return info;
}

static inline integer c_chbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, scomplex* ab, integer ldab,
                    const scomplex* bb, integer ldbb, scomplex* x,
                    integer ldx)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(chbgst,CHBGST)(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, dcomplex* ab, integer ldab,
                    const dcomplex* bb, integer ldbb, dcomplex* x,
                    integer ldx)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zhbgst,ZHBGST)(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_spbstf( char uplo, integer n, integer kb, float* bb,
                    integer ldbb )
{
    integer info;
    FFUNC(spbstf,SPBSTF)(&uplo, &n, &kb, bb, &ldbb, &info);
    return info;
}

static inline integer c_dpbstf( char uplo, integer n, integer kb, double* bb,
                    integer ldbb )
{
    integer info;
    FFUNC(dpbstf,DPBSTF)(&uplo, &n, &kb, bb, &ldbb, &info);
    return info;
}

static inline integer c_cpbstf( char uplo, integer n, integer kb, scomplex* bb,
                    integer ldbb )
{
    integer info;
    FFUNC(cpbstf,CPBSTF)(&uplo, &n, &kb, bb, &ldbb, &info);
    return info;
}

static inline integer c_zpbstf( char uplo, integer n, integer kb, dcomplex* bb,
                    integer ldbb )
{
    integer info;
    FFUNC(zpbstf,ZPBSTF)(&uplo, &n, &kb, bb, &ldbb, &info);
    return info;
}

static inline integer c_sgehrd( integer n, integer ilo, integer ihi, float* a,
                    integer lda, float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgehrd,SGEHRD)(&n, &ilo, &ihi, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgehrd,SGEHRD)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgehrd( integer n, integer ilo, integer ihi, double* a,
                    integer lda, double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgehrd,DGEHRD)(&n, &ilo, &ihi, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgehrd,DGEHRD)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgehrd( integer n, integer ilo, integer ihi, scomplex* a,
                    integer lda, scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cgehrd,CGEHRD)(&n, &ilo, &ihi, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgehrd,CGEHRD)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgehrd( integer n, integer ilo, integer ihi, dcomplex* a,
                    integer lda, dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zgehrd,ZGEHRD)(&n, &ilo, &ihi, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgehrd,ZGEHRD)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sorghr( integer n, integer ilo, integer ihi, float* a,
                    integer lda, const float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sorghr,SORGHR)(&n, &ilo, &ihi, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sorghr,SORGHR)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dorghr( integer n, integer ilo, integer ihi, double* a,
                    integer lda, const double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dorghr,DORGHR)(&n, &ilo, &ihi, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dorghr,DORGHR)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunghr( integer n, integer ilo, integer ihi, scomplex* a,
                    integer lda, const scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunghr,CUNGHR)(&n, &ilo, &ihi, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunghr,CUNGHR)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunghr( integer n, integer ilo, integer ihi, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunghr,ZUNGHR)(&n, &ilo, &ihi, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunghr,ZUNGHR)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sormhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const float* a,
                    integer lda, const float* tau, float* c,
                    integer ldc)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sormhr,SORMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sormhr,SORMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dormhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dormhr,DORMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dormhr,DORMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cunmhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const scomplex* a,
                    integer lda, const scomplex* tau, scomplex* c,
                    integer ldc)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cunmhr,CUNMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cunmhr,CUNMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zunmhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const dcomplex* a,
                    integer lda, const dcomplex* tau, dcomplex* c,
                    integer ldc)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zunmhr,ZUNMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zunmhr,ZUNMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sgebal( char job, integer n, float* a, integer lda,
                    integer* ilo, integer* ihi, float* scale )
{
    integer info;
    FFUNC(sgebal,SGEBAL)(&job, &n, a, &lda, ilo, ihi, scale, &info);
    return info;
}

static inline integer c_dgebal( char job, integer n, double* a, integer lda,
                    integer* ilo, integer* ihi, double* scale )
{
    integer info;
    FFUNC(dgebal,DGEBAL)(&job, &n, a, &lda, ilo, ihi, scale, &info);
    return info;
}

static inline integer c_cgebal( char job, integer n, scomplex* a, integer lda,
                    integer* ilo, integer* ihi, float* scale )
{
    integer info;
    FFUNC(cgebal,CGEBAL)(&job, &n, a, &lda, ilo, ihi, scale, &info);
    return info;
}

static inline integer c_zgebal( char job, integer n, dcomplex* a, integer lda,
                    integer* ilo, integer* ihi, double* scale )
{
    integer info;
    FFUNC(zgebal,ZGEBAL)(&job, &n, a, &lda, ilo, ihi, scale, &info);
    return info;
}

static inline integer c_sgebak( char job, char side, integer n, integer ilo,
                    integer ihi, const float* scale, integer m,
                    float* v, integer ldv )
{
    integer info;
    FFUNC(sgebak,SGEBAK)(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info);
    return info;
}

static inline integer c_dgebak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* scale, integer m,
                    double* v, integer ldv )
{
    integer info;
    FFUNC(dgebak,DGEBAK)(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info);
    return info;
}

static inline integer c_cgebak( char job, char side, integer n, integer ilo,
                    integer ihi, const float* scale, integer m,
                    scomplex* v, integer ldv )
{
    integer info;
    FFUNC(cgebak,CGEBAK)(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info);
    return info;
}

static inline integer c_zgebak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* scale, integer m,
                    dcomplex* v, integer ldv )
{
    integer info;
    FFUNC(zgebak,ZGEBAK)(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info);
    return info;
}

static inline integer c_shseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, float* h, integer ldh, scomplex* w,
                    float* z, integer ldz)
{
    integer info, lwork = -1;
    float *work, query;
    float *wr = SAFE_MALLOC(float,n);
    float *wi = SAFE_MALLOC(float,n);
    FFUNC(shseqr,SHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(shseqr,SHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, work, &lwork, &info);
    scopy(n,wr,1,(float*)w  ,2);
    scopy(n,wi,1,(float*)w+1,2);
    FREE(work);
    FREE(wr);
    FREE(wi);
    return info;
}

static inline integer c_dhseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, double* h, integer ldh, dcomplex* w,
                    double* z, integer ldz)
{
    integer info, lwork = -1;
    double *work, query;
    double *wr = SAFE_MALLOC(double,n);
    double *wi = SAFE_MALLOC(double,n);
    FFUNC(dhseqr,DHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dhseqr,DHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, work, &lwork, &info);
    dcopy(n,wr,1,(double*)w  ,2);
    dcopy(n,wi,1,(double*)w+1,2);
    FREE(work);
    FREE(wr);
    FREE(wi);
    return info;
}

static inline integer c_chseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, scomplex* h, integer ldh, scomplex* w,
                    scomplex* z, integer ldz)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(chseqr,CHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(chseqr,CHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zhseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, dcomplex* h, integer ldh, dcomplex* w,
                    dcomplex* z, integer ldz)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zhseqr,ZHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zhseqr,ZHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_shsein( char job, char eigsrc, char initv,
                    logical* select, integer n, const float* h,
                    integer ldh, const scomplex* w, float* vl,
                    integer ldvl, float* vr, integer ldvr,
                    integer mm, integer* m,
                    integer* ifaill, integer* ifailr )
{
    integer info;
    float *work = SAFE_MALLOC(float,(2+n)*n);
    float *wr = SAFE_MALLOC(float,n);
    float *wi = SAFE_MALLOC(float,n);
    scopy(n,(float*)w  ,2,wr,1);
    scopy(n,(float*)w+1,2,wi,1);
    FFUNC(shsein,SHSEIN)(&job, &eigsrc, &initv, select, &n, h, &ldh, wr, wi, vl, &ldvl, vr, &ldvr, &mm, m, work,
                            ifaill, ifailr, &info);
    FREE(wr);
    FREE(wi);
    FREE(work);
    return info;
}

static inline integer c_dhsein( char job, char eigsrc, char initv,
                    logical* select, integer n, const double* h,
                    integer ldh, const dcomplex* w, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer* m,
                    integer* ifaill, integer* ifailr )
{
    integer info;
    double *work = SAFE_MALLOC(double,(2+n)*n);
    double *wr = SAFE_MALLOC(double,n);
    double *wi = SAFE_MALLOC(double,n);
    dcopy(n,(double*)w  ,2,wr,1);
    dcopy(n,(double*)w+1,2,wi,1);
    FFUNC(dhsein,DHSEIN)(&job, &eigsrc, &initv, select, &n, h, &ldh, wr, wi, vl, &ldvl, vr, &ldvr, &mm, m, work,
                            ifaill, ifailr, &info);
    FREE(wr);
    FREE(wi);
    FREE(work);
    return info;
}

static inline integer c_chsein( char job, char eigsrc, char initv,
                    logical* select, integer n, const scomplex* h,
                    integer ldh, scomplex* w, scomplex* vl,
                    integer ldvl, scomplex* vr, integer ldvr,
                    integer mm, integer* m,
                    integer* ifaill, integer* ifailr )
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(chsein,CHSEIN)(&job, &eigsrc, &initv, select, &n, h, &ldh, w, vl, &ldvl, vr, &ldvr, &mm, m, work,
                            rwork, ifaill, ifailr, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhsein( char job, char eigsrc, char initv,
                    logical* select, integer n, const dcomplex* h,
                    integer ldh, dcomplex* w, dcomplex* vl,
                    integer ldvl, dcomplex* vr, integer ldvr,
                    integer mm, integer* m,
                    integer* ifaill, integer* ifailr )
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zhsein,ZHSEIN)(&job, &eigsrc, &initv, select, &n, h, &ldh, w, vl, &ldvl, vr, &ldvr, &mm, m, work,
                            rwork, ifaill, ifailr, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_strevc( char side, char howmny, const logical* select,
                    integer n, const float* t, integer ldt, float* vl,
                    integer ldvl, float* vr, integer ldvr,
                    integer mm, integer* m)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    FFUNC(strevc,STREVC)(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, m, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dtrevc( char side, char howmny, const logical* select,
                    integer n, const double* t, integer ldt, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer* m)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    FFUNC(dtrevc,DTREVC)(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, m, work, &info);
    FREE(work);
    return info;
}

static inline integer c_ctrevc( char side, char howmny, const logical* select,
                    integer n, const scomplex* t, integer ldt, scomplex* vl,
                    integer ldvl, scomplex* vr, integer ldvr,
                    integer mm, integer* m)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(ctrevc,CTREVC)(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, m, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ztrevc( char side, char howmny, const logical* select,
                    integer n, const dcomplex* t, integer ldt, dcomplex* vl,
                    integer ldvl, dcomplex* vr, integer ldvr,
                    integer mm, integer* m)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(ztrevc,ZTREVC)(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, m, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_strsna( char job, char howmny, const logical* select,
                    integer n, const float* t, integer ldt,
                    const float* vl, integer ldvl, const float* vr,
                    integer ldvr, float* s, float* sep, integer mm,
                    integer* m)
{
    integer info;
    float *work = SAFE_MALLOC(float,n*(n+6));
    integer *iwork = SAFE_MALLOC(integer,2*n);
    FFUNC(strsna,STRSNA)(&job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, m, work, &n, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtrsna( char job, char howmny, const logical* select,
                    integer n, const double* t, integer ldt,
                    const double* vl, integer ldvl, const double* vr,
                    integer ldvr, double* s, double* sep, integer mm,
                    integer* m)
{
    integer info;
    double *work = SAFE_MALLOC(double,n*(n+6));
    integer *iwork = SAFE_MALLOC(integer,2*n);
    FFUNC(dtrsna,DTRSNA)(&job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, m, work, &n, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctrsna( char job, char howmny, const logical* select,
                    integer n, const scomplex* t, integer ldt,
                    const scomplex* vl, integer ldvl, const scomplex* vr,
                    integer ldvr, float* s, float* sep, integer mm,
                    integer* m)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n*(n+6));
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(ctrsna,CTRSNA)(&job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, m, work, &n, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ztrsna( char job, char howmny, const logical* select,
                    integer n, const dcomplex* t, integer ldt,
                    const dcomplex* vl, integer ldvl, const dcomplex* vr,
                    integer ldvr, double* s, double* sep, integer mm,
                    integer* m)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n*(n+6));
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(ztrsna,ZTRSNA)(&job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, m, work, &n, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_strexc( char compq, integer n, float* t, integer ldt,
                    float* q, integer ldq, integer* ifst,
                    integer* ilst)
{
    integer info;
    float *work = SAFE_MALLOC(float,n);
    FFUNC(strexc,STREXC)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dtrexc( char compq, integer n, double* t, integer ldt,
                    double* q, integer ldq, integer* ifst,
                    integer* ilst)
{
    integer info;
    double *work = SAFE_MALLOC(double,n);
    FFUNC(dtrexc,DTREXC)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, work, &info);
    FREE(work);
    return info;
}

static inline integer c_ctrexc( char compq, integer n, scomplex* t, integer ldt,
                    scomplex* q, integer ldq, integer* ifst,
                    integer* ilst)
{
    integer info;
    FFUNC(ctrexc,CTREXC)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, &info);
    return info;
}

static inline integer c_ztrexc( char compq, integer n, dcomplex* t, integer ldt,
                    dcomplex* q, integer ldq, integer* ifst,
                    integer* ilst)
{
    integer info;
    FFUNC(ztrexc,ZTREXC)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, &info);
    return info;
}

static inline integer c_strsen( char job, char compq, const logical* select,
                    integer n, float* t, integer ldt, float* q,
                    integer ldq, scomplex* w, integer* m,
                    float* s, float* sep)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    float *wr = SAFE_MALLOC(float,n);
    float *wi = SAFE_MALLOC(float,n);
    FFUNC(strsen,STRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, m, s, sep, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(strsen,STRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, m, s, sep, work, &lwork, iwork, &liwork, &info);
    scopy(n,wr,1,(float*)w  ,2);
    scopy(n,wi,1,(float*)w+1,2);
    FREE(wr);
    FREE(wi);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtrsen( char job, char compq, const logical* select,
                    integer n, double* t, integer ldt, double* q,
                    integer ldq, dcomplex* w, integer* m,
                    double* s, double* sep)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    double *wr = SAFE_MALLOC(double,n);
    double *wi = SAFE_MALLOC(double,n);
    FFUNC(dtrsen,DTRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, m, s, sep, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dtrsen,DTRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, m, s, sep, work, &lwork, iwork, &liwork, &info);
    dcopy(n,wr,1,(double*)w  ,2);
    dcopy(n,wi,1,(double*)w+1,2);
    FREE(wr);
    FREE(wi);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctrsen( char job, char compq, const logical* select,
                    integer n, scomplex* t, integer ldt, scomplex* q,
                    integer ldq, scomplex* w, integer* m,
                    float* s, float* sep)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(ctrsen,CTRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, w, m, s, sep, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(ctrsen,CTRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, w, m, s, sep, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_ztrsen( char job, char compq, const logical* select,
                    integer n, dcomplex* t, integer ldt, dcomplex* q,
                    integer ldq, dcomplex* w, integer* m,
                    double* s, double* sep)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(ztrsen,ZTRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, w, m, s, sep, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(ztrsen,ZTRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, w, m, s, sep, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_strsyl( char trana, char tranb, integer isgn, integer m,
                    integer n, const float* a, integer lda,
                    const float* b, integer ldb, float* c,
                    integer ldc, float* scale )
{
    integer info;
    FFUNC(strsyl,STRSYL)(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, &info);
    return info;
}

static inline integer c_dtrsyl( char trana, char tranb, integer isgn, integer m,
                    integer n, const double* a, integer lda,
                    const double* b, integer ldb, double* c,
                    integer ldc, double* scale )
{
    integer info;
    FFUNC(dtrsyl,DTRSYL)(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, &info);
    return info;
}

static inline integer c_ctrsyl( char trana, char tranb, integer isgn, integer m,
                    integer n, const scomplex* a, integer lda,
                    const scomplex* b, integer ldb, scomplex* c,
                    integer ldc, float* scale )
{
    integer info;
    FFUNC(ctrsyl,CTRSYL)(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, &info);
    return info;
}

static inline integer c_ztrsyl( char trana, char tranb, integer isgn, integer m,
                    integer n, const dcomplex* a, integer lda,
                    const dcomplex* b, integer ldb, dcomplex* c,
                    integer ldc, double* scale )
{
    integer info;
    FFUNC(ztrsyl,ZTRSYL)(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, &info);
    return info;
}

static inline integer c_sgghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, float* a, integer lda, float* b,
                    integer ldb, float* q, integer ldq, float* z,
                    integer ldz )
{
    integer info;
    FFUNC(sgghrd,SGGHRD)(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info);
    return info;
}

static inline integer c_dgghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, double* a, integer lda, double* b,
                    integer ldb, double* q, integer ldq, double* z,
                    integer ldz )
{
    integer info;
    FFUNC(dgghrd,DGGHRD)(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info);
    return info;
}

static inline integer c_cgghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, scomplex* a, integer lda, scomplex* b,
                    integer ldb, scomplex* q, integer ldq, scomplex* z,
                    integer ldz )
{
    integer info;
    FFUNC(cgghrd,CGGHRD)(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info);
    return info;
}

static inline integer c_zgghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, dcomplex* a, integer lda, dcomplex* b,
                    integer ldb, dcomplex* q, integer ldq, dcomplex* z,
                    integer ldz )
{
    integer info;
    FFUNC(zgghrd,ZGGHRD)(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info);
    return info;
}

static inline integer c_sggbal( char job, integer n, float* a, integer lda,
                    float* b, integer ldb, integer* ilo,
                    integer* ihi, float* lscale, float* rscale)
{
    integer info;
    float *work = SAFE_MALLOC(float,6*n);
    FFUNC(sggbal,SGGBAL)(&job, &n, a, &lda, b, &ldb, ilo, ihi, lscale, rscale, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dggbal( char job, integer n, double* a, integer lda,
                    double* b, integer ldb, integer* ilo,
                    integer* ihi, double* lscale, double* rscale)
{
    integer info;
    double *work = SAFE_MALLOC(double,6*n);
    FFUNC(dggbal,DGGBAL)(&job, &n, a, &lda, b, &ldb, ilo, ihi, lscale, rscale, work, &info);
    FREE(work);
    return info;
}

static inline integer c_cggbal( char job, integer n, scomplex* a, integer lda,
                    scomplex* b, integer ldb, integer* ilo,
                    integer* ihi, float* lscale, float* rscale)
{
    integer info;
    float *work = SAFE_MALLOC(float,6*n);
    FFUNC(cggbal,CGGBAL)(&job, &n, a, &lda, b, &ldb, ilo, ihi, lscale, rscale, work, &info);
    FREE(work);
    return info;
}

static inline integer c_zggbal( char job, integer n, dcomplex* a, integer lda,
                    dcomplex* b, integer ldb, integer* ilo,
                    integer* ihi, double* lscale, double* rscale)
{
    integer info;
    double *work = SAFE_MALLOC(double,6*n);
    FFUNC(zggbal,ZGGBAL)(&job, &n, a, &lda, b, &ldb, ilo, ihi, lscale, rscale, work, &info);
    FREE(work);
    return info;
}

static inline integer c_sggbak( char job, char side, integer n, integer ilo,
                    integer ihi, const float* lscale, const float* rscale,
                    integer m, float* v, integer ldv )
{
    integer info;
    FFUNC(sggbak,SGGBAK)(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info);
    return info;
}

static inline integer c_dggbak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* lscale, const double* rscale,
                    integer m, double* v, integer ldv )
{
    integer info;
    FFUNC(dggbak,DGGBAK)(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info);
    return info;
}

static inline integer c_cggbak( char job, char side, integer n, integer ilo,
                    integer ihi, const float* lscale, const float* rscale,
                    integer m, scomplex* v, integer ldv )
{
    integer info;
    FFUNC(cggbak,CGGBAK)(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info);
    return info;
}

static inline integer c_zggbak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* lscale, const double* rscale,
                    integer m, dcomplex* v, integer ldv )
{
    integer info;
    FFUNC(zggbak,ZGGBAK)(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info);
    return info;
}

static inline integer c_shgeqz( char job, char compq, char compz, integer n,
                    integer ilo, integer ihi, float* h,
                    integer ldh, float* t, integer ldt, scomplex* alpha,
                    float* beta, float* q, integer ldq,
                    float* z, integer ldz)
{
    integer info, lwork = -1;
    float *work, query;
    float *alphar = SAFE_MALLOC(float,n);
    float *alphai = SAFE_MALLOC(float,n);
    FFUNC(shgeqz,SHGEQZ)(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alphar, alphai, beta, q, &ldq, z, &ldz,
                            &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(shgeqz,SHGEQZ)(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alphar, alphai, beta, q, &ldq, z, &ldz,
                            work, &lwork, &info);
    scopy(n,alphar,1,(float*)alpha  ,2);
    scopy(n,alphai,1,(float*)alpha+1,2);
    FREE(alphar);
    FREE(alphai);
    FREE(work);
    return info;
}

static inline integer c_dhgeqz( char job, char compq, char compz, integer n,
                    integer ilo, integer ihi, double* h,
                    integer ldh, double* t, integer ldt, dcomplex* alpha,
                    double* beta, double* q, integer ldq,
                    double* z, integer ldz)
{
    integer info, lwork = -1;
    double *work, query;
    double *alphar = SAFE_MALLOC(double,n);
    double *alphai = SAFE_MALLOC(double,n);
    FFUNC(dhgeqz,DHGEQZ)(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alphar, alphai, beta, q, &ldq, z, &ldz,
                            &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dhgeqz,DHGEQZ)(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alphar, alphai, beta, q, &ldq, z, &ldz,
                            work, &lwork, &info);
    dcopy(n,alphar,1,(double*)alpha  ,2);
    dcopy(n,alphai,1,(double*)alpha+1,2);
    FREE(alphar);
    FREE(alphai);
    FREE(work);
    return info;
}

static inline integer c_chgeqz( char job, char compq, char compz, integer n,
                    integer ilo, integer ihi, scomplex* h,
                    integer ldh, scomplex* t, integer ldt, scomplex* alpha,
                    scomplex* beta, scomplex* q, integer ldq,
                    scomplex* z, integer ldz)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,n);
    FFUNC(chgeqz,CHGEQZ)(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alpha, beta, q, &ldq, z, &ldz,
                            &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(chgeqz,CHGEQZ)(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alpha, beta, q, &ldq, z, &ldz,
                            work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhgeqz( char job, char compq, char compz, integer n,
                    integer ilo, integer ihi, dcomplex* h,
                    integer ldh, dcomplex* t, integer ldt, dcomplex* alpha,
                    dcomplex* beta, dcomplex* q, integer ldq,
                    dcomplex* z, integer ldz)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,n);
    FFUNC(zhgeqz,ZHGEQZ)(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alpha, beta, q, &ldq, z, &ldz,
                            &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zhgeqz,ZHGEQZ)(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alpha, beta, q, &ldq, z, &ldz,
                            work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_stgevc( char side, char howmny, const logical* select,
                    integer n, const float* s, integer lds,
                    const float* p, integer ldp, float* vl,
                    integer ldvl, float* vr, integer ldvr,
                    integer mm, integer* m)
{
    integer info;
    float *work = SAFE_MALLOC(float,6*n);
    FFUNC(stgevc,STGEVC)(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, m, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dtgevc( char side, char howmny, const logical* select,
                    integer n, const double* s, integer lds,
                    const double* p, integer ldp, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer* m)
{
    integer info;
    double *work = SAFE_MALLOC(double,6*n);
    FFUNC(dtgevc,DTGEVC)(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, m, work, &info);
    FREE(work);
    return info;
}

static inline integer c_ctgevc( char side, char howmny, const logical* select,
                    integer n, const scomplex* s, integer lds,
                    const scomplex* p, integer ldp, scomplex* vl,
                    integer ldvl, scomplex* vr, integer ldvr,
                    integer mm, integer* m)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(ctgevc,CTGEVC)(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, m, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ztgevc( char side, char howmny, const logical* select,
                    integer n, const dcomplex* s, integer lds,
                    const dcomplex* p, integer ldp, dcomplex* vl,
                    integer ldvl, dcomplex* vr, integer ldvr,
                    integer mm, integer* m)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(ztgevc,ZTGEVC)(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, m, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_stgexc( logical wantq, logical wantz, integer n,
                    float* a, integer lda, float* b, integer ldb,
                    float* q, integer ldq, float* z, integer ldz,
                    integer* ifst, integer* ilst)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(stgexc,STGEXC)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(stgexc,STGEXC)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dtgexc( logical wantq, logical wantz, integer n,
                    double* a, integer lda, double* b, integer ldb,
                    double* q, integer ldq, double* z, integer ldz,
                    integer* ifst, integer* ilst)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dtgexc,DTGEXC)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dtgexc,DTGEXC)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_ctgexc( logical wantq, logical wantz, integer n,
                    scomplex* a, integer lda, scomplex* b, integer ldb,
                    scomplex* q, integer ldq, scomplex* z, integer ldz,
                    integer* ifst, integer* ilst)
{
    integer info;
    FFUNC(ctgexc,CTGEXC)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, &info);
    return info;
}

static inline integer c_ztgexc( logical wantq, logical wantz, integer n,
                    dcomplex* a, integer lda, dcomplex* b, integer ldb,
                    dcomplex* q, integer ldq, dcomplex* z, integer ldz,
                    integer* ifst, integer* ilst)
{
    integer info;
    FFUNC(ztgexc,ZTGEXC)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, &info);
    return info;
}

static inline integer c_stgsen( integer ijob, logical wantq,
                    logical wantz, const logical* select,
                    integer n, float* a, integer lda, float* b,
                    integer ldb, scomplex* alpha,
                    float* beta, float* q, integer ldq, float* z,
                    integer ldz, integer* m, float* pl, float* pr,
                    float* dif)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    float *alphar = SAFE_MALLOC(float,n);
    float *alphai = SAFE_MALLOC(float,n);
    FFUNC(stgsen,STGSEN)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alphar, alphai, beta, q, &ldq, z, &ldz,
                            m, pl, pr, dif, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(stgsen,STGSEN)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alphar, alphai, beta, q, &ldq, z, &ldz,
                            m, pl, pr, dif, work, &lwork, iwork, &liwork, &info);
    scopy(n,alphar,1,(float*)alpha  ,2);
    scopy(n,alphai,1,(float*)alpha+1,2);
    FREE(alphar);
    FREE(alphai);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtgsen( integer ijob, logical wantq,
                    logical wantz, const logical* select,
                    integer n, double* a, integer lda, double* b,
                    integer ldb, dcomplex* alpha,
                    double* beta, double* q, integer ldq, double* z,
                    integer ldz, integer* m, double* pl, double* pr,
                    double* dif)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    double *alphar = SAFE_MALLOC(double,n);
    double *alphai = SAFE_MALLOC(double,n);
    FFUNC(dtgsen,DTGSEN)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alphar, alphai, beta, q, &ldq, z, &ldz,
                            m, pl, pr, dif, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dtgsen,DTGSEN)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alphar, alphai, beta, q, &ldq, z, &ldz,
                            m, pl, pr, dif, work, &lwork, iwork, &liwork, &info);
    dcopy(n,alphar,1,(double*)alpha  ,2);
    dcopy(n,alphai,1,(double*)alpha+1,2);
    FREE(alphar);
    FREE(alphai);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctgsen( integer ijob, logical wantq,
                    logical wantz, const logical* select,
                    integer n, scomplex* a, integer lda, scomplex* b,
                    integer ldb, scomplex* alpha,
                    scomplex* beta, scomplex* q, integer ldq, scomplex* z,
                    integer ldz, integer* m, float* pl, float* pr,
                    float* dif)
{
    integer info, lwork = -1, liwork = -1;
    scomplex *work, query;
    integer *iwork, iquery;
    FFUNC(ctgsen,CTGSEN)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alpha, beta, q, &ldq, z, &ldz,
                            m, pl, pr, dif, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)crealf(query);
    liwork = iquery;
    work = SAFE_MALLOC(scomplex,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(ctgsen,CTGSEN)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alpha, beta, q, &ldq, z, &ldz,
                            m, pl, pr, dif, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ztgsen( integer ijob, logical wantq,
                    logical wantz, const logical* select,
                    integer n, dcomplex* a, integer lda, dcomplex* b,
                    integer ldb, dcomplex* alpha,
                    dcomplex* beta, dcomplex* q, integer ldq, dcomplex* z,
                    integer ldz, integer* m, double* pl, double* pr,
                    double* dif)
{
    integer info, lwork = -1, liwork = -1;
    dcomplex *work, query;
    integer *iwork, iquery;
    FFUNC(ztgsen,ZTGSEN)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alpha, beta, q, &ldq, z, &ldz,
                            m, pl, pr, dif, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)creal(query);
    liwork = iquery;
    work = SAFE_MALLOC(dcomplex,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(ztgsen,ZTGSEN)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alpha, beta, q, &ldq, z, &ldz,
                            m, pl, pr, dif, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_stgsyl( char trans, integer ijob, integer m, integer n,
                    const float* a, integer lda, const float* b,
                    integer ldb, float* c, integer ldc,
                    const float* d, integer ldd, const float* e,
                    integer lde, float* f, integer ldf, float* scale,
                    float* dif)
{
    integer info, lwork = -1;
    float *work, query;
    integer *iwork = SAFE_MALLOC(integer,m+n+6);
    FFUNC(stgsyl,STGSYL)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, &query,
                            &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(stgsyl,STGSYL)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work,
                            &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtgsyl( char trans, integer ijob, integer m, integer n,
                    const double* a, integer lda, const double* b,
                    integer ldb, double* c, integer ldc,
                    const double* d, integer ldd, const double* e,
                    integer lde, double* f, integer ldf, double* scale,
                    double* dif)
{
    integer info, lwork = -1;
    double *work, query;
    integer *iwork = SAFE_MALLOC(integer,m+n+6);
    FFUNC(dtgsyl,DTGSYL)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, &query,
                            &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dtgsyl,DTGSYL)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work,
                            &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctgsyl( char trans, integer ijob, integer m, integer n,
                    const scomplex* a, integer lda, const scomplex* b,
                    integer ldb, scomplex* c, integer ldc,
                    const scomplex* d, integer ldd, const scomplex* e,
                    integer lde, scomplex* f, integer ldf, float* scale,
                    float* dif)
{
    integer info, lwork = -1;
    scomplex *work, query;
    integer *iwork = SAFE_MALLOC(integer,m+n+6);
    FFUNC(ctgsyl,CTGSYL)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, &query,
                            &lwork, iwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(ctgsyl,CTGSYL)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work,
                            &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ztgsyl( char trans, integer ijob, integer m, integer n,
                    const dcomplex* a, integer lda, const dcomplex* b,
                    integer ldb, dcomplex* c, integer ldc,
                    const dcomplex* d, integer ldd, const dcomplex* e,
                    integer lde, dcomplex* f, integer ldf, double* scale,
                    double* dif)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    integer *iwork = SAFE_MALLOC(integer,m+n+6);
    FFUNC(ztgsyl,ZTGSYL)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, &query,
                            &lwork, iwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(ztgsyl,ZTGSYL)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work,
                            &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_stgsna( char job, char howmny, const logical* select,
                    integer n, const float* a, integer lda,
                    const float* b, integer ldb, const float* vl,
                    integer ldvl, const float* vr, integer ldvr,
                    float* s, float* dif, integer mm, integer* m)
{
    integer info, lwork = -1;
    float *work, query;
    integer *iwork = SAFE_MALLOC(integer,n+6);
    FFUNC(stgsna,STGSNA)(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, m, &query,
                            &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(stgsna,STGSNA)(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, m, work,
                            &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dtgsna( char job, char howmny, const logical* select,
                    integer n, const double* a, integer lda,
                    const double* b, integer ldb, const double* vl,
                    integer ldvl, const double* vr, integer ldvr,
                    double* s, double* dif, integer mm, integer* m)
{
    integer info, lwork = -1;
    double *work, query;
    integer *iwork = SAFE_MALLOC(integer,n+6);
    FFUNC(dtgsna,DTGSNA)(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, m, &query,
                            &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dtgsna,DTGSNA)(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, m, work,
                            &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ctgsna( char job, char howmny, const logical* select,
                    integer n, const scomplex* a, integer lda,
                    const scomplex* b, integer ldb, const scomplex* vl,
                    integer ldvl, const scomplex* vr, integer ldvr,
                    float* s, float* dif, integer mm, integer* m)
{
    integer info, lwork = -1;
    scomplex *work, query;
    integer *iwork = SAFE_MALLOC(integer,n+6);
    FFUNC(ctgsna,CTGSNA)(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, m, &query,
                            &lwork, iwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(ctgsna,CTGSNA)(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, m, work,
                            &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_ztgsna( char job, char howmny, const logical* select,
                    integer n, const dcomplex* a, integer lda,
                    const dcomplex* b, integer ldb, const dcomplex* vl,
                    integer ldvl, const dcomplex* vr, integer ldvr,
                    double* s, double* dif, integer mm, integer* m)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    integer *iwork = SAFE_MALLOC(integer,n+6);
    FFUNC(ztgsna,ZTGSNA)(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, m, &query,
                            &lwork, iwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(ztgsna,ZTGSNA)(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, m, work,
                            &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_sggsvp( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, float* a, integer lda,
                    float* b, integer ldb, float tola, float tolb,
                    integer* k, integer* l, float* u, integer ldu,
                    float* v, integer ldv, float* q, integer ldq)
{
    integer info;
    integer *iwork = SAFE_MALLOC(integer,n);
    float *tau = SAFE_MALLOC(float,n);
    float *work = SAFE_MALLOC(float,MAX(MAX(3*n,m),p));
    FFUNC(sggsvp,SGGSVP)(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, &tola, &tolb, k, l, u, &ldu, v, &ldv, q, &ldq,
                            iwork, tau, work, &info);
    FREE(iwork);
    FREE(tau);
    FREE(work);
    return info;
}

static inline integer c_dggsvp( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, double* a, integer lda,
                    double* b, integer ldb, double tola, double tolb,
                    integer* k, integer* l, double* u, integer ldu,
                    double* v, integer ldv, double* q, integer ldq)
{
    integer info;
    integer *iwork = SAFE_MALLOC(integer,n);
    double *tau = SAFE_MALLOC(double,n);
    double *work = SAFE_MALLOC(double,MAX(MAX(3*n,m),p));
    FFUNC(dggsvp,DGGSVP)(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, &tola, &tolb, k, l, u, &ldu, v, &ldv, q, &ldq,
                            iwork, tau, work, &info);
    FREE(iwork);
    FREE(tau);
    FREE(work);
    return info;
}

static inline integer c_cggsvp( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, scomplex* a, integer lda,
                    scomplex* b, integer ldb, float tola, float tolb,
                    integer* k, integer* l, scomplex* u, integer ldu,
                    scomplex* v, integer ldv, scomplex* q, integer ldq)
{
    integer info;
    integer *iwork = SAFE_MALLOC(integer,n);
    float *rwork = SAFE_MALLOC(float,2*n);
    scomplex *tau = SAFE_MALLOC(scomplex,n);
    scomplex *work = SAFE_MALLOC(scomplex,MAX(MAX(3*n,m),p));
    FFUNC(cggsvp,CGGSVP)(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, &tola, &tolb, k, l, u, &ldu, v, &ldv, q, &ldq,
                            iwork, rwork, tau, work, &info);
    FREE(iwork);
    FREE(rwork);
    FREE(tau);
    FREE(work);
    return info;
}

static inline integer c_zggsvp( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, dcomplex* a, integer lda,
                    dcomplex* b, integer ldb, double tola, double tolb,
                    integer* k, integer* l, dcomplex* u, integer ldu,
                    dcomplex* v, integer ldv, dcomplex* q, integer ldq)
{
    integer info;
    integer *iwork = SAFE_MALLOC(integer,n);
    double *rwork = SAFE_MALLOC(double,2*n);
    dcomplex *tau = SAFE_MALLOC(dcomplex,n);
    dcomplex *work = SAFE_MALLOC(dcomplex,MAX(MAX(3*n,m),p));
    FFUNC(zggsvp,ZGGSVP)(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, &tola, &tolb, k, l, u, &ldu, v, &ldv, q, &ldq,
                            iwork, rwork, tau, work, &info);
    FREE(iwork);
    FREE(rwork);
    FREE(tau);
    FREE(work);
    return info;
}

static inline integer c_stgsja( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, integer k, integer l,
                    float* a, integer lda, float* b, integer ldb,
                    float tola, float tolb, float* alpha, float* beta,
                    float* u, integer ldu, float* v, integer ldv,
                    float* q, integer ldq, integer* ncycle )
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    FFUNC(stgsja,STGSJA)(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, &tola, &tolb, alpha, beta, u, &ldu,
                            v, &ldv, q, &ldq, work, ncycle, &info);
    FREE(work);
    return info;
}

static inline integer c_dtgsja( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, integer k, integer l,
                    double* a, integer lda, double* b, integer ldb,
                    double tola, double tolb, double* alpha, double* beta,
                    double* u, integer ldu, double* v, integer ldv,
                    double* q, integer ldq, integer* ncycle )
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    FFUNC(dtgsja,DTGSJA)(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, &tola, &tolb, alpha, beta, u, &ldu,
                            v, &ldv, q, &ldq, work, ncycle, &info);
    FREE(work);
    return info;
}

static inline integer c_ctgsja( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, integer k, integer l,
                    scomplex* a, integer lda, scomplex* b, integer ldb,
                    float tola, float tolb, float* alpha, float* beta,
                    scomplex* u, integer ldu, scomplex* v, integer ldv,
                    scomplex* q, integer ldq, integer* ncycle )
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    FFUNC(ctgsja,CTGSJA)(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, &tola, &tolb, alpha, beta, u, &ldu,
                            v, &ldv, q, &ldq, work, ncycle, &info);
    FREE(work);
    return info;
}

static inline integer c_ztgsja( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, integer k, integer l,
                    dcomplex* a, integer lda, dcomplex* b, integer ldb,
                    double tola, double tolb, double* alpha, double* beta,
                    dcomplex* u, integer ldu, dcomplex* v, integer ldv,
                    dcomplex* q, integer ldq, integer* ncycle )
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    FFUNC(ztgsja,ZTGSJA)(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, &tola, &tolb, alpha, beta, u, &ldu,
                            v, &ldv, q, &ldq, work, ncycle, &info);
    FREE(work);
    return info;
}

static inline integer c_sgels( char trans, integer m, integer n, integer nrhs,
                   float* a, integer lda, float* b, integer ldb)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgels,SGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgels,SGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgels( char trans, integer m, integer n, integer nrhs,
                   double* a, integer lda, double* b, integer ldb)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgels,DGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgels,DGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgels( char trans, integer m, integer n, integer nrhs,
                   scomplex* a, integer lda, scomplex* b, integer ldb)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cgels,CGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgels,CGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgels( char trans, integer m, integer n, integer nrhs,
                   dcomplex* a, integer lda, dcomplex* b, integer ldb)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zgels,ZGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgels,ZGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sgelsy( integer m, integer n, integer nrhs, float* a,
                    integer lda, float* b, integer ldb,
                    integer* jpvt, float rcond, integer* rank)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgelsy,SGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgelsy,SGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgelsy( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb,
                    integer* jpvt, double rcond, integer* rank)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgelsy,DGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgelsy,DGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgelsy( integer m, integer n, integer nrhs, scomplex* a,
                    integer lda, scomplex* b, integer ldb,
                    integer* jpvt, float rcond, integer* rank)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgelsy,CGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgelsy,CGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgelsy( integer m, integer n, integer nrhs, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb,
                    integer* jpvt, double rcond, integer* rank)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgelsy,ZGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgelsy,ZGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgelss( integer m, integer n, integer nrhs, float* a,
                    integer lda, float* b, integer ldb, float* s,
                    float rcond, integer* rank)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgelss,SGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgelss,SGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgelss( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* s,
                    double rcond, integer* rank)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgelss,DGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgelss,DGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgelss( integer m, integer n, integer nrhs, scomplex* a,
                    integer lda, scomplex* b, integer ldb, float* s,
                    float rcond, integer* rank)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgelss,CGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgelss,CGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgelss( integer m, integer n, integer nrhs, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb, double* s,
                    double rcond, integer* rank)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgelss,ZGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgelss,ZGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgglse( integer m, integer n, integer p, float* a,
                    integer lda, float* b, integer ldb, float* c,
                    float* d, float* x)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgglse,SGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgglse,SGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgglse( integer m, integer n, integer p, double* a,
                    integer lda, double* b, integer ldb, double* c,
                    double* d, double* x)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgglse,DGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgglse,DGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgglse( integer m, integer n, integer p, scomplex* a,
                    integer lda, scomplex* b, integer ldb, scomplex* c,
                    scomplex* d, scomplex* x)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cgglse,CGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgglse,CGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgglse( integer m, integer n, integer p, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb, dcomplex* c,
                    dcomplex* d, dcomplex* x)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zgglse,ZGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgglse,ZGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sggglm( integer n, integer m, integer p, float* a,
                    integer lda, float* b, integer ldb, float* d,
                    float* x, float* y)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sggglm,SGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sggglm,SGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dggglm( integer n, integer m, integer p, double* a,
                    integer lda, double* b, integer ldb, double* d,
                    double* x, double* y)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dggglm,DGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dggglm,DGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cggglm( integer n, integer m, integer p, scomplex* a,
                    integer lda, scomplex* b, integer ldb, scomplex* d,
                    scomplex* x, scomplex* y)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cggglm,CGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cggglm,CGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zggglm( integer n, integer m, integer p, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb, dcomplex* d,
                    dcomplex* x, dcomplex* y)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zggglm,ZGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zggglm,ZGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_ssyev( char jobz, char uplo, integer n, float* a,
                   integer lda, float* w)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(ssyev,SSYEV)(&jobz, &uplo, &n, a, &lda, w, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(ssyev,SSYEV)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dsyev( char jobz, char uplo, integer n, double* a,
                   integer lda, double* w)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dsyev,DSYEV)(&jobz, &uplo, &n, a, &lda, w, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dsyev,DSYEV)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cheev( char jobz, char uplo, integer n, scomplex* a,
                   integer lda, float* w)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,3*n);
    FFUNC(cheev,CHEEV)(&jobz, &uplo, &n, a, &lda, w, &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cheev,CHEEV)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zheev( char jobz, char uplo, integer n, dcomplex* a,
                   integer lda, double* w)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,3*n);
    FFUNC(zheev,ZHEEV)(&jobz, &uplo, &n, a, &lda, w, &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zheev,ZHEEV)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ssyevd( char jobz, char uplo, integer n, float* a,
                    integer lda, float* w)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(ssyevd,SSYEVD)(&jobz, &uplo, &n, a, &lda, w, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(ssyevd,SSYEVD)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsyevd( char jobz, char uplo, integer n, double* a,
                    integer lda, double* w)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dsyevd,DSYEVD)(&jobz, &uplo, &n, a, &lda, w, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dsyevd,DSYEVD)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cheevd( char jobz, char uplo, integer n, scomplex* a,
                    integer lda, float* w)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    scomplex *work, query;
    float *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(cheevd,CHEEVD)(&jobz, &uplo, &n, a, &lda, w, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)crealf(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(scomplex,lwork);
    rwork = SAFE_MALLOC(float,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(cheevd,CHEEVD)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zheevd( char jobz, char uplo, integer n, dcomplex* a,
                    integer lda, double* w)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    dcomplex *work, query;
    double *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(zheevd,ZHEEVD)(&jobz, &uplo, &n, a, &lda, w, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)creal(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(dcomplex,lwork);
    rwork = SAFE_MALLOC(double,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zheevd,ZHEEVD)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_ssyevx( char jobz, char range, char uplo, integer n,
                    float* a, integer lda, float vl, float vu,
                    integer il, integer iu, float abstol,
                    integer* m, float* w, float* z, integer ldz,
                    integer* ifail )
{
    integer info, lwork = -1;
    float *work, query;
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(ssyevx,SSYEVX)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, &query, &lwork,
                            iwork, ifail, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(ssyevx,SSYEVX)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork,
                            iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsyevx( char jobz, char range, char uplo, integer n,
                    double* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer* m, double* w, double* z, integer ldz,
                    integer* ifail )
{
    integer info, lwork = -1;
    double *work, query;
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(dsyevx,DSYEVX)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, &query, &lwork,
                            iwork, ifail, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dsyevx,DSYEVX)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork,
                            iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cheevx( char jobz, char range, char uplo, integer n,
                    scomplex* a, integer lda, float vl, float vu,
                    integer il, integer iu, float abstol,
                    integer* m, float* w, scomplex* z, integer ldz,
                    integer* ifail )
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(cheevx,CHEEVX)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, &query, &lwork,
                            rwork, iwork, ifail, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cheevx,CHEEVX)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork,
                            rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zheevx( char jobz, char range, char uplo, integer n,
                    dcomplex* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer* m, double* w, dcomplex* z, integer ldz,
                    integer* ifail )
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(zheevx,ZHEEVX)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, &query, &lwork,
                            rwork, iwork, ifail, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zheevx,ZHEEVX)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork,
                            rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_ssyevr( char jobz, char range, char uplo, integer n,
                    float* a, integer lda, float vl, float vu,
                    integer il, integer iu, float abstol,
                    integer* m, float* w, float* z, integer ldz,
                    integer* isuppz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(ssyevr,SSYEVR)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(ssyevr,SSYEVR)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsyevr( char jobz, char range, char uplo, integer n,
                    double* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer* m, double* w, double* z, integer ldz,
                    integer* isuppz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dsyevr,DSYEVR)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dsyevr,DSYEVR)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cheevr( char jobz, char range, char uplo, integer n,
                    scomplex* a, integer lda, float vl, float vu,
                    integer il, integer iu, float abstol,
                    integer* m, float* w, scomplex* z, integer ldz,
                    integer* isuppz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    scomplex *work, query;
    float *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(cheevr,CHEEVR)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)crealf(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(scomplex,lwork);
    rwork = SAFE_MALLOC(float,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(cheevr,CHEEVR)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zheevr( char jobz, char range, char uplo, integer n,
                    dcomplex* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer* m, double* w, dcomplex* z, integer ldz,
                    integer* isuppz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    dcomplex *work, query;
    double *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(zheevr,ZHEEVR)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)creal(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(dcomplex,lwork);
    rwork = SAFE_MALLOC(double,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zheevr,ZHEEVR)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_sspev( char jobz, char uplo, integer n, float* ap, float* w,
                   float* z, integer ldz)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    FFUNC(sspev,SSPEV)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dspev( char jobz, char uplo, integer n, double* ap, double* w,
                   double* z, integer ldz)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    FFUNC(dspev,DSPEV)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_chpev( char jobz, char uplo, integer n, scomplex* ap, float* w,
                   scomplex* z, integer ldz)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,3*n);
    FFUNC(chpev,CHPEV)(&jobz, &uplo, &n, ap, w, z, &ldz, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhpev( char jobz, char uplo, integer n, dcomplex* ap, double* w,
                   dcomplex* z, integer ldz)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,3*n);
    FFUNC(zhpev,ZHPEV)(&jobz, &uplo, &n, ap, w, z, &ldz, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sspevd( char jobz, char uplo, integer n, float* ap,
                    float* w, float* z, integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(sspevd,SSPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(sspevd,SSPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dspevd( char jobz, char uplo, integer n, double* ap,
                    double* w, double* z, integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dspevd,DSPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dspevd,DSPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chpevd( char jobz, char uplo, integer n, scomplex* ap,
                    float* w, scomplex* z, integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    scomplex *work, query;
    float *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(chpevd,CHPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)crealf(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(scomplex,lwork);
    rwork = SAFE_MALLOC(float,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(chpevd,CHPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhpevd( char jobz, char uplo, integer n, dcomplex* ap,
                    double* w, dcomplex* z, integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    dcomplex *work, query;
    double *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(zhpevd,ZHPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)creal(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(dcomplex,lwork);
    rwork = SAFE_MALLOC(double,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zhpevd,ZHPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_sspevx( char jobz, char range, char uplo, integer n,
                    float* ap, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    integer info;
    float *work = SAFE_MALLOC(float,8*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(sspevx,SSPEVX)(&jobz, &range, &uplo, &n, ap, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dspevx( char jobz, char range, char uplo, integer n,
                    double* ap, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    integer info;
    double *work = SAFE_MALLOC(double,8*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(dspevx,DSPEVX)(&jobz, &range, &uplo, &n, ap, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chpevx( char jobz, char range, char uplo, integer n,
                    scomplex* ap, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    scomplex* z, integer ldz, integer* ifail )
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(chpevx,CHPEVX)(&jobz, &range, &uplo, &n, ap, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhpevx( char jobz, char range, char uplo, integer n,
                    dcomplex* ap, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    dcomplex* z, integer ldz, integer* ifail )
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(zhpevx,ZHPEVX)(&jobz, &range, &uplo, &n, ap, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_ssbev( char jobz, char uplo, integer n, integer kd,
                   float* ab, integer ldab, float* w, float* z,
                   integer ldz)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    FFUNC(ssbev,SSBEV)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dsbev( char jobz, char uplo, integer n, integer kd,
                   double* ab, integer ldab, double* w, double* z,
                   integer ldz)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    FFUNC(dsbev,DSBEV)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_chbev( char jobz, char uplo, integer n, integer kd,
                   scomplex* ab, integer ldab, float* w, scomplex* z,
                   integer ldz)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    float *rwork = SAFE_MALLOC(float,3*n);
    FFUNC(chbev,CHBEV)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhbev( char jobz, char uplo, integer n, integer kd,
                   dcomplex* ab, integer ldab, double* w, dcomplex* z,
                   integer ldz)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    double *rwork = SAFE_MALLOC(double,3*n);
    FFUNC(zhbev,ZHBEV)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ssbevd( char jobz, char uplo, integer n, integer kd,
                    float* ab, integer ldab, float* w, float* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(ssbevd,SSBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(ssbevd,SSBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsbevd( char jobz, char uplo, integer n, integer kd,
                    double* ab, integer ldab, double* w, double* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dsbevd,DSBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dsbevd,DSBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chbevd( char jobz, char uplo, integer n, integer kd,
                    scomplex* ab, integer ldab, float* w, scomplex* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    scomplex *work, query;
    float *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(chbevd,CHBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)crealf(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(scomplex,lwork);
    rwork = SAFE_MALLOC(float,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(chbevd,CHBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhbevd( char jobz, char uplo, integer n, integer kd,
                    dcomplex* ab, integer ldab, double* w, dcomplex* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    dcomplex *work, query;
    double *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(zhbevd,ZHBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)creal(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(dcomplex,lwork);
    rwork = SAFE_MALLOC(double,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zhbevd,ZHBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_ssbevx( char jobz, char range, char uplo, integer n,
                    integer kd, float* ab, integer ldab, float* q,
                    integer ldq, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    integer info;
    float *work = SAFE_MALLOC(float,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(ssbevx,SSBEVX)(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                           work, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsbevx( char jobz, char range, char uplo, integer n,
                    integer kd, double* ab, integer ldab, double* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    integer info;
    double *work = SAFE_MALLOC(double,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(dsbevx,DSBEVX)(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                           work, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chbevx( char jobz, char range, char uplo, integer n,
                    integer kd, scomplex* ab, integer ldab, scomplex* q,
                    integer ldq, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    scomplex* z, integer ldz, integer* ifail )
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    float *rwork = SAFE_MALLOC(float,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(chbevx,CHBEVX)(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                           work, rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhbevx( char jobz, char range, char uplo, integer n,
                    integer kd, dcomplex* ab, integer ldab, dcomplex* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    dcomplex* z, integer ldz, integer* ifail )
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    double *rwork = SAFE_MALLOC(double,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(zhbevx,ZHBEVX)(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                           work, rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_sstev( char jobz, integer n, float* d, float* e, float* z,
                   integer ldz)
{
    integer info;
    float *work = SAFE_MALLOC(float,2*n);
    FFUNC(sstev,SSTEV)(&jobz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dstev( char jobz, integer n, double* d, double* e, double* z,
                   integer ldz)
{
    integer info;
    double *work = SAFE_MALLOC(double,2*n);
    FFUNC(dstev,DSTEV)(&jobz, &n, d, e, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_sstevd( char jobz, integer n, float* d, float* e, float* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(sstevd,SSTEVD)(&jobz, &n, d, e, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(sstevd,SSTEVD)(&jobz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dstevd( char jobz, integer n, double* d, double* e, double* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dstevd,DSTEVD)(&jobz, &n, d, e, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dstevd,DSTEVD)(&jobz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_sstevx( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    integer info;
    float *work = SAFE_MALLOC(float,5*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(sstevx,SSTEVX)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dstevx( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    integer info;
    double *work = SAFE_MALLOC(double,5*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(dstevx,DSTEVX)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_sstevr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    float* z, integer ldz, integer* isuppz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(sstevr,SSTEVR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(sstevr,SSTEVR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dstevr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, integer* isuppz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dstevr,DSTEVR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, &query, &lwork,
                            &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dstevr,DSTEVR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_sgeev( char jobvl, char jobvr, integer n, float* a,
                   integer lda, scomplex* w, float* vl,
                   integer ldvl, float* vr, integer ldvr)
{
    integer info, lwork = -1;
    float *work, query;
    float *wr = SAFE_MALLOC(float,n);
    float *wi = SAFE_MALLOC(float,n);
    FFUNC(sgeev,SGEEV)(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgeev,SGEEV)(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    scopy(n,wr,1,(float*)w  ,2);
    scopy(n,wi,1,(float*)w+1,2);
    FREE(wr);
    FREE(wi);
    FREE(work);
    return info;
}

static inline integer c_dgeev( char jobvl, char jobvr, integer n, double* a,
                   integer lda, dcomplex* w, double* vl,
                   integer ldvl, double* vr, integer ldvr)
{
    integer info, lwork = -1;
    double *work, query;
    double *wr = SAFE_MALLOC(double,n);
    double *wi = SAFE_MALLOC(double,n);
    FFUNC(dgeev,DGEEV)(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgeev,DGEEV)(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    dcopy(n,wr,1,(double*)w  ,2);
    dcopy(n,wi,1,(double*)w+1,2);
    FREE(wr);
    FREE(wi);
    FREE(work);
    return info;
}

static inline integer c_cgeev( char jobvl, char jobvr, integer n, scomplex* a,
                   integer lda, scomplex* w, scomplex* vl,
                   integer ldvl, scomplex* vr, integer ldvr)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgeev,CGEEV)(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgeev,CGEEV)(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgeev( char jobvl, char jobvr, integer n, dcomplex* a,
                   integer lda, dcomplex* w, dcomplex* vl,
                   integer ldvl, dcomplex* vr, integer ldvr)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgeev,ZGEEV)(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgeev,ZGEEV)(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgeevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, float* a, integer lda, scomplex* w,
                    float* vl, integer ldvl, float* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    float* scale, float* abnrm, float* rconde,
                    float* rcondv)
{
    integer info, lwork = -1;
    float *work, query;
    float *wr = SAFE_MALLOC(float,n);
    float *wi = SAFE_MALLOC(float,n);
    integer *iwork = SAFE_MALLOC(integer,2*n);
    FFUNC(sgeevx,SGEEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, ilo, ihi, scale, abnrm,
                            rconde, rcondv, &query, &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgeevx,SGEEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, ilo, ihi, scale, abnrm,
                            rconde, rcondv, work, &lwork, iwork, &info);
    scopy(n,wr,1,(float*)w  ,2);
    scopy(n,wi,1,(float*)w+1,2);
    FREE(wr);
    FREE(wi);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgeevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, double* a, integer lda, dcomplex* w,
                    double* vl, integer ldvl, double* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    double* scale, double* abnrm, double* rconde,
                    double* rcondv)
{
    integer info, lwork = -1;
    double *work, query;
    double *wr = SAFE_MALLOC(double,n);
    double *wi = SAFE_MALLOC(double,n);
    integer *iwork = SAFE_MALLOC(integer,2*n);
    FFUNC(dgeevx,DGEEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, ilo, ihi, scale, abnrm,
                            rconde, rcondv, &query, &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgeevx,DGEEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, ilo, ihi, scale, abnrm,
                            rconde, rcondv, work, &lwork, iwork, &info);
    dcopy(n,wr,1,(double*)w  ,2);
    dcopy(n,wi,1,(double*)w+1,2);
    FREE(wr);
    FREE(wi);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgeevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, scomplex* a, integer lda, scomplex* w,
                    scomplex* vl, integer ldvl, scomplex* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    float* scale, float* abnrm, float* rconde,
                    float* rcondv)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,2*n);
    FFUNC(cgeevx,CGEEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, ilo, ihi, scale, abnrm,
                            rconde, rcondv, &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgeevx,CGEEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, ilo, ihi, scale, abnrm,
                            rconde, rcondv, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgeevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, dcomplex* a, integer lda, dcomplex* w,
                    dcomplex* vl, integer ldvl, dcomplex* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    double* scale, double* abnrm, double* rconde,
                    double* rcondv)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,2*n);
    FFUNC(zgeevx,ZGEEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, ilo, ihi, scale, abnrm,
                            rconde, rcondv, &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgeevx,ZGEEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, ilo, ihi, scale, abnrm,
                            rconde, rcondv, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgesvd( char jobu, char jobvt, integer m, integer n,
                    float* a, integer lda, float* s, float* u,
                    integer ldu, float* vt, integer ldvt)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgesvd,SGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgesvd,SGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgesvd( char jobu, char jobvt, integer m, integer n,
                    double* a, integer lda, double* s, double* u,
                    integer ldu, double* vt, integer ldvt)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgesvd,DGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgesvd,DGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgesvd( char jobu, char jobvt, integer m, integer n,
                    scomplex* a, integer lda, float* s, scomplex* u,
                    integer ldu, scomplex* vt, integer ldvt)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,5*MIN(m,n));
    FFUNC(cgesvd,CGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgesvd,CGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zgesvd( char jobu, char jobvt, integer m, integer n,
                    dcomplex* a, integer lda, double* s, dcomplex* u,
                    integer ldu, dcomplex* vt, integer ldvt)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,5*MIN(m,n));
    FFUNC(zgesvd,ZGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgesvd,ZGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sgesdd( char jobz, integer m, integer n, float* a,
                    integer lda, float* s, float* u, integer ldu,
                    float* vt, integer ldvt)
{
    integer info, lwork = -1;
    float *work, query;
    integer *iwork = SAFE_MALLOC(integer,8*MIN(m,n));
    FFUNC(sgesdd,SGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &query, &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgesdd,SGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgesdd( char jobz, integer m, integer n, double* a,
                    integer lda, double* s, double* u, integer ldu,
                    double* vt, integer ldvt)
{
    integer info, lwork = -1;
    double *work, query;
    integer *iwork = SAFE_MALLOC(integer,8*MIN(m,n));
    FFUNC(dgesdd,DGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &query, &lwork, iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgesdd,DGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cgesdd( char jobz, integer m, integer n, scomplex* a,
                    integer lda, float* s, scomplex* u, integer ldu,
                    scomplex* vt, integer ldvt)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,5*MIN(m,n)*MIN(m,n) + 7*MIN(m,n));
    integer *iwork = SAFE_MALLOC(integer,8*MIN(m,n));
    FFUNC(cgesdd,CGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &query, &lwork, rwork, iwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgesdd,CGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zgesdd( char jobz, integer m, integer n, dcomplex* a,
                    integer lda, double* s, dcomplex* u, integer ldu,
                    dcomplex* vt, integer ldvt)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,5*MIN(m,n)*MIN(m,n) + 7*MIN(m,n));
    integer *iwork = SAFE_MALLOC(integer,8*MIN(m,n));
    FFUNC(zgesdd,ZGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &query, &lwork, rwork, iwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgesdd,ZGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_sgejsv( char joba, char jobu, char jobv, char jobr, char jobt,
                    char jobp, integer m, integer n, float* a,
                    integer lda, float* sva, float* u, integer ldu,
                    float* v, integer ldv)
{
    integer info, lwork = -1;
    float *work, query;
    integer *iwork = SAFE_MALLOC(integer,m+3*n);
    FFUNC(sgejsv,SGEJSV)(&joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, a, &lda, sva, u, &ldu, v, &ldv, &query, &lwork,
                            iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgejsv,SGEJSV)(&joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, a, &lda, sva, u, &ldu, v, &ldv, work, &lwork,
                            iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dgejsv( char joba, char jobu, char jobv, char jobr, char jobt,
                    char jobp, integer m, integer n, double* a,
                    integer lda, double* sva, double* u, integer ldu,
                    double* v, integer ldv)
{
    integer info, lwork = -1;
    double *work, query;
    integer *iwork = SAFE_MALLOC(integer,m+3*n);
    FFUNC(dgejsv,DGEJSV)(&joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, a, &lda, sva, u, &ldu, v, &ldv, &query, &lwork,
                            iwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgejsv,DGEJSV)(&joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, a, &lda, sva, u, &ldu, v, &ldv, work, &lwork,
                            iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_sgesvj( char joba, char jobu, char jobv, integer m,
                    integer n, float* a, integer lda, float* sva,
                    integer mv, float* v, integer ldv)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgesvj,SGESVJ)(&joba, &jobu, &jobv, &m, &n, a, &lda, sva, &mv, v, &ldv, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgesvj,SGESVJ)(&joba, &jobu, &jobv, &m, &n, a, &lda, sva, &mv, v, &ldv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgesvj( char joba, char jobu, char jobv, integer m,
                    integer n, double* a, integer lda, double* sva,
                    integer mv, double* v, integer ldv)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgesvj,DGESVJ)(&joba, &jobu, &jobv, &m, &n, a, &lda, sva, &mv, v, &ldv, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgesvj,DGESVJ)(&joba, &jobu, &jobv, &m, &n, a, &lda, sva, &mv, v, &ldv, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_sggsvd( char jobu, char jobv, char jobq, integer m,
                    integer n, integer p, integer* k, integer* l,
                    float* a, integer lda, float* b, integer ldb,
                    float* alpha, float* beta, float* u, integer ldu,
                    float* v, integer ldv, float* q, integer ldq)
{
    integer info;
    float *work = SAFE_MALLOC(float,MAX(MAX(3*n,m),p)+n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(sggsvd,SGGSVD)(&jobu, &jobv, &jobq, &m, &n, &p, k, l, a, &lda, b, &ldb, alpha, beta, u, &ldu, v, &ldv, q, &ldq,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dggsvd( char jobu, char jobv, char jobq, integer m,
                    integer n, integer p, integer* k, integer* l,
                    double* a, integer lda, double* b, integer ldb,
                    double* alpha, double* beta, double* u, integer ldu,
                    double* v, integer ldv, double* q, integer ldq)
{
    integer info;
    double *work = SAFE_MALLOC(double,MAX(MAX(3*n,m),p)+n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(dggsvd,DGGSVD)(&jobu, &jobv, &jobq, &m, &n, &p, k, l, a, &lda, b, &ldb, alpha, beta, u, &ldu, v, &ldv, q, &ldq,
                            work, iwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_cggsvd( char jobu, char jobv, char jobq, integer m,
                    integer n, integer p, integer* k, integer* l,
                    scomplex* a, integer lda, scomplex* b, integer ldb,
                    float* alpha, float* beta, scomplex* u, integer ldu,
                    scomplex* v, integer ldv, scomplex* q, integer ldq)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,MAX(MAX(3*n,m),p)+n);
    float *rwork = SAFE_MALLOC(float,2*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(cggsvd,CGGSVD)(&jobu, &jobv, &jobq, &m, &n, &p, k, l, a, &lda, b, &ldb, alpha, beta, u, &ldu, v, &ldv, q, &ldq,
                            work, rwork, iwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zggsvd( char jobu, char jobv, char jobq, integer m,
                    integer n, integer p, integer* k, integer* l,
                    dcomplex* a, integer lda, dcomplex* b, integer ldb,
                    double* alpha, double* beta, dcomplex* u, integer ldu,
                    dcomplex* v, integer ldv, dcomplex* q, integer ldq)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,MAX(MAX(3*n,m),p)+n);
    double *rwork = SAFE_MALLOC(double,2*n);
    integer *iwork = SAFE_MALLOC(integer,n);
    FFUNC(zggsvd,ZGGSVD)(&jobu, &jobv, &jobq, &m, &n, &p, k, l, a, &lda, b, &ldb, alpha, beta, u, &ldu, v, &ldv, q, &ldq,
                            work, rwork, iwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_ssygv( integer itype, char jobz, char uplo, integer n,
                   float* a, integer lda, float* b, integer ldb,
                   float* w)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(ssygv,SSYGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(ssygv,SSYGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dsygv( integer itype, char jobz, char uplo, integer n,
                   double* a, integer lda, double* b, integer ldb,
                   double* w)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dsygv,DSYGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dsygv,DSYGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_chegv( integer itype, char jobz, char uplo, integer n,
                   scomplex* a, integer lda, scomplex* b, integer ldb,
                   float* w)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,3*n);
    FFUNC(chegv,CHEGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(chegv,CHEGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhegv( integer itype, char jobz, char uplo, integer n,
                   dcomplex* a, integer lda, dcomplex* b, integer ldb,
                   double* w)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,3*n);
    FFUNC(zhegv,ZHEGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zhegv,ZHEGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ssygvd( integer itype, char jobz, char uplo, integer n,
                    float* a, integer lda, float* b, integer ldb,
                    float* w)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(ssygvd,SSYGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w,
                        &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(ssygvd,SSYGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w,
                        work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsygvd( integer itype, char jobz, char uplo, integer n,
                    double* a, integer lda, double* b, integer ldb,
                    double* w)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dsygvd,DSYGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w,
                        &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dsygvd,DSYGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w,
                        work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chegvd( integer itype, char jobz, char uplo, integer n,
                    scomplex* a, integer lda, scomplex* b, integer ldb,
                    float* w)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    scomplex *work, query;
    float *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(chegvd,CHEGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w,
                        &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)crealf(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(scomplex,lwork);
    rwork = SAFE_MALLOC(float,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(chegvd,CHEGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w,
                        work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhegvd( integer itype, char jobz, char uplo, integer n,
                    dcomplex* a, integer lda, dcomplex* b, integer ldb,
                    double* w)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    dcomplex *work, query;
    double *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(zhegvd,ZHEGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w,
                        &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)creal(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(dcomplex,lwork);
    rwork = SAFE_MALLOC(double,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zhegvd,ZHEGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w,
                        work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_ssygvx( integer itype, char jobz, char range, char uplo,
                    integer n, float* a, integer lda, float* b,
                    integer ldb, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    integer info, lwork = -1;
    float *work, query;
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(ssygvx,SSYGVX)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                            &query, &lwork, iwork, ifail, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(ssygvx,SSYGVX)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                            work, &lwork, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsygvx( integer itype, char jobz, char range, char uplo,
                    integer n, double* a, integer lda, double* b,
                    integer ldb, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    integer info, lwork = -1;
    double *work, query;
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(dsygvx,DSYGVX)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                            &query, &lwork, iwork, ifail, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dsygvx,DSYGVX)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                            work, &lwork, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chegvx( integer itype, char jobz, char range, char uplo,
                    integer n, scomplex* a, integer lda, scomplex* b,
                    integer ldb, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    scomplex* z, integer ldz, integer* ifail )
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(chegvx,CHEGVX)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                            &query, &lwork, rwork, iwork, ifail, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(chegvx,CHEGVX)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                            work, &lwork, rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhegvx( integer itype, char jobz, char range, char uplo,
                    integer n, dcomplex* a, integer lda, dcomplex* b,
                    integer ldb, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    dcomplex* z, integer ldz, integer* ifail )
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(zhegvx,ZHEGVX)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                            &query, &lwork, rwork, iwork, ifail, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zhegvx,ZHEGVX)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                            work, &lwork, rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_sspgv( integer itype, char jobz, char uplo, integer n,
                    float* ap, float* bp, float* w, float* z,
                   integer ldz)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    FFUNC(sspgv,SSPGV)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dspgv( integer itype, char jobz, char uplo, integer n,
                   double* ap, double* bp, double* w, double* z,
                   integer ldz)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    FFUNC(dspgv,DSPGV)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_chpgv( integer itype, char jobz, char uplo, integer n,
                   scomplex* ap, scomplex* bp, float* w, scomplex* z,
                   integer ldz)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,3*n);
    float *rwork = SAFE_MALLOC(float,3*n);
    FFUNC(chpgv,CHPGV)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhpgv( integer itype, char jobz, char uplo, integer n,
                   dcomplex* ap, dcomplex* bp, double* w, dcomplex* z,
                   integer ldz)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,3*n);
    double *rwork = SAFE_MALLOC(double,3*n);
    FFUNC(zhpgv,ZHPGV)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sspgvd( integer itype, char jobz, char uplo, integer n,
                    float* ap, float* bp, float* w, float* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(sspgvd,SSPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(sspgvd,SSPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dspgvd( integer itype, char jobz, char uplo, integer n,
                    double* ap, double* bp, double* w, double* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dspgvd,DSPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dspgvd,DSPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chpgvd( integer itype, char jobz, char uplo, integer n,
                    scomplex* ap, scomplex* bp, float* w, scomplex* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    scomplex *work, query;
    float *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(chpgvd,CHPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)crealf(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(scomplex,lwork);
    rwork = SAFE_MALLOC(float,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(chpgvd,CHPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhpgvd( integer itype, char jobz, char uplo, integer n,
                    dcomplex* ap, dcomplex* bp, double* w, dcomplex* z,
                    integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    dcomplex *work, query;
    double *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(zhpgvd,ZHPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)creal(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(dcomplex,lwork);
    rwork = SAFE_MALLOC(double,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zhpgvd,ZHPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_sspgvx( integer itype, char jobz, char range, char uplo,
                    integer n, float* ap, float* bp, float vl,
                    float vu, integer il, integer iu, float abstol,
                    integer* m, float* w, float* z, integer ldz,
                    integer* ifail )
{
    integer info;
    float *work = SAFE_MALLOC(float,8*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(sspgvx,SSPGVX)(&itype, &jobz, &range, &uplo, &n, ap, bp, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork,
                            ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dspgvx( integer itype, char jobz, char range, char uplo,
                    integer n, double* ap, double* bp, double vl,
                    double vu, integer il, integer iu, double abstol,
                    integer* m, double* w, double* z, integer ldz,
                    integer* ifail )
{
    integer info;
    double *work = SAFE_MALLOC(double,8*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(dspgvx,DSPGVX)(&itype, &jobz, &range, &uplo, &n, ap, bp, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork,
                            ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chpgvx( integer itype, char jobz, char range, char uplo,
                    integer n, scomplex* ap, scomplex* bp, float vl,
                    float vu, integer il, integer iu, float abstol,
                    integer* m, float* w, scomplex* z, integer ldz,
                    integer* ifail )
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,2*n);
    float *rwork = SAFE_MALLOC(float,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(chpgvx,CHPGVX)(&itype, &jobz, &range, &uplo, &n, ap, bp, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, rwork, iwork,
                            ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhpgvx( integer itype, char jobz, char range, char uplo,
                    integer n, dcomplex* ap, dcomplex* bp, double vl,
                    double vu, integer il, integer iu, double abstol,
                    integer* m, double* w, dcomplex* z, integer ldz,
                    integer* ifail )
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,2*n);
    double *rwork = SAFE_MALLOC(double,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(zhpgvx,ZHPGVX)(&itype, &jobz, &range, &uplo, &n, ap, bp, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, rwork, iwork,
                            ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_ssbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, float* ab, integer ldab, float* bb,
                   integer ldbb, float* w, float* z, integer ldz)
{
    integer info;
    float *work = SAFE_MALLOC(float,3*n);
    FFUNC(ssbgv,SSBGV)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_dsbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, double* ab, integer ldab, double* bb,
                   integer ldbb, double* w, double* z, integer ldz)
{
    integer info;
    double *work = SAFE_MALLOC(double,3*n);
    FFUNC(dsbgv,DSBGV)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &info);
    FREE(work);
    return info;
}

static inline integer c_chbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, scomplex* ab, integer ldab, scomplex* bb,
                   integer ldbb, float* w, scomplex* z, integer ldz)
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    float *rwork = SAFE_MALLOC(float,3*n);
    FFUNC(chbgv,CHBGV)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zhbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, dcomplex* ab, integer ldab, dcomplex* bb,
                   integer ldbb, double* w, dcomplex* z, integer ldz)
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    double *rwork = SAFE_MALLOC(double,3*n);
    FFUNC(zhbgv,ZHBGV)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_ssbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, float* ab, integer ldab, float* bb,
                    integer ldbb, float* w, float* z, integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    float *work, query;
    integer *iwork, iquery;
    FFUNC(ssbgvd,SSBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(float,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(ssbgvd,SSBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, double* ab, integer ldab, double* bb,
                    integer ldbb, double* w, double* z, integer ldz)
{
    integer info, lwork = -1, liwork = -1;
    double *work, query;
    integer *iwork, iquery;
    FFUNC(dsbgvd,DSBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, &query, &lwork, &iquery, &liwork, &info);
    lwork = (integer)query;
    liwork = iquery;
    work = SAFE_MALLOC(double,lwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(dsbgvd,DSBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, scomplex* ab, integer ldab, scomplex* bb,
                    integer ldbb, float* w, scomplex* z, integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    scomplex *work, query;
    float *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(chbgvd,CHBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)crealf(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(scomplex,lwork);
    rwork = SAFE_MALLOC(float,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(chbgvd,CHBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, dcomplex* ab, integer ldab, dcomplex* bb,
                    integer ldbb, double* w, dcomplex* z, integer ldz)
{
    integer info, lwork = -1, liwork = -1, lrwork = -1;
    dcomplex *work, query;
    double *rwork, rquery;
    integer *iwork, iquery;
    FFUNC(zhbgvd,ZHBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, &query, &lwork, &rquery, &lrwork, &iquery, &liwork, &info);
    lwork = (integer)creal(query);
    lrwork = (integer)rquery;
    liwork = iquery;
    work = SAFE_MALLOC(dcomplex,lwork);
    rwork = SAFE_MALLOC(double,lrwork);
    iwork = SAFE_MALLOC(integer,liwork);
    FFUNC(zhbgvd,ZHBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_ssbgvx( char jobz, char range, char uplo, integer n,
                    integer ka, integer kb, float* ab,
                    integer ldab, float* bb, integer ldbb, float* q,
                    integer ldq, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    integer info;
    float *work = SAFE_MALLOC(float,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(ssbgvx,SSBGVX)(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, &vl, &vu, &il, &iu, &abstol,
                            m, w, z, &ldz, work, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_dsbgvx( char jobz, char range, char uplo, integer n,
                    integer ka, integer kb, double* ab,
                    integer ldab, double* bb, integer ldbb, double* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    integer info;
    double *work = SAFE_MALLOC(double,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(dsbgvx,DSBGVX)(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, &vl, &vu, &il, &iu, &abstol,
                            m, w, z, &ldz, work, iwork, ifail, &info);
    FREE(work);
    FREE(iwork);
    return info;
}

static inline integer c_chbgvx( char jobz, char range, char uplo, integer n,
                    integer ka, integer kb, scomplex* ab,
                    integer ldab, scomplex* bb, integer ldbb, scomplex* q,
                    integer ldq, float vl, float vu, integer il,
                    integer iu, float abstol, integer* m, float* w,
                    scomplex* z, integer ldz, integer* ifail )
{
    integer info;
    scomplex *work = SAFE_MALLOC(scomplex,n);
    float *rwork = SAFE_MALLOC(float,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(chbgvx,CHBGVX)(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, &vl, &vu, &il, &iu, &abstol,
                            m, w, z, &ldz, work, rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_zhbgvx( char jobz, char range, char uplo, integer n,
                    integer ka, integer kb, dcomplex* ab,
                    integer ldab, dcomplex* bb, integer ldbb, dcomplex* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    dcomplex* z, integer ldz, integer* ifail )
{
    integer info;
    dcomplex *work = SAFE_MALLOC(dcomplex,n);
    double *rwork = SAFE_MALLOC(double,7*n);
    integer *iwork = SAFE_MALLOC(integer,5*n);
    FFUNC(zhbgvx,ZHBGVX)(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, &vl, &vu, &il, &iu, &abstol,
                            m, w, z, &ldz, work, rwork, iwork, ifail, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    return info;
}

static inline integer c_sggev( char jobvl, char jobvr, integer n, float* a,
                   integer lda, float* b, integer ldb, scomplex* alpha,
                   float* beta, float* vl, integer ldvl,
                   float* vr, integer ldvr)
{
    integer info, lwork = -1;
    float *work, query;
    float *alphar = SAFE_MALLOC(float,n);
    float *alphai = SAFE_MALLOC(float,n);
    FFUNC(sggev,SGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sggev,SGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    scopy(n,alphar,1,(float*)alpha  ,2);
    scopy(n,alphai,1,(float*)alpha+1,2);
    FREE(alphar);
    FREE(alphai);
    FREE(work);
    return info;
}

static inline integer c_dggev( char jobvl, char jobvr, integer n, double* a,
                   integer lda, double* b, integer ldb, dcomplex* alpha,
                   double* beta, double* vl, integer ldvl,
                   double* vr, integer ldvr)
{
    integer info, lwork = -1;
    double *work, query;
    double *alphar = SAFE_MALLOC(double,n);
    double *alphai = SAFE_MALLOC(double,n);
    FFUNC(dggev,DGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dggev,DGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    dcopy(n,alphar,1,(double*)alpha  ,2);
    dcopy(n,alphai,1,(double*)alpha+1,2);
    FREE(alphar);
    FREE(alphai);
    FREE(work);
    return info;
}

static inline integer c_cggev( char jobvl, char jobvr, integer n, scomplex* a,
                   integer lda, scomplex* b, integer ldb, scomplex* alpha,
                   scomplex* beta, scomplex* vl, integer ldvl,
                   scomplex* vr, integer ldvr)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,8*n);
    FFUNC(cggev,CGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, &query, &lwork, rwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cggev,CGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_zggev( char jobvl, char jobvr, integer n, dcomplex* a,
                   integer lda, dcomplex* b, integer ldb, dcomplex* alpha,
                   dcomplex* beta, dcomplex* vl, integer ldvl,
                   dcomplex* vr, integer ldvr)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,8*n);
    FFUNC(zggev,ZGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, &query, &lwork, rwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zggev,ZGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
    FREE(work);
    FREE(rwork);
    return info;
}

static inline integer c_sggevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, float* a, integer lda, float* b,
                    integer ldb, scomplex* alpha,
                    float* beta, float* vl, integer ldvl, float* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    float* lscale, float* rscale, float* abnrm,
                    float* bbnrm, float* rconde, float* rcondv)
{
    integer info, lwork = -1;
    float *work, query;
    float *alphar = SAFE_MALLOC(float,n);
    float *alphai = SAFE_MALLOC(float,n);
    integer *iwork = SAFE_MALLOC(integer,n+6);
    logical *bwork = SAFE_MALLOC(logical,n);
    FFUNC(sggevx,SGGEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr,
                            ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, &query, &lwork, iwork, bwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sggevx,SGGEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr,
                            ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, iwork, bwork, &info);
    scopy(n,alphar,1,(float*)alpha  ,2);
    scopy(n,alphai,1,(float*)alpha+1,2);
    FREE(alphar);
    FREE(alphai);
    FREE(work);
    FREE(iwork);
    FREE(bwork);
    return info;
}

static inline integer c_dggevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, double* a, integer lda, double* b,
                    integer ldb, dcomplex* alpha,
                    double* beta, double* vl, integer ldvl, double* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    double* lscale, double* rscale, double* abnrm,
                    double* bbnrm, double* rconde, double* rcondv)
{
    integer info, lwork = -1;
    double *work, query;
    double *alphar = SAFE_MALLOC(double,n);
    double *alphai = SAFE_MALLOC(double,n);
    integer *iwork = SAFE_MALLOC(integer,n+6);
    logical *bwork = SAFE_MALLOC(logical,n);
    FFUNC(dggevx,DGGEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr,
                            ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, &query, &lwork, iwork, bwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dggevx,DGGEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr,
                            ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, iwork, bwork, &info);
    dcopy(n,alphar,1,(double*)alpha  ,2);
    dcopy(n,alphai,1,(double*)alpha+1,2);
    FREE(alphar);
    FREE(alphai);
    FREE(work);
    FREE(iwork);
    FREE(bwork);
    return info;
}

static inline integer c_cggevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, scomplex* a, integer lda, scomplex* b,
                    integer ldb, scomplex* alpha,
                    scomplex* beta, scomplex* vl, integer ldvl, scomplex* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    float* lscale, float* rscale, float* abnrm,
                    float* bbnrm, float* rconde, float* rcondv)
{
    integer info, lwork = -1;
    scomplex *work, query;
    float *rwork = SAFE_MALLOC(float,6*n);
    integer *iwork = SAFE_MALLOC(integer,n+6);
    logical *bwork = SAFE_MALLOC(logical,n);
    FFUNC(cggevx,CGGEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr,
                            ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, &query, &lwork, rwork, iwork, bwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cggevx,CGGEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr,
                            ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    FREE(bwork);
    return info;
}

static inline integer c_zggevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, dcomplex* a, integer lda, dcomplex* b,
                    integer ldb, dcomplex* alpha,
                    dcomplex* beta, dcomplex* vl, integer ldvl, dcomplex* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    double* lscale, double* rscale, double* abnrm,
                    double* bbnrm, double* rconde, double* rcondv)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    double *rwork = SAFE_MALLOC(double,6*n);
    integer *iwork = SAFE_MALLOC(integer,n+6);
    logical *bwork = SAFE_MALLOC(logical,n);
    FFUNC(zggevx,ZGGEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr,
                            ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, &query, &lwork, rwork, iwork, bwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zggevx,ZGGEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr,
                            ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, &info);
    FREE(work);
    FREE(rwork);
    FREE(iwork);
    FREE(bwork);
    return info;
}

static inline void c_ssfrk( char transr, char uplo, char trans, integer n,
                   integer k, float alpha, const float* a,
                   integer lda, float beta, float* c )
{
    FFUNC(ssfrk,SSFRK)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c);
}

static inline void c_dsfrk( char transr, char uplo, char trans, integer n,
                   integer k, double alpha, const double* a,
                   integer lda, double beta, double* c )
{
    FFUNC(dsfrk,DSFRK)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c);
}

static inline void c_chfrk( char transr, char uplo, char trans, integer n,
                   integer k, float alpha, const scomplex* a,
                   integer lda, float beta, scomplex* c )
{
    FFUNC(chfrk,CHFRK)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c);
}

static inline void c_zhfrk( char transr, char uplo, char trans, integer n,
                   integer k, double alpha, const dcomplex* a,
                   integer lda, double beta, dcomplex* c )
{
    FFUNC(zhfrk,ZHFRK)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c);
}

static inline void c_stfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, float alpha,
                   const float* a, float* b, integer ldb )
{
    FFUNC(stfsm,STFSM)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb);
}

static inline void c_dtfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, double alpha,
                   const double* a, double* b, integer ldb )
{
    FFUNC(dtfsm,DTFSM)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb);
}

static inline void c_ctfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, scomplex alpha,
                   const scomplex* a, scomplex* b, integer ldb )
{
    FFUNC(ctfsm,CTFSM)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb);
}

static inline void c_ztfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, dcomplex alpha,
                   const dcomplex* a, dcomplex* b, integer ldb )
{
    FFUNC(ztfsm,ZTFSM)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb);
}

static inline integer c_stfttp( char transr, char uplo, integer n, const float* arf,
                    float* ap )
{
    integer info;
    FFUNC(stfttp,STFTTP)(&transr, &uplo, &n, arf, ap, &info);
    return info;
}

static inline integer c_dtfttp( char transr, char uplo, integer n, const double* arf,
                    double* ap )
{
    integer info;
    FFUNC(dtfttp,DTFTTP)(&transr, &uplo, &n, arf, ap, &info);
    return info;
}

static inline integer c_ctfttp( char transr, char uplo, integer n, const scomplex* arf,
                    scomplex* ap )
{
    integer info;
    FFUNC(ctfttp,CTFTTP)(&transr, &uplo, &n, arf, ap, &info);
    return info;
}

static inline integer c_ztfttp( char transr, char uplo, integer n, const dcomplex* arf,
                    dcomplex* ap )
{
    integer info;
    FFUNC(ztfttp,ZTFTTP)(&transr, &uplo, &n, arf, ap, &info);
    return info;
}

static inline integer c_stfttr( char transr, char uplo, integer n, const float* arf,
                    float* a, integer lda )
{
    integer info;
    FFUNC(stfttr,STFTTR)(&transr, &uplo, &n, arf, a, &lda, &info);
    return info;
}

static inline integer c_dtfttr( char transr, char uplo, integer n, const double* arf,
                    double* a, integer lda )
{
    integer info;
    FFUNC(dtfttr,DTFTTR)(&transr, &uplo, &n, arf, a, &lda, &info);
    return info;
}

static inline integer c_ctfttr( char transr, char uplo, integer n, const scomplex* arf,
                    scomplex* a, integer lda )
{
    integer info;
    FFUNC(ctfttr,CTFTTR)(&transr, &uplo, &n, arf, a, &lda, &info);
    return info;
}

static inline integer c_ztfttr( char transr, char uplo, integer n, const dcomplex* arf,
                    dcomplex* a, integer lda )
{
    integer info;
    FFUNC(ztfttr,ZTFTTR)(&transr, &uplo, &n, arf, a, &lda, &info);
    return info;
}

static inline integer c_stpttf( char transr, char uplo, integer n, const float* ap,
                    float* arf )
{
    integer info;
    FFUNC(stpttf,STPTTF)(&transr, &uplo, &n, ap, arf, &info);
    return info;
}

static inline integer c_dtpttf( char transr, char uplo, integer n, const double* ap,
                    double* arf )
{
    integer info;
    FFUNC(dtpttf,DTPTTF)(&transr, &uplo, &n, ap, arf, &info);
    return info;
}

static inline integer c_ctpttf( char transr, char uplo, integer n, const scomplex* ap,
                    scomplex* arf )
{
    integer info;
    FFUNC(ctpttf,CTPTTF)(&transr, &uplo, &n, ap, arf, &info);
    return info;
}

static inline integer c_ztpttf( char transr, char uplo, integer n, const dcomplex* ap,
                    dcomplex* arf )
{
    integer info;
    FFUNC(ztpttf,ZTPTTF)(&transr, &uplo, &n, ap, arf, &info);
    return info;
}

static inline integer c_stpttr( char uplo, integer n, const float* ap, float* a,
                    integer lda )
{
    integer info;
    FFUNC(stpttr,STPTTR)(&uplo, &n, ap, a, &lda, &info);
    return info;
}

static inline integer c_dtpttr( char uplo, integer n, const double* ap, double* a,
                    integer lda )
{
    integer info;
    FFUNC(dtpttr,DTPTTR)(&uplo, &n, ap, a, &lda, &info);
    return info;
}

static inline integer c_ctpttr( char uplo, integer n, const scomplex* ap, scomplex* a,
                    integer lda )
{
    integer info;
    FFUNC(ctpttr,CTPTTR)(&uplo, &n, ap, a, &lda, &info);
    return info;
}

static inline integer c_ztpttr( char uplo, integer n, const dcomplex* ap, dcomplex* a,
                    integer lda )
{
    integer info;
    FFUNC(ztpttr,ZTPTTR)(&uplo, &n, ap, a, &lda, &info);
    return info;
}

static inline integer c_strttf( char transr, char uplo, integer n, const float* a,
                    integer lda, float* arf )
{
    integer info;
    FFUNC(strttf,STRTTF)(&transr, &uplo, &n, a, &lda, arf, &info);
    return info;
}

static inline integer c_dtrttf( char transr, char uplo, integer n, const double* a,
                    integer lda, double* arf )
{
    integer info;
    FFUNC(dtrttf,DTRTTF)(&transr, &uplo, &n, a, &lda, arf, &info);
    return info;
}

static inline integer c_ctrttf( char transr, char uplo, integer n, const scomplex* a,
                    integer lda, scomplex* arf )
{
    integer info;
    FFUNC(ctrttf,CTRTTF)(&transr, &uplo, &n, a, &lda, arf, &info);
    return info;
}

static inline integer c_ztrttf( char transr, char uplo, integer n, const dcomplex* a,
                    integer lda, dcomplex* arf )
{
    integer info;
    FFUNC(ztrttf,ZTRTTF)(&transr, &uplo, &n, a, &lda, arf, &info);
    return info;
}

static inline integer c_strttp( char uplo, integer n, const float* a, integer lda,
                    float* ap )
{
    integer info;
    FFUNC(strttp,STRTTP)(&uplo, &n, a, &lda, ap, &info);
    return info;
}

static inline integer c_dtrttp( char uplo, integer n, const double* a, integer lda,
                    double* ap )
{
    integer info;
    FFUNC(dtrttp,DTRTTP)(&uplo, &n, a, &lda, ap, &info);
    return info;
}

static inline integer c_ctrttp( char uplo, integer n, const scomplex* a, integer lda,
                    scomplex* ap )
{
    integer info;
    FFUNC(ctrttp,CTRTTP)(&uplo, &n, a, &lda, ap, &info);
    return info;
}

static inline integer c_ztrttp( char uplo, integer n, const dcomplex* a, integer lda,
                    dcomplex* ap )
{
    integer info;
    FFUNC(ztrttp,ZTRTTP)(&uplo, &n, a, &lda, ap, &info);
    return info;
}

static inline integer c_sgeqrfp( integer m, integer n, float* a, integer lda,
                     float* tau)
{
    integer info, lwork = -1;
    float *work, query;
    FFUNC(sgeqrfp,SGEQRFP)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(float,lwork);
    FFUNC(sgeqrfp,SGEQRFP)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_dgeqrfp( integer m, integer n, double* a, integer lda,
                     double* tau)
{
    integer info, lwork = -1;
    double *work, query;
    FFUNC(dgeqrfp,DGEQRFP)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)query;
    work = SAFE_MALLOC(double,lwork);
    FFUNC(dgeqrfp,DGEQRFP)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_cgeqrfp( integer m, integer n, scomplex* a, integer lda,
                     scomplex* tau)
{
    integer info, lwork = -1;
    scomplex *work, query;
    FFUNC(cgeqrfp,CGEQRFP)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)crealf(query);
    work = SAFE_MALLOC(scomplex,lwork);
    FFUNC(cgeqrfp,CGEQRFP)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

static inline integer c_zgeqrfp( integer m, integer n, dcomplex* a, integer lda,
                     dcomplex* tau)
{
    integer info, lwork = -1;
    dcomplex *work, query;
    FFUNC(zgeqrfp,ZGEQRFP)(&m, &n, a, &lda, tau, &query, &lwork, &info);
    lwork = (integer)creal(query);
    work = SAFE_MALLOC(dcomplex,lwork);
    FFUNC(zgeqrfp,ZGEQRFP)(&m, &n, a, &lda, tau, work, &lwork, &info);
    FREE(work);
    return info;
}

#ifdef __cplusplus
}
#endif

/*
 * #define more familiar names for the C versions
 */
#define sgetrf c_sgetrf
#define dgetrf c_dgetrf
#define cgetrf c_cgetrf
#define zgetrf c_zgetrf
#define sgbtrf c_sgbtrf
#define dgbtrf c_dgbtrf
#define cgbtrf c_cgbtrf
#define zgbtrf c_zgbtrf
#define sgttrf c_sgttrf
#define dgttrf c_dgttrf
#define cgttrf c_cgttrf
#define zgttrf c_zgttrf
#define spotrf c_spotrf
#define dpotrf c_dpotrf
#define cpotrf c_cpotrf
#define zpotrf c_zpotrf
#define spstrf c_spstrf
#define dpstrf c_dpstrf
#define cpstrf c_cpstrf
#define zpstrf c_zpstrf
#define spftrf c_spftrf
#define dpftrf c_dpftrf
#define cpftrf c_cpftrf
#define zpftrf c_zpftrf
#define spptrf c_spptrf
#define dpptrf c_dpptrf
#define cpptrf c_cpptrf
#define zpptrf c_zpptrf
#define spbtrf c_spbtrf
#define dpbtrf c_dpbtrf
#define cpbtrf c_cpbtrf
#define zpbtrf c_zpbtrf
#define spttrf c_spttrf
#define dpttrf c_dpttrf
#define cpttrf c_cpttrf
#define zpttrf c_zpttrf
#define ssytrf c_ssytrf
#define dsytrf c_dsytrf
#define csytrf c_csytrf
#define zsytrf c_zsytrf
#define chetrf c_chetrf
#define zhetrf c_zhetrf
#define ssptrf c_ssptrf
#define dsptrf c_dsptrf
#define csptrf c_csptrf
#define zsptrf c_zsptrf
#define chptrf c_chptrf
#define zhptrf c_zhptrf
#define sgetrs c_sgetrs
#define dgetrs c_dgetrs
#define cgetrs c_cgetrs
#define zgetrs c_zgetrs
#define sgbtrs c_sgbtrs
#define dgbtrs c_dgbtrs
#define cgbtrs c_cgbtrs
#define zgbtrs c_zgbtrs
#define sgttrs c_sgttrs
#define dgttrs c_dgttrs
#define cgttrs c_cgttrs
#define zgttrs c_zgttrs
#define spotrs c_spotrs
#define dpotrs c_dpotrs
#define cpotrs c_cpotrs
#define zpotrs c_zpotrs
#define spftrs c_spftrs
#define dpftrs c_dpftrs
#define cpftrs c_cpftrs
#define zpftrs c_zpftrs
#define spptrs c_spptrs
#define dpptrs c_dpptrs
#define cpptrs c_cpptrs
#define zpptrs c_zpptrs
#define spbtrs c_spbtrs
#define dpbtrs c_dpbtrs
#define cpbtrs c_cpbtrs
#define zpbtrs c_zpbtrs
#define spttrs c_spttrs
#define dpttrs c_dpttrs
#define cpttrs c_cpttrs
#define zpttrs c_zpttrs
#define ssytrs c_ssytrs
#define dsytrs c_dsytrs
#define csytrs c_csytrs
#define zsytrs c_zsytrs
#define chetrs c_chetrs
#define zhetrs c_zhetrs
#define ssptrs c_ssptrs
#define dsptrs c_dsptrs
#define csptrs c_csptrs
#define zsptrs c_zsptrs
#define chptrs c_chptrs
#define zhptrs c_zhptrs
#define strtrs c_strtrs
#define dtrtrs c_dtrtrs
#define ctrtrs c_ctrtrs
#define ztrtrs c_ztrtrs
#define stptrs c_stptrs
#define dtptrs c_dtptrs
#define ctptrs c_ctptrs
#define ztptrs c_ztptrs
#define stbtrs c_stbtrs
#define dtbtrs c_dtbtrs
#define ctbtrs c_ctbtrs
#define ztbtrs c_ztbtrs
#define sgecon c_sgecon
#define dgecon c_dgecon
#define cgecon c_cgecon
#define zgecon c_zgecon
#define sgbcon c_sgbcon
#define dgbcon c_dgbcon
#define cgbcon c_cgbcon
#define zgbcon c_zgbcon
#define sgtcon c_sgtcon
#define dgtcon c_dgtcon
#define cgtcon c_cgtcon
#define zgtcon c_zgtcon
#define spocon c_spocon
#define dpocon c_dpocon
#define cpocon c_cpocon
#define zpocon c_zpocon
#define sppcon c_sppcon
#define dppcon c_dppcon
#define cppcon c_cppcon
#define zppcon c_zppcon
#define spbcon c_spbcon
#define dpbcon c_dpbcon
#define cpbcon c_cpbcon
#define zpbcon c_zpbcon
#define sptcon c_sptcon
#define dptcon c_dptcon
#define cptcon c_cptcon
#define zptcon c_zptcon
#define ssycon c_ssycon
#define dsycon c_dsycon
#define csycon c_csycon
#define zsycon c_zsycon
#define checon c_checon
#define zhecon c_zhecon
#define sspcon c_sspcon
#define dspcon c_dspcon
#define cspcon c_cspcon
#define zspcon c_zspcon
#define chpcon c_chpcon
#define zhpcon c_zhpcon
#define strcon c_strcon
#define dtrcon c_dtrcon
#define ctrcon c_ctrcon
#define ztrcon c_ztrcon
#define stpcon c_stpcon
#define dtpcon c_dtpcon
#define ctpcon c_ctpcon
#define ztpcon c_ztpcon
#define stbcon c_stbcon
#define dtbcon c_dtbcon
#define ctbcon c_ctbcon
#define ztbcon c_ztbcon
#define sgerfs c_sgerfs
#define dgerfs c_dgerfs
#define cgerfs c_cgerfs
#define zgerfs c_zgerfs
#define sgerfsx c_sgerfsx
#define dgerfsx c_dgerfsx
#define cgerfsx c_cgerfsx
#define zgerfsx c_zgerfsx
#define sgbrfs c_sgbrfs
#define dgbrfs c_dgbrfs
#define cgbrfs c_cgbrfs
#define zgbrfs c_zgbrfs
#define sgbrfsx c_sgbrfsx
#define dgbrfsx c_dgbrfsx
#define cgbrfsx c_cgbrfsx
#define zgbrfsx c_zgbrfsx
#define sgtrfs c_sgtrfs
#define dgtrfs c_dgtrfs
#define cgtrfs c_cgtrfs
#define zgtrfs c_zgtrfs
#define sporfs c_sporfs
#define dporfs c_dporfs
#define cporfs c_cporfs
#define zporfs c_zporfs
#define sporfsx c_sporfsx
#define dporfsx c_dporfsx
#define cporfsx c_cporfsx
#define zporfsx c_zporfsx
#define spprfs c_spprfs
#define dpprfs c_dpprfs
#define cpprfs c_cpprfs
#define zpprfs c_zpprfs
#define spbrfs c_spbrfs
#define dpbrfs c_dpbrfs
#define cpbrfs c_cpbrfs
#define zpbrfs c_zpbrfs
#define sptrfs c_sptrfs
#define dptrfs c_dptrfs
#define cptrfs c_cptrfs
#define zptrfs c_zptrfs
#define ssyrfs c_ssyrfs
#define dsyrfs c_dsyrfs
#define csyrfs c_csyrfs
#define zsyrfs c_zsyrfs
#define cherfs c_cherfs
#define zherfs c_zherfs
#define ssyrfsx c_ssyrfsx
#define dsyrfsx c_dsyrfsx
#define csyrfsx c_csyrfsx
#define zsyrfsx c_zsyrfsx
#define cherfsx c_cherfsx
#define zherfsx c_zherfsx
#define ssprfs c_ssprfs
#define dsprfs c_dsprfs
#define csprfs c_csprfs
#define zsprfs c_zsprfs
#define chprfs c_chprfs
#define zhprfs c_zhprfs
#define strrfs c_strrfs
#define dtrrfs c_dtrrfs
#define ctrrfs c_ctrrfs
#define ztrrfs c_ztrrfs
#define stprfs c_stprfs
#define dtprfs c_dtprfs
#define ctprfs c_ctprfs
#define ztprfs c_ztprfs
#define stbrfs c_stbrfs
#define dtbrfs c_dtbrfs
#define ctbrfs c_ctbrfs
#define ztbrfs c_ztbrfs
#define sgetri c_sgetri
#define dgetri c_dgetri
#define cgetri c_cgetri
#define zgetri c_zgetri
#define spotri c_spotri
#define dpotri c_dpotri
#define cpotri c_cpotri
#define zpotri c_zpotri
#define spftri c_spftri
#define dpftri c_dpftri
#define cpftri c_cpftri
#define zpftri c_zpftri
#define spptri c_spptri
#define dpptri c_dpptri
#define cpptri c_cpptri
#define zpptri c_zpptri
#define ssytri c_ssytri
#define dsytri c_dsytri
#define csytri c_csytri
#define zsytri c_zsytri
#define chetri c_chetri
#define zhetri c_zhetri
#define ssptri c_ssptri
#define dsptri c_dsptri
#define csptri c_csptri
#define zsptri c_zsptri
#define chptri c_chptri
#define zhptri c_zhptri
#define strtri c_strtri
#define dtrtri c_dtrtri
#define ctrtri c_ctrtri
#define ztrtri c_ztrtri
#define stftri c_stftri
#define dtftri c_dtftri
#define ctftri c_ctftri
#define ztftri c_ztftri
#define stptri c_stptri
#define dtptri c_dtptri
#define ctptri c_ctptri
#define ztptri c_ztptri
#define sgeequ c_sgeequ
#define dgeequ c_dgeequ
#define cgeequ c_cgeequ
#define zgeequ c_zgeequ
#define sgeequb c_sgeequb
#define dgeequb c_dgeequb
#define cgeequb c_cgeequb
#define zgeequb c_zgeequb
#define sgbequ c_sgbequ
#define dgbequ c_dgbequ
#define cgbequ c_cgbequ
#define zgbequ c_zgbequ
#define sgbequb c_sgbequb
#define dgbequb c_dgbequb
#define cgbequb c_cgbequb
#define zgbequb c_zgbequb
#define spoequ c_spoequ
#define dpoequ c_dpoequ
#define cpoequ c_cpoequ
#define zpoequ c_zpoequ
#define spoequb c_spoequb
#define dpoequb c_dpoequb
#define cpoequb c_cpoequb
#define zpoequb c_zpoequb
#define sppequ c_sppequ
#define dppequ c_dppequ
#define cppequ c_cppequ
#define zppequ c_zppequ
#define spbequ c_spbequ
#define dpbequ c_dpbequ
#define cpbequ c_cpbequ
#define zpbequ c_zpbequ
#define ssyequb c_ssyequb
#define dsyequb c_dsyequb
#define csyequb c_csyequb
#define zsyequb c_zsyequb
#define cheequb c_cheequb
#define zheequb c_zheequb
#define sgesv c_sgesv
#define dgesv c_dgesv
#define cgesv c_cgesv
#define zgesv c_zgesv
#define dsgesv c_dsgesv
#define zcgesv c_zcgesv
#define sgesvx c_sgesvx
#define dgesvx c_dgesvx
#define cgesvx c_cgesvx
#define zgesvx c_zgesvx
#define sgesvxx c_sgesvxx
#define dgesvxx c_dgesvxx
#define cgesvxx c_cgesvxx
#define zgesvxx c_zgesvxx
#define sgbsv c_sgbsv
#define dgbsv c_dgbsv
#define cgbsv c_cgbsv
#define zgbsv c_zgbsv
#define sgbsvx c_sgbsvx
#define dgbsvx c_dgbsvx
#define cgbsvx c_cgbsvx
#define zgbsvx c_zgbsvx
#define sgbsvxx c_sgbsvxx
#define dgbsvxx c_dgbsvxx
#define cgbsvxx c_cgbsvxx
#define zgbsvxx c_zgbsvxx
#define sgtsv c_sgtsv
#define dgtsv c_dgtsv
#define cgtsv c_cgtsv
#define zgtsv c_zgtsv
#define sgtsvx c_sgtsvx
#define dgtsvx c_dgtsvx
#define cgtsvx c_cgtsvx
#define zgtsvx c_zgtsvx
#define sposv c_sposv
#define dposv c_dposv
#define cposv c_cposv
#define zposv c_zposv
#define dsposv c_dsposv
#define zcposv c_zcposv
#define sposvx c_sposvx
#define dposvx c_dposvx
#define cposvx c_cposvx
#define zposvx c_zposvx
#define sposvxx c_sposvxx
#define dposvxx c_dposvxx
#define cposvxx c_cposvxx
#define zposvxx c_zposvxx
#define sppsv c_sppsv
#define dppsv c_dppsv
#define cppsv c_cppsv
#define zppsv c_zppsv
#define sppsvx c_sppsvx
#define dppsvx c_dppsvx
#define cppsvx c_cppsvx
#define zppsvx c_zppsvx
#define spbsv c_spbsv
#define dpbsv c_dpbsv
#define cpbsv c_cpbsv
#define zpbsv c_zpbsv
#define spbsvx c_spbsvx
#define dpbsvx c_dpbsvx
#define cpbsvx c_cpbsvx
#define zpbsvx c_zpbsvx
#define sptsv c_sptsv
#define dptsv c_dptsv
#define cptsv c_cptsv
#define zptsv c_zptsv
#define sptsvx c_sptsvx
#define dptsvx c_dptsvx
#define cptsvx c_cptsvx
#define zptsvx c_zptsvx
#define ssysv c_ssysv
#define dsysv c_dsysv
#define csysv c_csysv
#define zsysv c_zsysv
#define chesv c_chesv
#define zhesv c_zhesv
#define ssysvx c_ssysvx
#define dsysvx c_dsysvx
#define csysvx c_csysvx
#define zsysvx c_zsysvx
#define chesvx c_chesvx
#define zhesvx c_zhesvx
#define ssysvxx c_ssysvxx
#define dsysvxx c_dsysvxx
#define csysvxx c_csysvxx
#define zsysvxx c_zsysvxx
#define chesvxx c_chesvxx
#define zhesvxx c_zhesvxx
#define sspsv c_sspsv
#define dspsv c_dspsv
#define cspsv c_cspsv
#define zspsv c_zspsv
#define chpsv c_chpsv
#define zhpsv c_zhpsv
#define sspsvx c_sspsvx
#define dspsvx c_dspsvx
#define cspsvx c_cspsvx
#define zspsvx c_zspsvx
#define chpsvx c_chpsvx
#define zhpsvx c_zhpsvx
#define sgeqrf c_sgeqrf
#define dgeqrf c_dgeqrf
#define cgeqrf c_cgeqrf
#define zgeqrf c_zgeqrf
#define sgeqpf c_sgeqpf
#define dgeqpf c_dgeqpf
#define cgeqpf c_cgeqpf
#define zgeqpf c_zgeqpf
#define sgeqp3 c_sgeqp3
#define dgeqp3 c_dgeqp3
#define cgeqp3 c_cgeqp3
#define zgeqp3 c_zgeqp3
#define sorgqr c_sorgqr
#define dorgqr c_dorgqr
#define cungqr c_cungqr
#define zungqr c_zungqr
#define sormqr c_sormqr
#define dormqr c_dormqr
#define cunmqr c_cunmqr
#define zunmqr c_zunmqr
#define sgelqf c_sgelqf
#define dgelqf c_dgelqf
#define cgelqf c_cgelqf
#define zgelqf c_zgelqf
#define sorglq c_sorglq
#define dorglq c_dorglq
#define cunglq c_cunglq
#define zunglq c_zunglq
#define sormlq c_sormlq
#define dormlq c_dormlq
#define cunmlq c_cunmlq
#define zunmlq c_zunmlq
#define sgeqlf c_sgeqlf
#define dgeqlf c_dgeqlf
#define cgeqlf c_cgeqlf
#define zgeqlf c_zgeqlf
#define sorgql c_sorgql
#define dorgql c_dorgql
#define cungql c_cungql
#define zungql c_zungql
#define sormql c_sormql
#define dormql c_dormql
#define cunmql c_cunmql
#define zunmql c_zunmql
#define sgerqf c_sgerqf
#define dgerqf c_dgerqf
#define cgerqf c_cgerqf
#define zgerqf c_zgerqf
#define sorgrq c_sorgrq
#define dorgrq c_dorgrq
#define cungrq c_cungrq
#define zungrq c_zungrq
#define sormrq c_sormrq
#define dormrq c_dormrq
#define cunmrq c_cunmrq
#define zunmrq c_zunmrq
#define stzrzf c_stzrzf
#define dtzrzf c_dtzrzf
#define ctzrzf c_ctzrzf
#define ztzrzf c_ztzrzf
#define sormrz c_sormrz
#define dormrz c_dormrz
#define cunmrz c_cunmrz
#define zunmrz c_zunmrz
#define sggqrf c_sggqrf
#define dggqrf c_dggqrf
#define cggqrf c_cggqrf
#define zggqrf c_zggqrf
#define sggrqf c_sggrqf
#define dggrqf c_dggrqf
#define cggrqf c_cggrqf
#define zggrqf c_zggrqf
#define sgebrd c_sgebrd
#define dgebrd c_dgebrd
#define cgebrd c_cgebrd
#define zgebrd c_zgebrd
#define sgbbrd c_sgbbrd
#define dgbbrd c_dgbbrd
#define cgbbrd c_cgbbrd
#define zgbbrd c_zgbbrd
#define sorgbr c_sorgbr
#define dorgbr c_dorgbr
#define cungbr c_cungbr
#define zungbr c_zungbr
#define sormbr c_sormbr
#define dormbr c_dormbr
#define cunmbr c_cunmbr
#define zunmbr c_zunmbr
#define sbdsqr c_sbdsqr
#define dbdsqr c_dbdsqr
#define cbdsqr c_cbdsqr
#define zbdsqr c_zbdsqr
#define sbdsdc c_sbdsdc
#define dbdsdc c_dbdsdc
#define ssytrd c_ssytrd
#define dsytrd c_dsytrd
#define chetrd c_chetrd
#define zhetrd c_zhetrd
#define sorgtr c_sorgtr
#define dorgtr c_dorgtr
#define cungtr c_cungtr
#define zungtr c_zungtr
#define sormtr c_sormtr
#define dormtr c_dormtr
#define cunmtr c_cunmtr
#define zunmtr c_zunmtr
#define ssptrd c_ssptrd
#define dsptrd c_dsptrd
#define chptrd c_chptrd
#define zhptrd c_zhptrd
#define sopgtr c_sopgtr
#define dopgtr c_dopgtr
#define cupgtr c_cupgtr
#define zupgtr c_zupgtr
#define sopmtr c_sopmtr
#define dopmtr c_dopmtr
#define cupmtr c_cupmtr
#define zupmtr c_zupmtr
#define ssbtrd c_ssbtrd
#define dsbtrd c_dsbtrd
#define chbtrd c_chbtrd
#define zhbtrd c_zhbtrd
#define ssterf c_ssterf
#define dsterf c_dsterf
#define ssteqr c_ssteqr
#define dsteqr c_dsteqr
#define csteqr c_csteqr
#define zsteqr c_zsteqr
#define sstemr c_sstemr
#define dstemr c_dstemr
#define cstemr c_cstemr
#define zstemr c_zstemr
#define sstedc c_sstedc
#define dstedc c_dstedc
#define cstedc c_cstedc
#define zstedc c_zstedc
#define sstegr c_sstegr
#define dstegr c_dstegr
#define cstegr c_cstegr
#define zstegr c_zstegr
#define spteqr c_spteqr
#define dpteqr c_dpteqr
#define cpteqr c_cpteqr
#define zpteqr c_zpteqr
#define sstebz c_sstebz
#define dstebz c_dstebz
#define sstein c_sstein
#define dstein c_dstein
#define cstein c_cstein
#define zstein c_zstein
#define sdisna c_sdisna
#define ddisna c_ddisna
#define ssygst c_ssygst
#define dsygst c_dsygst
#define chegst c_chegst
#define zhegst c_zhegst
#define sspgst c_sspgst
#define dspgst c_dspgst
#define chpgst c_chpgst
#define zhpgst c_zhpgst
#define ssbgst c_ssbgst
#define dsbgst c_dsbgst
#define chbgst c_chbgst
#define zhbgst c_zhbgst
#define spbstf c_spbstf
#define dpbstf c_dpbstf
#define cpbstf c_cpbstf
#define zpbstf c_zpbstf
#define sgehrd c_sgehrd
#define dgehrd c_dgehrd
#define cgehrd c_cgehrd
#define zgehrd c_zgehrd
#define sorghr c_sorghr
#define dorghr c_dorghr
#define cunghr c_cunghr
#define zunghr c_zunghr
#define sormhr c_sormhr
#define dormhr c_dormhr
#define cunmhr c_cunmhr
#define zunmhr c_zunmhr
#define sgebal c_sgebal
#define dgebal c_dgebal
#define cgebal c_cgebal
#define zgebal c_zgebal
#define sgebak c_sgebak
#define dgebak c_dgebak
#define cgebak c_cgebak
#define zgebak c_zgebak
#define shseqr c_shseqr
#define dhseqr c_dhseqr
#define chseqr c_chseqr
#define zhseqr c_zhseqr
#define shsein c_shsein
#define dhsein c_dhsein
#define chsein c_chsein
#define zhsein c_zhsein
#define strevc c_strevc
#define dtrevc c_dtrevc
#define ctrevc c_ctrevc
#define ztrevc c_ztrevc
#define strsna c_strsna
#define dtrsna c_dtrsna
#define ctrsna c_ctrsna
#define ztrsna c_ztrsna
#define strexc c_strexc
#define dtrexc c_dtrexc
#define ctrexc c_ctrexc
#define ztrexc c_ztrexc
#define strsen c_strsen
#define dtrsen c_dtrsen
#define ctrsen c_ctrsen
#define ztrsen c_ztrsen
#define strsyl c_strsyl
#define dtrsyl c_dtrsyl
#define ctrsyl c_ctrsyl
#define ztrsyl c_ztrsyl
#define sgghrd c_sgghrd
#define dgghrd c_dgghrd
#define cgghrd c_cgghrd
#define zgghrd c_zgghrd
#define sggbal c_sggbal
#define dggbal c_dggbal
#define cggbal c_cggbal
#define zggbal c_zggbal
#define sggbak c_sggbak
#define dggbak c_dggbak
#define cggbak c_cggbak
#define zggbak c_zggbak
#define shgeqz c_shgeqz
#define dhgeqz c_dhgeqz
#define chgeqz c_chgeqz
#define zhgeqz c_zhgeqz
#define stgevc c_stgevc
#define dtgevc c_dtgevc
#define ctgevc c_ctgevc
#define ztgevc c_ztgevc
#define stgexc c_stgexc
#define dtgexc c_dtgexc
#define ctgexc c_ctgexc
#define ztgexc c_ztgexc
#define stgsen c_stgsen
#define dtgsen c_dtgsen
#define ctgsen c_ctgsen
#define ztgsen c_ztgsen
#define stgsyl c_stgsyl
#define dtgsyl c_dtgsyl
#define ctgsyl c_ctgsyl
#define ztgsyl c_ztgsyl
#define stgsna c_stgsna
#define dtgsna c_dtgsna
#define ctgsna c_ctgsna
#define ztgsna c_ztgsna
#define sggsvp c_sggsvp
#define dggsvp c_dggsvp
#define cggsvp c_cggsvp
#define zggsvp c_zggsvp
#define stgsja c_stgsja
#define dtgsja c_dtgsja
#define ctgsja c_ctgsja
#define ztgsja c_ztgsja
#define sgels c_sgels
#define dgels c_dgels
#define cgels c_cgels
#define zgels c_zgels
#define sgelsy c_sgelsy
#define dgelsy c_dgelsy
#define cgelsy c_cgelsy
#define zgelsy c_zgelsy
#define sgelss c_sgelss
#define dgelss c_dgelss
#define cgelss c_cgelss
#define zgelss c_zgelss
#define sgglse c_sgglse
#define dgglse c_dgglse
#define cgglse c_cgglse
#define zgglse c_zgglse
#define sggglm c_sggglm
#define dggglm c_dggglm
#define cggglm c_cggglm
#define zggglm c_zggglm
#define ssyev c_ssyev
#define dsyev c_dsyev
#define cheev c_cheev
#define zheev c_zheev
#define ssyevd c_ssyevd
#define dsyevd c_dsyevd
#define cheevd c_cheevd
#define zheevd c_zheevd
#define ssyevx c_ssyevx
#define dsyevx c_dsyevx
#define cheevx c_cheevx
#define zheevx c_zheevx
#define ssyevr c_ssyevr
#define dsyevr c_dsyevr
#define cheevr c_cheevr
#define zheevr c_zheevr
#define sspev c_sspev
#define dspev c_dspev
#define chpev c_chpev
#define zhpev c_zhpev
#define sspevd c_sspevd
#define dspevd c_dspevd
#define chpevd c_chpevd
#define zhpevd c_zhpevd
#define sspevx c_sspevx
#define dspevx c_dspevx
#define chpevx c_chpevx
#define zhpevx c_zhpevx
#define ssbev c_ssbev
#define dsbev c_dsbev
#define chbev c_chbev
#define zhbev c_zhbev
#define ssbevd c_ssbevd
#define dsbevd c_dsbevd
#define chbevd c_chbevd
#define zhbevd c_zhbevd
#define ssbevx c_ssbevx
#define dsbevx c_dsbevx
#define chbevx c_chbevx
#define zhbevx c_zhbevx
#define sstev c_sstev
#define dstev c_dstev
#define sstevd c_sstevd
#define dstevd c_dstevd
#define sstevx c_sstevx
#define dstevx c_dstevx
#define sstevr c_sstevr
#define dstevr c_dstevr
#define sgeev c_sgeev
#define dgeev c_dgeev
#define cgeev c_cgeev
#define zgeev c_zgeev
#define sgeevx c_sgeevx
#define dgeevx c_dgeevx
#define cgeevx c_cgeevx
#define zgeevx c_zgeevx
#define sgesvd c_sgesvd
#define dgesvd c_dgesvd
#define cgesvd c_cgesvd
#define zgesvd c_zgesvd
#define sgesdd c_sgesdd
#define dgesdd c_dgesdd
#define cgesdd c_cgesdd
#define zgesdd c_zgesdd
#define sgejsv c_sgejsv
#define dgejsv c_dgejsv
#define sgesvj c_sgesvj
#define dgesvj c_dgesvj
#define sggsvd c_sggsvd
#define dggsvd c_dggsvd
#define cggsvd c_cggsvd
#define zggsvd c_zggsvd
#define ssygv c_ssygv
#define dsygv c_dsygv
#define chegv c_chegv
#define zhegv c_zhegv
#define ssygvd c_ssygvd
#define dsygvd c_dsygvd
#define chegvd c_chegvd
#define zhegvd c_zhegvd
#define ssygvx c_ssygvx
#define dsygvx c_dsygvx
#define chegvx c_chegvx
#define zhegvx c_zhegvx
#define sspgv c_sspgv
#define dspgv c_dspgv
#define chpgv c_chpgv
#define zhpgv c_zhpgv
#define sspgvd c_sspgvd
#define dspgvd c_dspgvd
#define chpgvd c_chpgvd
#define zhpgvd c_zhpgvd
#define sspgvx c_sspgvx
#define dspgvx c_dspgvx
#define chpgvx c_chpgvx
#define zhpgvx c_zhpgvx
#define ssbgv c_ssbgv
#define dsbgv c_dsbgv
#define chbgv c_chbgv
#define zhbgv c_zhbgv
#define ssbgvd c_ssbgvd
#define dsbgvd c_dsbgvd
#define chbgvd c_chbgvd
#define zhbgvd c_zhbgvd
#define ssbgvx c_ssbgvx
#define dsbgvx c_dsbgvx
#define chbgvx c_chbgvx
#define zhbgvx c_zhbgvx
#define sggev c_sggev
#define dggev c_dggev
#define cggev c_cggev
#define zggev c_zggev
#define sggevx c_sggevx
#define dggevx c_dggevx
#define cggevx c_cggevx
#define zggevx c_zggevx
#define ssfrk c_ssfrk
#define dsfrk c_dsfrk
#define chfrk c_chfrk
#define zhfrk c_zhfrk
#define stfsm c_stfsm
#define dtfsm c_dtfsm
#define ctfsm c_ctfsm
#define ztfsm c_ztfsm
#define stfttp c_stfttp
#define dtfttp c_dtfttp
#define ctfttp c_ctfttp
#define ztfttp c_ztfttp
#define stfttr c_stfttr
#define dtfttr c_dtfttr
#define ctfttr c_ctfttr
#define ztfttr c_ztfttr
#define stpttf c_stpttf
#define dtpttf c_dtpttf
#define ctpttf c_ctpttf
#define ztpttf c_ztpttf
#define stpttr c_stpttr
#define dtpttr c_dtpttr
#define ctpttr c_ctpttr
#define ztpttr c_ztpttr
#define strttf c_strttf
#define dtrttf c_dtrttf
#define ctrttf c_ctrttf
#define ztrttf c_ztrttf
#define strttp c_strttp
#define dtrttp c_dtrttp
#define ctrttp c_ctrttp
#define ztrttp c_ztrttp
#define sgeqrfp c_sgeqrfp
#define dgeqrfp c_dgeqrfp
#define cgeqrfp c_cgeqrfp
#define zgeqrfp c_zgeqrfp

#ifdef __cplusplus

/*
 * C++ overloads
 */
inline integer getrf( integer m, integer n, float* a, integer lda,
                    integer* ipiv )
{
    return sgetrf(m, n, a, lda, ipiv);
}

inline integer getrf( integer m, integer n, double* a, integer lda,
                    integer* ipiv )
{
    return dgetrf(m, n, a, lda, ipiv);
}

inline integer getrf( integer m, integer n, scomplex* a, integer lda,
                    integer* ipiv )
{
    return cgetrf(m, n, a, lda, ipiv);
}

inline integer getrf( integer m, integer n, dcomplex* a, integer lda,
                    integer* ipiv )
{
    return zgetrf(m, n, a, lda, ipiv);
}

inline integer gbtrf( integer m, integer n, integer kl,
                    integer ku, float* ab, integer ldab,
                    integer* ipiv )
{
    return sgbtrf(m, n, kl, ku, ab, ldab, ipiv);
}

inline integer gbtrf( integer m, integer n, integer kl,
                    integer ku, double* ab, integer ldab,
                    integer* ipiv )
{
    return dgbtrf(m, n, kl, ku, ab, ldab, ipiv);
}

inline integer gbtrf( integer m, integer n, integer kl,
                    integer ku, scomplex* ab, integer ldab,
                    integer* ipiv )
{
    return cgbtrf(m, n, kl, ku, ab, ldab, ipiv);
}

inline integer gbtrf( integer m, integer n, integer kl,
                    integer ku, dcomplex* ab, integer ldab,
                    integer* ipiv )
{
    return zgbtrf(m, n, kl, ku, ab, ldab, ipiv);
}

inline integer gttrf( integer n, float* dl, float* d, float* du,
                    float* du2, integer* ipiv )
{
    return sgttrf(n, dl, d, du, du2, ipiv);
}

inline integer gttrf( integer n, double* dl, double* d, double* du,
                    double* du2, integer* ipiv )
{
    return dgttrf(n, dl, d, du, du2, ipiv);
}

inline integer gttrf( integer n, scomplex* dl, scomplex* d, scomplex* du,
                    scomplex* du2, integer* ipiv )
{
    return cgttrf(n, dl, d, du, du2, ipiv);
}

inline integer gttrf( integer n, dcomplex* dl, dcomplex* d, dcomplex* du,
                    dcomplex* du2, integer* ipiv )
{
    return zgttrf(n, dl, d, du, du2, ipiv);
}

inline integer potrf( char uplo, integer n, float* a, integer lda )
{
    return spotrf(uplo, n, a, lda);
}

inline integer potrf( char uplo, integer n, double* a, integer lda )
{
    return dpotrf(uplo, n, a, lda);
}

inline integer potrf( char uplo, integer n, scomplex* a, integer lda )
{
    return cpotrf(uplo, n, a, lda);
}

inline integer potrf( char uplo, integer n, dcomplex* a, integer lda )
{
    return zpotrf(uplo, n, a, lda);
}

inline integer pstrf( char uplo, integer n, float* a, integer lda,
                    integer* piv, integer& rank, float tol)
{
    return spstrf(uplo, n, a, lda, piv, &rank, tol);
}

inline integer pstrf( char uplo, integer n, double* a, integer lda,
                    integer* piv, integer& rank, double tol)
{
    return dpstrf(uplo, n, a, lda, piv, &rank, tol);
}

inline integer pstrf( char uplo, integer n, scomplex* a, integer lda,
                    integer* piv, integer& rank, float tol)
{
    return cpstrf(uplo, n, a, lda, piv, &rank, tol);
}

inline integer pstrf( char uplo, integer n, dcomplex* a, integer lda,
                    integer* piv, integer& rank, double tol)
{
    return zpstrf(uplo, n, a, lda, piv, &rank, tol);
}

inline integer pftrf( char transr, char uplo, integer n, float* a )
{
    return spftrf(transr, uplo, n, a);
}

inline integer pftrf( char transr, char uplo, integer n, double* a )
{
    return dpftrf(transr, uplo, n, a);
}

inline integer pftrf( char transr, char uplo, integer n, scomplex* a )
{
    return cpftrf(transr, uplo, n, a);
}

inline integer pftrf( char transr, char uplo, integer n, dcomplex* a )
{
    return zpftrf(transr, uplo, n, a);
}

inline integer pptrf( char uplo, integer n, float* ap )
{
    return spptrf(uplo, n, ap);
}

inline integer pptrf( char uplo, integer n, double* ap )
{
    return dpptrf(uplo, n, ap);
}

inline integer pptrf( char uplo, integer n, scomplex* ap )
{
    return cpptrf(uplo, n, ap);
}

inline integer pptrf( char uplo, integer n, dcomplex* ap )
{
    return zpptrf(uplo, n, ap);
}

inline integer pbtrf( char uplo, integer n, integer kd, float* ab,
                    integer ldab )
{
    return spbtrf(uplo, n, kd, ab, ldab);
}

inline integer pbtrf( char uplo, integer n, integer kd, double* ab,
                    integer ldab )
{
    return dpbtrf(uplo, n, kd, ab, ldab);
}

inline integer pbtrf( char uplo, integer n, integer kd, scomplex* ab,
                    integer ldab )
{
    return cpbtrf(uplo, n, kd, ab, ldab);
}

inline integer pbtrf( char uplo, integer n, integer kd, dcomplex* ab,
                    integer ldab )
{
    return zpbtrf(uplo, n, kd, ab, ldab);
}

inline integer pttrf( integer n, float* d, float* e )
{
    return spttrf(n, d, e);
}

inline integer pttrf( integer n, double* d, double* e )
{
    return dpttrf(n, d, e);
}

inline integer pttrf( integer n, float* d, scomplex* e )
{
    return cpttrf(n, d, e);
}

inline integer pttrf( integer n, double* d, dcomplex* e )
{
    return zpttrf(n, d, e);
}

inline integer sytrf( char uplo, integer n, float* a, integer lda,
                    integer* ipiv)
{
    return ssytrf(uplo, n, a, lda, ipiv);
}

inline integer sytrf( char uplo, integer n, double* a, integer lda,
                    integer* ipiv)
{
    return dsytrf(uplo, n, a, lda, ipiv);
}

inline integer sytrf( char uplo, integer n, scomplex* a, integer lda,
                    integer* ipiv)
{
    return csytrf(uplo, n, a, lda, ipiv);
}

inline integer sytrf( char uplo, integer n, dcomplex* a, integer lda,
                    integer* ipiv)
{
    return zsytrf(uplo, n, a, lda, ipiv);
}

inline integer hetrf( char uplo, integer n, float* a, integer lda,
                    integer* ipiv)
{
    return ssytrf(uplo, n, a, lda, ipiv);
}

inline integer hetrf( char uplo, integer n, double* a, integer lda,
                    integer* ipiv)
{
    return dsytrf(uplo, n, a, lda, ipiv);
}

inline integer hetrf( char uplo, integer n, scomplex* a, integer lda,
                    integer* ipiv)
{
    return chetrf(uplo, n, a, lda, ipiv);
}

inline integer hetrf( char uplo, integer n, dcomplex* a, integer lda,
                    integer* ipiv)
{
    return zhetrf(uplo, n, a, lda, ipiv);
}

inline integer sptrf( char uplo, integer n, float* ap, integer* ipiv )
{
    return ssptrf(uplo, n, ap, ipiv);
}

inline integer sptrf( char uplo, integer n, double* ap, integer* ipiv )
{
    return dsptrf(uplo, n, ap, ipiv);
}

inline integer sptrf( char uplo, integer n, scomplex* ap, integer* ipiv )
{
    return csptrf(uplo, n, ap, ipiv);
}

inline integer sptrf( char uplo, integer n, dcomplex* ap, integer* ipiv )
{
    return zsptrf(uplo, n, ap, ipiv);
}

inline integer hptrf( char uplo, integer n, float* ap, integer* ipiv )
{
    return ssptrf(uplo, n, ap, ipiv);
}

inline integer hptrf( char uplo, integer n, double* ap, integer* ipiv )
{
    return dsptrf(uplo, n, ap, ipiv);
}

inline integer hptrf( char uplo, integer n, scomplex* ap, integer* ipiv )
{
    return chptrf(uplo, n, ap, ipiv);
}

inline integer hptrf( char uplo, integer n, dcomplex* ap, integer* ipiv )
{
    return zhptrf(uplo, n, ap, ipiv);
}

inline integer getrs( char trans, integer n, integer nrhs,
                    const float* a, integer lda, const integer* ipiv,
                    float* b, integer ldb )
{
    return sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer getrs( char trans, integer n, integer nrhs,
                    const double* a, integer lda, const integer* ipiv,
                    double* b, integer ldb )
{
    return dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer getrs( char trans, integer n, integer nrhs,
                    const scomplex* a, integer lda, const integer* ipiv,
                    scomplex* b, integer ldb )
{
    return cgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer getrs( char trans, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const integer* ipiv,
                    dcomplex* b, integer ldb )
{
    return zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer gbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const float* ab, integer ldab,
                    const integer* ipiv, float* b, integer ldb )
{
    return sgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

inline integer gbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const double* ab, integer ldab,
                    const integer* ipiv, double* b, integer ldb )
{
    return dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

inline integer gbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const scomplex* ab, integer ldab,
                    const integer* ipiv, scomplex* b, integer ldb )
{
    return cgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

inline integer gbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const dcomplex* ab, integer ldab,
                    const integer* ipiv, dcomplex* b, integer ldb )
{
    return zgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

inline integer gttrs( char trans, integer n, integer nrhs,
                    const float* dl, const float* d, const float* du,
                    const float* du2, const integer* ipiv, float* b,
                    integer ldb )
{
    return sgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}

inline integer gttrs( char trans, integer n, integer nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* du2, const integer* ipiv, double* b,
                    integer ldb )
{
    return dgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}

inline integer gttrs( char trans, integer n, integer nrhs,
                    const scomplex* dl, const scomplex* d, const scomplex* du,
                    const scomplex* du2, const integer* ipiv, scomplex* b,
                    integer ldb )
{
    return cgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}

inline integer gttrs( char trans, integer n, integer nrhs,
                    const dcomplex* dl, const dcomplex* d, const dcomplex* du,
                    const dcomplex* du2, const integer* ipiv, dcomplex* b,
                    integer ldb )
{
    return zgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}

inline integer potrs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, float* b,
                    integer ldb )
{
    return spotrs(uplo, n, nrhs, a, lda, b, ldb);
}

inline integer potrs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, double* b,
                    integer ldb )
{
    return dpotrs(uplo, n, nrhs, a, lda, b, ldb);
}

inline integer potrs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, scomplex* b,
                    integer ldb )
{
    return cpotrs(uplo, n, nrhs, a, lda, b, ldb);
}

inline integer potrs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, dcomplex* b,
                    integer ldb )
{
    return zpotrs(uplo, n, nrhs, a, lda, b, ldb);
}

inline integer pftrs( char transr, char uplo, integer n, integer nrhs,
                    const float* a, float* b, integer ldb )
{
    return spftrs(transr, uplo, n, nrhs, a, b, ldb);
}

inline integer pftrs( char transr, char uplo, integer n, integer nrhs,
                    const double* a, double* b, integer ldb )
{
    return dpftrs(transr, uplo, n, nrhs, a, b, ldb);
}

inline integer pftrs( char transr, char uplo, integer n, integer nrhs,
                    const scomplex* a, scomplex* b, integer ldb )
{
    return cpftrs(transr, uplo, n, nrhs, a, b, ldb);
}

inline integer pftrs( char transr, char uplo, integer n, integer nrhs,
                    const dcomplex* a, dcomplex* b, integer ldb )
{
    return zpftrs(transr, uplo, n, nrhs, a, b, ldb);
}

inline integer pptrs( char uplo, integer n, integer nrhs,
                    const float* ap, float* b, integer ldb )
{
    return spptrs(uplo, n, nrhs, ap, b, ldb);
}

inline integer pptrs( char uplo, integer n, integer nrhs,
                    const double* ap, double* b, integer ldb )
{
    return dpptrs(uplo, n, nrhs, ap, b, ldb);
}

inline integer pptrs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, scomplex* b, integer ldb )
{
    return cpptrs(uplo, n, nrhs, ap, b, ldb);
}

inline integer pptrs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, dcomplex* b, integer ldb )
{
    return zpptrs(uplo, n, nrhs, ap, b, ldb);
}

inline integer pbtrs( char uplo, integer n, integer kd, integer nrhs,
                    const float* ab, integer ldab, float* b,
                    integer ldb )
{
    return spbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer pbtrs( char uplo, integer n, integer kd, integer nrhs,
                    const double* ab, integer ldab, double* b,
                    integer ldb )
{
    return dpbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer pbtrs( char uplo, integer n, integer kd, integer nrhs,
                    const scomplex* ab, integer ldab, scomplex* b,
                    integer ldb )
{
    return cpbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer pbtrs( char uplo, integer n, integer kd, integer nrhs,
                    const dcomplex* ab, integer ldab, dcomplex* b,
                    integer ldb )
{
    return zpbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer pttrs( integer n, integer nrhs, const float* d,
                    const float* e, float* b, integer ldb )
{
    return spttrs(n, nrhs, d, e, b, ldb);
}

inline integer pttrs( integer n, integer nrhs, const double* d,
                    const double* e, double* b, integer ldb )
{
    return dpttrs(n, nrhs, d, e, b, ldb);
}

inline integer pttrs( integer n, integer nrhs, const float* d,
                    const scomplex* e, scomplex* b, integer ldb )
{
    return cpttrs(n, nrhs, d, e, b, ldb);
}

inline integer pttrs( integer n, integer nrhs, const double* d,
                    const dcomplex* e, dcomplex* b, integer ldb )
{
    return zpttrs(n, nrhs, d, e, b, ldb);
}

inline integer sytrs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, const integer* ipiv,
                    float* b, integer ldb )
{
    return ssytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer sytrs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const integer* ipiv,
                    double* b, integer ldb )
{
    return dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer sytrs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const integer* ipiv,
                    scomplex* b, integer ldb )
{
    return csytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer sytrs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const integer* ipiv,
                    dcomplex* b, integer ldb )
{
    return zsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer hetrs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, const integer* ipiv,
                    float* b, integer ldb )
{
    return ssytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer hetrs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const integer* ipiv,
                    double* b, integer ldb )
{
    return dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer hetrs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const integer* ipiv,
                    scomplex* b, integer ldb )
{
    return chetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer hetrs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const integer* ipiv,
                    dcomplex* b, integer ldb )
{
    return zhetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer sptrs( char uplo, integer n, integer nrhs,
                    const float* ap, const integer* ipiv, float* b,
                    integer ldb )
{
    return ssptrs(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer sptrs( char uplo, integer n, integer nrhs,
                    const double* ap, const integer* ipiv, double* b,
                    integer ldb )
{
    return dsptrs(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer sptrs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const integer* ipiv, scomplex* b,
                    integer ldb )
{
    return csptrs(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer sptrs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const integer* ipiv, dcomplex* b,
                    integer ldb )
{
    return zsptrs(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer hptrs( char uplo, integer n, integer nrhs,
                    const float* ap, const integer* ipiv, float* b,
                    integer ldb )
{
    return ssptrs(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer hptrs( char uplo, integer n, integer nrhs,
                    const double* ap, const integer* ipiv, double* b,
                    integer ldb )
{
    return dsptrs(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer hptrs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const integer* ipiv, scomplex* b,
                    integer ldb )
{
    return chptrs(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer hptrs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const integer* ipiv, dcomplex* b,
                    integer ldb )
{
    return zhptrs(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer trtrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const float* a, integer lda,
                    float* b, integer ldb )
{
    return strtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb);
}

inline integer trtrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* a, integer lda,
                    double* b, integer ldb )
{
    return dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb);
}

inline integer trtrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const scomplex* a, integer lda,
                    scomplex* b, integer ldb )
{
    return ctrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb);
}

inline integer trtrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const dcomplex* a, integer lda,
                    dcomplex* b, integer ldb )
{
    return ztrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb);
}

inline integer tptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const float* ap, float* b,
                    integer ldb )
{
    return stptrs(uplo, trans, diag, n, nrhs, ap, b, ldb);
}

inline integer tptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* ap, double* b,
                    integer ldb )
{
    return dtptrs(uplo, trans, diag, n, nrhs, ap, b, ldb);
}

inline integer tptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const scomplex* ap, scomplex* b,
                    integer ldb )
{
    return ctptrs(uplo, trans, diag, n, nrhs, ap, b, ldb);
}

inline integer tptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const dcomplex* ap, dcomplex* b,
                    integer ldb )
{
    return ztptrs(uplo, trans, diag, n, nrhs, ap, b, ldb);
}

inline integer tbtrs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const float* ab,
                    integer ldab, float* b, integer ldb )
{
    return stbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer tbtrs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const double* ab,
                    integer ldab, double* b, integer ldb )
{
    return dtbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer tbtrs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const scomplex* ab,
                    integer ldab, scomplex* b, integer ldb )
{
    return ctbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer tbtrs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const dcomplex* ab,
                    integer ldab, dcomplex* b, integer ldb )
{
    return ztbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer gecon( char norm, integer n, const float* a, integer lda,
                    float anorm, float& rcond)
{
    return sgecon(norm, n, a, lda, anorm, &rcond);
}

inline integer gecon( char norm, integer n, const double* a, integer lda,
                    double anorm, double& rcond)
{
    return dgecon(norm, n, a, lda, anorm, &rcond);
}

inline integer gecon( char norm, integer n, const scomplex* a, integer lda,
                    float anorm, float& rcond)
{
    return cgecon(norm, n, a, lda, anorm, &rcond);
}

inline integer gecon( char norm, integer n, const dcomplex* a, integer lda,
                    double anorm, double& rcond)
{
    return zgecon(norm, n, a, lda, anorm, &rcond);
}

inline integer gbcon( char norm, integer n, integer kl, integer ku,
                    const float* ab, integer ldab, const integer* ipiv,
                    float anorm, float& rcond)
{
    return sgbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, &rcond);
}

inline integer gbcon( char norm, integer n, integer kl, integer ku,
                    const double* ab, integer ldab, const integer* ipiv,
                    double anorm, double& rcond)
{
    return dgbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, &rcond);
}

inline integer gbcon( char norm, integer n, integer kl, integer ku,
                    const scomplex* ab, integer ldab, const integer* ipiv,
                    float anorm, float& rcond)
{
    return cgbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, &rcond);
}

inline integer gbcon( char norm, integer n, integer kl, integer ku,
                    const dcomplex* ab, integer ldab, const integer* ipiv,
                    double anorm, double& rcond)
{
    return zgbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, &rcond);
}

inline integer gtcon( char norm, integer n, const float* dl,
                    const float* d, const float* du, const float* du2,
                    const integer* ipiv, float anorm, float& rcond)
{
    return sgtcon(norm, n, dl, d, du, du2, ipiv, anorm, &rcond);
}

inline integer gtcon( char norm, integer n, const double* dl,
                    const double* d, const double* du, const double* du2,
                    const integer* ipiv, double anorm, double& rcond)
{
    return dgtcon(norm, n, dl, d, du, du2, ipiv, anorm, &rcond);
}

inline integer gtcon( char norm, integer n, const scomplex* dl,
                    const scomplex* d, const scomplex* du, const scomplex* du2,
                    const integer* ipiv, float anorm, float& rcond)
{
    return cgtcon(norm, n, dl, d, du, du2, ipiv, anorm, &rcond);
}

inline integer gtcon( char norm, integer n, const dcomplex* dl,
                    const dcomplex* d, const dcomplex* du, const dcomplex* du2,
                    const integer* ipiv, double anorm, double& rcond)
{
    return zgtcon(norm, n, dl, d, du, du2, ipiv, anorm, &rcond);
}

inline integer pocon( char uplo, integer n, const float* a, integer lda,
                    float anorm, float& rcond)
{
    return spocon(uplo, n, a, lda, anorm, &rcond);
}

inline integer pocon( char uplo, integer n, const double* a, integer lda,
                    double anorm, double& rcond)
{
    return dpocon(uplo, n, a, lda, anorm, &rcond);
}

inline integer pocon( char uplo, integer n, const scomplex* a, integer lda,
                    float anorm, float& rcond)
{
    return cpocon(uplo, n, a, lda, anorm, &rcond);
}

inline integer pocon( char uplo, integer n, const dcomplex* a, integer lda,
                    double anorm, double& rcond)
{
    return zpocon(uplo, n, a, lda, anorm, &rcond);
}

inline integer ppcon( char uplo, integer n, const float* ap, float anorm, float& rcond)
{
    return sppcon(uplo, n, ap, anorm, &rcond);
}

inline integer ppcon( char uplo, integer n, const double* ap, double anorm, double& rcond)
{
    return dppcon(uplo, n, ap, anorm, &rcond);
}

inline integer ppcon( char uplo, integer n, const scomplex* ap, float anorm, float& rcond)
{
    return cppcon(uplo, n, ap, anorm, &rcond);
}

inline integer ppcon( char uplo, integer n, const dcomplex* ap, double anorm, double& rcond)
{
    return zppcon(uplo, n, ap, anorm, &rcond);
}

inline integer pbcon( char uplo, integer n, integer kd, const float* ab,
                    integer ldab, float anorm, float& rcond)
{
    return spbcon(uplo, n, kd, ab, ldab, anorm, &rcond);
}

inline integer pbcon( char uplo, integer n, integer kd, const double* ab,
                    integer ldab, double anorm, double& rcond)
{
    return dpbcon(uplo, n, kd, ab, ldab, anorm, &rcond);
}

inline integer pbcon( char uplo, integer n, integer kd, const scomplex* ab,
                    integer ldab, float anorm, float& rcond)
{
    return cpbcon(uplo, n, kd, ab, ldab, anorm, &rcond);
}

inline integer pbcon( char uplo, integer n, integer kd, const dcomplex* ab,
                    integer ldab, double anorm, double& rcond)
{
    return zpbcon(uplo, n, kd, ab, ldab, anorm, &rcond);
}

inline integer ptcon( integer n, const float* d, const float* e,
                    float anorm, float& rcond)
{
    return sptcon(n, d, e, anorm, &rcond);
}

inline integer ptcon( integer n, const double* d, const double* e,
                    double anorm, double& rcond)
{
    return dptcon(n, d, e, anorm, &rcond);
}

inline integer ptcon( integer n, const float* d, const scomplex* e,
                    float anorm, float& rcond)
{
    return cptcon(n, d, e, anorm, &rcond);
}

inline integer ptcon( integer n, const double* d, const dcomplex* e,
                    double anorm, double& rcond)
{
    return zptcon(n, d, e, anorm, &rcond);
}

inline integer sycon( char uplo, integer n, const float* a, integer lda,
                    const integer* ipiv, float anorm, float& rcond)
{
    return ssycon(uplo, n, a, lda, ipiv, anorm, &rcond);
}

inline integer sycon( char uplo, integer n, const double* a, integer lda,
                    const integer* ipiv, double anorm, double& rcond)
{
    return dsycon(uplo, n, a, lda, ipiv, anorm, &rcond);
}

inline integer sycon( char uplo, integer n, const scomplex* a, integer lda,
                    const integer* ipiv, float anorm, float& rcond)
{
    return csycon(uplo, n, a, lda, ipiv, anorm, &rcond);
}

inline integer sycon( char uplo, integer n, const dcomplex* a, integer lda,
                    const integer* ipiv, double anorm, double& rcond)
{
    return zsycon(uplo, n, a, lda, ipiv, anorm, &rcond);
}

inline integer hecon( char uplo, integer n, const float* a, integer lda,
                    const integer* ipiv, float anorm, float& rcond)
{
    return ssycon(uplo, n, a, lda, ipiv, anorm, &rcond);
}

inline integer hecon( char uplo, integer n, const double* a, integer lda,
                    const integer* ipiv, double anorm, double& rcond)
{
    return dsycon(uplo, n, a, lda, ipiv, anorm, &rcond);
}

inline integer hecon( char uplo, integer n, const scomplex* a, integer lda,
                    const integer* ipiv, float anorm, float& rcond)
{
    return checon(uplo, n, a, lda, ipiv, anorm, &rcond);
}

inline integer hecon( char uplo, integer n, const dcomplex* a, integer lda,
                    const integer* ipiv, double anorm, double& rcond)
{
    return zhecon(uplo, n, a, lda, ipiv, anorm, &rcond);
}

inline integer spcon( char uplo, integer n, const float* ap,
                    const integer* ipiv, float anorm, float& rcond)
{
    return sspcon(uplo, n, ap, ipiv, anorm, &rcond);
}

inline integer spcon( char uplo, integer n, const double* ap,
                    const integer* ipiv, double anorm, double& rcond)
{
    return dspcon(uplo, n, ap, ipiv, anorm, &rcond);
}

inline integer spcon( char uplo, integer n, const scomplex* ap,
                    const integer* ipiv, float anorm, float& rcond)
{
    return cspcon(uplo, n, ap, ipiv, anorm, &rcond);
}

inline integer spcon( char uplo, integer n, const dcomplex* ap,
                    const integer* ipiv, double anorm, double& rcond)
{
    return zspcon(uplo, n, ap, ipiv, anorm, &rcond);
}

inline integer hpcon( char uplo, integer n, const float* ap,
                    const integer* ipiv, float anorm, float& rcond)
{
    return sspcon(uplo, n, ap, ipiv, anorm, &rcond);
}

inline integer hpcon( char uplo, integer n, const double* ap,
                    const integer* ipiv, double anorm, double& rcond)
{
    return dspcon(uplo, n, ap, ipiv, anorm, &rcond);
}

inline integer hpcon( char uplo, integer n, const scomplex* ap,
                    const integer* ipiv, float anorm, float& rcond)
{
    return chpcon(uplo, n, ap, ipiv, anorm, &rcond);
}

inline integer hpcon( char uplo, integer n, const dcomplex* ap,
                    const integer* ipiv, double anorm, double& rcond)
{
    return zhpcon(uplo, n, ap, ipiv, anorm, &rcond);
}

inline integer trcon( char norm, char uplo, char diag, integer n,
                    const float* a, integer lda, float& rcond)
{
    return strcon(norm, uplo, diag, n, a, lda, &rcond);
}

inline integer trcon( char norm, char uplo, char diag, integer n,
                    const double* a, integer lda, double& rcond)
{
    return dtrcon(norm, uplo, diag, n, a, lda, &rcond);
}

inline integer trcon( char norm, char uplo, char diag, integer n,
                    const scomplex* a, integer lda, float& rcond)
{
    return ctrcon(norm, uplo, diag, n, a, lda, &rcond);
}

inline integer trcon( char norm, char uplo, char diag, integer n,
                    const dcomplex* a, integer lda, double& rcond)
{
    return ztrcon(norm, uplo, diag, n, a, lda, &rcond);
}

inline integer tpcon( char norm, char uplo, char diag, integer n,
                    const float* ap, float& rcond)
{
    return stpcon(norm, uplo, diag, n, ap, &rcond);
}

inline integer tpcon( char norm, char uplo, char diag, integer n,
                    const double* ap, double& rcond)
{
    return dtpcon(norm, uplo, diag, n, ap, &rcond);
}

inline integer tpcon( char norm, char uplo, char diag, integer n,
                    const scomplex* ap, float& rcond)
{
    return ctpcon(norm, uplo, diag, n, ap, &rcond);
}

inline integer tpcon( char norm, char uplo, char diag, integer n,
                    const dcomplex* ap, double& rcond)
{
    return ztpcon(norm, uplo, diag, n, ap, &rcond);
}

inline integer tbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const float* ab, integer ldab, float& rcond)
{
    return stbcon(norm, uplo, diag, n, kd, ab, ldab, &rcond);
}

inline integer tbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const double* ab, integer ldab, double& rcond)
{
    return dtbcon(norm, uplo, diag, n, kd, ab, ldab, &rcond);
}

inline integer tbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const scomplex* ab, integer ldab, float& rcond)
{
    return ctbcon(norm, uplo, diag, n, kd, ab, ldab, &rcond);
}

inline integer tbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const dcomplex* ab, integer ldab, double& rcond)
{
    return ztbcon(norm, uplo, diag, n, kd, ab, ldab, &rcond);
}

inline integer gerfs( char trans, integer n, integer nrhs,
                    const float* a, integer lda, const float* af,
                    integer ldaf, const integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float* ferr,
                    float* berr)
{
    return sgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gerfs( char trans, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr)
{
    return dgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gerfs( char trans, integer n, integer nrhs,
                    const scomplex* a, integer lda, const scomplex* af,
                    integer ldaf, const integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    return cgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gerfs( char trans, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const dcomplex* af,
                    integer ldaf, const integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    return zgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gerfsx( char trans, char equed, integer n, integer nrhs,
                     const float* a, integer lda, const float* af,
                     integer ldaf, const integer* ipiv, const float* r,
                     const float* c, const float* b, integer ldb,
                     float* x, integer ldx, float& rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return sgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gerfsx( char trans, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const integer* ipiv, const double* r,
                     const double* c, const double* b, integer ldb,
                     double* x, integer ldx, double& rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return dgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gerfsx( char trans, char equed, integer n, integer nrhs,
                     const scomplex* a, integer lda, const scomplex* af,
                     integer ldaf, const integer* ipiv, const float* r,
                     const float* c, const scomplex* b, integer ldb,
                     scomplex* x, integer ldx, float& rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return cgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gerfsx( char trans, char equed, integer n, integer nrhs,
                     const dcomplex* a, integer lda, const dcomplex* af,
                     integer ldaf, const integer* ipiv, const double* r,
                     const double* c, const dcomplex* b, integer ldb,
                     dcomplex* x, integer ldx, double& rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return zgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gbrfs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const float* ab, integer ldab,
                    const float* afb, integer ldafb,
                    const integer* ipiv, const float* b, integer ldb,
                    float* x, integer ldx, float* ferr, float* berr)
{
    return sgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gbrfs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const double* ab, integer ldab,
                    const double* afb, integer ldafb,
                    const integer* ipiv, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr)
{
    return dgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gbrfs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const scomplex* ab, integer ldab,
                    const scomplex* afb, integer ldafb,
                    const integer* ipiv, const scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float* ferr, float* berr)
{
    return cgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gbrfs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const dcomplex* ab, integer ldab,
                    const dcomplex* afb, integer ldafb,
                    const integer* ipiv, const dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double* ferr, double* berr)
{
    return zgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gbrfsx( char trans, char equed, integer n, integer kl,
                     integer ku, integer nrhs, const float* ab,
                     integer ldab, const float* afb, integer ldafb,
                     const integer* ipiv, const float* r, const float* c,
                     const float* b, integer ldb, float* x,
                     integer ldx, float& rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return sgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gbrfsx( char trans, char equed, integer n, integer kl,
                     integer ku, integer nrhs, const double* ab,
                     integer ldab, const double* afb, integer ldafb,
                     const integer* ipiv, const double* r, const double* c,
                     const double* b, integer ldb, double* x,
                     integer ldx, double& rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return dgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gbrfsx( char trans, char equed, integer n, integer kl,
                     integer ku, integer nrhs, const scomplex* ab,
                     integer ldab, const scomplex* afb, integer ldafb,
                     const integer* ipiv, const float* r, const float* c,
                     const scomplex* b, integer ldb, scomplex* x,
                     integer ldx, float& rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return cgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gbrfsx( char trans, char equed, integer n, integer kl,
                     integer ku, integer nrhs, const dcomplex* ab,
                     integer ldab, const dcomplex* afb, integer ldafb,
                     const integer* ipiv, const double* r, const double* c,
                     const dcomplex* b, integer ldb, dcomplex* x,
                     integer ldx, double& rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return zgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gtrfs( char trans, integer n, integer nrhs,
                    const float* dl, const float* d, const float* du,
                    const float* dlf, const float* df, const float* duf,
                    const float* du2, const integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float* ferr,
                    float* berr)
{
    return sgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gtrfs( char trans, integer n, integer nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* dlf, const double* df, const double* duf,
                    const double* du2, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr)
{
    return dgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gtrfs( char trans, integer n, integer nrhs,
                    const scomplex* dl, const scomplex* d, const scomplex* du,
                    const scomplex* dlf, const scomplex* df, const scomplex* duf,
                    const scomplex* du2, const integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    return cgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer gtrfs( char trans, integer n, integer nrhs,
                    const dcomplex* dl, const dcomplex* d, const dcomplex* du,
                    const dcomplex* dlf, const dcomplex* df, const dcomplex* duf,
                    const dcomplex* du2, const integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    return zgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer porfs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, const float* af,
                    integer ldaf, const float* b, integer ldb,
                    float* x, integer ldx, float* ferr, float* berr)
{
    return sporfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr);
}

inline integer porfs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr)
{
    return dporfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr);
}

inline integer porfs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const scomplex* af,
                    integer ldaf, const scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float* ferr, float* berr)
{
    return cporfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr);
}

inline integer porfs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const dcomplex* af,
                    integer ldaf, const dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double* ferr, double* berr)
{
    return zporfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr);
}

inline integer porfsx( char uplo, char equed, integer n, integer nrhs,
                     const float* a, integer lda, const float* af,
                     integer ldaf, const float* s, const float* b,
                     integer ldb, float* x, integer ldx, float& rcond,
                     float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    return sporfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldx, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer porfsx( char uplo, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const double* s, const double* b,
                     integer ldb, double* x, integer ldx, double& rcond,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    return dporfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldx, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer porfsx( char uplo, char equed, integer n, integer nrhs,
                     const scomplex* a, integer lda, const scomplex* af,
                     integer ldaf, const float* s, const scomplex* b,
                     integer ldb, scomplex* x, integer ldx, float& rcond,
                     float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    return cporfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldx, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer porfsx( char uplo, char equed, integer n, integer nrhs,
                     const dcomplex* a, integer lda, const dcomplex* af,
                     integer ldaf, const double* s, const dcomplex* b,
                     integer ldb, dcomplex* x, integer ldx, double& rcond,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    return zporfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldx, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer pprfs( char uplo, integer n, integer nrhs,
                    const float* ap, const float* afp, const float* b,
                    integer ldb, float* x, integer ldx, float* ferr,
                    float* berr)
{
    return spprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr);
}

inline integer pprfs( char uplo, integer n, integer nrhs,
                    const double* ap, const double* afp, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr)
{
    return dpprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr);
}

inline integer pprfs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const scomplex* afp, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    return cpprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr);
}

inline integer pprfs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const dcomplex* afp, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    return zpprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr);
}

inline integer pbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const float* ab, integer ldab, const float* afb,
                    integer ldafb, const float* b, integer ldb,
                    float* x, integer ldx, float* ferr, float* berr)
{
    return spbrfs(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr);
}

inline integer pbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const double* ab, integer ldab, const double* afb,
                    integer ldafb, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr)
{
    return dpbrfs(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr);
}

inline integer pbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const scomplex* ab, integer ldab, const scomplex* afb,
                    integer ldafb, const scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float* ferr, float* berr)
{
    return cpbrfs(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr);
}

inline integer pbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const dcomplex* ab, integer ldab, const dcomplex* afb,
                    integer ldafb, const dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double* ferr, double* berr)
{
    return zpbrfs(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr);
}

inline integer ptrfs( integer n, integer nrhs, const float* d,
                    const float* e, const float* df, const float* ef,
                    const float* b, integer ldb, float* x,
                    integer ldx, float* ferr, float* berr)
{
    return sptrfs(n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr);
}

inline integer ptrfs( integer n, integer nrhs, const double* d,
                    const double* e, const double* df, const double* ef,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* ferr, double* berr)
{
    return dptrfs(n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr);
}

inline integer ptrfs( integer n, integer nrhs, const float* d,
                    const scomplex* e, const float* df, const scomplex* ef,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float* ferr, float* berr)
{
    return cptrfs(n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr);
}

inline integer ptrfs( integer n, integer nrhs, const double* d,
                    const dcomplex* e, const double* df, const dcomplex* ef,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double* ferr, double* berr)
{
    return zptrfs(n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr);
}

inline integer syrfs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, const float* af,
                    integer ldaf, const integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float* ferr,
                    float* berr)
{
    return ssyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer syrfs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr)
{
    return dsyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer syrfs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const scomplex* af,
                    integer ldaf, const integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    return csyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer syrfs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const dcomplex* af,
                    integer ldaf, const integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    return zsyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer herfs( char uplo, integer n, integer nrhs,
                    const float* a, integer lda, const float* af,
                    integer ldaf, const integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float* ferr,
                    float* berr)
{
    return ssyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer herfs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr)
{
    return dsyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer herfs( char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, const scomplex* af,
                    integer ldaf, const integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    return cherfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer herfs( char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, const dcomplex* af,
                    integer ldaf, const integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    return zherfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer syrfsx( char uplo, char equed, integer n, integer nrhs,
                     const float* a, integer lda, const float* af,
                     integer ldaf, const integer* ipiv, const float* s,
                     const float* b, integer ldb, float* x,
                     integer ldx, float& rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return ssyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer syrfsx( char uplo, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const integer* ipiv, const double* s,
                     const double* b, integer ldb, double* x,
                     integer ldx, double& rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return dsyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer syrfsx( char uplo, char equed, integer n, integer nrhs,
                     const scomplex* a, integer lda, const scomplex* af,
                     integer ldaf, const integer* ipiv, const float* s,
                     const scomplex* b, integer ldb, scomplex* x,
                     integer ldx, float& rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return csyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer syrfsx( char uplo, char equed, integer n, integer nrhs,
                     const dcomplex* a, integer lda, const dcomplex* af,
                     integer ldaf, const integer* ipiv, const double* s,
                     const dcomplex* b, integer ldb, dcomplex* x,
                     integer ldx, double& rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return zsyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer herfsx( char uplo, char equed, integer n, integer nrhs,
                     const float* a, integer lda, const float* af,
                     integer ldaf, const integer* ipiv, const float* s,
                     const float* b, integer ldb, float* x,
                     integer ldx, float& rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return ssyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer herfsx( char uplo, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const integer* ipiv, const double* s,
                     const double* b, integer ldb, double* x,
                     integer ldx, double& rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return dsyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer herfsx( char uplo, char equed, integer n, integer nrhs,
                     const scomplex* a, integer lda, const scomplex* af,
                     integer ldaf, const integer* ipiv, const float* s,
                     const scomplex* b, integer ldb, scomplex* x,
                     integer ldx, float& rcond, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return cherfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer herfsx( char uplo, char equed, integer n, integer nrhs,
                     const dcomplex* a, integer lda, const dcomplex* af,
                     integer ldaf, const integer* ipiv, const double* s,
                     const dcomplex* b, integer ldb, dcomplex* x,
                     integer ldx, double& rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return zherfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx,
                   &rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer sprfs( char uplo, integer n, integer nrhs,
                    const float* ap, const float* afp, const integer* ipiv,
                    const float* b, integer ldb, float* x,
                    integer ldx, float* ferr, float* berr)
{
    return ssprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer sprfs( char uplo, integer n, integer nrhs,
                    const double* ap, const double* afp, const integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* ferr, double* berr)
{
    return dsprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer sprfs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const scomplex* afp, const integer* ipiv,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float* ferr, float* berr)
{
    return csprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer sprfs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const dcomplex* afp, const integer* ipiv,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double* ferr, double* berr)
{
    return zsprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer hprfs( char uplo, integer n, integer nrhs,
                    const float* ap, const float* afp, const integer* ipiv,
                    const float* b, integer ldb, float* x,
                    integer ldx, float* ferr, float* berr)
{
    return ssprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer hprfs( char uplo, integer n, integer nrhs,
                    const double* ap, const double* afp, const integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* ferr, double* berr)
{
    return dsprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer hprfs( char uplo, integer n, integer nrhs,
                    const scomplex* ap, const scomplex* afp, const integer* ipiv,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float* ferr, float* berr)
{
    return chprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer hprfs( char uplo, integer n, integer nrhs,
                    const dcomplex* ap, const dcomplex* afp, const integer* ipiv,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double* ferr, double* berr)
{
    return zhprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr);
}

inline integer trrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const float* a, integer lda,
                    const float* b, integer ldb, const float* x,
                    integer ldx, float* ferr, float* berr)
{
    return strrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr);
}

inline integer trrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* a, integer lda,
                    const double* b, integer ldb, const double* x,
                    integer ldx, double* ferr, double* berr)
{
    return dtrrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr);
}

inline integer trrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const scomplex* a, integer lda,
                    const scomplex* b, integer ldb, const scomplex* x,
                    integer ldx, float* ferr, float* berr)
{
    return ctrrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr);
}

inline integer trrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const dcomplex* a, integer lda,
                    const dcomplex* b, integer ldb, const dcomplex* x,
                    integer ldx, double* ferr, double* berr)
{
    return ztrrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr);
}

inline integer tprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const float* ap, const float* b,
                    integer ldb, const float* x, integer ldx,
                    float* ferr, float* berr)
{
    return stprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr);
}

inline integer tprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* ap, const double* b,
                    integer ldb, const double* x, integer ldx,
                    double* ferr, double* berr)
{
    return dtprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr);
}

inline integer tprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const scomplex* ap, const scomplex* b,
                    integer ldb, const scomplex* x, integer ldx,
                    float* ferr, float* berr)
{
    return ctprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr);
}

inline integer tprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const dcomplex* ap, const dcomplex* b,
                    integer ldb, const dcomplex* x, integer ldx,
                    double* ferr, double* berr)
{
    return ztprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr);
}

inline integer tbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const float* ab,
                    integer ldab, const float* b, integer ldb,
                    const float* x, integer ldx, float* ferr,
                    float* berr)
{
    return stbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr);
}

inline integer tbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const double* ab,
                    integer ldab, const double* b, integer ldb,
                    const double* x, integer ldx, double* ferr,
                    double* berr)
{
    return dtbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr);
}

inline integer tbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const scomplex* ab,
                    integer ldab, const scomplex* b, integer ldb,
                    const scomplex* x, integer ldx, float* ferr,
                    float* berr)
{
    return ctbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr);
}

inline integer tbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const dcomplex* ab,
                    integer ldab, const dcomplex* b, integer ldb,
                    const dcomplex* x, integer ldx, double* ferr,
                    double* berr)
{
    return ztbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr);
}

inline integer getri( integer n, float* a, integer lda,
                    const integer* ipiv)
{
    return sgetri(n, a, lda, ipiv);
}

inline integer getri( integer n, double* a, integer lda,
                    const integer* ipiv)
{
    return dgetri(n, a, lda, ipiv);
}

inline integer getri( integer n, scomplex* a, integer lda,
                    const integer* ipiv)
{
    return cgetri(n, a, lda, ipiv);
}

inline integer getri( integer n, dcomplex* a, integer lda,
                    const integer* ipiv)
{
    return zgetri(n, a, lda, ipiv);
}

inline integer potri( char uplo, integer n, float* a, integer lda )
{
    return spotri(uplo, n, a, lda);
}

inline integer potri( char uplo, integer n, double* a, integer lda )
{
    return dpotri(uplo, n, a, lda);
}

inline integer potri( char uplo, integer n, scomplex* a, integer lda )
{
    return cpotri(uplo, n, a, lda);
}

inline integer potri( char uplo, integer n, dcomplex* a, integer lda )
{
    return zpotri(uplo, n, a, lda);
}

inline integer pftri( char transr, char uplo, integer n, float* a )
{
    return spftri(transr, uplo, n, a);
}

inline integer pftri( char transr, char uplo, integer n, double* a )
{
    return dpftri(transr, uplo, n, a);
}

inline integer pftri( char transr, char uplo, integer n, scomplex* a )
{
    return cpftri(transr, uplo, n, a);
}

inline integer pftri( char transr, char uplo, integer n, dcomplex* a )
{
    return zpftri(transr, uplo, n, a);
}

inline integer pptri( char uplo, integer n, float* ap )
{
    return spptri(uplo, n, ap);
}

inline integer pptri( char uplo, integer n, double* ap )
{
    return dpptri(uplo, n, ap);
}

inline integer pptri( char uplo, integer n, scomplex* ap )
{
    return cpptri(uplo, n, ap);
}

inline integer pptri( char uplo, integer n, dcomplex* ap )
{
    return zpptri(uplo, n, ap);
}

inline integer sytri( char uplo, integer n, float* a, integer lda,
                    const integer* ipiv)
{
    return ssytri(uplo, n, a, lda, ipiv);
}

inline integer sytri( char uplo, integer n, double* a, integer lda,
                    const integer* ipiv)
{
    return dsytri(uplo, n, a, lda, ipiv);
}

inline integer sytri( char uplo, integer n, scomplex* a, integer lda,
                    const integer* ipiv)
{
    return csytri(uplo, n, a, lda, ipiv);
}

inline integer sytri( char uplo, integer n, dcomplex* a, integer lda,
                    const integer* ipiv)
{
    return zsytri(uplo, n, a, lda, ipiv);
}

inline integer hetri( char uplo, integer n, float* a, integer lda,
                    const integer* ipiv)
{
    return ssytri(uplo, n, a, lda, ipiv);
}

inline integer hetri( char uplo, integer n, double* a, integer lda,
                    const integer* ipiv)
{
    return dsytri(uplo, n, a, lda, ipiv);
}

inline integer hetri( char uplo, integer n, scomplex* a, integer lda,
                    const integer* ipiv)
{
    return chetri(uplo, n, a, lda, ipiv);
}

inline integer hetri( char uplo, integer n, dcomplex* a, integer lda,
                    const integer* ipiv)
{
    return zhetri(uplo, n, a, lda, ipiv);
}

inline integer sptri( char uplo, integer n, float* ap,
                    const integer* ipiv)
{
    return ssptri(uplo, n, ap, ipiv);
}

inline integer sptri( char uplo, integer n, double* ap,
                    const integer* ipiv)
{
    return dsptri(uplo, n, ap, ipiv);
}

inline integer sptri( char uplo, integer n, scomplex* ap,
                    const integer* ipiv)
{
    return csptri(uplo, n, ap, ipiv);
}

inline integer sptri( char uplo, integer n, dcomplex* ap,
                    const integer* ipiv)
{
    return zsptri(uplo, n, ap, ipiv);
}

inline integer hptri( char uplo, integer n, float* ap,
                    const integer* ipiv)
{
    return ssptri(uplo, n, ap, ipiv);
}

inline integer hptri( char uplo, integer n, double* ap,
                    const integer* ipiv)
{
    return dsptri(uplo, n, ap, ipiv);
}

inline integer hptri( char uplo, integer n, scomplex* ap,
                    const integer* ipiv)
{
    return chptri(uplo, n, ap, ipiv);
}

inline integer hptri( char uplo, integer n, dcomplex* ap,
                    const integer* ipiv)
{
    return zhptri(uplo, n, ap, ipiv);
}

inline integer trtri( char uplo, char diag, integer n, float* a,
                    integer lda )
{
    return strtri(uplo, diag, n, a, lda);
}

inline integer trtri( char uplo, char diag, integer n, double* a,
                    integer lda )
{
    return dtrtri(uplo, diag, n, a, lda);
}

inline integer trtri( char uplo, char diag, integer n, scomplex* a,
                    integer lda )
{
    return ctrtri(uplo, diag, n, a, lda);
}

inline integer trtri( char uplo, char diag, integer n, dcomplex* a,
                    integer lda )
{
    return ztrtri(uplo, diag, n, a, lda);
}

inline integer tftri( char transr, char uplo, char diag, integer n,
                    float* a )
{
    return stftri(transr, uplo, diag, n, a);
}

inline integer tftri( char transr, char uplo, char diag, integer n,
                    double* a )
{
    return dtftri(transr, uplo, diag, n, a);
}

inline integer tftri( char transr, char uplo, char diag, integer n,
                    scomplex* a )
{
    return ctftri(transr, uplo, diag, n, a);
}

inline integer tftri( char transr, char uplo, char diag, integer n,
                    dcomplex* a )
{
    return ztftri(transr, uplo, diag, n, a);
}

inline integer tptri( char uplo, char diag, integer n, float* ap )
{
    return stptri(uplo, diag, n, ap);
}

inline integer tptri( char uplo, char diag, integer n, double* ap )
{
    return dtptri(uplo, diag, n, ap);
}

inline integer tptri( char uplo, char diag, integer n, scomplex* ap )
{
    return ctptri(uplo, diag, n, ap);
}

inline integer tptri( char uplo, char diag, integer n, dcomplex* ap )
{
    return ztptri(uplo, diag, n, ap);
}

inline integer geequ( integer m, integer n, const float* a,
                    integer lda, float* r, float* c, float& rowcnd,
                    float& colcnd, float& amax )
{
    return sgeequ(m, n, a, lda, r, c, &rowcnd, &colcnd, &amax);
}

inline integer geequ( integer m, integer n, const double* a,
                    integer lda, double* r, double* c, double& rowcnd,
                    double& colcnd, double& amax )
{
    return dgeequ(m, n, a, lda, r, c, &rowcnd, &colcnd, &amax);
}

inline integer geequ( integer m, integer n, const scomplex* a,
                    integer lda, float* r, float* c, float& rowcnd,
                    float& colcnd, float& amax )
{
    return cgeequ(m, n, a, lda, r, c, &rowcnd, &colcnd, &amax);
}

inline integer geequ( integer m, integer n, const dcomplex* a,
                    integer lda, double* r, double* c, double& rowcnd,
                    double& colcnd, double& amax )
{
    return zgeequ(m, n, a, lda, r, c, &rowcnd, &colcnd, &amax);
}

inline integer geequb( integer m, integer n, const float* a,
                     integer lda, float* r, float* c, float& rowcnd,
                     float& colcnd, float& amax )
{
    return sgeequb(m, n, a, lda, r, c, &rowcnd, &colcnd, &amax);
}

inline integer geequb( integer m, integer n, const double* a,
                     integer lda, double* r, double* c, double& rowcnd,
                     double& colcnd, double& amax )
{
    return dgeequb(m, n, a, lda, r, c, &rowcnd, &colcnd, &amax);
}

inline integer geequb( integer m, integer n, const scomplex* a,
                     integer lda, float* r, float* c, float& rowcnd,
                     float& colcnd, float& amax )
{
    return cgeequb(m, n, a, lda, r, c, &rowcnd, &colcnd, &amax);
}

inline integer geequb( integer m, integer n, const dcomplex* a,
                     integer lda, double* r, double* c, double& rowcnd,
                     double& colcnd, double& amax )
{
    return zgeequb(m, n, a, lda, r, c, &rowcnd, &colcnd, &amax);
}

inline integer gbequ( integer m, integer n, integer kl,
                    integer ku, const float* ab, integer ldab,
                    float* r, float* c, float& rowcnd, float& colcnd,
                    float& amax )
{
    return sgbequ(m, n, kl, ku, ab, ldab, r, c, &rowcnd, &colcnd, &amax);
}

inline integer gbequ( integer m, integer n, integer kl,
                    integer ku, const double* ab, integer ldab,
                    double* r, double* c, double& rowcnd, double& colcnd,
                    double& amax )
{
    return dgbequ(m, n, kl, ku, ab, ldab, r, c, &rowcnd, &colcnd, &amax);
}

inline integer gbequ( integer m, integer n, integer kl,
                    integer ku, const scomplex* ab, integer ldab,
                    float* r, float* c, float& rowcnd, float& colcnd,
                    float& amax )
{
    return cgbequ(m, n, kl, ku, ab, ldab, r, c, &rowcnd, &colcnd, &amax);
}

inline integer gbequ( integer m, integer n, integer kl,
                    integer ku, const dcomplex* ab, integer ldab,
                    double* r, double* c, double& rowcnd, double& colcnd,
                    double& amax )
{
    return zgbequ(m, n, kl, ku, ab, ldab, r, c, &rowcnd, &colcnd, &amax);
}

inline integer gbequb( integer m, integer n, integer kl,
                     integer ku, const float* ab, integer ldab,
                     float* r, float* c, float& rowcnd, float& colcnd,
                     float& amax )
{
    return sgbequb(m, n, kl, ku, ab, ldab, r, c, &rowcnd, &colcnd, &amax);
}

inline integer gbequb( integer m, integer n, integer kl,
                     integer ku, const double* ab, integer ldab,
                     double* r, double* c, double& rowcnd, double& colcnd,
                     double& amax )
{
    return dgbequb(m, n, kl, ku, ab, ldab, r, c, &rowcnd, &colcnd, &amax);
}

inline integer gbequb( integer m, integer n, integer kl,
                     integer ku, const scomplex* ab, integer ldab,
                     float* r, float* c, float& rowcnd, float& colcnd,
                     float& amax )
{
    return cgbequb(m, n, kl, ku, ab, ldab, r, c, &rowcnd, &colcnd, &amax);
}

inline integer gbequb( integer m, integer n, integer kl,
                     integer ku, const dcomplex* ab, integer ldab,
                     double* r, double* c, double& rowcnd, double& colcnd,
                     double& amax )
{
    return zgbequb(m, n, kl, ku, ab, ldab, r, c, &rowcnd, &colcnd, &amax);
}

inline integer poequ( integer n, const float* a, integer lda, float* s,
                    float& scond, float& amax )
{
    return spoequ(n, a, lda, s, &scond, &amax);
}

inline integer poequ( integer n, const double* a, integer lda, double* s,
                    double& scond, double& amax )
{
    return dpoequ(n, a, lda, s, &scond, &amax);
}

inline integer poequ( integer n, const scomplex* a, integer lda, float* s,
                    float& scond, float& amax )
{
    return cpoequ(n, a, lda, s, &scond, &amax);
}

inline integer poequ( integer n, const dcomplex* a, integer lda, double* s,
                    double& scond, double& amax )
{
    return zpoequ(n, a, lda, s, &scond, &amax);
}

inline integer poequb( integer n, const float* a, integer lda, float* s,
                     float& scond, float& amax )
{
    return spoequb(n, a, lda, s, &scond, &amax);
}

inline integer poequb( integer n, const double* a, integer lda, double* s,
                     double& scond, double& amax )
{
    return dpoequb(n, a, lda, s, &scond, &amax);
}

inline integer poequb( integer n, const scomplex* a, integer lda, float* s,
                    float& scond, float& amax )
{
    return cpoequb(n, a, lda, s, &scond, &amax);
}

inline integer poequb( integer n, const dcomplex* a, integer lda, double* s,
                     double& scond, double& amax )
{
    return zpoequb(n, a, lda, s, &scond, &amax);
}

inline integer ppequ( char uplo, integer n, const float* ap, float* s,
                    float& scond, float& amax )
{
    return sppequ(uplo, n, ap, s, &scond, &amax);
}

inline integer ppequ( char uplo, integer n, const double* ap, double* s,
                    double& scond, double& amax )
{
    return dppequ(uplo, n, ap, s, &scond, &amax);
}

inline integer ppequ( char uplo, integer n, const scomplex* ap, float* s,
                    float& scond, float& amax )
{
    return cppequ(uplo, n, ap, s, &scond, &amax);
}

inline integer ppequ( char uplo, integer n, const dcomplex* ap, double* s,
                    double& scond, double& amax )
{
    return zppequ(uplo, n, ap, s, &scond, &amax);
}

inline integer pbequ( char uplo, integer n, integer kd, const float* ab,
                    integer ldab, float* s, float& scond, float& amax )
{
    return spbequ(uplo, n, kd, ab, ldab, s, &scond, &amax);
}

inline integer pbequ( char uplo, integer n, integer kd, const double* ab,
                    integer ldab, double* s, double& scond, double& amax )
{
    return dpbequ(uplo, n, kd, ab, ldab, s, &scond, &amax);
}

inline integer pbequ( char uplo, integer n, integer kd, const scomplex* ab,
                    integer ldab, float* s, float& scond, float& amax )
{
    return cpbequ(uplo, n, kd, ab, ldab, s, &scond, &amax);
}

inline integer pbequ( char uplo, integer n, integer kd, const dcomplex* ab,
                    integer ldab, double* s, double& scond, double& amax )
{
    return zpbequ(uplo, n, kd, ab, ldab, s, &scond, &amax);
}

inline integer syequb( char uplo, integer n, const float* a,
                     integer lda, float* s, float& scond, float& amax)
{
    return ssyequb(uplo, n, a, lda, s, &scond, &amax);
}

inline integer syequb( char uplo, integer n, const double* a,
                     integer lda, double* s, double& scond, double& amax)
{
    return dsyequb(uplo, n, a, lda, s, &scond, &amax);
}

inline integer syequb( char uplo, integer n, const scomplex* a,
                     integer lda, float* s, float& scond, float& amax)
{
    return csyequb(uplo, n, a, lda, s, &scond, &amax);
}

inline integer syequb( char uplo, integer n, const dcomplex* a,
                     integer lda, double* s, double& scond, double& amax)
{
    return zsyequb(uplo, n, a, lda, s, &scond, &amax);
}

inline integer heequb( char uplo, integer n, const float* a,
                     integer lda, float* s, float& scond, float& amax)
{
    return ssyequb(uplo, n, a, lda, s, &scond, &amax);
}

inline integer heequb( char uplo, integer n, const double* a,
                     integer lda, double* s, double& scond, double& amax)
{
    return dsyequb(uplo, n, a, lda, s, &scond, &amax);
}

inline integer heequb( char uplo, integer n, const scomplex* a,
                     integer lda, float* s, float& scond, float& amax)
{
    return cheequb(uplo, n, a, lda, s, &scond, &amax);
}

inline integer heequb( char uplo, integer n, const dcomplex* a,
                     integer lda, double* s, double& scond, double& amax)
{
    return zheequb(uplo, n, a, lda, s, &scond, &amax);
}

inline integer gesv( integer n, integer nrhs, float* a, integer lda,
                   integer* ipiv, float* b, integer ldb )
{
    return sgesv(n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer gesv( integer n, integer nrhs, double* a, integer lda,
                   integer* ipiv, double* b, integer ldb )
{
    return dgesv(n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer gesv( integer n, integer nrhs, scomplex* a, integer lda,
                   integer* ipiv, scomplex* b, integer ldb )
{
    return cgesv(n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer gesv( integer n, integer nrhs, dcomplex* a, integer lda,
                   integer* ipiv, dcomplex* b, integer ldb )
{
    return zgesv(n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer gesv( integer n, integer nrhs, double* a, integer lda,
                    integer* ipiv, double* b, integer ldb, double* x,
                    integer ldx, integer& iter )
{
    return dsgesv(n, nrhs, a, lda, ipiv, b, ldb, x, ldx, &iter);
}

inline integer gesv( integer n, integer nrhs, dcomplex* a, integer lda,
                    integer* ipiv, dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, integer& iter )
{
    return zcgesv(n, nrhs, a, lda, ipiv, b, ldb, x, ldx, &iter);
}

inline integer gesvx( char fact, char trans, integer n, integer nrhs,
                    float* a, integer lda, float* af, integer ldaf,
                    integer* ipiv, char* equed, float* r, float* c,
                    float* b, integer ldb, float* x, integer ldx,
                    float& rcond, float* ferr, float* berr)
{
    return sgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gesvx( char fact, char trans, integer n, integer nrhs,
                    double* a, integer lda, double* af, integer ldaf,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer ldb, double* x, integer ldx,
                    double& rcond, double* ferr, double* berr)
{
    return dgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gesvx( char fact, char trans, integer n, integer nrhs,
                    scomplex* a, integer lda, scomplex* af, integer ldaf,
                    integer* ipiv, char* equed, float* r, float* c,
                    scomplex* b, integer ldb, scomplex* x, integer ldx,
                    float& rcond, float* ferr, float* berr)
{
    return cgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gesvx( char fact, char trans, integer n, integer nrhs,
                    dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                    integer* ipiv, char* equed, double* r, double* c,
                    dcomplex* b, integer ldb, dcomplex* x, integer ldx,
                    double& rcond, double* ferr, double* berr)
{
    return zgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gesvxx( char fact, char trans, integer n, integer nrhs,
                    float* a, integer lda, float* af, integer ldaf,
                     integer* ipiv, char* equed, float* r, float* c,
                     float* b, integer ldb, float* x, integer ldx,
                     float& rcond, float& rpvgrw, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return sgesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gesvxx( char fact, char trans, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer ldb, double* x, integer ldx,
                     double& rcond, double& rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return dgesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gesvxx( char fact, char trans, integer n, integer nrhs,
                     scomplex* a, integer lda, scomplex* af, integer ldaf,
                     integer* ipiv, char* equed, float* r, float* c,
                     scomplex* b, integer ldb, scomplex* x, integer ldx,
                     float& rcond, float& rpvgrw, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return cgesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gesvxx( char fact, char trans, integer n, integer nrhs,
                    dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                     integer* ipiv, char* equed, double* r, double* c,
                     dcomplex* b, integer ldb, dcomplex* x, integer ldx,
                     double& rcond, double& rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return zgesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gbsv( integer n, integer kl, integer ku,
                   integer nrhs, float* ab, integer ldab,
                   integer* ipiv, float* b, integer ldb )
{
    return sgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

inline integer gbsv( integer n, integer kl, integer ku,
                   integer nrhs, double* ab, integer ldab,
                   integer* ipiv, double* b, integer ldb )
{
    return dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

inline integer gbsv( integer n, integer kl, integer ku,
                   integer nrhs, scomplex* ab, integer ldab,
                   integer* ipiv, scomplex* b, integer ldb )
{
    return cgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

inline integer gbsv( integer n, integer kl, integer ku,
                   integer nrhs, dcomplex* ab, integer ldab,
                   integer* ipiv, dcomplex* b, integer ldb )
{
    return zgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
}

inline integer gbsvx( char fact, char trans, integer n, integer kl,
                    integer ku, integer nrhs, float* ab,
                    integer ldab, float* afb, integer ldafb,
                    integer* ipiv, char* equed, float* r, float* c,
                    float* b, integer ldb, float* x, integer ldx,
                    float& rcond, float* ferr, float* berr)
{
    return sgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gbsvx( char fact, char trans, integer n, integer kl,
                    integer ku, integer nrhs, double* ab,
                    integer ldab, double* afb, integer ldafb,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer ldb, double* x, integer ldx,
                    double& rcond, double* ferr, double* berr)
{
    return dgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gbsvx( char fact, char trans, integer n, integer kl,
                    integer ku, integer nrhs, scomplex* ab,
                    integer ldab, scomplex* afb, integer ldafb,
                    integer* ipiv, char* equed, float* r, float* c,
                    scomplex* b, integer ldb, scomplex* x, integer ldx,
                    float& rcond, float* ferr, float* berr)
{
    return cgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gbsvx( char fact, char trans, integer n, integer kl,
                    integer ku, integer nrhs, dcomplex* ab,
                    integer ldab, dcomplex* afb, integer ldafb,
                    integer* ipiv, char* equed, double* r, double* c,
                    dcomplex* b, integer ldb, dcomplex* x, integer ldx,
                    double& rcond, double* ferr, double* berr)
{
    return zgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gbsvxx( char fact, char trans, integer n, integer kl,
                     integer ku, integer nrhs, float* ab,
                     integer ldab, float* afb, integer ldafb,
                     integer* ipiv, char* equed, float* r, float* c,
                     float* b, integer ldb, float* x, integer ldx,
                     float& rcond, float& rpvgrw, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return sgbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gbsvxx( char fact, char trans, integer n, integer kl,
                     integer ku, integer nrhs, double* ab,
                     integer ldab, double* afb, integer ldafb,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer ldb, double* x, integer ldx,
                     double& rcond, double& rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return dgbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gbsvxx( char fact, char trans, integer n, integer kl,
                     integer ku, integer nrhs, scomplex* ab,
                     integer ldab, scomplex* afb, integer ldafb,
                     integer* ipiv, char* equed, float* r, float* c,
                     scomplex* b, integer ldb, scomplex* x, integer ldx,
                     float& rcond, float& rpvgrw, float* berr,
                     integer n_err_bnds, float* err_bnds_norm,
                     float* err_bnds_comp, integer nparams, float* params)
{
    return cgbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gbsvxx( char fact, char trans, integer n, integer kl,
                     integer ku, integer nrhs, dcomplex* ab,
                     integer ldab, dcomplex* afb, integer ldafb,
                     integer* ipiv, char* equed, double* r, double* c,
                     dcomplex* b, integer ldb, dcomplex* x, integer ldx,
                     double& rcond, double& rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params)
{
    return zgbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer gtsv( integer n, integer nrhs, float* dl, float* d,
                   float* du, float* b, integer ldb )
{
    return sgtsv(n, nrhs, dl, d, du, b, ldb);
}

inline integer gtsv( integer n, integer nrhs, double* dl, double* d,
                   double* du, double* b, integer ldb )
{
    return dgtsv(n, nrhs, dl, d, du, b, ldb);
}

inline integer gtsv( integer n, integer nrhs, scomplex* dl, scomplex* d,
                   scomplex* du, scomplex* b, integer ldb )
{
    return cgtsv(n, nrhs, dl, d, du, b, ldb);
}

inline integer gtsv( integer n, integer nrhs, dcomplex* dl, dcomplex* d,
                   dcomplex* du, dcomplex* b, integer ldb )
{
    return zgtsv(n, nrhs, dl, d, du, b, ldb);
}

inline integer gtsvx( char fact, char trans, integer n, integer nrhs,
                    const float* dl, const float* d, const float* du,
                    float* dlf, float* df, float* duf, float* du2,
                    integer* ipiv, const float* b, integer ldb,
                    float* x, integer ldx, float& rcond, float* ferr,
                    float* berr)
{
    return sgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gtsvx( char fact, char trans, integer n, integer nrhs,
                    const double* dl, const double* d, const double* du,
                    double* dlf, double* df, double* duf, double* du2,
                    integer* ipiv, const double* b, integer ldb,
                    double* x, integer ldx, double& rcond, double* ferr,
                    double* berr)
{
    return dgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gtsvx( char fact, char trans, integer n, integer nrhs,
                    const scomplex* dl, const scomplex* d, const scomplex* du,
                    scomplex* dlf, scomplex* df, scomplex* duf, scomplex* du2,
                    integer* ipiv, const scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float& rcond, float* ferr,
                    float* berr)
{
    return cgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer gtsvx( char fact, char trans, integer n, integer nrhs,
                    const dcomplex* dl, const dcomplex* d, const dcomplex* du,
                    dcomplex* dlf, dcomplex* df, dcomplex* duf, dcomplex* du2,
                    integer* ipiv, const dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double& rcond, double* ferr,
                    double* berr)
{
    return zgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv,
                  b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer posv( char uplo, integer n, integer nrhs, float* a,
                   integer lda, float* b, integer ldb )
{
    return sposv(uplo, n, nrhs, a, lda, b, ldb);
}

inline integer posv( char uplo, integer n, integer nrhs, double* a,
                   integer lda, double* b, integer ldb )
{
    return dposv(uplo, n, nrhs, a, lda, b, ldb);
}

inline integer posv( char uplo, integer n, integer nrhs, scomplex* a,
                   integer lda, scomplex* b, integer ldb )
{
    return cposv(uplo, n, nrhs, a, lda, b, ldb);
}

inline integer posv( char uplo, integer n, integer nrhs, dcomplex* a,
                   integer lda, dcomplex* b, integer ldb )
{
    return zposv(uplo, n, nrhs, a, lda, b, ldb);
}

inline integer posv( char uplo, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* x,
                    integer ldx, integer& iter )
{
    return dsposv(uplo, n, nrhs, a, lda, b, ldb, x, ldx, &iter);
}

inline integer posv( char uplo, integer n, integer nrhs, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, integer& iter )
{
    return zcposv(uplo, n, nrhs, a, lda, b, ldb, x, ldx, &iter);
}

inline integer posvx( char fact, char uplo, integer n, integer nrhs,
                    float* a, integer lda, float* af, integer ldaf,
                    char* equed, float* s, float* b, integer ldb,
                    float* x, integer ldx, float& rcond, float* ferr,
                    float* berr)
{
    return sposvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb,
                  x, ldx, &rcond, ferr, berr);
}

inline integer posvx( char fact, char uplo, integer n, integer nrhs,
                    double* a, integer lda, double* af, integer ldaf,
                    char* equed, double* s, double* b, integer ldb,
                    double* x, integer ldx, double& rcond, double* ferr,
                    double* berr)
{
    return dposvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb,
                  x, ldx, &rcond, ferr, berr);
}

inline integer posvx( char fact, char uplo, integer n, integer nrhs,
                    scomplex* a, integer lda, scomplex* af, integer ldaf,
                    char* equed, float* s, scomplex* b, integer ldb,
                    scomplex* x, integer ldx, float& rcond, float* ferr,
                    float* berr)
{
    return cposvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb,
                  x, ldx, &rcond, ferr, berr);
}

inline integer posvx( char fact, char uplo, integer n, integer nrhs,
                    dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                    char* equed, double* s, dcomplex* b, integer ldb,
                    dcomplex* x, integer ldx, double& rcond, double* ferr,
                    double* berr)
{
    return zposvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb,
                  x, ldx, &rcond, ferr, berr);
}

inline integer posvxx( char fact, char uplo, integer n, integer nrhs,
                    float* a, integer lda, float* af, integer ldaf,
                     char* equed, float* s, float* b, integer ldb,
                     float* x, integer ldx, float& rcond, float& rpvgrw,
                     float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    return sposvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer posvxx( char fact, char uplo, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     char* equed, double* s, double* b, integer ldb,
                     double* x, integer ldx, double& rcond, double& rpvgrw,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    return dposvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer posvxx( char fact, char uplo, integer n, integer nrhs,
                     scomplex* a, integer lda, scomplex* af, integer ldaf,
                     char* equed, float* s, scomplex* b, integer ldb,
                     scomplex* x, integer ldx, float& rcond, float& rpvgrw,
                     float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    return cposvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer posvxx( char fact, char uplo, integer n, integer nrhs,
                     dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                     char* equed, double* s, dcomplex* b, integer ldb,
                     dcomplex* x, integer ldx, double& rcond, double& rpvgrw,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    return zposvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer ppsv( char uplo, integer n, integer nrhs, float* ap,
                    float* b, integer ldb )
{
    return sppsv(uplo, n, nrhs, ap, b, ldb);
}

inline integer ppsv( char uplo, integer n, integer nrhs, double* ap,
                   double* b, integer ldb )
{
    return dppsv(uplo, n, nrhs, ap, b, ldb);
}

inline integer ppsv( char uplo, integer n, integer nrhs, scomplex* ap,
                   scomplex* b, integer ldb )
{
    return cppsv(uplo, n, nrhs, ap, b, ldb);
}

inline integer ppsv( char uplo, integer n, integer nrhs, dcomplex* ap,
                   dcomplex* b, integer ldb )
{
    return zppsv(uplo, n, nrhs, ap, b, ldb);
}

inline integer ppsvx( char fact, char uplo, integer n, integer nrhs,
                    float* ap, float* afp, char* equed, float* s, float* b,
                    integer ldb, float* x, integer ldx, float& rcond,
                    float* ferr, float* berr)
{
    return sppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx,
                  &rcond, ferr, berr);
}

inline integer ppsvx( char fact, char uplo, integer n, integer nrhs,
                    double* ap, double* afp, char* equed, double* s, double* b,
                    integer ldb, double* x, integer ldx, double& rcond,
                    double* ferr, double* berr)
{
    return dppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx,
                  &rcond, ferr, berr);
}

inline integer ppsvx( char fact, char uplo, integer n, integer nrhs,
                    scomplex* ap, scomplex* afp, char* equed, float* s, scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float& rcond,
                    float* ferr, float* berr)
{
    return cppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx,
                  &rcond, ferr, berr);
}

inline integer ppsvx( char fact, char uplo, integer n, integer nrhs,
                    dcomplex* ap, dcomplex* afp, char* equed, double* s, dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double& rcond,
                    double* ferr, double* berr)
{
    return zppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx,
                  &rcond, ferr, berr);
}

inline integer pbsv( char uplo, integer n, integer kd, integer nrhs,
                   float* ab, integer ldab, float* b, integer ldb )
{
    return spbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer pbsv( char uplo, integer n, integer kd, integer nrhs,
                   double* ab, integer ldab, double* b, integer ldb )
{
    return dpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer pbsv( char uplo, integer n, integer kd, integer nrhs,
                   scomplex* ab, integer ldab, scomplex* b, integer ldb )
{
    return cpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer pbsv( char uplo, integer n, integer kd, integer nrhs,
                   dcomplex* ab, integer ldab, dcomplex* b, integer ldb )
{
    return zpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb);
}

inline integer pbsvx( char fact, char uplo, integer n, integer kd,
                    integer nrhs, float* ab, integer ldab, float* afb,
                    integer ldafb, char* equed, float* s, float* b,
                    integer ldb, float* x, integer ldx, float& rcond,
                    float* ferr, float* berr)
{
    return spbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb,
                  x, ldx, &rcond, ferr, berr);
}

inline integer pbsvx( char fact, char uplo, integer n, integer kd,
                    integer nrhs, double* ab, integer ldab, double* afb,
                    integer ldafb, char* equed, double* s, double* b,
                    integer ldb, double* x, integer ldx, double& rcond,
                    double* ferr, double* berr)
{
    return dpbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb,
                  x, ldx, &rcond, ferr, berr);
}

inline integer pbsvx( char fact, char uplo, integer n, integer kd,
                    integer nrhs, scomplex* ab, integer ldab, scomplex* afb,
                    integer ldafb, char* equed, float* s, scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float& rcond,
                    float* ferr, float* berr)
{
    return cpbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb,
                  x, ldx, &rcond, ferr, berr);
}

inline integer pbsvx( char fact, char uplo, integer n, integer kd,
                    integer nrhs, dcomplex* ab, integer ldab, dcomplex* afb,
                    integer ldafb, char* equed, double* s, dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double& rcond,
                    double* ferr, double* berr)
{
    return zpbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb,
                  x, ldx, &rcond, ferr, berr);
}

inline integer ptsv( integer n, integer nrhs, float* d, float* e,
                   float* b, integer ldb )
{
    return sptsv(n, nrhs, d, e, b, ldb);
}

inline integer ptsv( integer n, integer nrhs, double* d, double* e,
                   double* b, integer ldb )
{
    return dptsv(n, nrhs, d, e, b, ldb);
}

inline integer ptsv( integer n, integer nrhs, float* d, scomplex* e,
                   scomplex* b, integer ldb )
{
    return cptsv(n, nrhs, d, e, b, ldb);
}

inline integer ptsv( integer n, integer nrhs, double* d, dcomplex* e,
                   dcomplex* b, integer ldb )
{
    return zptsv(n, nrhs, d, e, b, ldb);
}

inline integer ptsvx( char fact, integer n, integer nrhs,
                    const float* d, const float* e, float* df, float* ef,
                    const float* b, integer ldb, float* x,
                    integer ldx, float& rcond, float* ferr, float* berr)
{
    return sptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer ptsvx( char fact, integer n, integer nrhs,
                    const double* d, const double* e, double* df, double* ef,
                    const double* b, integer ldb, double* x,
                    integer ldx, double& rcond, double* ferr, double* berr)
{
    return dptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer ptsvx( char fact, integer n, integer nrhs,
                    const float* d, const scomplex* e, float* df, scomplex* ef,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float& rcond, float* ferr, float* berr)
{
    return cptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer ptsvx( char fact, integer n, integer nrhs,
                    const double* d, const dcomplex* e, double* df, dcomplex* ef,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double& rcond, double* ferr, double* berr)
{
    return zptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer sysv( char uplo, integer n, integer nrhs, float* a,
                   integer lda, integer* ipiv, float* b,
                   integer ldb)
{
    return ssysv(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer sysv( char uplo, integer n, integer nrhs, double* a,
                   integer lda, integer* ipiv, double* b,
                   integer ldb)
{
    return dsysv(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer sysv( char uplo, integer n, integer nrhs, scomplex* a,
                   integer lda, integer* ipiv, scomplex* b,
                   integer ldb)
{
    return csysv(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer sysv( char uplo, integer n, integer nrhs, dcomplex* a,
                   integer lda, integer* ipiv, dcomplex* b,
                   integer ldb)
{
    return zsysv(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer hesv( char uplo, integer n, integer nrhs, float* a,
                   integer lda, integer* ipiv, float* b,
                   integer ldb)
{
    return ssysv(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer hesv( char uplo, integer n, integer nrhs, double* a,
                   integer lda, integer* ipiv, double* b,
                   integer ldb)
{
    return dsysv(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer hesv( char uplo, integer n, integer nrhs, scomplex* a,
                   integer lda, integer* ipiv, scomplex* b,
                   integer ldb)
{
    return chesv(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer hesv( char uplo, integer n, integer nrhs, dcomplex* a,
                   integer lda, integer* ipiv, dcomplex* b,
                   integer ldb)
{
    return zhesv(uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

inline integer sysvx( char fact, char uplo, integer n, integer nrhs,
                    const float* a, integer lda, float* af,
                    integer ldaf, integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float& rcond,
                    float* ferr, float* berr)
{
    return ssysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer sysvx( char fact, char uplo, integer n, integer nrhs,
                    const double* a, integer lda, double* af,
                    integer ldaf, integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double& rcond,
                    double* ferr, double* berr)
{
    return dsysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer sysvx( char fact, char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, scomplex* af,
                    integer ldaf, integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float& rcond,
                    float* ferr, float* berr)
{
    return csysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer sysvx( char fact, char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, dcomplex* af,
                    integer ldaf, integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double& rcond,
                    double* ferr, double* berr)
{
    return zsysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer hesvx( char fact, char uplo, integer n, integer nrhs,
                    const float* a, integer lda, float* af,
                    integer ldaf, integer* ipiv, const float* b,
                    integer ldb, float* x, integer ldx, float& rcond,
                    float* ferr, float* berr)
{
    return ssysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer hesvx( char fact, char uplo, integer n, integer nrhs,
                    const double* a, integer lda, double* af,
                    integer ldaf, integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double& rcond,
                    double* ferr, double* berr)
{
    return dsysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer hesvx( char fact, char uplo, integer n, integer nrhs,
                    const scomplex* a, integer lda, scomplex* af,
                    integer ldaf, integer* ipiv, const scomplex* b,
                    integer ldb, scomplex* x, integer ldx, float& rcond,
                    float* ferr, float* berr)
{
    return chesvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer hesvx( char fact, char uplo, integer n, integer nrhs,
                    const dcomplex* a, integer lda, dcomplex* af,
                    integer ldaf, integer* ipiv, const dcomplex* b,
                    integer ldb, dcomplex* x, integer ldx, double& rcond,
                    double* ferr, double* berr)
{
    return zhesvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer sysvxx( char fact, char uplo, integer n, integer nrhs,
                     float* a, integer lda, float* af, integer ldaf,
                     integer* ipiv, char* equed, float* s, float* b,
                     integer ldb, float* x, integer ldx, float& rcond,
                     float& rpvgrw, float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    return ssysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer sysvxx( char fact, char uplo, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     integer* ipiv, char* equed, double* s, double* b,
                     integer ldb, double* x, integer ldx, double& rcond,
                     double& rpvgrw, double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    return dsysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer sysvxx( char fact, char uplo, integer n, integer nrhs,
                     scomplex* a, integer lda, scomplex* af, integer ldaf,
                     integer* ipiv, char* equed, float* s, scomplex* b,
                     integer ldb, scomplex* x, integer ldx, float& rcond,
                     float& rpvgrw, float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    return csysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer sysvxx( char fact, char uplo, integer n, integer nrhs,
                     dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                     integer* ipiv, char* equed, double* s, dcomplex* b,
                     integer ldb, dcomplex* x, integer ldx, double& rcond,
                     double& rpvgrw, double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    return zsysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer hesvxx( char fact, char uplo, integer n, integer nrhs,
                     float* a, integer lda, float* af, integer ldaf,
                     integer* ipiv, char* equed, float* s, float* b,
                     integer ldb, float* x, integer ldx, float& rcond,
                     float& rpvgrw, float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    return ssysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer hesvxx( char fact, char uplo, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     integer* ipiv, char* equed, double* s, double* b,
                     integer ldb, double* x, integer ldx, double& rcond,
                     double& rpvgrw, double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    return dsysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer hesvxx( char fact, char uplo, integer n, integer nrhs,
                     scomplex* a, integer lda, scomplex* af, integer ldaf,
                     integer* ipiv, char* equed, float* s, scomplex* b,
                     integer ldb, scomplex* x, integer ldx, float& rcond,
                     float& rpvgrw, float* berr, integer n_err_bnds,
                     float* err_bnds_norm, float* err_bnds_comp,
                     integer nparams, float* params)
{
    return chesvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer hesvxx( char fact, char uplo, integer n, integer nrhs,
                     dcomplex* a, integer lda, dcomplex* af, integer ldaf,
                     integer* ipiv, char* equed, double* s, dcomplex* b,
                     integer ldb, dcomplex* x, integer ldx, double& rcond,
                     double& rpvgrw, double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params)
{
    return zhesvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                   &rcond, &rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params);
}

inline integer spsv( char uplo, integer n, integer nrhs, float* ap,
                   integer* ipiv, float* b, integer ldb )
{
    return sspsv(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer spsv( char uplo, integer n, integer nrhs, double* ap,
                   integer* ipiv, double* b, integer ldb )
{
    return dspsv(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer spsv( char uplo, integer n, integer nrhs, scomplex* ap,
                   integer* ipiv, scomplex* b, integer ldb )
{
    return cspsv(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer spsv( char uplo, integer n, integer nrhs, dcomplex* ap,
                   integer* ipiv, dcomplex* b, integer ldb )
{
    return zspsv(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer hpsv( char uplo, integer n, integer nrhs, float* ap,
                   integer* ipiv, float* b, integer ldb )
{
    return sspsv(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer hpsv( char uplo, integer n, integer nrhs, double* ap,
                   integer* ipiv, double* b, integer ldb )
{
    return dspsv(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer hpsv( char uplo, integer n, integer nrhs, scomplex* ap,
                   integer* ipiv, scomplex* b, integer ldb )
{
    return chpsv(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer hpsv( char uplo, integer n, integer nrhs, dcomplex* ap,
                   integer* ipiv, dcomplex* b, integer ldb )
{
    return zhpsv(uplo, n, nrhs, ap, ipiv, b, ldb);
}

inline integer spsvx( char fact, char uplo, integer n, integer nrhs,
                    const float* ap, float* afp, integer* ipiv,
                    const float* b, integer ldb, float* x,
                    integer ldx, float& rcond, float* ferr, float* berr)
{
    return sspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer spsvx( char fact, char uplo, integer n, integer nrhs,
                    const double* ap, double* afp, integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double& rcond, double* ferr, double* berr)
{
    return dspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer spsvx( char fact, char uplo, integer n, integer nrhs,
                    const scomplex* ap, scomplex* afp, integer* ipiv,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float& rcond, float* ferr, float* berr)
{
    return cspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer spsvx( char fact, char uplo, integer n, integer nrhs,
                    const dcomplex* ap, dcomplex* afp, integer* ipiv,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double& rcond, double* ferr, double* berr)
{
    return zspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer hpsvx( char fact, char uplo, integer n, integer nrhs,
                    const float* ap, float* afp, integer* ipiv,
                    const float* b, integer ldb, float* x,
                    integer ldx, float& rcond, float* ferr, float* berr)
{
    return sspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer hpsvx( char fact, char uplo, integer n, integer nrhs,
                    const double* ap, double* afp, integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double& rcond, double* ferr, double* berr)
{
    return dspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer hpsvx( char fact, char uplo, integer n, integer nrhs,
                    const scomplex* ap, scomplex* afp, integer* ipiv,
                    const scomplex* b, integer ldb, scomplex* x,
                    integer ldx, float& rcond, float* ferr, float* berr)
{
    return chpsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer hpsvx( char fact, char uplo, integer n, integer nrhs,
                    const dcomplex* ap, dcomplex* afp, integer* ipiv,
                    const dcomplex* b, integer ldb, dcomplex* x,
                    integer ldx, double& rcond, double* ferr, double* berr)
{
    return zhpsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, &rcond, ferr, berr);
}

inline integer geqrf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    return sgeqrf(m, n, a, lda, tau);
}

inline integer geqrf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    return dgeqrf(m, n, a, lda, tau);
}

inline integer geqrf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    return cgeqrf(m, n, a, lda, tau);
}

inline integer geqrf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    return zgeqrf(m, n, a, lda, tau);
}

inline integer geqpf( integer m, integer n, float* a, integer lda,
                    integer* jpvt, float* tau)
{
    return sgeqpf(m, n, a, lda, jpvt, tau);
}

inline integer geqpf( integer m, integer n, double* a, integer lda,
                    integer* jpvt, double* tau)
{
    return dgeqpf(m, n, a, lda, jpvt, tau);
}

inline integer geqpf( integer m, integer n, scomplex* a, integer lda,
                    integer* jpvt, scomplex* tau)
{
    return cgeqpf(m, n, a, lda, jpvt, tau);
}

inline integer geqpf( integer m, integer n, dcomplex* a, integer lda,
                    integer* jpvt, dcomplex* tau)
{
    return zgeqpf(m, n, a, lda, jpvt, tau);
}

inline integer geqp3( integer m, integer n, float* a, integer lda,
                    integer* jpvt, float* tau)
{
    return sgeqp3(m, n, a, lda, jpvt, tau);
}

inline integer geqp3( integer m, integer n, double* a, integer lda,
                    integer* jpvt, double* tau)
{
    return dgeqp3(m, n, a, lda, jpvt, tau);
}

inline integer geqp3( integer m, integer n, scomplex* a, integer lda,
                    integer* jpvt, scomplex* tau)
{
    return cgeqp3(m, n, a, lda, jpvt, tau);
}

inline integer geqp3( integer m, integer n, dcomplex* a, integer lda,
                    integer* jpvt, dcomplex* tau)
{
    return zgeqp3(m, n, a, lda, jpvt, tau);
}

inline integer ungqr( integer m, integer n, integer k, float* a,
                    integer lda, const float* tau)
{
    return sorgqr(m, n, k, a, lda, tau);
}

inline integer ungqr( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau)
{
    return dorgqr(m, n, k, a, lda, tau);
}

inline integer ungqr( integer m, integer n, integer k, scomplex* a,
                    integer lda, const scomplex* tau)
{
    return cungqr(m, n, k, a, lda, tau);
}

inline integer ungqr( integer m, integer n, integer k, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    return zungqr(m, n, k, a, lda, tau);
}

inline integer unmqr( char side, char trans, integer m, integer n,
                    integer k, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    return sormqr(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmqr( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    return dormqr(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmqr( char side, char trans, integer m, integer n,
                    integer k, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    return cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmqr( char side, char trans, integer m, integer n,
                    integer k, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    return zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer gelqf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    return sgelqf(m, n, a, lda, tau);
}

inline integer gelqf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    return dgelqf(m, n, a, lda, tau);
}

inline integer gelqf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    return cgelqf(m, n, a, lda, tau);
}

inline integer gelqf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    return zgelqf(m, n, a, lda, tau);
}

inline integer unglq( integer m, integer n, integer k, float* a,
                    integer lda, const float* tau)
{
    return sorglq(m, n, k, a, lda, tau);
}

inline integer unglq( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau)
{
    return dorglq(m, n, k, a, lda, tau);
}

inline integer unglq( integer m, integer n, integer k, scomplex* a,
                    integer lda, const scomplex* tau)
{
    return cunglq(m, n, k, a, lda, tau);
}

inline integer unglq( integer m, integer n, integer k, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    return zunglq(m, n, k, a, lda, tau);
}

inline integer unmlq( char side, char trans, integer m, integer n,
                    integer k, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    return sormlq(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmlq( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    return dormlq(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmlq( char side, char trans, integer m, integer n,
                    integer k, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    return cunmlq(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmlq( char side, char trans, integer m, integer n,
                    integer k, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    return zunmlq(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer geqlf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    return sgeqlf(m, n, a, lda, tau);
}

inline integer geqlf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    return dgeqlf(m, n, a, lda, tau);
}

inline integer geqlf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    return cgeqlf(m, n, a, lda, tau);
}

inline integer geqlf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    return zgeqlf(m, n, a, lda, tau);
}

inline integer ungql( integer m, integer n, integer k, float* a,
                    integer lda, const float* tau)
{
    return sorgql(m, n, k, a, lda, tau);
}

inline integer ungql( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau)
{
    return dorgql(m, n, k, a, lda, tau);
}

inline integer ungql( integer m, integer n, integer k, scomplex* a,
                    integer lda, const scomplex* tau)
{
    return cungql(m, n, k, a, lda, tau);
}

inline integer ungql( integer m, integer n, integer k, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    return zungql(m, n, k, a, lda, tau);
}

inline integer unmql( char side, char trans, integer m, integer n,
                    integer k, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    return sormql(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmql( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    return dormql(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmql( char side, char trans, integer m, integer n,
                    integer k, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    return cunmql(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmql( char side, char trans, integer m, integer n,
                    integer k, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    return zunmql(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer gerqf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    return sgerqf(m, n, a, lda, tau);
}

inline integer gerqf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    return dgerqf(m, n, a, lda, tau);
}

inline integer gerqf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    return cgerqf(m, n, a, lda, tau);
}

inline integer gerqf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    return zgerqf(m, n, a, lda, tau);
}

inline integer ungrq( integer m, integer n, integer k, float* a,
                    integer lda, const float* tau)
{
    return sorgrq(m, n, k, a, lda, tau);
}

inline integer ungrq( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau)
{
    return dorgrq(m, n, k, a, lda, tau);
}

inline integer ungrq( integer m, integer n, integer k, scomplex* a,
                    integer lda, const scomplex* tau)
{
    return cungrq(m, n, k, a, lda, tau);
}

inline integer ungrq( integer m, integer n, integer k, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    return zungrq(m, n, k, a, lda, tau);
}

inline integer unmrq( char side, char trans, integer m, integer n,
                    integer k, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    return sormrq(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmrq( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    return dormrq(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmrq( char side, char trans, integer m, integer n,
                    integer k, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    return cunmrq(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmrq( char side, char trans, integer m, integer n,
                    integer k, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    return zunmrq(side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer tzrzf( integer m, integer n, float* a, integer lda,
                    float* tau)
{
    return stzrzf(m, n, a, lda, tau);
}

inline integer tzrzf( integer m, integer n, double* a, integer lda,
                    double* tau)
{
    return dtzrzf(m, n, a, lda, tau);
}

inline integer tzrzf( integer m, integer n, scomplex* a, integer lda,
                    scomplex* tau)
{
    return ctzrzf(m, n, a, lda, tau);
}

inline integer tzrzf( integer m, integer n, dcomplex* a, integer lda,
                    dcomplex* tau)
{
    return ztzrzf(m, n, a, lda, tau);
}

inline integer unmrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const float* a,
                    integer lda, const float* tau, float* c,
                    integer ldc)
{
    return sormrz(side, trans, m, n, k, l, a, lda, tau, c, ldc);
}

inline integer unmrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc)
{
    return dormrz(side, trans, m, n, k, l, a, lda, tau, c, ldc);
}

inline integer unmrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const scomplex* a,
                    integer lda, const scomplex* tau, scomplex* c,
                    integer ldc)
{
    return cunmrz(side, trans, m, n, k, l, a, lda, tau, c, ldc);
}

inline integer unmrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const dcomplex* a,
                    integer lda, const dcomplex* tau, dcomplex* c,
                    integer ldc)
{
    return zunmrz(side, trans, m, n, k, l, a, lda, tau, c, ldc);
}

inline integer ggqrf( integer n, integer m, integer p, float* a,
                    integer lda, float* taua, float* b, integer ldb,
                    float* taub)
{
    return sggqrf(n, m, p, a, lda, taua, b, ldb, taub);
}

inline integer ggqrf( integer n, integer m, integer p, double* a,
                    integer lda, double* taua, double* b, integer ldb,
                    double* taub)
{
    return dggqrf(n, m, p, a, lda, taua, b, ldb, taub);
}

inline integer ggqrf( integer n, integer m, integer p, scomplex* a,
                    integer lda, scomplex* taua, scomplex* b, integer ldb,
                    scomplex* taub)
{
    return cggqrf(n, m, p, a, lda, taua, b, ldb, taub);
}

inline integer ggqrf( integer n, integer m, integer p, dcomplex* a,
                    integer lda, dcomplex* taua, dcomplex* b, integer ldb,
                    dcomplex* taub)
{
    return zggqrf(n, m, p, a, lda, taua, b, ldb, taub);
}

inline integer ggrqf( integer m, integer p, integer n, float* a,
                    integer lda, float* taua, float* b, integer ldb,
                    float* taub)
{
    return sggrqf(m, p, n, a, lda, taua, b, ldb, taub);
}

inline integer ggrqf( integer m, integer p, integer n, double* a,
                    integer lda, double* taua, double* b, integer ldb,
                    double* taub)
{
    return dggrqf(m, p, n, a, lda, taua, b, ldb, taub);
}

inline integer ggrqf( integer m, integer p, integer n, scomplex* a,
                    integer lda, scomplex* taua, scomplex* b, integer ldb,
                    scomplex* taub)
{
    return cggrqf(m, p, n, a, lda, taua, b, ldb, taub);
}

inline integer ggrqf( integer m, integer p, integer n, dcomplex* a,
                    integer lda, dcomplex* taua, dcomplex* b, integer ldb,
                    dcomplex* taub)
{
    return zggrqf(m, p, n, a, lda, taua, b, ldb, taub);
}

inline integer gebrd( integer m, integer n, float* a, integer lda,
                    float* d, float* e, float* tauq, float* taup)
{
    return sgebrd(m, n, a, lda, d, e, tauq, taup);
}

inline integer gebrd( integer m, integer n, double* a, integer lda,
                    double* d, double* e, double* tauq, double* taup)
{
    return dgebrd(m, n, a, lda, d, e, tauq, taup);
}

inline integer gebrd( integer m, integer n, scomplex* a, integer lda,
                    float* d, float* e, scomplex* tauq, scomplex* taup)
{
    return cgebrd(m, n, a, lda, d, e, tauq, taup);
}

inline integer gebrd( integer m, integer n, dcomplex* a, integer lda,
                    double* d, double* e, dcomplex* tauq, dcomplex* taup)
{
    return zgebrd(m, n, a, lda, d, e, tauq, taup);
}

inline integer gbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, float* ab,
                    integer ldab, float* d, float* e, float* q,
                    integer ldq, float* pt, integer ldpt, float* c,
                    integer ldc)
{
    return sgbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc);
}

inline integer gbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, double* ab,
                    integer ldab, double* d, double* e, double* q,
                    integer ldq, double* pt, integer ldpt, double* c,
                    integer ldc)
{
    return dgbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc);
}

inline integer gbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, scomplex* ab,
                    integer ldab, float* d, float* e, scomplex* q,
                    integer ldq, scomplex* pt, integer ldpt, scomplex* c,
                    integer ldc)
{
    return cgbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc);
}

inline integer gbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, dcomplex* ab,
                    integer ldab, double* d, double* e, dcomplex* q,
                    integer ldq, dcomplex* pt, integer ldpt, dcomplex* c,
                    integer ldc)
{
    return zgbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc);
}

inline integer ungbr( char vect, integer m, integer n, integer k,
                    float* a, integer lda, const float* tau)
{
    return sorgbr(vect, m, n, k, a, lda, tau);
}

inline integer ungbr( char vect, integer m, integer n, integer k,
                    double* a, integer lda, const double* tau)
{
    return dorgbr(vect, m, n, k, a, lda, tau);
}

inline integer ungbr( char vect, integer m, integer n, integer k,
                    scomplex* a, integer lda, const scomplex* tau)
{
    return cungbr(vect, m, n, k, a, lda, tau);
}

inline integer ungbr( char vect, integer m, integer n, integer k,
                    dcomplex* a, integer lda, const dcomplex* tau)
{
    return zungbr(vect, m, n, k, a, lda, tau);
}

inline integer unmbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const float* a,
                    integer lda, const float* tau, float* c,
                    integer ldc)
{
    return sormbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc)
{
    return dormbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const scomplex* a,
                    integer lda, const scomplex* tau, scomplex* c,
                    integer ldc)
{
    return cunmbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer unmbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const dcomplex* a,
                    integer lda, const dcomplex* tau, dcomplex* c,
                    integer ldc)
{
    return zunmbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc);
}

inline integer bdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, float* d, float* e,
                    float* vt, integer ldvt, float* u, integer ldu,
                    float* c, integer ldc)
{
    return sbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc);
}

inline integer bdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, double* d, double* e,
                    double* vt, integer ldvt, double* u, integer ldu,
                    double* c, integer ldc)
{
    return dbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc);
}

inline integer bdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, float* d, float* e,
                    scomplex* vt, integer ldvt, scomplex* u, integer ldu,
                    scomplex* c, integer ldc)
{
    return cbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc);
}

inline integer bdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, double* d, double* e,
                    dcomplex* vt, integer ldvt, dcomplex* u, integer ldu,
                    dcomplex* c, integer ldc)
{
    return zbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc);
}

inline integer bdsdc( char uplo, char compq, integer n, float* d,
                    float* e, float* u, integer ldu, float* vt,
                    integer ldvt, float* q, integer* iq)
{
    return sbdsdc(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq);
}

inline integer bdsdc( char uplo, char compq, integer n, double* d,
                    double* e, double* u, integer ldu, double* vt,
                    integer ldvt, double* q, integer* iq)
{
    return dbdsdc(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq);
}

inline integer hetrd( char uplo, integer n, float* a, integer lda,
                    float* d, float* e, float* tau)
{
    return ssytrd(uplo, n, a, lda, d, e, tau);
}

inline integer hetrd( char uplo, integer n, double* a, integer lda,
                    double* d, double* e, double* tau)
{
    return dsytrd(uplo, n, a, lda, d, e, tau);
}

inline integer hetrd( char uplo, integer n, scomplex* a, integer lda,
                    float* d, float* e, scomplex* tau)
{
    return chetrd(uplo, n, a, lda, d, e, tau);
}

inline integer hetrd( char uplo, integer n, dcomplex* a, integer lda,
                    double* d, double* e, dcomplex* tau)
{
    return zhetrd(uplo, n, a, lda, d, e, tau);
}

inline integer ungtr( char uplo, integer n, float* a, integer lda,
                    const float* tau)
{
    return sorgtr(uplo, n, a, lda, tau);
}

inline integer ungtr( char uplo, integer n, double* a, integer lda,
                    const double* tau)
{
    return dorgtr(uplo, n, a, lda, tau);
}

inline integer ungtr( char uplo, integer n, scomplex* a, integer lda,
                    const scomplex* tau)
{
    return cungtr(uplo, n, a, lda, tau);
}

inline integer ungtr( char uplo, integer n, dcomplex* a, integer lda,
                    const dcomplex* tau)
{
    return zungtr(uplo, n, a, lda, tau);
}

inline integer unmtr( char side, char uplo, char trans, integer m,
                    integer n, const float* a, integer lda,
                    const float* tau, float* c, integer ldc)
{
    return sormtr(side, uplo, trans, m, n, a, lda, tau, c, ldc);
}

inline integer unmtr( char side, char uplo, char trans, integer m,
                    integer n, const double* a, integer lda,
                    const double* tau, double* c, integer ldc)
{
    return dormtr(side, uplo, trans, m, n, a, lda, tau, c, ldc);
}

inline integer unmtr( char side, char uplo, char trans, integer m,
                    integer n, const scomplex* a, integer lda,
                    const scomplex* tau, scomplex* c, integer ldc)
{
    return cunmtr(side, uplo, trans, m, n, a, lda, tau, c, ldc);
}

inline integer unmtr( char side, char uplo, char trans, integer m,
                    integer n, const dcomplex* a, integer lda,
                    const dcomplex* tau, dcomplex* c, integer ldc)
{
    return zunmtr(side, uplo, trans, m, n, a, lda, tau, c, ldc);
}

inline integer hptrd( char uplo, integer n, float* ap, float* d, float* e,
                    float* tau )
{
    return ssptrd(uplo, n, ap, d, e, tau);
}

inline integer hptrd( char uplo, integer n, double* ap, double* d, double* e,
                    double* tau )
{
    return dsptrd(uplo, n, ap, d, e, tau);
}

inline integer hptrd( char uplo, integer n, scomplex* ap, float* d, float* e,
                    scomplex* tau )
{
    return chptrd(uplo, n, ap, d, e, tau);
}

inline integer hptrd( char uplo, integer n, dcomplex* ap, double* d, double* e,
                    dcomplex* tau )
{
    return zhptrd(uplo, n, ap, d, e, tau);
}

inline integer upgtr( char uplo, integer n, const float* ap,
                    const float* tau, float* q, integer ldq)
{
    return sopgtr(uplo, n, ap, tau, q, ldq);
}

inline integer upgtr( char uplo, integer n, const double* ap,
                    const double* tau, double* q, integer ldq)
{
    return dopgtr(uplo, n, ap, tau, q, ldq);
}

inline integer upgtr( char uplo, integer n, const scomplex* ap,
                    const scomplex* tau, scomplex* q, integer ldq)
{
    return cupgtr(uplo, n, ap, tau, q, ldq);
}

inline integer upgtr( char uplo, integer n, const dcomplex* ap,
                    const dcomplex* tau, dcomplex* q, integer ldq)
{
    return zupgtr(uplo, n, ap, tau, q, ldq);
}

inline integer upmtr( char side, char uplo, char trans, integer m,
                    integer n, const float* ap, const float* tau,
                    float* c, integer ldc)
{
    return sopmtr(side, uplo, trans, m, n, ap, tau, c, ldc);
}

inline integer upmtr( char side, char uplo, char trans, integer m,
                    integer n, const double* ap, const double* tau,
                    double* c, integer ldc)
{
    return dopmtr(side, uplo, trans, m, n, ap, tau, c, ldc);
}

inline integer upmtr( char side, char uplo, char trans, integer m,
                    integer n, const scomplex* ap, const scomplex* tau,
                    scomplex* c, integer ldc)
{
    return cupmtr(side, uplo, trans, m, n, ap, tau, c, ldc);
}

inline integer upmtr( char side, char uplo, char trans, integer m,
                    integer n, const dcomplex* ap, const dcomplex* tau,
                    dcomplex* c, integer ldc)
{
    return zupmtr(side, uplo, trans, m, n, ap, tau, c, ldc);
}

inline integer hbtrd( char vect, char uplo, integer n, integer kd,
                    float* ab, integer ldab, float* d, float* e,
                    float* q, integer ldq)
{
    return ssbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq);
}

inline integer hbtrd( char vect, char uplo, integer n, integer kd,
                    double* ab, integer ldab, double* d, double* e,
                    double* q, integer ldq)
{
    return dsbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq);
}

inline integer hbtrd( char vect, char uplo, integer n, integer kd,
                    scomplex* ab, integer ldab, float* d, float* e,
                    scomplex* q, integer ldq)
{
    return chbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq);
}

inline integer hbtrd( char vect, char uplo, integer n, integer kd,
                    dcomplex* ab, integer ldab, double* d, double* e,
                    dcomplex* q, integer ldq)
{
    return zhbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq);
}

inline integer sterf( integer n, float* d, float* e )
{
    return ssterf(n, d, e);
}

inline integer sterf( integer n, double* d, double* e )
{
    return dsterf(n, d, e);
}

inline integer steqr( char compz, integer n, float* d, float* e, float* z,
                    integer ldz)
{
    return ssteqr(compz, n, d, e, z, ldz);
}

inline integer steqr( char compz, integer n, double* d, double* e, double* z,
                    integer ldz)
{
    return dsteqr(compz, n, d, e, z, ldz);
}

inline integer steqr( char compz, integer n, float* d, float* e, scomplex* z,
                    integer ldz)
{
    return csteqr(compz, n, d, e, z, ldz);
}

inline integer steqr( char compz, integer n, double* d, double* e, dcomplex* z,
                    integer ldz)
{
    return zsteqr(compz, n, d, e, z, ldz);
}

inline integer stemr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, integer& m, float* w, float* z,
                    integer ldz, integer nzc, integer* isuppz,
                    logical& tryrac)
{
    return sstemr(jobz, range, n, d, e, vl, vu, il, iu, &m, w, z, ldz, nzc, isuppz, &tryrac);
}

inline integer stemr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, integer& m, double* w, double* z,
                    integer ldz, integer nzc, integer* isuppz,
                    logical& tryrac)
{
    return dstemr(jobz, range, n, d, e, vl, vu, il, iu, &m, w, z, ldz, nzc, isuppz, &tryrac);
}

inline integer stemr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, integer& m, float* w, scomplex* z,
                    integer ldz, integer nzc, integer* isuppz,
                    logical& tryrac)
{
    return cstemr(jobz, range, n, d, e, vl, vu, il, iu, &m, w, z, ldz, nzc, isuppz, &tryrac);
}

inline integer stemr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, integer& m, double* w, dcomplex* z,
                    integer ldz, integer nzc, integer* isuppz,
                    logical& tryrac)
{
    return zstemr(jobz, range, n, d, e, vl, vu, il, iu, &m, w, z, ldz, nzc, isuppz, &tryrac);
}

inline integer stedc( char compz, integer n, float* d, float* e, float* z,
                    integer ldz)
{
    return sstedc(compz, n, d, e, z, ldz);
}

inline integer stedc( char compz, integer n, double* d, double* e, double* z,
                    integer ldz)
{
    return dstedc(compz, n, d, e, z, ldz);
}

inline integer stedc( char compz, integer n, float* d, float* e, scomplex* z,
                    integer ldz)
{
    return cstedc(compz, n, d, e, z, ldz);
}

inline integer stedc( char compz, integer n, double* d, double* e, dcomplex* z,
                    integer ldz)
{
    return zstedc(compz, n, d, e, z, ldz);
}

inline integer stegr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    float* z, integer ldz, integer* isuppz)
{
    return sstegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer stegr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    double* z, integer ldz, integer* isuppz)
{
    return dstegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer stegr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    scomplex* z, integer ldz, integer* isuppz)
{
    return cstegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer stegr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    dcomplex* z, integer ldz, integer* isuppz)
{
    return zstegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer pteqr( char compz, integer n, float* d, float* e, float* z,
                    integer ldz)
{
    return spteqr(compz, n, d, e, z, ldz);
}

inline integer pteqr( char compz, integer n, double* d, double* e, double* z,
                    integer ldz)
{
    return dpteqr(compz, n, d, e, z, ldz);
}

inline integer pteqr( char compz, integer n, float* d, float* e, scomplex* z,
                    integer ldz)
{
    return cpteqr(compz, n, d, e, z, ldz);
}

inline integer pteqr( char compz, integer n, double* d, double* e, dcomplex* z,
                    integer ldz)
{
    return zpteqr(compz, n, d, e, z, ldz);
}

inline integer stebz( char range, char order, integer n, float vl,
                    float vu, integer il, integer iu, float abstol,
                    const float* d, const float* e, integer& m,
                    integer& nsplit, float* w, integer* iblock,
                    integer* isplit)
{
    return sstebz(range, order, n, vl, vu, il, iu, abstol, d, e, &m, &nsplit, w, iblock, isplit);
}

inline integer stebz( char range, char order, integer n, double vl,
                    double vu, integer il, integer iu, double abstol,
                    const double* d, const double* e, integer& m,
                    integer& nsplit, double* w, integer* iblock,
                    integer* isplit)
{
    return dstebz(range, order, n, vl, vu, il, iu, abstol, d, e, &m, &nsplit, w, iblock, isplit);
}

inline integer stein( integer n, const float* d, const float* e,
                    integer m, const float* w, const integer* iblock,
                    const integer* isplit, float* z, integer ldz,
                    integer* ifailv )
{
    return sstein(n, d, e, m, w, iblock, isplit, z, ldz, ifailv);
}

inline integer stein( integer n, const double* d, const double* e,
                    integer m, const double* w, const integer* iblock,
                    const integer* isplit, double* z, integer ldz,
                    integer* ifailv )
{
    return dstein(n, d, e, m, w, iblock, isplit, z, ldz, ifailv);
}

inline integer stein( integer n, const float* d, const float* e,
                    integer m, const float* w, const integer* iblock,
                    const integer* isplit, scomplex* z, integer ldz,
                    integer* ifailv )
{
    return cstein(n, d, e, m, w, iblock, isplit, z, ldz, ifailv);
}

inline integer stein( integer n, const double* d, const double* e,
                    integer m, const double* w, const integer* iblock,
                    const integer* isplit, dcomplex* z, integer ldz,
                    integer* ifailv )
{
    return zstein(n, d, e, m, w, iblock, isplit, z, ldz, ifailv);
}

inline integer disna( char job, integer m, integer n, const float* d,
                    float* sep )
{
    return sdisna(job, m, n, d, sep);
}

inline integer disna( char job, integer m, integer n, const double* d,
                    double* sep )
{
    return ddisna(job, m, n, d, sep);
}

inline integer hegst( integer itype, char uplo, integer n, float* a,
                    integer lda, const float* b, integer ldb )
{
    return ssygst(itype, uplo, n, a, lda, b, ldb);
}

inline integer hegst( integer itype, char uplo, integer n, double* a,
                    integer lda, const double* b, integer ldb )
{
    return dsygst(itype, uplo, n, a, lda, b, ldb);
}

inline integer hegst( integer itype, char uplo, integer n, scomplex* a,
                    integer lda, const scomplex* b, integer ldb )
{
    return chegst(itype, uplo, n, a, lda, b, ldb);
}

inline integer hegst( integer itype, char uplo, integer n, dcomplex* a,
                    integer lda, const dcomplex* b, integer ldb )
{
    return zhegst(itype, uplo, n, a, lda, b, ldb);
}

inline integer hpgst( integer itype, char uplo, integer n, float* ap,
                    const float* bp )
{
    return sspgst(itype, uplo, n, ap, bp);
}

inline integer hpgst( integer itype, char uplo, integer n, double* ap,
                    const double* bp )
{
    return dspgst(itype, uplo, n, ap, bp);
}

inline integer hpgst( integer itype, char uplo, integer n, scomplex* ap,
                    const scomplex* bp )
{
    return chpgst(itype, uplo, n, ap, bp);
}

inline integer hpgst( integer itype, char uplo, integer n, dcomplex* ap,
                    const dcomplex* bp )
{
    return zhpgst(itype, uplo, n, ap, bp);
}

inline integer hbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, float* ab, integer ldab,
                    const float* bb, integer ldbb, float* x,
                    integer ldx)
{
    return ssbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx);
}

inline integer hbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, double* ab, integer ldab,
                    const double* bb, integer ldbb, double* x,
                    integer ldx)
{
    return dsbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx);
}

inline integer hbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, scomplex* ab, integer ldab,
                    const scomplex* bb, integer ldbb, scomplex* x,
                    integer ldx)
{
    return chbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx);
}

inline integer hbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, dcomplex* ab, integer ldab,
                    const dcomplex* bb, integer ldbb, dcomplex* x,
                    integer ldx)
{
    return zhbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx);
}

inline integer pbstf( char uplo, integer n, integer kb, float* bb,
                    integer ldbb )
{
    return spbstf(uplo, n, kb, bb, ldbb);
}

inline integer pbstf( char uplo, integer n, integer kb, double* bb,
                    integer ldbb )
{
    return dpbstf(uplo, n, kb, bb, ldbb);
}

inline integer pbstf( char uplo, integer n, integer kb, scomplex* bb,
                    integer ldbb )
{
    return cpbstf(uplo, n, kb, bb, ldbb);
}

inline integer pbstf( char uplo, integer n, integer kb, dcomplex* bb,
                    integer ldbb )
{
    return zpbstf(uplo, n, kb, bb, ldbb);
}

inline integer gehrd( integer n, integer ilo, integer ihi, float* a,
                    integer lda, float* tau)
{
    return sgehrd(n, ilo, ihi, a, lda, tau);
}

inline integer gehrd( integer n, integer ilo, integer ihi, double* a,
                    integer lda, double* tau)
{
    return dgehrd(n, ilo, ihi, a, lda, tau);
}

inline integer gehrd( integer n, integer ilo, integer ihi, scomplex* a,
                    integer lda, scomplex* tau)
{
    return cgehrd(n, ilo, ihi, a, lda, tau);
}

inline integer gehrd( integer n, integer ilo, integer ihi, dcomplex* a,
                    integer lda, dcomplex* tau)
{
    return zgehrd(n, ilo, ihi, a, lda, tau);
}

inline integer unghr( integer n, integer ilo, integer ihi, float* a,
                    integer lda, const float* tau)
{
    return sorghr(n, ilo, ihi, a, lda, tau);
}

inline integer unghr( integer n, integer ilo, integer ihi, double* a,
                    integer lda, const double* tau)
{
    return dorghr(n, ilo, ihi, a, lda, tau);
}

inline integer unghr( integer n, integer ilo, integer ihi, scomplex* a,
                    integer lda, const scomplex* tau)
{
    return cunghr(n, ilo, ihi, a, lda, tau);
}

inline integer unghr( integer n, integer ilo, integer ihi, dcomplex* a,
                    integer lda, const dcomplex* tau)
{
    return zunghr(n, ilo, ihi, a, lda, tau);
}

inline integer unmhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const float* a,
                    integer lda, const float* tau, float* c,
                    integer ldc)
{
    return sormhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc);
}

inline integer unmhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc)
{
    return dormhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc);
}

inline integer unmhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const scomplex* a,
                    integer lda, const scomplex* tau, scomplex* c,
                    integer ldc)
{
    return cunmhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc);
}

inline integer unmhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const dcomplex* a,
                    integer lda, const dcomplex* tau, dcomplex* c,
                    integer ldc)
{
    return zunmhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc);
}

inline integer gebal( char job, integer n, float* a, integer lda,
                    integer& ilo, integer& ihi, float* scale )
{
    return sgebal(job, n, a, lda, &ilo, &ihi, scale);
}

inline integer gebal( char job, integer n, double* a, integer lda,
                    integer& ilo, integer& ihi, double* scale )
{
    return dgebal(job, n, a, lda, &ilo, &ihi, scale);
}

inline integer gebal( char job, integer n, scomplex* a, integer lda,
                    integer& ilo, integer& ihi, float* scale )
{
    return cgebal(job, n, a, lda, &ilo, &ihi, scale);
}

inline integer gebal( char job, integer n, dcomplex* a, integer lda,
                    integer& ilo, integer& ihi, double* scale )
{
    return zgebal(job, n, a, lda, &ilo, &ihi, scale);
}

inline integer gebak( char job, char side, integer n, integer ilo,
                    integer ihi, const float* scale, integer m,
                    float* v, integer ldv )
{
    return sgebak(job, side, n, ilo, ihi, scale, m, v, ldv);
}

inline integer gebak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* scale, integer m,
                    double* v, integer ldv )
{
    return dgebak(job, side, n, ilo, ihi, scale, m, v, ldv);
}

inline integer gebak( char job, char side, integer n, integer ilo,
                    integer ihi, const float* scale, integer m,
                    scomplex* v, integer ldv )
{
    return cgebak(job, side, n, ilo, ihi, scale, m, v, ldv);
}

inline integer gebak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* scale, integer m,
                    dcomplex* v, integer ldv )
{
    return zgebak(job, side, n, ilo, ihi, scale, m, v, ldv);
}

inline integer hseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, float* h, integer ldh, scomplex* w,
                    float* z, integer ldz)
{
    return shseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz);
}

inline integer hseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, double* h, integer ldh, dcomplex* w,
                    double* z, integer ldz)
{
    return dhseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz);
}

inline integer hseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, scomplex* h, integer ldh, scomplex* w,
                    scomplex* z, integer ldz)
{
    return chseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz);
}

inline integer hseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, dcomplex* h, integer ldh, dcomplex* w,
                    dcomplex* z, integer ldz)
{
    return zhseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz);
}

inline integer hsein( char job, char eigsrc, char initv,
                    logical* select, integer n, const float* h,
                    integer ldh, const scomplex* w, float* vl,
                    integer ldvl, float* vr, integer ldvr,
                    integer mm, integer& m,
                    integer* ifaill, integer* ifailr )
{
    return shsein(job, eigsrc, initv, select, n, h, ldh, w, vl, ldvl, vr, ldvr, mm, &m, ifaill, ifailr);
}

inline integer hsein( char job, char eigsrc, char initv,
                    logical* select, integer n, const double* h,
                    integer ldh, const dcomplex* w, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer& m,
                    integer* ifaill, integer* ifailr )
{
    return dhsein(job, eigsrc, initv, select, n, h, ldh, w, vl, ldvl, vr, ldvr, mm, &m, ifaill, ifailr);
}

inline integer hsein( char job, char eigsrc, char initv,
                    logical* select, integer n, const scomplex* h,
                    integer ldh, scomplex* w, scomplex* vl,
                    integer ldvl, scomplex* vr, integer ldvr,
                    integer mm, integer& m,
                    integer* ifaill, integer* ifailr )
{
    return chsein(job, eigsrc, initv, select, n, h, ldh, w, vl, ldvl, vr, ldvr, mm, &m, ifaill, ifailr);
}

inline integer hsein( char job, char eigsrc, char initv,
                    logical* select, integer n, const dcomplex* h,
                    integer ldh, dcomplex* w, dcomplex* vl,
                    integer ldvl, dcomplex* vr, integer ldvr,
                    integer mm, integer& m,
                    integer* ifaill, integer* ifailr )
{
    return zhsein(job, eigsrc, initv, select, n, h, ldh, w, vl, ldvl, vr, ldvr, mm, &m, ifaill, ifailr);
}

inline integer trevc( char side, char howmny, const logical* select,
                    integer n, const float* t, integer ldt, float* vl,
                    integer ldvl, float* vr, integer ldvr,
                    integer mm, integer& m)
{
    return strevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, &m);
}

inline integer trevc( char side, char howmny, const logical* select,
                    integer n, const double* t, integer ldt, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer& m)
{
    return dtrevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, &m);
}

inline integer trevc( char side, char howmny, const logical* select,
                    integer n, const scomplex* t, integer ldt, scomplex* vl,
                    integer ldvl, scomplex* vr, integer ldvr,
                    integer mm, integer& m)
{
    return ctrevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, &m);
}

inline integer trevc( char side, char howmny, const logical* select,
                    integer n, const dcomplex* t, integer ldt, dcomplex* vl,
                    integer ldvl, dcomplex* vr, integer ldvr,
                    integer mm, integer& m)
{
    return ztrevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, &m);
}

inline integer trsna( char job, char howmny, const logical* select,
                    integer n, const float* t, integer ldt,
                    const float* vl, integer ldvl, const float* vr,
                    integer ldvr, float* s, float* sep, integer mm,
                    integer& m)
{
    return strsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, &m);
}

inline integer trsna( char job, char howmny, const logical* select,
                    integer n, const double* t, integer ldt,
                    const double* vl, integer ldvl, const double* vr,
                    integer ldvr, double* s, double* sep, integer mm,
                    integer& m)
{
    return dtrsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, &m);
}

inline integer trsna( char job, char howmny, const logical* select,
                    integer n, const scomplex* t, integer ldt,
                    const scomplex* vl, integer ldvl, const scomplex* vr,
                    integer ldvr, float* s, float* sep, integer mm,
                    integer& m)
{
    return ctrsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, &m);
}

inline integer trsna( char job, char howmny, const logical* select,
                    integer n, const dcomplex* t, integer ldt,
                    const dcomplex* vl, integer ldvl, const dcomplex* vr,
                    integer ldvr, double* s, double* sep, integer mm,
                    integer& m)
{
    return ztrsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, &m);
}

inline integer trexc( char compq, integer n, float* t, integer ldt,
                    float* q, integer ldq, integer& ifst,
                    integer& ilst)
{
    return strexc(compq, n, t, ldt, q, ldq, &ifst, &ilst);
}

inline integer trexc( char compq, integer n, double* t, integer ldt,
                    double* q, integer ldq, integer& ifst,
                    integer& ilst)
{
    return dtrexc(compq, n, t, ldt, q, ldq, &ifst, &ilst);
}

inline integer trexc( char compq, integer n, scomplex* t, integer ldt,
                    scomplex* q, integer ldq, integer& ifst,
                    integer& ilst)
{
    return ctrexc(compq, n, t, ldt, q, ldq, &ifst, &ilst);
}

inline integer trexc( char compq, integer n, dcomplex* t, integer ldt,
                    dcomplex* q, integer ldq, integer& ifst,
                    integer& ilst)
{
    return ztrexc(compq, n, t, ldt, q, ldq, &ifst, &ilst);
}

inline integer trsen( char job, char compq, const logical* select,
                    integer n, float* t, integer ldt, float* q,
                    integer ldq, scomplex* w, integer& m,
                    float& s, float& sep)
{
    return strsen(job, compq, select, n, t, ldt, q, ldq, w, &m, &s, &sep);
}

inline integer trsen( char job, char compq, const logical* select,
                    integer n, double* t, integer ldt, double* q,
                    integer ldq, dcomplex* w, integer& m,
                    double& s, double& sep)
{
    return dtrsen(job, compq, select, n, t, ldt, q, ldq, w, &m, &s, &sep);
}

inline integer trsen( char job, char compq, const logical* select,
                    integer n, scomplex* t, integer ldt, scomplex* q,
                    integer ldq, scomplex* w, integer& m,
                    float& s, float& sep)
{
    return ctrsen(job, compq, select, n, t, ldt, q, ldq, w, &m, &s, &sep);
}

inline integer trsen( char job, char compq, const logical* select,
                    integer n, dcomplex* t, integer ldt, dcomplex* q,
                    integer ldq, dcomplex* w, integer& m,
                    double& s, double& sep)
{
    return ztrsen(job, compq, select, n, t, ldt, q, ldq, w, &m, &s, &sep);
}

inline integer trsyl( char trana, char tranb, integer isgn, integer m,
                    integer n, const float* a, integer lda,
                    const float* b, integer ldb, float* c,
                    integer ldc, float& scale )
{
    return strsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, &scale);
}

inline integer trsyl( char trana, char tranb, integer isgn, integer m,
                    integer n, const double* a, integer lda,
                    const double* b, integer ldb, double* c,
                    integer ldc, double& scale )
{
    return dtrsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, &scale);
}

inline integer trsyl( char trana, char tranb, integer isgn, integer m,
                    integer n, const scomplex* a, integer lda,
                    const scomplex* b, integer ldb, scomplex* c,
                    integer ldc, float& scale )
{
    return ctrsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, &scale);
}

inline integer trsyl( char trana, char tranb, integer isgn, integer m,
                    integer n, const dcomplex* a, integer lda,
                    const dcomplex* b, integer ldb, dcomplex* c,
                    integer ldc, double& scale )
{
    return ztrsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, &scale);
}

inline integer gghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, float* a, integer lda, float* b,
                    integer ldb, float* q, integer ldq, float* z,
                    integer ldz )
{
    return sgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
}

inline integer gghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, double* a, integer lda, double* b,
                    integer ldb, double* q, integer ldq, double* z,
                    integer ldz )
{
    return dgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
}

inline integer gghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, scomplex* a, integer lda, scomplex* b,
                    integer ldb, scomplex* q, integer ldq, scomplex* z,
                    integer ldz )
{
    return cgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
}

inline integer gghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, dcomplex* a, integer lda, dcomplex* b,
                    integer ldb, dcomplex* q, integer ldq, dcomplex* z,
                    integer ldz )
{
    return zgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
}

inline integer ggbal( char job, integer n, float* a, integer lda,
                    float* b, integer ldb, integer& ilo,
                    integer& ihi, float* lscale, float* rscale)
{
    return sggbal(job, n, a, lda, b, ldb, &ilo, &ihi, lscale, rscale);
}

inline integer ggbal( char job, integer n, double* a, integer lda,
                    double* b, integer ldb, integer& ilo,
                    integer& ihi, double* lscale, double* rscale)
{
    return dggbal(job, n, a, lda, b, ldb, &ilo, &ihi, lscale, rscale);
}

inline integer ggbal( char job, integer n, scomplex* a, integer lda,
                    scomplex* b, integer ldb, integer& ilo,
                    integer& ihi, float* lscale, float* rscale)
{
    return cggbal(job, n, a, lda, b, ldb, &ilo, &ihi, lscale, rscale);
}

inline integer ggbal( char job, integer n, dcomplex* a, integer lda,
                    dcomplex* b, integer ldb, integer& ilo,
                    integer& ihi, double* lscale, double* rscale)
{
    return zggbal(job, n, a, lda, b, ldb, &ilo, &ihi, lscale, rscale);
}

inline integer ggbak( char job, char side, integer n, integer ilo,
                    integer ihi, const float* lscale, const float* rscale,
                    integer m, float* v, integer ldv )
{
    return sggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv);
}

inline integer ggbak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* lscale, const double* rscale,
                    integer m, double* v, integer ldv )
{
    return dggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv);
}

inline integer ggbak( char job, char side, integer n, integer ilo,
                    integer ihi, const float* lscale, const float* rscale,
                    integer m, scomplex* v, integer ldv )
{
    return cggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv);
}

inline integer ggbak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* lscale, const double* rscale,
                    integer m, dcomplex* v, integer ldv )
{
    return zggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv);
}

inline integer hgeqz( char job, char compq, char compz, integer n,
                    integer ilo, integer ihi, float* h,
                    integer ldh, float* t, integer ldt, scomplex* alpha,
                    float* beta, float* q, integer ldq,
                    float* z, integer ldz)
{
    return shgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz);
}

inline integer hgeqz( char job, char compq, char compz, integer n,
                    integer ilo, integer ihi, double* h,
                    integer ldh, double* t, integer ldt, dcomplex* alpha,
                    double* beta, double* q, integer ldq,
                    double* z, integer ldz)
{
    return dhgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz);
}

inline integer hgeqz( char job, char compq, char compz, integer n,
                    integer ilo, integer ihi, scomplex* h,
                    integer ldh, scomplex* t, integer ldt, scomplex* alpha,
                    scomplex* beta, scomplex* q, integer ldq,
                    scomplex* z, integer ldz)
{
    return chgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz);
}

inline integer hgeqz( char job, char compq, char compz, integer n,
                    integer ilo, integer ihi, dcomplex* h,
                    integer ldh, dcomplex* t, integer ldt, dcomplex* alpha,
                    dcomplex* beta, dcomplex* q, integer ldq,
                    dcomplex* z, integer ldz)
{
    return zhgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz);
}

inline integer tgevc( char side, char howmny, const logical* select,
                    integer n, const float* s, integer lds,
                    const float* p, integer ldp, float* vl,
                    integer ldvl, float* vr, integer ldvr,
                    integer mm, integer& m)
{
    return stgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, &m);
}

inline integer tgevc( char side, char howmny, const logical* select,
                    integer n, const double* s, integer lds,
                    const double* p, integer ldp, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer& m)
{
    return dtgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, &m);
}

inline integer tgevc( char side, char howmny, const logical* select,
                    integer n, const scomplex* s, integer lds,
                    const scomplex* p, integer ldp, scomplex* vl,
                    integer ldvl, scomplex* vr, integer ldvr,
                    integer mm, integer& m)
{
    return ctgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, &m);
}

inline integer tgevc( char side, char howmny, const logical* select,
                    integer n, const dcomplex* s, integer lds,
                    const dcomplex* p, integer ldp, dcomplex* vl,
                    integer ldvl, dcomplex* vr, integer ldvr,
                    integer mm, integer& m)
{
    return ztgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, &m);
}

inline integer tgexc( logical wantq, logical wantz, integer n,
                    float* a, integer lda, float* b, integer ldb,
                    float* q, integer ldq, float* z, integer ldz,
                    integer& ifst, integer& ilst)
{
    return stgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, &ifst, &ilst);
}

inline integer tgexc( logical wantq, logical wantz, integer n,
                    double* a, integer lda, double* b, integer ldb,
                    double* q, integer ldq, double* z, integer ldz,
                    integer& ifst, integer& ilst)
{
    return dtgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, &ifst, &ilst);
}

inline integer tgexc( logical wantq, logical wantz, integer n,
                    scomplex* a, integer lda, scomplex* b, integer ldb,
                    scomplex* q, integer ldq, scomplex* z, integer ldz,
                    integer& ifst, integer& ilst)
{
    return ctgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, &ifst, &ilst);
}

inline integer tgexc( logical wantq, logical wantz, integer n,
                    dcomplex* a, integer lda, dcomplex* b, integer ldb,
                    dcomplex* q, integer ldq, dcomplex* z, integer ldz,
                    integer& ifst, integer& ilst)
{
    return ztgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, &ifst, &ilst);
}

inline integer tgsen( integer ijob, logical wantq,
                    logical wantz, const logical* select,
                    integer n, float* a, integer lda, float* b,
                    integer ldb, scomplex* alpha,
                    float* beta, float* q, integer ldq, float* z,
                    integer ldz, integer& m, float& pl, float& pr,
                    float* dif)
{
    return stgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, &m, &pl, &pr, dif);
}

inline integer tgsen( integer ijob, logical wantq,
                    logical wantz, const logical* select,
                    integer n, double* a, integer lda, double* b,
                    integer ldb, dcomplex* alpha,
                    double* beta, double* q, integer ldq, double* z,
                    integer ldz, integer& m, double& pl, double& pr,
                    double* dif)
{
    return dtgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, &m, &pl, &pr, dif);
}

inline integer tgsen( integer ijob, logical wantq,
                    logical wantz, const logical* select,
                    integer n, scomplex* a, integer lda, scomplex* b,
                    integer ldb, scomplex* alpha,
                    scomplex* beta, scomplex* q, integer ldq, scomplex* z,
                    integer ldz, integer& m, float& pl, float& pr,
                    float* dif)
{
    return ctgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, &m, &pl, &pr, dif);
}

inline integer tgsen( integer ijob, logical wantq,
                    logical wantz, const logical* select,
                    integer n, dcomplex* a, integer lda, dcomplex* b,
                    integer ldb, dcomplex* alpha,
                    dcomplex* beta, dcomplex* q, integer ldq, dcomplex* z,
                    integer ldz, integer& m, double& pl, double& pr,
                    double* dif)
{
    return ztgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, &m, &pl, &pr, dif);
}

inline integer tgsyl( char trans, integer ijob, integer m, integer n,
                    const float* a, integer lda, const float* b,
                    integer ldb, float* c, integer ldc,
                    const float* d, integer ldd, const float* e,
                    integer lde, float* f, integer ldf, float& scale,
                    float& dif)
{
    return stgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, &scale, &dif);
}

inline integer tgsyl( char trans, integer ijob, integer m, integer n,
                    const double* a, integer lda, const double* b,
                    integer ldb, double* c, integer ldc,
                    const double* d, integer ldd, const double* e,
                    integer lde, double* f, integer ldf, double& scale,
                    double& dif)
{
    return dtgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, &scale, &dif);
}

inline integer tgsyl( char trans, integer ijob, integer m, integer n,
                    const scomplex* a, integer lda, const scomplex* b,
                    integer ldb, scomplex* c, integer ldc,
                    const scomplex* d, integer ldd, const scomplex* e,
                    integer lde, scomplex* f, integer ldf, float& scale,
                    float& dif)
{
    return ctgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, &scale, &dif);
}

inline integer tgsyl( char trans, integer ijob, integer m, integer n,
                    const dcomplex* a, integer lda, const dcomplex* b,
                    integer ldb, dcomplex* c, integer ldc,
                    const dcomplex* d, integer ldd, const dcomplex* e,
                    integer lde, dcomplex* f, integer ldf, double& scale,
                    double& dif)
{
    return ztgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, &scale, &dif);
}

inline integer tgsna( char job, char howmny, const logical* select,
                    integer n, const float* a, integer lda,
                    const float* b, integer ldb, const float* vl,
                    integer ldvl, const float* vr, integer ldvr,
                    float* s, float* dif, integer mm, integer& m)
{
    return stgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, &m);
}

inline integer tgsna( char job, char howmny, const logical* select,
                    integer n, const double* a, integer lda,
                    const double* b, integer ldb, const double* vl,
                    integer ldvl, const double* vr, integer ldvr,
                    double* s, double* dif, integer mm, integer& m)
{
    return dtgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, &m);
}

inline integer tgsna( char job, char howmny, const logical* select,
                    integer n, const scomplex* a, integer lda,
                    const scomplex* b, integer ldb, const scomplex* vl,
                    integer ldvl, const scomplex* vr, integer ldvr,
                    float* s, float* dif, integer mm, integer& m)
{
    return ctgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, &m);
}

inline integer tgsna( char job, char howmny, const logical* select,
                    integer n, const dcomplex* a, integer lda,
                    const dcomplex* b, integer ldb, const dcomplex* vl,
                    integer ldvl, const dcomplex* vr, integer ldvr,
                    double* s, double* dif, integer mm, integer& m)
{
    return ztgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, &m);
}

inline integer ggsvp( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, float* a, integer lda,
                    float* b, integer ldb, float tola, float tolb,
                    integer* k, integer* l, float* u, integer ldu,
                    float* v, integer ldv, float* q, integer ldq)
{
    return sggsvp(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq);
}

inline integer ggsvp( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, double* a, integer lda,
                    double* b, integer ldb, double tola, double tolb,
                    integer* k, integer* l, double* u, integer ldu,
                    double* v, integer ldv, double* q, integer ldq)
{
    return dggsvp(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq);
}

inline integer ggsvp( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, scomplex* a, integer lda,
                    scomplex* b, integer ldb, float tola, float tolb,
                    integer* k, integer* l, scomplex* u, integer ldu,
                    scomplex* v, integer ldv, scomplex* q, integer ldq)
{
    return cggsvp(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq);
}

inline integer ggsvp( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, dcomplex* a, integer lda,
                    dcomplex* b, integer ldb, double tola, double tolb,
                    integer* k, integer* l, dcomplex* u, integer ldu,
                    dcomplex* v, integer ldv, dcomplex* q, integer ldq)
{
    return zggsvp(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq);
}

inline integer tgsja( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, integer k, integer l,
                    float* a, integer lda, float* b, integer ldb,
                    float tola, float tolb, float* alpha, float* beta,
                    float* u, integer ldu, float* v, integer ldv,
                    float* q, integer ldq, integer& ncycle )
{
    return stgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, &ncycle);
}

inline integer tgsja( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, integer k, integer l,
                    double* a, integer lda, double* b, integer ldb,
                    double tola, double tolb, double* alpha, double* beta,
                    double* u, integer ldu, double* v, integer ldv,
                    double* q, integer ldq, integer& ncycle )
{
    return dtgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, &ncycle);
}

inline integer tgsja( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, integer k, integer l,
                    scomplex* a, integer lda, scomplex* b, integer ldb,
                    float tola, float tolb, float* alpha, float* beta,
                    scomplex* u, integer ldu, scomplex* v, integer ldv,
                    scomplex* q, integer ldq, integer& ncycle )
{
    return ctgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, &ncycle);
}

inline integer tgsja( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, integer k, integer l,
                    dcomplex* a, integer lda, dcomplex* b, integer ldb,
                    double tola, double tolb, double* alpha, double* beta,
                    dcomplex* u, integer ldu, dcomplex* v, integer ldv,
                    dcomplex* q, integer ldq, integer& ncycle )
{
    return ztgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, &ncycle);
}

inline integer gels( char trans, integer m, integer n, integer nrhs,
                   float* a, integer lda, float* b, integer ldb)
{
    return sgels(trans, m, n, nrhs, a, lda, b, ldb);
}

inline integer gels( char trans, integer m, integer n, integer nrhs,
                   double* a, integer lda, double* b, integer ldb)
{
    return dgels(trans, m, n, nrhs, a, lda, b, ldb);
}

inline integer gels( char trans, integer m, integer n, integer nrhs,
                   scomplex* a, integer lda, scomplex* b, integer ldb)
{
    return cgels(trans, m, n, nrhs, a, lda, b, ldb);
}

inline integer gels( char trans, integer m, integer n, integer nrhs,
                   dcomplex* a, integer lda, dcomplex* b, integer ldb)
{
    return zgels(trans, m, n, nrhs, a, lda, b, ldb);
}

inline integer gelsy( integer m, integer n, integer nrhs, float* a,
                    integer lda, float* b, integer ldb,
                    integer* jpvt, float rcond, integer& rank)
{
    return sgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, &rank);
}

inline integer gelsy( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb,
                    integer* jpvt, double rcond, integer& rank)
{
    return dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, &rank);
}

inline integer gelsy( integer m, integer n, integer nrhs, scomplex* a,
                    integer lda, scomplex* b, integer ldb,
                    integer* jpvt, float rcond, integer& rank)
{
    return cgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, &rank);
}

inline integer gelsy( integer m, integer n, integer nrhs, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb,
                    integer* jpvt, double rcond, integer& rank)
{
    return zgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, &rank);
}

inline integer gelss( integer m, integer n, integer nrhs, float* a,
                    integer lda, float* b, integer ldb, float* s,
                    float rcond, integer& rank)
{
    return sgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, &rank);
}

inline integer gelss( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* s,
                    double rcond, integer& rank)
{
    return dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, &rank);
}

inline integer gelss( integer m, integer n, integer nrhs, scomplex* a,
                    integer lda, scomplex* b, integer ldb, float* s,
                    float rcond, integer& rank)
{
    return cgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, &rank);
}

inline integer gelss( integer m, integer n, integer nrhs, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb, double* s,
                    double rcond, integer& rank)
{
    return zgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, &rank);
}

inline integer gglse( integer m, integer n, integer p, float* a,
                    integer lda, float* b, integer ldb, float* c,
                    float* d, float* x)
{
    return sgglse(m, n, p, a, lda, b, ldb, c, d, x);
}

inline integer gglse( integer m, integer n, integer p, double* a,
                    integer lda, double* b, integer ldb, double* c,
                    double* d, double* x)
{
    return dgglse(m, n, p, a, lda, b, ldb, c, d, x);
}

inline integer gglse( integer m, integer n, integer p, scomplex* a,
                    integer lda, scomplex* b, integer ldb, scomplex* c,
                    scomplex* d, scomplex* x)
{
    return cgglse(m, n, p, a, lda, b, ldb, c, d, x);
}

inline integer gglse( integer m, integer n, integer p, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb, dcomplex* c,
                    dcomplex* d, dcomplex* x)
{
    return zgglse(m, n, p, a, lda, b, ldb, c, d, x);
}

inline integer ggglm( integer n, integer m, integer p, float* a,
                    integer lda, float* b, integer ldb, float* d,
                    float* x, float* y)
{
    return sggglm(n, m, p, a, lda, b, ldb, d, x, y);
}

inline integer ggglm( integer n, integer m, integer p, double* a,
                    integer lda, double* b, integer ldb, double* d,
                    double* x, double* y)
{
    return dggglm(n, m, p, a, lda, b, ldb, d, x, y);
}

inline integer ggglm( integer n, integer m, integer p, scomplex* a,
                    integer lda, scomplex* b, integer ldb, scomplex* d,
                    scomplex* x, scomplex* y)
{
    return cggglm(n, m, p, a, lda, b, ldb, d, x, y);
}

inline integer ggglm( integer n, integer m, integer p, dcomplex* a,
                    integer lda, dcomplex* b, integer ldb, dcomplex* d,
                    dcomplex* x, dcomplex* y)
{
    return zggglm(n, m, p, a, lda, b, ldb, d, x, y);
}

inline integer heev( char jobz, char uplo, integer n, float* a,
                   integer lda, float* w)
{
    return ssyev(jobz, uplo, n, a, lda, w);
}

inline integer heev( char jobz, char uplo, integer n, double* a,
                   integer lda, double* w)
{
    return dsyev(jobz, uplo, n, a, lda, w);
}

inline integer heev( char jobz, char uplo, integer n, scomplex* a,
                   integer lda, float* w)
{
    return cheev(jobz, uplo, n, a, lda, w);
}

inline integer heev( char jobz, char uplo, integer n, dcomplex* a,
                   integer lda, double* w)
{
    return zheev(jobz, uplo, n, a, lda, w);
}

inline integer heevd( char jobz, char uplo, integer n, float* a,
                    integer lda, float* w)
{
    return ssyevd(jobz, uplo, n, a, lda, w);
}

inline integer heevd( char jobz, char uplo, integer n, double* a,
                    integer lda, double* w)
{
    return dsyevd(jobz, uplo, n, a, lda, w);
}

inline integer heevd( char jobz, char uplo, integer n, scomplex* a,
                    integer lda, float* w)
{
    return cheevd(jobz, uplo, n, a, lda, w);
}

inline integer heevd( char jobz, char uplo, integer n, dcomplex* a,
                    integer lda, double* w)
{
    return zheevd(jobz, uplo, n, a, lda, w);
}

inline integer heevx( char jobz, char range, char uplo, integer n,
                    float* a, integer lda, float vl, float vu,
                    integer il, integer iu, float abstol,
                    integer& m, float* w, float* z, integer ldz,
                    integer* ifail )
{
    return ssyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer heevx( char jobz, char range, char uplo, integer n,
                    double* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer& m, double* w, double* z, integer ldz,
                    integer* ifail )
{
    return dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer heevx( char jobz, char range, char uplo, integer n,
                    scomplex* a, integer lda, float vl, float vu,
                    integer il, integer iu, float abstol,
                    integer& m, float* w, scomplex* z, integer ldz,
                    integer* ifail )
{
    return cheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer heevx( char jobz, char range, char uplo, integer n,
                    dcomplex* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer& m, double* w, dcomplex* z, integer ldz,
                    integer* ifail )
{
    return zheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer heevr( char jobz, char range, char uplo, integer n,
                    float* a, integer lda, float vl, float vu,
                    integer il, integer iu, float abstol,
                    integer& m, float* w, float* z, integer ldz,
                    integer* isuppz)
{
    return ssyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer heevr( char jobz, char range, char uplo, integer n,
                    double* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer& m, double* w, double* z, integer ldz,
                    integer* isuppz)
{
    return dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer heevr( char jobz, char range, char uplo, integer n,
                    scomplex* a, integer lda, float vl, float vu,
                    integer il, integer iu, float abstol,
                    integer& m, float* w, scomplex* z, integer ldz,
                    integer* isuppz)
{
    return cheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer heevr( char jobz, char range, char uplo, integer n,
                    dcomplex* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer& m, double* w, dcomplex* z, integer ldz,
                    integer* isuppz)
{
    return zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer hpev( char jobz, char uplo, integer n, float* ap, float* w,
                   float* z, integer ldz)
{
    return sspev(jobz, uplo, n, ap, w, z, ldz);
}

inline integer hpev( char jobz, char uplo, integer n, double* ap, double* w,
                   double* z, integer ldz)
{
    return dspev(jobz, uplo, n, ap, w, z, ldz);
}

inline integer hpev( char jobz, char uplo, integer n, scomplex* ap, float* w,
                   scomplex* z, integer ldz)
{
    return chpev(jobz, uplo, n, ap, w, z, ldz);
}

inline integer hpev( char jobz, char uplo, integer n, dcomplex* ap, double* w,
                   dcomplex* z, integer ldz)
{
    return zhpev(jobz, uplo, n, ap, w, z, ldz);
}

inline integer hpevd( char jobz, char uplo, integer n, float* ap,
                    float* w, float* z, integer ldz)
{
    return sspevd(jobz, uplo, n, ap, w, z, ldz);
}

inline integer hpevd( char jobz, char uplo, integer n, double* ap,
                    double* w, double* z, integer ldz)
{
    return dspevd(jobz, uplo, n, ap, w, z, ldz);
}

inline integer hpevd( char jobz, char uplo, integer n, scomplex* ap,
                    float* w, scomplex* z, integer ldz)
{
    return chpevd(jobz, uplo, n, ap, w, z, ldz);
}

inline integer hpevd( char jobz, char uplo, integer n, dcomplex* ap,
                    double* w, dcomplex* z, integer ldz)
{
    return zhpevd(jobz, uplo, n, ap, w, z, ldz);
}

inline integer hpevx( char jobz, char range, char uplo, integer n,
                    float* ap, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    return sspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hpevx( char jobz, char range, char uplo, integer n,
                    double* ap, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    return dspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hpevx( char jobz, char range, char uplo, integer n,
                    scomplex* ap, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    scomplex* z, integer ldz, integer* ifail )
{
    return chpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hpevx( char jobz, char range, char uplo, integer n,
                    dcomplex* ap, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    dcomplex* z, integer ldz, integer* ifail )
{
    return zhpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hbev( char jobz, char uplo, integer n, integer kd,
                   float* ab, integer ldab, float* w, float* z,
                   integer ldz)
{
    return ssbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz);
}

inline integer hbev( char jobz, char uplo, integer n, integer kd,
                   double* ab, integer ldab, double* w, double* z,
                   integer ldz)
{
    return dsbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz);
}

inline integer hbev( char jobz, char uplo, integer n, integer kd,
                   scomplex* ab, integer ldab, float* w, scomplex* z,
                   integer ldz)
{
    return chbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz);
}

inline integer hbev( char jobz, char uplo, integer n, integer kd,
                   dcomplex* ab, integer ldab, double* w, dcomplex* z,
                   integer ldz)
{
    return zhbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz);
}

inline integer hbevd( char jobz, char uplo, integer n, integer kd,
                    float* ab, integer ldab, float* w, float* z,
                    integer ldz)
{
    return ssbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz);
}

inline integer hbevd( char jobz, char uplo, integer n, integer kd,
                    double* ab, integer ldab, double* w, double* z,
                    integer ldz)
{
    return dsbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz);
}

inline integer hbevd( char jobz, char uplo, integer n, integer kd,
                    scomplex* ab, integer ldab, float* w, scomplex* z,
                    integer ldz)
{
    return chbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz);
}

inline integer hbevd( char jobz, char uplo, integer n, integer kd,
                    dcomplex* ab, integer ldab, double* w, dcomplex* z,
                    integer ldz)
{
    return zhbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz);
}

inline integer hbevx( char jobz, char range, char uplo, integer n,
                    integer kd, float* ab, integer ldab, float* q,
                    integer ldq, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    return ssbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hbevx( char jobz, char range, char uplo, integer n,
                    integer kd, double* ab, integer ldab, double* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    return dsbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hbevx( char jobz, char range, char uplo, integer n,
                    integer kd, scomplex* ab, integer ldab, scomplex* q,
                    integer ldq, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    scomplex* z, integer ldz, integer* ifail )
{
    return chbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hbevx( char jobz, char range, char uplo, integer n,
                    integer kd, dcomplex* ab, integer ldab, dcomplex* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    dcomplex* z, integer ldz, integer* ifail )
{
    return zhbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer stev( char jobz, integer n, float* d, float* e, float* z,
                   integer ldz)
{
    return sstev(jobz, n, d, e, z, ldz);
}

inline integer stev( char jobz, integer n, double* d, double* e, double* z,
                   integer ldz)
{
    return dstev(jobz, n, d, e, z, ldz);
}

inline integer stevd( char jobz, integer n, float* d, float* e, float* z,
                    integer ldz)
{
    return sstevd(jobz, n, d, e, z, ldz);
}

inline integer stevd( char jobz, integer n, double* d, double* e, double* z,
                    integer ldz)
{
    return dstevd(jobz, n, d, e, z, ldz);
}

inline integer stevx( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    return sstevx(jobz, range, n, d, e, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer stevx( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    return dstevx(jobz, range, n, d, e, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer stevr( char jobz, char range, integer n, float* d,
                    float* e, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    float* z, integer ldz, integer* isuppz)
{
    return sstevr(jobz, range, n, d, e, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer stevr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    double* z, integer ldz, integer* isuppz)
{
    return dstevr(jobz, range, n, d, e, vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
}

inline integer geev( char jobvl, char jobvr, integer n, float* a,
                   integer lda, scomplex* w, float* vl,
                   integer ldvl, float* vr, integer ldvr)
{
    return sgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
}

inline integer geev( char jobvl, char jobvr, integer n, double* a,
                   integer lda, dcomplex* w, double* vl,
                   integer ldvl, double* vr, integer ldvr)
{
    return dgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
}

inline integer geev( char jobvl, char jobvr, integer n, scomplex* a,
                   integer lda, scomplex* w, scomplex* vl,
                   integer ldvl, scomplex* vr, integer ldvr)
{
    return cgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
}

inline integer geev( char jobvl, char jobvr, integer n, dcomplex* a,
                   integer lda, dcomplex* w, dcomplex* vl,
                   integer ldvl, dcomplex* vr, integer ldvr)
{
    return zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr);
}

inline integer geevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, float* a, integer lda, scomplex* w,
                    float* vl, integer ldvl, float* vr,
                    integer ldvr, integer& ilo, integer& ihi,
                    float* scale, float& abnrm, float* rconde,
                    float* rcondv)
{
    return sgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv);
}

inline integer geevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, double* a, integer lda, dcomplex* w,
                    double* vl, integer ldvl, double* vr,
                    integer ldvr, integer& ilo, integer& ihi,
                    double* scale, double& abnrm, double* rconde,
                    double* rcondv)
{
    return dgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv);
}

inline integer geevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, scomplex* a, integer lda, scomplex* w,
                    scomplex* vl, integer ldvl, scomplex* vr,
                    integer ldvr, integer& ilo, integer& ihi,
                    float* scale, float& abnrm, float* rconde,
                    float* rcondv)
{
    return cgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv);
}

inline integer geevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, dcomplex* a, integer lda, dcomplex* w,
                    dcomplex* vl, integer ldvl, dcomplex* vr,
                    integer ldvr, integer& ilo, integer& ihi,
                    double* scale, double& abnrm, double* rconde,
                    double* rcondv)
{
    return zgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv);
}

inline integer gesvd( char jobu, char jobvt, integer m, integer n,
                    float* a, integer lda, float* s, float* u,
                    integer ldu, float* vt, integer ldvt)
{
    return sgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt);
}

inline integer gesvd( char jobu, char jobvt, integer m, integer n,
                    double* a, integer lda, double* s, double* u,
                    integer ldu, double* vt, integer ldvt)
{
    return dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt);
}

inline integer gesvd( char jobu, char jobvt, integer m, integer n,
                    scomplex* a, integer lda, float* s, scomplex* u,
                    integer ldu, scomplex* vt, integer ldvt)
{
    return cgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt);
}

inline integer gesvd( char jobu, char jobvt, integer m, integer n,
                    dcomplex* a, integer lda, double* s, dcomplex* u,
                    integer ldu, dcomplex* vt, integer ldvt)
{
    return zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt);
}

inline integer gesdd( char jobz, integer m, integer n, float* a,
                    integer lda, float* s, float* u, integer ldu,
                    float* vt, integer ldvt)
{
    return sgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}

inline integer gesdd( char jobz, integer m, integer n, double* a,
                    integer lda, double* s, double* u, integer ldu,
                    double* vt, integer ldvt)
{
    return dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}

inline integer gesdd( char jobz, integer m, integer n, scomplex* a,
                    integer lda, float* s, scomplex* u, integer ldu,
                    scomplex* vt, integer ldvt)
{
    return cgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}

inline integer gesdd( char jobz, integer m, integer n, dcomplex* a,
                    integer lda, double* s, dcomplex* u, integer ldu,
                    dcomplex* vt, integer ldvt)
{
    return zgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
}

inline integer gejsv( char joba, char jobu, char jobv, char jobr, char jobt,
                    char jobp, integer m, integer n, float* a,
                    integer lda, float* sva, float* u, integer ldu,
                    float* v, integer ldv)
{
    return sgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv);
}

inline integer gejsv( char joba, char jobu, char jobv, char jobr, char jobt,
                    char jobp, integer m, integer n, double* a,
                    integer lda, double* sva, double* u, integer ldu,
                    double* v, integer ldv)
{
    return dgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv);
}

inline integer gesvj( char joba, char jobu, char jobv, integer m,
                    integer n, float* a, integer lda, float* sva,
                    integer mv, float* v, integer ldv)
{
    return sgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv);
}

inline integer gesvj( char joba, char jobu, char jobv, integer m,
                    integer n, double* a, integer lda, double* sva,
                    integer mv, double* v, integer ldv)
{
    return dgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv);
}

inline integer ggsvd( char jobu, char jobv, char jobq, integer m,
                    integer n, integer p, integer& k, integer& l,
                    float* a, integer lda, float* b, integer ldb,
                    float* alpha, float* beta, float* u, integer ldu,
                    float* v, integer ldv, float* q, integer ldq)
{
    return sggsvd(jobu, jobv, jobq, m, n, p, &k, &l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq);
}

inline integer ggsvd( char jobu, char jobv, char jobq, integer m,
                    integer n, integer p, integer& k, integer& l,
                    double* a, integer lda, double* b, integer ldb,
                    double* alpha, double* beta, double* u, integer ldu,
                    double* v, integer ldv, double* q, integer ldq)
{
    return dggsvd(jobu, jobv, jobq, m, n, p, &k, &l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq);
}

inline integer ggsvd( char jobu, char jobv, char jobq, integer m,
                    integer n, integer p, integer& k, integer& l,
                    scomplex* a, integer lda, scomplex* b, integer ldb,
                    float* alpha, float* beta, scomplex* u, integer ldu,
                    scomplex* v, integer ldv, scomplex* q, integer ldq)
{
    return cggsvd(jobu, jobv, jobq, m, n, p, &k, &l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq);
}

inline integer ggsvd( char jobu, char jobv, char jobq, integer m,
                    integer n, integer p, integer& k, integer& l,
                    dcomplex* a, integer lda, dcomplex* b, integer ldb,
                    double* alpha, double* beta, dcomplex* u, integer ldu,
                    dcomplex* v, integer ldv, dcomplex* q, integer ldq)
{
    return zggsvd(jobu, jobv, jobq, m, n, p, &k, &l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq);
}

inline integer hegv( integer itype, char jobz, char uplo, integer n,
                   float* a, integer lda, float* b, integer ldb,
                   float* w)
{
    return ssygv(itype, jobz, uplo, n, a, lda, b, ldb, w);
}

inline integer hegv( integer itype, char jobz, char uplo, integer n,
                   double* a, integer lda, double* b, integer ldb,
                   double* w)
{
    return dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w);
}

inline integer hegv( integer itype, char jobz, char uplo, integer n,
                   scomplex* a, integer lda, scomplex* b, integer ldb,
                   float* w)
{
    return chegv(itype, jobz, uplo, n, a, lda, b, ldb, w);
}

inline integer hegv( integer itype, char jobz, char uplo, integer n,
                   dcomplex* a, integer lda, dcomplex* b, integer ldb,
                   double* w)
{
    return zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w);
}

inline integer hegvd( integer itype, char jobz, char uplo, integer n,
                    float* a, integer lda, float* b, integer ldb,
                    float* w)
{
    return ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w);
}

inline integer hegvd( integer itype, char jobz, char uplo, integer n,
                    double* a, integer lda, double* b, integer ldb,
                    double* w)
{
    return dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w);
}

inline integer hegvd( integer itype, char jobz, char uplo, integer n,
                    scomplex* a, integer lda, scomplex* b, integer ldb,
                    float* w)
{
    return chegvd(itype, jobz, uplo, n, a, lda, b, ldb, w);
}

inline integer hegvd( integer itype, char jobz, char uplo, integer n,
                    dcomplex* a, integer lda, dcomplex* b, integer ldb,
                    double* w)
{
    return zhegvd(itype, jobz, uplo, n, a, lda, b, ldb, w);
}

inline integer hegvx( integer itype, char jobz, char range, char uplo,
                    integer n, float* a, integer lda, float* b,
                    integer ldb, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    return ssygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hegvx( integer itype, char jobz, char range, char uplo,
                    integer n, double* a, integer lda, double* b,
                    integer ldb, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    return dsygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hegvx( integer itype, char jobz, char range, char uplo,
                    integer n, scomplex* a, integer lda, scomplex* b,
                    integer ldb, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    scomplex* z, integer ldz, integer* ifail )
{
    return chegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hegvx( integer itype, char jobz, char range, char uplo,
                    integer n, dcomplex* a, integer lda, dcomplex* b,
                    integer ldb, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    dcomplex* z, integer ldz, integer* ifail )
{
    return zhegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hpgv( integer itype, char jobz, char uplo, integer n,
                    float* ap, float* bp, float* w, float* z,
                   integer ldz)
{
    return sspgv(itype, jobz, uplo, n, ap, bp, w, z, ldz);
}

inline integer hpgv( integer itype, char jobz, char uplo, integer n,
                   double* ap, double* bp, double* w, double* z,
                   integer ldz)
{
    return dspgv(itype, jobz, uplo, n, ap, bp, w, z, ldz);
}

inline integer hpgv( integer itype, char jobz, char uplo, integer n,
                   scomplex* ap, scomplex* bp, float* w, scomplex* z,
                   integer ldz)
{
    return chpgv(itype, jobz, uplo, n, ap, bp, w, z, ldz);
}

inline integer hpgv( integer itype, char jobz, char uplo, integer n,
                   dcomplex* ap, dcomplex* bp, double* w, dcomplex* z,
                   integer ldz)
{
    return zhpgv(itype, jobz, uplo, n, ap, bp, w, z, ldz);
}

inline integer hpgvd( integer itype, char jobz, char uplo, integer n,
                    float* ap, float* bp, float* w, float* z,
                    integer ldz)
{
    return sspgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz);
}

inline integer hpgvd( integer itype, char jobz, char uplo, integer n,
                    double* ap, double* bp, double* w, double* z,
                    integer ldz)
{
    return dspgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz);
}

inline integer hpgvd( integer itype, char jobz, char uplo, integer n,
                    scomplex* ap, scomplex* bp, float* w, scomplex* z,
                    integer ldz)
{
    return chpgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz);
}

inline integer hpgvd( integer itype, char jobz, char uplo, integer n,
                    dcomplex* ap, dcomplex* bp, double* w, dcomplex* z,
                    integer ldz)
{
    return zhpgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz);
}

inline integer hpgvx( integer itype, char jobz, char range, char uplo,
                    integer n, float* ap, float* bp, float vl,
                    float vu, integer il, integer iu, float abstol,
                    integer& m, float* w, float* z, integer ldz,
                    integer* ifail )
{
    return sspgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hpgvx( integer itype, char jobz, char range, char uplo,
                    integer n, double* ap, double* bp, double vl,
                    double vu, integer il, integer iu, double abstol,
                    integer& m, double* w, double* z, integer ldz,
                    integer* ifail )
{
    return dspgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hpgvx( integer itype, char jobz, char range, char uplo,
                    integer n, scomplex* ap, scomplex* bp, float vl,
                    float vu, integer il, integer iu, float abstol,
                    integer& m, float* w, scomplex* z, integer ldz,
                    integer* ifail )
{
    return chpgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hpgvx( integer itype, char jobz, char range, char uplo,
                    integer n, dcomplex* ap, dcomplex* bp, double vl,
                    double vu, integer il, integer iu, double abstol,
                    integer& m, double* w, dcomplex* z, integer ldz,
                    integer* ifail )
{
    return zhpgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, float* ab, integer ldab, float* bb,
                   integer ldbb, float* w, float* z, integer ldz)
{
    return ssbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz);
}

inline integer hbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, double* ab, integer ldab, double* bb,
                   integer ldbb, double* w, double* z, integer ldz)
{
    return dsbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz);
}

inline integer hbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, scomplex* ab, integer ldab, scomplex* bb,
                   integer ldbb, float* w, scomplex* z, integer ldz)
{
    return chbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz);
}

inline integer hbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, dcomplex* ab, integer ldab, dcomplex* bb,
                   integer ldbb, double* w, dcomplex* z, integer ldz)
{
    return zhbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz);
}

inline integer hbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, float* ab, integer ldab, float* bb,
                    integer ldbb, float* w, float* z, integer ldz)
{
    return ssbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz);
}

inline integer hbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, double* ab, integer ldab, double* bb,
                    integer ldbb, double* w, double* z, integer ldz)
{
    return dsbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz);
}

inline integer hbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, scomplex* ab, integer ldab, scomplex* bb,
                    integer ldbb, float* w, scomplex* z, integer ldz)
{
    return chbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz);
}

inline integer hbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, dcomplex* ab, integer ldab, dcomplex* bb,
                    integer ldbb, double* w, dcomplex* z, integer ldz)
{
    return zhbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz);
}

inline integer hbgvx( char jobz, char range, char uplo, integer n,
                    integer ka, integer kb, float* ab,
                    integer ldab, float* bb, integer ldbb, float* q,
                    integer ldq, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    float* z, integer ldz, integer* ifail )
{
    return ssbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hbgvx( char jobz, char range, char uplo, integer n,
                    integer ka, integer kb, double* ab,
                    integer ldab, double* bb, integer ldbb, double* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    double* z, integer ldz, integer* ifail )
{
    return dsbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hbgvx( char jobz, char range, char uplo, integer n,
                    integer ka, integer kb, scomplex* ab,
                    integer ldab, scomplex* bb, integer ldbb, scomplex* q,
                    integer ldq, float vl, float vu, integer il,
                    integer iu, float abstol, integer& m, float* w,
                    scomplex* z, integer ldz, integer* ifail )
{
    return chbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer hbgvx( char jobz, char range, char uplo, integer n,
                    integer ka, integer kb, dcomplex* ab,
                    integer ldab, dcomplex* bb, integer ldbb, dcomplex* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer& m, double* w,
                    dcomplex* z, integer ldz, integer* ifail )
{
    return zhbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
}

inline integer ggev( char jobvl, char jobvr, integer n, float* a,
                   integer lda, float* b, integer ldb, scomplex* alpha,
                   float* beta, float* vl, integer ldvl,
                   float* vr, integer ldvr)
{
    return sggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr);
}

inline integer ggev( char jobvl, char jobvr, integer n, double* a,
                   integer lda, double* b, integer ldb, dcomplex* alpha,
                   double* beta, double* vl, integer ldvl,
                   double* vr, integer ldvr)
{
    return dggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr);
}

inline integer ggev( char jobvl, char jobvr, integer n, scomplex* a,
                   integer lda, scomplex* b, integer ldb, scomplex* alpha,
                   scomplex* beta, scomplex* vl, integer ldvl,
                   scomplex* vr, integer ldvr)
{
    return cggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr);
}

inline integer ggev( char jobvl, char jobvr, integer n, dcomplex* a,
                   integer lda, dcomplex* b, integer ldb, dcomplex* alpha,
                   dcomplex* beta, dcomplex* vl, integer ldvl,
                   dcomplex* vr, integer ldvr)
{
    return zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr);
}

inline integer ggevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, float* a, integer lda, float* b,
                    integer ldb, scomplex* alpha,
                    float* beta, float* vl, integer ldvl, float* vr,
                    integer ldvr, integer& ilo, integer& ihi,
                    float* lscale, float* rscale, float& abnrm,
                    float& bbnrm, float* rconde, float* rcondv)
{
    return sggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr,
                  &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm, rconde, rcondv);
}

inline integer ggevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, double* a, integer lda, double* b,
                    integer ldb, dcomplex* alpha,
                    double* beta, double* vl, integer ldvl, double* vr,
                    integer ldvr, integer& ilo, integer& ihi,
                    double* lscale, double* rscale, double& abnrm,
                    double& bbnrm, double* rconde, double* rcondv)
{
    return dggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr,
                  &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm, rconde, rcondv);
}

inline integer ggevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, scomplex* a, integer lda, scomplex* b,
                    integer ldb, scomplex* alpha,
                    scomplex* beta, scomplex* vl, integer ldvl, scomplex* vr,
                    integer ldvr, integer& ilo, integer& ihi,
                    float* lscale, float* rscale, float& abnrm,
                    float& bbnrm, float* rconde, float* rcondv)
{
    return cggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr,
                  &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm, rconde, rcondv);
}

inline integer ggevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, dcomplex* a, integer lda, dcomplex* b,
                    integer ldb, dcomplex* alpha,
                    dcomplex* beta, dcomplex* vl, integer ldvl, dcomplex* vr,
                    integer ldvr, integer& ilo, integer& ihi,
                    double* lscale, double* rscale, double& abnrm,
                    double& bbnrm, double* rconde, double* rcondv)
{
    return zggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr,
                  &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm, rconde, rcondv);
}

inline void hfrk( char transr, char uplo, char trans, integer n,
                   integer k, float alpha, const float* a,
                   integer lda, float beta, float* c )
{
    ssfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}

inline void hfrk( char transr, char uplo, char trans, integer n,
                   integer k, double alpha, const double* a,
                   integer lda, double beta, double* c )
{
    dsfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}

inline void hfrk( char transr, char uplo, char trans, integer n,
                   integer k, float alpha, const scomplex* a,
                   integer lda, float beta, scomplex* c )
{
    chfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}

inline void hfrk( char transr, char uplo, char trans, integer n,
                   integer k, double alpha, const dcomplex* a,
                   integer lda, double beta, dcomplex* c )
{
    zhfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}

inline void tfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, float alpha,
                   const float* a, float* b, integer ldb )
{
    stfsm(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}

inline void tfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, double alpha,
                   const double* a, double* b, integer ldb )
{
    dtfsm(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}

inline void tfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, scomplex alpha,
                   const scomplex* a, scomplex* b, integer ldb )
{
    ctfsm(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}

inline void tfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, dcomplex alpha,
                   const dcomplex* a, dcomplex* b, integer ldb )
{
    ztfsm(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}

inline integer tfttp( char transr, char uplo, integer n, const float* arf,
                    float* ap )
{
    return stfttp(transr, uplo, n, arf, ap);
}

inline integer tfttp( char transr, char uplo, integer n, const double* arf,
                    double* ap )
{
    return dtfttp(transr, uplo, n, arf, ap);
}

inline integer tfttp( char transr, char uplo, integer n, const scomplex* arf,
                    scomplex* ap )
{
    return ctfttp(transr, uplo, n, arf, ap);
}

inline integer tfttp( char transr, char uplo, integer n, const dcomplex* arf,
                    dcomplex* ap )
{
    return ztfttp(transr, uplo, n, arf, ap);
}

inline integer tfttr( char transr, char uplo, integer n, const float* arf,
                    float* a, integer lda )
{
    return stfttr(transr, uplo, n, arf, a, lda);
}

inline integer tfttr( char transr, char uplo, integer n, const double* arf,
                    double* a, integer lda )
{
    return dtfttr(transr, uplo, n, arf, a, lda);
}

inline integer tfttr( char transr, char uplo, integer n, const scomplex* arf,
                    scomplex* a, integer lda )
{
    return ctfttr(transr, uplo, n, arf, a, lda);
}

inline integer tfttr( char transr, char uplo, integer n, const dcomplex* arf,
                    dcomplex* a, integer lda )
{
    return ztfttr(transr, uplo, n, arf, a, lda);
}

inline integer tpttf( char transr, char uplo, integer n, const float* ap,
                    float* arf )
{
    return stpttf(transr, uplo, n, ap, arf);
}

inline integer tpttf( char transr, char uplo, integer n, const double* ap,
                    double* arf )
{
    return dtpttf(transr, uplo, n, ap, arf);
}

inline integer tpttf( char transr, char uplo, integer n, const scomplex* ap,
                    scomplex* arf )
{
    return ctpttf(transr, uplo, n, ap, arf);
}

inline integer tpttf( char transr, char uplo, integer n, const dcomplex* ap,
                    dcomplex* arf )
{
    return ztpttf(transr, uplo, n, ap, arf);
}

inline integer tpttr( char uplo, integer n, const float* ap, float* a,
                    integer lda )
{
    return stpttr(uplo, n, ap, a, lda);
}

inline integer tpttr( char uplo, integer n, const double* ap, double* a,
                    integer lda )
{
    return dtpttr(uplo, n, ap, a, lda);
}

inline integer tpttr( char uplo, integer n, const scomplex* ap, scomplex* a,
                    integer lda )
{
    return ctpttr(uplo, n, ap, a, lda);
}

inline integer tpttr( char uplo, integer n, const dcomplex* ap, dcomplex* a,
                    integer lda )
{
    return ztpttr(uplo, n, ap, a, lda);
}

inline integer trttf( char transr, char uplo, integer n, const float* a,
                    integer lda, float* arf )
{
    return strttf(transr, uplo, n, a, lda, arf);
}

inline integer trttf( char transr, char uplo, integer n, const double* a,
                    integer lda, double* arf )
{
    return dtrttf(transr, uplo, n, a, lda, arf);
}

inline integer trttf( char transr, char uplo, integer n, const scomplex* a,
                    integer lda, scomplex* arf )
{
    return ctrttf(transr, uplo, n, a, lda, arf);
}

inline integer trttf( char transr, char uplo, integer n, const dcomplex* a,
                    integer lda, dcomplex* arf )
{
    return ztrttf(transr, uplo, n, a, lda, arf);
}

inline integer trttp( char uplo, integer n, const float* a, integer lda,
                    float* ap )
{
    return strttp(uplo, n, a, lda, ap);
}

inline integer trttp( char uplo, integer n, const double* a, integer lda,
                    double* ap )
{
    return dtrttp(uplo, n, a, lda, ap);
}

inline integer trttp( char uplo, integer n, const scomplex* a, integer lda,
                    scomplex* ap )
{
    return ctrttp(uplo, n, a, lda, ap);
}

inline integer trttp( char uplo, integer n, const dcomplex* a, integer lda,
                    dcomplex* ap )
{
    return ztrttp(uplo, n, a, lda, ap);
}

inline integer geqrfp( integer m, integer n, float* a, integer lda,
                     float* tau)
{
    return sgeqrfp(m, n, a, lda, tau);
}

inline integer geqrfp( integer m, integer n, double* a, integer lda,
                     double* tau)
{
    return dgeqrfp(m, n, a, lda, tau);
}

inline integer geqrfp( integer m, integer n, scomplex* a, integer lda,
                     scomplex* tau)
{
    return cgeqrfp(m, n, a, lda, tau);
}

inline integer geqrfp( integer m, integer n, dcomplex* a, integer lda,
                     dcomplex* tau)
{
    return zgeqrfp(m, n, a, lda, tau);
}

#endif

#endif
