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

#include "fortran.h"

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

void FFUNC(zhetrf,ZHETRF)( char* uplo, integer* n, dcomplex* a, integer* lda,
                    integer* ipiv, dcomplex* work, integer* lwork,
                    integer *info );

void FFUNC(csytrf,CSYTRF)( char* uplo, integer* n, scomplex* a, integer* lda,
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
                    double* anorm, double* rcond, scomplex* work,
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
                    double* rcond, scomplex* work, float* rwork,
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
                     const integer* ipiv, const double* r, const double* c,
                     const scomplex* b, integer* ldb, scomplex* x,
                     integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
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
                    integer* nrhs, const dcomplex* ap, const scomplex* b,
                    integer* ldb, const scomplex* x, integer* ldx,
                    double* ferr, double* berr, scomplex* work, double* rwork,
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
                    integer* ldx, dcomplex* work, scomplex* swork, float* rwork,
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
                    integer* ldx, dcomplex* work, scomplex* swork, float* rwork,
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
                    integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* afb,
                    integer* ldafb, char* equed, float* s, dcomplex* b,
                    integer* ldb, dcomplex* x, integer* ldx, float* rcond,
                    float* ferr, float* berr, dcomplex* work, float* rwork,
                    integer *info );

void FFUNC(zpbsvx,ZPBSVX)( char* fact, char* uplo, integer* n, integer* kd,
                    integer* nrhs, scomplex* ab, integer* ldab, scomplex* afb,
                    integer* ldafb, char* equed, double* s, scomplex* b,
                    integer* ldb, scomplex* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, scomplex* work, double* rwork,
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
                    integer* jpvt, scomplex* tau, scomplex* work, float* rwork,
                    integer* lwork, integer *info );

void FFUNC(zgeqp3,ZGEQP3)( integer* m, integer* n, dcomplex* a, integer* lda,
                    integer* jpvt, dcomplex* tau, dcomplex* work, double* rwork,
                    integer* lwork, integer *info );

void FFUNC(sorgqr,sORGQR)( integer* m, integer* n, integer* k, float* a,
                    integer* lda, const float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dorgqr,DORGQR)( integer* m, integer* n, integer* k, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(sormqr,SORMQR)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormqr,DORMQR)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
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

void FFUNC(sormlq,SORMLQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormlq,DORMLQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
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

void FFUNC(sormql,SORMQL)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormql,DORMQL)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
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

void FFUNC(sormrq,SORMRQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormrq,DORMRQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
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

void FFUNC(sormtr,SORMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const float* a, integer* lda,
                    const float* tau, float* c, integer* ldc, float* work,
                    integer* lwork, integer *info );

void FFUNC(dormtr,DORMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
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

void FFUNC(sopmtr,SOPMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const float* ap, const float* tau,
                    float* c, integer* ldc, float* work,
                    integer *info );

void FFUNC(dopmtr,DOPMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const double* ap, const double* tau,
                    double* c, integer* ldc, double* work,
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

void FFUNC(csbgst,CSBGST)( char* vect, char* uplo, integer* n, integer* ka,
                    integer* kb, scomplex* ab, integer* ldab,
                    const scomplex* bb, integer* ldbb, scomplex* x,
                    integer* ldx, scomplex* work, float* rwork, integer *info );

void FFUNC(zsbgst,ZSBGST)( char* vect, char* uplo, integer* n, integer* ka,
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

void FFUNC(sgehrd,SGEHRD)( integer* n, integer* ilo, integer* ihi, double* a,
                    integer* lda, float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dgehrd,DGEHRD)( integer* n, integer* ilo, integer* ihi, double* a,
                    integer* lda, double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(cgehrd,CGEHRD)( integer* n, integer* ilo, integer* ihi, double* a,
                    integer* lda, scomplex* tau, scomplex* work,
                    integer* lwork, integer *info );

void FFUNC(zgehrd,ZGEHRD)( integer* n, integer* ilo, integer* ihi, double* a,
                    integer* lda, dcomplex* tau, dcomplex* work,
                    integer* lwork, integer *info );

void FFUNC(sorghr,SORGHR)( integer* n, integer* ilo, integer* ihi, float* a,
                    integer* lda, const float* tau, float* work,
                    integer* lwork, integer *info );

void FFUNC(dorghr,DORGHR)( integer* n, integer* ilo, integer* ihi, double* a,
                    integer* lda, const double* tau, double* work,
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
                    integer* ihi, const double* scale, integer* m,
                    scomplex* v, integer* ldv, integer *info );

void FFUNC(zgebak,ZGEBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const float* scale, integer* m,
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

void FFUNC(strevc,STREVC)( char* side, char* howmny, logical* select,
                    integer* n, const float* t, integer* ldt, float* vl,
                    integer* ldvl, float* vr, integer* ldvr,
                    integer* mm, integer* m, float* work,
                    integer *info );

void FFUNC(dtrevc,DTREVC)( char* side, char* howmny, logical* select,
                    integer* n, const double* t, integer* ldt, double* vl,
                    integer* ldvl, double* vr, integer* ldvr,
                    integer* mm, integer* m, double* work,
                    integer *info );

void FFUNC(ctrevc,CTREVC)( char* side, char* howmny, logical* select,
                    integer* n, const scomplex* t, integer* ldt, scomplex* vl,
                    integer* ldvl, scomplex* vr, integer* ldvr,
                    integer* mm, integer* m, scomplex* work, float* rwork,
                    integer *info );

void FFUNC(ztrevc,ZTREVC)( char* side, char* howmny, logical* select,
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
                    float* s, float* sep, double* work, integer* lwork, integer *info );

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

void FFUNC(sgelss,SGELSS)( integer* m, integer* n, integer* nrhs, double* a,
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

void FFUNC(sgelsd,SGELSD)( integer* m, integer* n, integer* nrhs, float* a,
                    integer* lda, float* b, integer* ldb, float* s,
                    float* rcond, integer* rank, float* work,
                    integer* lwork, integer* iwork, integer *info );

void FFUNC(dgelsd,DGELSD)( integer* m, integer* n, integer* nrhs, double* a,
                    integer* lda, double* b, integer* ldb, double* s,
                    double* rcond, integer* rank, double* work,
                    integer* lwork, integer* iwork, integer *info );

void FFUNC(cgelsd,CGELSD)( integer* m, integer* n, integer* nrhs, scomplex* a,
                    integer* lda, scomplex* b, integer* ldb, float* s,
                    float* rcond, integer* rank, scomplex* work,
                    integer* lwork, float* rwork, integer* iwork, integer *info );

void FFUNC(zgelsd,ZGELSD)( integer* m, integer* n, integer* nrhs, dcomplex* a,
                    integer* lda, dcomplex* b, integer* ldb, double* s,
                    double* rcond, integer* rank, dcomplex* work,
                    integer* lwork, double* rwork, integer* iwork, integer *info );

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

void FFUNC(csyev,CSYEV)( char* jobz, char* uplo, integer* n, scomplex* a,
                   integer* lda, float* w, scomplex* work, integer* lwork,
                   float* rwork, integer *info );

void FFUNC(zsyev,ZSYEV)( char* jobz, char* uplo, integer* n, scomplex* a,
                   integer* lda, double* w, dcomplex* work, integer* lwork,
                   double* rwork, integer *info );

void FFUNC(ssyevd,SSYEVD)( char* jobz, char* uplo, integer* n, float* a,
                    integer* lda, float* w, float* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dsyevd,DSYEVD)( char* jobz, char* uplo, integer* n, double* a,
                    integer* lda, double* w, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(csyevd,CSYEVD)( char* jobz, char* uplo, integer* n, scomplex* a,
                    integer* lda, float* w, scomplex* work, integer* lwork,
                    float* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(zsyevd,ZSYEVD)( char* jobz, char* uplo, integer* n, dcomplex* a,
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

void FFUNC(csyevx,CSYEVX)( char* jobz, char* range, char* uplo, integer* n,
                    scomplex* a, integer* lda, float* vl, float* vu,
                    integer* il, integer* iu, float* abstol,
                    integer* m, float* w, scomplex* z, integer* ldz,
                    scomplex* work, integer* lwork, float* rwork, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(zsyevx,ZSYEVX)( char* jobz, char* range, char* uplo, integer* n,
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

void FFUNC(csyevr,CSYEVR)( char* jobz, char* range, char* uplo, integer* n,
                    scomplex* a, integer* lda, float* vl, float* vu,
                    integer* il, integer* iu, float* abstol,
                    integer* m, float* w, scomplex* z, integer* ldz,
                    integer* isuppz, scomplex* work, integer* lwork,
                    float* rwork, integer* lrwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(zsyevr,ZSYEVR)( char* jobz, char* range, char* uplo, integer* n,
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

void FFUNC(cspev,CSPEV)( char* jobz, char* uplo, integer* n, scomplex* ap, float* w,
                   scomplex* z, integer* ldz, scomplex* work, float* rwork, integer *info );

void FFUNC(zspev,ZSPEV)( char* jobz, char* uplo, integer* n, dcomplex* ap, double* w,
                   dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer *info );

void FFUNC(dspevd,DSPEVD)( char* jobz, char* uplo, integer* n, double* ap,
                    double* w, double* z, integer* ldz, double* work,
                    integer* lwork, integer* iwork, integer* liwork,
                    integer *info );

void FFUNC(dspevx,DSPEVX)( char* jobz, char* range, char* uplo, integer* n,
                    double* ap, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dsbev,DSBEV)( char* jobz, char* uplo, integer* n, integer* kd,
                   double* ab, integer* ldab, double* w, double* z,
                   integer* ldz, double* work, integer *info );

void FFUNC(dsbevd,DSBEVD)( char* jobz, char* uplo, integer* n, integer* kd,
                    double* ab, integer* ldab, double* w, double* z,
                    integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dsbevx,DSBEVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* kd, double* ab, integer* ldab, double* q,
                    integer* ldq, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dstev,DSTEV)( char* jobz, integer* n, double* d, double* e, double* z,
                   integer* ldz, double* work, integer *info );

void FFUNC(dstevd,DSTEVD)( char* jobz, integer* n, double* d, double* e, double* z,
                    integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dstevx,DSTEVX)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dstevr,DSTEVR)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, integer* isuppz,
                    double* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(dgeev,DGEEV)( char* jobvl, char* jobvr, integer* n, double* a,
                   integer* lda, double* wr, double* wi, double* vl,
                   integer* ldvl, double* vr, integer* ldvr, double* work,
                   integer* lwork, integer *info );

void FFUNC(dgeevx,DGEEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, double* a, integer* lda, double* wr,
                    double* wi, double* vl, integer* ldvl, double* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    double* scale, double* abnrm, double* rconde,
                    double* rcondv, double* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(dgesvd,DGESVD)( char* jobu, char* jobvt, integer* m, integer* n,
                    double* a, integer* lda, double* s, double* u,
                    integer* ldu, double* vt, integer* ldvt, double* work,
                    integer* lwork, integer *info );

void FFUNC(dgesdd,DGESDD)( char* jobz, integer* m, integer* n, double* a,
                    integer* lda, double* s, double* u, integer* ldu,
                    double* vt, integer* ldvt, double* work,
                    integer* lwork, integer* iwork, integer *info );

void FFUNC(dgejsv,DGEJSV)( char* joba, char* jobu, char* jobv, char* jobr, char* jobt,
                    char* jobp, integer* m, integer* n, double* a,
                    integer* lda, double* sva, double* u, integer* ldu,
                    double* v, integer* ldv, double* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(dgesvj,DGESVJ)( char* joba, char* jobu, char* jobv, integer* m,
                    integer* n, double* a, integer* lda, double* sva,
                    integer* mv, double* v, integer* ldv, double* work,
                    integer* lwork, integer *info );

void FFUNC(dggsvd,DGGSVD)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* n, integer* p, integer* k, integer* l,
                    double* a, integer* lda, double* b, integer* ldb,
                    double* alpha, double* beta, double* u, integer* ldu,
                    double* v, integer* ldv, double* q, integer* ldq,
                    double* work, integer* iwork, integer *info );

void FFUNC(dsygv,DSYGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   double* a, integer* lda, double* b, integer* ldb,
                   double* w, double* work, integer* lwork,
                   integer *info );

void FFUNC(dsygvd,DSYGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    double* a, integer* lda, double* b, integer* ldb,
                    double* w, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dsygvx,DSYGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, double* a, integer* lda, double* b,
                    integer* ldb, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* ifail, integer *info );

void FFUNC(dspgv,DSPGV)( integer* itype, char* jobz, char* uplo, integer* n,
                   double* ap, double* bp, double* w, double* z,
                   integer* ldz, double* work, integer *info );

void FFUNC(dspgvd,DSPGVD)( integer* itype, char* jobz, char* uplo, integer* n,
                    double* ap, double* bp, double* w, double* z,
                    integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dspgvx,DSPGVX)( integer* itype, char* jobz, char* range, char* uplo,
                    integer* n, double* ap, double* bp, double* vl,
                    double* vu, integer* il, integer* iu, double* abstol,
                    integer* m, double* w, double* z, integer* ldz,
                    double* work, integer* iwork, integer* ifail,
                    integer *info );

void FFUNC(dsbgv,DSBGV)( char* jobz, char* uplo, integer* n, integer* ka,
                   integer* kb, double* ab, integer* ldab, double* bb,
                   integer* ldbb, double* w, double* z, integer* ldz,
                   double* work, integer *info );

void FFUNC(dsbgvd,DSBGVD)( char* jobz, char* uplo, integer* n, integer* ka,
                    integer* kb, double* ab, integer* ldab, double* bb,
                    integer* ldbb, double* w, double* z, integer* ldz,
                    double* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(dsbgvx,DSBGVX)( char* jobz, char* range, char* uplo, integer* n,
                    integer* ka, integer* kb, double* ab,
                    integer* ldab, double* bb, integer* ldbb, double* q,
                    integer* ldq, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, double* work, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dggev,DGGEV)( char* jobvl, char* jobvr, integer* n, double* a,
                   integer* lda, double* b, integer* ldb, double* alphar,
                   double* alphai, double* beta, double* vl, integer* ldvl,
                   double* vr, integer* ldvr, double* work,
                   integer* lwork, integer *info );

void FFUNC(dggevx,DGGEVX)( char* balanc, char* jobvl, char* jobvr, char* sense,
                    integer* n, double* a, integer* lda, double* b,
                    integer* ldb, double* alphar, double* alphai,
                    double* beta, double* vl, integer* ldvl, double* vr,
                    integer* ldvr, integer* ilo, integer* ihi,
                    double* lscale, double* rscale, double* abnrm,
                    double* bbnrm, double* rconde, double* rcondv, double* work,
                    integer* lwork, integer* iwork, logical* bwork,
                    integer *info );

void FFUNC(dsfrk,DSFRK)( char* transr, char* uplo, char* trans, integer* n,
                   integer* k, double* alpha, const double* a,
                   integer* lda, double* beta, double* c );

void FFUNC(dtfsm,DTFSM)( char* transr, char* side, char* uplo, char* trans,
                   char* diag, integer* m, integer* n, double* alpha,
                   const double* a, double* b, integer* ldb );

void FFUNC(dtfttp,DTFTTP)( char* transr, char* uplo, integer* n, const double* arf,
                    double* ap, integer *info );

void FFUNC(dtfttr,DTFTTR)( char* transr, char* uplo, integer* n, const double* arf,
                    double* a, integer* lda, integer *info );

void FFUNC(dtpttf,DTPTTF)( char* transr, char* uplo, integer* n, const double* ap,
                    double* arf, integer *info );

void FFUNC(dtpttr,DTPTTR)( char* uplo, integer* n, const double* ap, double* a,
                    integer* lda, integer *info );

void FFUNC(dtrttf,DTRTTF)( char* transr, char* uplo, integer* n, const double* a,
                    integer* lda, double* arf, integer *info );

void FFUNC(dtrttp,DTRTTP)( char* uplo, integer* n, const double* a, integer* lda,
                    double* ap, integer *info );

void FFUNC(dgeqrfp,DGEQRFP)( integer* m, integer* n, double* a, integer* lda,
                     double* tau, double* work, integer* lwork,
                     integer *info );

/*
 * C wrappers
 */
static inline integer c_dgetrf( integer m, integer n, double* a, integer lda,
                    integer* ipiv )
{
    integer info;
    FFUNC(dgetrf,DGETRF)(&m, &n, a, &lda, ipiv, &info);
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

static inline integer c_dgttrf( integer n, double* dl, double* d, double* du,
                    double* du2, integer* ipiv )
{
    integer info;
    FFUNC(dgttrf,DGTTRF)(&n, dl, d, du, du2, ipiv, &info);
    return info;
}

static inline integer c_dpotrf( char uplo, integer n, double* a, integer lda )
{
    integer info;
    FFUNC(dpotrf,DPOTRF)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_dpstrf( char uplo, integer n, double* a, integer lda,
                    integer* piv, integer* rank, double tol,
                    double* work )
{
    integer info;
    FFUNC(dpstrf,DPSTRF)(&uplo, &n, a, &lda, piv, rank, &tol, work, &info);
    return info;
}

static inline integer c_dpftrf( char transr, char uplo, integer n, double* a )
{
    integer info;
    FFUNC(dpftrf,DPFTRF)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_dpptrf( char uplo, integer n, double* ap )
{
    integer info;
    FFUNC(dpptrf,DPPTRF)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_dpbtrf( char uplo, integer n, integer kd, double* ab,
                    integer ldab )
{
    integer info;
    FFUNC(dpbtrf,DPBTRF)(&uplo, &n, &kd, ab, &ldab, &info);
    return info;
}

static inline integer c_dpttrf( integer n, double* d, double* e )
{
    integer info;
    FFUNC(dpttrf,DPTTRF)(&n, d, e, &info);
    return info;
}

static inline integer c_dsytrf( char uplo, integer n, double* a, integer lda,
                    integer* ipiv, double* work, integer* lwork )
{
    integer info;
    FFUNC(dsytrf,DSYTRF)(&uplo, &n, a, &lda, ipiv, work, lwork, &info);
    return info;
}

static inline integer c_dsptrf( char uplo, integer n, double* ap, integer* ipiv )
{
    integer info;
    FFUNC(dsptrf,DSPTRF)(&uplo, &n, ap, ipiv, &info);
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

static inline integer c_dgbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const double* ab, integer ldab,
                    const integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dgbtrs,DGBTRS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
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

static inline integer c_dpotrs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dpotrs,DPOTRS)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_dpftrs( char transr, char uplo, integer n, integer nrhs,
                    const double* a, double* b, integer ldb )
{
    integer info;
    FFUNC(dpftrs,DPFTRS)(&transr, &uplo, &n, &nrhs, a, b, &ldb, &info);
    return info;
}

static inline integer c_dpptrs( char uplo, integer n, integer nrhs,
                    const double* ap, double* b, integer ldb )
{
    integer info;
    FFUNC(dpptrs,DPPTRS)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
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

static inline integer c_dpttrs( integer n, integer nrhs, const double* d,
                    const double* e, double* b, integer ldb )
{
    integer info;
    FFUNC(dpttrs,DPTTRS)(&n, &nrhs, d, e, b, &ldb, &info);
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

static inline integer c_dsptrs( char uplo, integer n, integer nrhs,
                    const double* ap, const integer* ipiv, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dsptrs,DSPTRS)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
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

static inline integer c_dtptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* ap, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dtptrs,DTPTRS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info);
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

static inline integer c_dgecon( char norm, integer n, const double* a, integer lda,
                    double anorm, double* rcond, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dgecon,DGECON)(&norm, &n, a, &lda, &anorm, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dgbcon( char norm, integer n, integer kl, integer ku,
                    const double* ab, integer ldab, const integer* ipiv,
                    double anorm, double* rcond, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dgbcon,DGBCON)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dgtcon( char norm, integer n, const double* dl,
                    const double* d, const double* du, const double* du2,
                    const integer* ipiv, double anorm, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dgtcon,DGTCON)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dpocon( char uplo, integer n, const double* a, integer lda,
                    double anorm, double* rcond, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dpocon,DPOCON)(&uplo, &n, a, &lda, &anorm, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dppcon( char uplo, integer n, const double* ap, double anorm,
                    double* rcond, double* work, integer* iwork )
{
    integer info;
    FFUNC(dppcon,DPPCON)(&uplo, &n, ap, &anorm, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dpbcon( char uplo, integer n, integer kd, const double* ab,
                    integer ldab, double anorm, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dpbcon,DPBCON)(&uplo, &n, &kd, ab, &ldab, &anorm, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dptcon( integer n, const double* d, const double* e,
                    double anorm, double* rcond, double* work )
{
    integer info;
    FFUNC(dptcon,DPTCON)(&n, d, e, &anorm, rcond, work, &info);
    return info;
}

static inline integer c_dsycon( char uplo, integer n, const double* a, integer lda,
                    const integer* ipiv, double anorm, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dsycon,DSYCON)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dspcon( char uplo, integer n, const double* ap,
                    const integer* ipiv, double anorm, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dspcon,DSPCON)(&uplo, &n, ap, ipiv, &anorm, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dtrcon( char norm, char uplo, char diag, integer n,
                    const double* a, integer lda, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dtrcon,DTRCON)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dtpcon( char norm, char uplo, char diag, integer n,
                    const double* ap, double* rcond, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dtpcon,DTPCON)(&norm, &uplo, &diag, &n, ap, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dtbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const double* ab, integer ldab,
                    double* rcond, double* work, integer* iwork )
{
    integer info;
    FFUNC(dtbcon,DTBCON)(&norm, &uplo, &diag, &n, &kd, ab, &ldab, rcond, work, iwork, &info);
    return info;
}

static inline integer c_dgerfs( char trans, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dgerfs,DGERFS)(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dgerfsx( char trans, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const integer* ipiv, const double* r,
                     const double* c, const double* b, integer ldb,
                     double* x, integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params,
                     double* work, integer* iwork )
{
    integer info;
    FFUNC(dgerfsx,DGERFSX)(&trans, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    return info;
}

static inline integer c_dgbrfs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const double* ab, integer ldab,
                    const double* afb, integer ldafb,
                    const integer* ipiv, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dgbrfs,DGBRFS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb, x, &ldx,
                            ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dgbrfsx( char trans, char equed, integer n, integer kl,
                     integer ku, integer nrhs, const double* ab,
                     integer ldab, const double* afb, integer ldafb,
                     const integer* ipiv, const double* r, const double* c,
                     const double* b, integer ldb, double* x,
                     integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params,
                     double* work, integer* iwork )
{
    integer info;
    FFUNC(dgbrfsx,DGBRFSX)(&trans, &equed, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, r, c, b, &ldb, x, &ldx,
                              rcond, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    return info;
}

static inline integer c_dgtrfs( char trans, integer n, integer nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* dlf, const double* df, const double* duf,
                    const double* du2, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dgtrfs,DGTRFS)(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr,
                            work, iwork, &info);
    return info;
}

static inline integer c_dporfs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dporfs,DPORFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dporfsx( char uplo, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const double* s, const double* b,
                     integer ldb, double* x, integer ldx, double* rcond,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params, double* work,
                     integer* iwork )
{
    integer info;
    FFUNC(dporfsx,DPORFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, s, b, &ldb, x, &ldx, rcond, berr, &n_err_bnds,
                              err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    return info;
}

static inline integer c_dpprfs( char uplo, integer n, integer nrhs,
                    const double* ap, const double* afp, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dpprfs,DPPRFS)(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dpbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const double* ab, integer ldab, const double* afb,
                    integer ldafb, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dpbrfs,DPBRFS)(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dptrfs( integer n, integer nrhs, const double* d,
                    const double* e, const double* df, const double* ef,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* ferr, double* berr, double* work )
{
    integer info;
    FFUNC(dptrfs,DPTRFS)(&n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, &info);
    return info;
}

static inline integer c_dsyrfs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dsyrfs,DSYRFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dsyrfsx( char uplo, char equed, integer n, integer nrhs,
                     const double* a, integer lda, const double* af,
                     integer ldaf, const integer* ipiv, const double* s,
                     const double* b, integer ldb, double* x,
                     integer ldx, double* rcond, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params,
                     double* work, integer* iwork )
{
    integer info;
    FFUNC(dsyrfsx,DSYRFSX)(&uplo, &equed, &n, &nrhs, a, &lda, af, &ldaf, ipiv, s, b, &ldb, x, &ldx, rcond, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    return info;
}

static inline integer c_dsprfs( char uplo, integer n, integer nrhs,
                    const double* ap, const double* afp, const integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* ferr, double* berr, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dsprfs,DSPRFS)(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dtrrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* a, integer lda,
                    const double* b, integer ldb, const double* x,
                    integer ldx, double* ferr, double* berr, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dtrrfs,DTRRFS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dtprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* ap, const double* b,
                    integer ldb, const double* x, integer ldx,
                    double* ferr, double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dtprfs,DTPRFS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dtbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const double* ab,
                    integer ldab, const double* b, integer ldb,
                    const double* x, integer ldx, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dtbrfs,DTBRFS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dgetri( integer n, double* a, integer lda,
                    const integer* ipiv, double* work, integer lwork )
{
    integer info;
    FFUNC(dgetri,DGETRI)(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

static inline integer c_dpotri( char uplo, integer n, double* a, integer lda )
{
    integer info;
    FFUNC(dpotri,DPOTRI)(&uplo, &n, a, &lda, &info);
    return info;
}

static inline integer c_dpftri( char transr, char uplo, integer n, double* a )
{
    integer info;
    FFUNC(dpftri,DPFTRI)(&transr, &uplo, &n, a, &info);
    return info;
}

static inline integer c_dpptri( char uplo, integer n, double* ap )
{
    integer info;
    FFUNC(dpptri,DPPTRI)(&uplo, &n, ap, &info);
    return info;
}

static inline integer c_dsytri( char uplo, integer n, double* a, integer lda,
                    const integer* ipiv, double* work )
{
    integer info;
    FFUNC(dsytri,DSYTRI)(&uplo, &n, a, &lda, ipiv, work, &info);
    return info;
}

static inline integer c_dsptri( char uplo, integer n, double* ap,
                    const integer* ipiv, double* work )
{
    integer info;
    FFUNC(dsptri,DSPTRI)(&uplo, &n, ap, ipiv, work, &info);
    return info;
}

static inline integer c_dtrtri( char uplo, char diag, integer n, double* a,
                    integer lda )
{
    integer info;
    FFUNC(dtrtri,DTRTRI)(&uplo, &diag, &n, a, &lda, &info);
    return info;
}

static inline integer c_dtftri( char transr, char uplo, char diag, integer n,
                    double* a )
{
    integer info;
    FFUNC(dtftri,DTFTRI)(&transr, &uplo, &diag, &n, a, &info);
    return info;
}

static inline integer c_dtptri( char uplo, char diag, integer n, double* ap )
{
    integer info;
    FFUNC(dtptri,DTPTRI)(&uplo, &diag, &n, ap, &info);
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

static inline integer c_dgeequb( integer m, integer n, const double* a,
                     integer lda, double* r, double* c, double* rowcnd,
                     double* colcnd, double* amax )
{
    integer info;
    FFUNC(dgeequb,DGEEQUB)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
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

static inline integer c_dgbequb( integer m, integer n, integer kl,
                     integer ku, const double* ab, integer ldab,
                     double* r, double* c, double* rowcnd, double* colcnd,
                     double* amax )
{
    integer info;
    FFUNC(dgbequb,DGBEQUB)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

static inline integer c_dpoequ( integer n, const double* a, integer lda, double* s,
                    double* scond, double* amax )
{
    integer info;
    FFUNC(dpoequ,DPOEQU)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_dpoequb( integer n, const double* a, integer lda, double* s,
                     double* scond, double* amax )
{
    integer info;
    FFUNC(dpoequb,DPOEQUB)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

static inline integer c_dppequ( char uplo, integer n, const double* ap, double* s,
                    double* scond, double* amax )
{
    integer info;
    FFUNC(dppequ,DPPEQU)(&uplo, &n, ap, s, scond, amax, &info);
    return info;
}

static inline integer c_dpbequ( char uplo, integer n, integer kd, const double* ab,
                    integer ldab, double* s, double* scond, double* amax )
{
    integer info;
    FFUNC(dpbequ,DPBEQU)(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info);
    return info;
}

static inline integer c_dsyequb( char uplo, integer n, const double* a,
                     integer lda, double* s, double* scond, double* amax,
                     double* work )
{
    integer info;
    FFUNC(dsyequb,DSYEQUB)(&uplo, &n, a, &lda, s, scond, amax, work, &info);
    return info;
}

static inline integer c_dgesv( integer n, integer nrhs, double* a, integer lda,
                   integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dgesv,DGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dsgesv( integer n, integer nrhs, double* a, integer lda,
                    integer* ipiv, double* b, integer ldb, double* x,
                    integer ldx, double* work, float* swork,
                    integer* iter )
{
    integer info;
    FFUNC(dsgesv,DSGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, x, &ldx, work, swork, iter, &info);
    return info;
}

static inline integer c_dgesvx( char fact, char trans, integer n, integer nrhs,
                    double* a, integer lda, double* af, integer ldaf,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer ldb, double* x, integer ldx,
                    double* rcond, double* ferr, double* berr, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dgesvx,DGESVX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldx, x, &ldx, rcond,
                           ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dgesvxx( char fact, char trans, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer ldb, double* x, integer ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params,
                     double* work, integer* iwork )
{
    integer info;
    FFUNC(dgesvxx,DGESVXX)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldb, x, &ldx, rcond, rpvgrw,
                              berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
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

static inline integer c_dgbsvx( char fact, char trans, integer n, integer kl,
                    integer ku, integer nrhs, double* ab,
                    integer ldab, double* afb, integer ldafb,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer ldb, double* x, integer ldx,
                    double* rcond, double* ferr, double* berr, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dgbsvx,DGBSVX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                           rcond, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dgbsvxx( char fact, char trans, integer n, integer kl,
                     integer ku, integer nrhs, double* ab,
                     integer ldab, double* afb, integer ldafb,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer ldb, double* x, integer ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer nparams, double* params,
                     double* work, integer* iwork )
{
    integer info;
    FFUNC(dgbsvxx,DGBSVXX)(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, equed, r, c, b, &ldb, x, &ldx,
                             rcond, rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    return info;
}

static inline integer c_dgtsv( integer n, integer nrhs, double* dl, double* d,
                   double* du, double* b, integer ldb )
{
    integer info;
    FFUNC(dgtsv,DGTSV)(&n, &nrhs, dl, d, du, b, &ldb, &info);
    return info;
}

static inline integer c_dgtsvx( char fact, char trans, integer n, integer nrhs,
                    const double* dl, const double* d, const double* du,
                    double* dlf, double* df, double* duf, double* du2,
                    integer* ipiv, const double* b, integer ldb,
                    double* x, integer ldx, double* rcond, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dgtsvx,DGTSVX)(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, iwork, &info);
    return info;
}

static inline integer c_dposv( char uplo, integer n, integer nrhs, double* a,
                   integer lda, double* b, integer ldb )
{
    integer info;
    FFUNC(dposv,DPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_dsposv( char uplo, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* x,
                    integer ldx, double* work, float* swork,
                    integer* iter )
{
    integer info;
    FFUNC(dsposv,DSPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, work, swork, iter, &info);
    return info;
}

static inline integer c_dposvx( char fact, char uplo, integer n, integer nrhs,
                    double* a, integer lda, double* af, integer ldaf,
                    char* equed, double* s, double* b, integer ldb,
                    double* x, integer ldx, double* rcond, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dposvx,DPOSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, iwork, &info);
    return info;
}

static inline integer c_dposvxx( char fact, char uplo, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     char* equed, double* s, double* b, integer ldb,
                     double* x, integer ldx, double* rcond, double* rpvgrw,
                     double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params, double* work,
                     integer* iwork )
{
    integer info;
    FFUNC(dposvxx,DPOSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, rpvgrw, berr,
                              &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    return info;
}

static inline integer c_dppsv( char uplo, integer n, integer nrhs, double* ap,
                   double* b, integer ldb )
{
    integer info;
    FFUNC(dppsv,DPPSV)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

static inline integer c_dppsvx( char fact, char uplo, integer n, integer nrhs,
                    double* ap, double* afp, char* equed, double* s, double* b,
                    integer ldb, double* x, integer ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dppsvx,DPPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dpbsv( char uplo, integer n, integer kd, integer nrhs,
                   double* ab, integer ldab, double* b, integer ldb )
{
    integer info;
    FFUNC(dpbsv,DPBSV)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

static inline integer c_dpbsvx( char fact, char uplo, integer n, integer kd,
                    integer nrhs, double* ab, integer ldab, double* afb,
                    integer ldafb, char* equed, double* s, double* b,
                    integer ldb, double* x, integer ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dpbsvx,DPBSVX)(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, iwork, &info);
    return info;
}

static inline integer c_dptsv( integer n, integer nrhs, double* d, double* e,
                   double* b, integer ldb )
{
    integer info;
    FFUNC(dptsv,DPTSV)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

static inline integer c_dptsvx( char fact, integer n, integer nrhs,
                    const double* d, const double* e, double* df, double* ef,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* rcond, double* ferr, double* berr,
                    double* work )
{
    integer info;
    FFUNC(dptsvx,DPTSVX)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, &info);
    return info;
}

static inline integer c_dsysv( char uplo, integer n, integer nrhs, double* a,
                   integer lda, integer* ipiv, double* b,
                   integer ldb, double* work, integer lwork )
{
    integer info;
    FFUNC(dsysv,DSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);
    return info;
}

static inline integer c_dsysvx( char fact, char uplo, integer n, integer nrhs,
                    const double* a, integer lda, double* af,
                    integer ldaf, integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer lwork,
                    integer* iwork )
{
    integer info;
    FFUNC(dsysvx,DSYSVX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr,
                            work, &lwork, iwork, &info);
    return info;
}

static inline integer c_dsysvxx( char fact, char uplo, integer n, integer nrhs,
                     double* a, integer lda, double* af, integer ldaf,
                     integer* ipiv, char* equed, double* s, double* b,
                     integer ldb, double* x, integer ldx, double* rcond,
                     double* rpvgrw, double* berr, integer n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer nparams, double* params, double* work,
                     integer* iwork )
{
    integer info;
    FFUNC(dsysvxx,DSYSVXX)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, s, b, &ldb, x, &ldx, rcond,
                              rpvgrw, berr, &n_err_bnds, err_bnds_norm, err_bnds_comp, &nparams, params, work, iwork, &info);
    return info;
}

static inline integer c_dspsv( char uplo, integer n, integer nrhs, double* ap,
                   integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dspsv,DSPSV)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

static inline integer c_dspsvx( char fact, char uplo, integer n, integer nrhs,
                    const double* ap, double* afp, integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* rcond, double* ferr, double* berr,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dspsvx,DSPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, &info);
    return info;
}

static inline integer c_dgeqrf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dgeqrf,DGEQRF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dgeqpf( integer m, integer n, double* a, integer lda,
                    integer* jpvt, double* tau, double* work )
{
    integer info;
    FFUNC(dgeqpf,DGEQPF)(&m, &n, a, &lda, jpvt, tau, work, &info);
    return info;
}

static inline integer c_dgeqp3( integer m, integer n, double* a, integer lda,
                    integer* jpvt, double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgeqp3,DGEQP3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dorgqr( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorgqr,DORGQR)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dormqr( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormqr,DORMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

static inline integer c_dgelqf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dgelqf,DGELQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dorglq( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorglq,DORGLQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dormlq( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormlq,DORMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

static inline integer c_dgeqlf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dgeqlf,DGEQLF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dorgql( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorgql,DORGQL)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dormql( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormql,DORMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

static inline integer c_dgerqf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dgerqf,DGERQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dorgrq( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorgrq,DORGRQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dormrq( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormrq,DORMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

static inline integer c_dtzrzf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dtzrzf,DTZRZF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dormrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc, double* work, integer lwork )
{
    integer info;
    FFUNC(dormrz,DORMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

static inline integer c_dggqrf( integer n, integer m, integer p, double* a,
                    integer lda, double* taua, double* b, integer ldb,
                    double* taub, double* work, integer lwork )
{
    integer info;
    FFUNC(dggqrf,DGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    return info;
}

static inline integer c_dggrqf( integer m, integer p, integer n, double* a,
                    integer lda, double* taua, double* b, integer ldb,
                    double* taub, double* work, integer lwork )
{
    integer info;
    FFUNC(dggrqf,DGGRQF)(&m, &p, &n, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    return info;
}

static inline integer c_dgebrd( integer m, integer n, double* a, integer lda,
                    double* d, double* e, double* tauq, double* taup,
                    double* work, integer lwork )
{
    integer info;
    FFUNC(dgebrd,DGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
    return info;
}

static inline integer c_dgbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, double* ab,
                    integer ldab, double* d, double* e, double* q,
                    integer ldq, double* pt, integer ldpt, double* c,
                    integer ldc, double* work )
{
    integer info;
    FFUNC(dgbbrd,DGBBRD)(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info);
    return info;
}

static inline integer c_dorgbr( char vect, integer m, integer n, integer k,
                    double* a, integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorgbr,DORGBR)(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dormbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc, double* work, integer lwork )
{
    integer info;
    FFUNC(dormbr,DORMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

static inline integer c_dbdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, double* d, double* e,
                    double* vt, integer ldvt, double* u, integer ldu,
                    double* c, integer ldc, double* work )
{
    integer info;
    FFUNC(dbdsqr,DBDSQR)(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);
    return info;
}

static inline integer c_dbdsdc( char uplo, char compq, integer n, double* d,
                    double* e, double* u, integer ldu, double* vt,
                    integer ldvt, double* q, integer* iq, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dbdsdc,DBDSDC)(&uplo, &compq, &n, d, e, u, &ldu, vt, &ldvt, q, iq, work, iwork, &info);
    return info;
}

static inline integer c_dsytrd( char uplo, integer n, double* a, integer lda,
                    double* d, double* e, double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dsytrd,DSYTRD)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dorgtr( char uplo, integer n, double* a, integer lda,
                    const double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dorgtr,DORGTR)(&uplo, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dormtr( char side, char uplo, char trans, integer m,
                    integer n, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormtr,DORMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

static inline integer c_dsptrd( char uplo, integer n, double* ap, double* d, double* e,
                    double* tau )
{
    integer info;
    FFUNC(dsptrd,DSPTRD)(&uplo, &n, ap, d, e, tau, &info);
    return info;
}

static inline integer c_dopgtr( char uplo, integer n, const double* ap,
                    const double* tau, double* q, integer ldq, double* work )
{
    integer info;
    FFUNC(dopgtr,DOPGTR)(&uplo, &n, ap, tau, q, &ldq, work, &info);
    return info;
}

static inline integer c_dopmtr( char side, char uplo, char trans, integer m,
                    integer n, const double* ap, const double* tau,
                    double* c, integer ldc, double* work )
{
    integer info;
    FFUNC(dopmtr,DOPMTR)(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info);
    return info;
}

static inline integer c_dsbtrd( char vect, char uplo, integer n, integer kd,
                    double* ab, integer ldab, double* d, double* e,
                    double* q, integer ldq, double* work )
{
    integer info;
    FFUNC(dsbtrd,DSBTRD)(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info);
    return info;
}

static inline integer c_dsterf( integer n, double* d, double* e )
{
    integer info;
    FFUNC(dsterf,DSTERF)(&n, d, e, &info);
    return info;
}

static inline integer c_dsteqr( char compz, integer n, double* d, double* e, double* z,
                    integer ldz, double* work )
{
    integer info;
    FFUNC(dsteqr,DSTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    return info;
}

static inline integer c_dstemr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, integer* m, double* w, double* z,
                    integer ldz, integer nzc, integer* isuppz,
                    logical* tryrac, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dstemr,DSTEMR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac, work, &lwork,
                            iwork, &liwork, &info);
    return info;
}

static inline integer c_dstedc( char compz, integer n, double* d, double* e, double* z,
                    integer ldz, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dstedc,DSTEDC)(&compz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

static inline integer c_dstegr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, integer* isuppz,
                    double* work, integer lwork, integer* iwork,
                    integer liwork )
{
    integer info;
    FFUNC(dstegr,DSTEGR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    return info;
}

static inline integer c_dpteqr( char compz, integer n, double* d, double* e, double* z,
                    integer ldz, double* work )
{
    integer info;
    FFUNC(dpteqr,DPTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    return info;
}

static inline integer c_dstebz( char range, char order, integer n, double vl,
                    double vu, integer il, integer iu, double abstol,
                    const double* d, const double* e, integer* m,
                    integer* nsplit, double* w, integer* iblock,
                    integer* isplit, double* work, integer* iwork )
{
    integer info;
    FFUNC(dstebz,DSTEBZ)(&range, &order, &n, &vl, &vu, &il, &iu, &abstol, d, e, m, nsplit, w, iblock, isplit,
                            work, iwork, &info);
    return info;
}

static inline integer c_dstein( integer n, const double* d, const double* e,
                    integer m, const double* w, const integer* iblock,
                    const integer* isplit, double* z, integer ldz,
                    double* work, integer* iwork, integer* ifailv )
{
    integer info;
    FFUNC(dstein,DSTEIN)(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifailv, &info);
    return info;
}

static inline integer c_ddisna( char job, integer m, integer n, const double* d,
                    double* sep )
{
    integer info;
    FFUNC(ddisna,DDISNA)(&job, &m, &n, d, sep, &info);
    return info;
}

static inline integer c_dsygst( integer itype, char uplo, integer n, double* a,
                    integer lda, const double* b, integer ldb )
{
    integer info;
    FFUNC(dsygst,DSYGST)(&itype, &uplo, &n, a, &lda, b, &ldb, &info);
    return info;
}

static inline integer c_dspgst( integer itype, char uplo, integer n, double* ap,
                    const double* bp )
{
    integer info;
    FFUNC(dspgst,DSPGST)(&itype, &uplo, &n, ap, bp, &info);
    return info;
}

static inline integer c_dsbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, double* ab, integer ldab,
                    const double* bb, integer ldbb, double* x,
                    integer ldx, double* work )
{
    integer info;
    FFUNC(dsbgst,DSBGST)(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, &info);
    return info;
}

static inline integer c_dpbstf( char uplo, integer n, integer kb, double* bb,
                    integer ldbb )
{
    integer info;
    FFUNC(dpbstf,DPBSTF)(&uplo, &n, &kb, bb, &ldbb, &info);
    return info;
}

static inline integer c_dgehrd( integer n, integer ilo, integer ihi, double* a,
                    integer lda, double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgehrd,DGEHRD)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dorghr( integer n, integer ilo, integer ihi, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorghr,DORGHR)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    return info;
}

static inline integer c_dormhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc, double* work, integer lwork )
{
    integer info;
    FFUNC(dormhr,DORMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

static inline integer c_dgebal( char job, integer n, double* a, integer lda,
                    integer* ilo, integer* ihi, double* scale )
{
    integer info;
    FFUNC(dgebal,DGEBAL)(&job, &n, a, &lda, ilo, ihi, scale, &info);
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

static inline integer c_dhseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, double* h, integer ldh, double* wr,
                    double* wi, double* z, integer ldz, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dhseqr,DHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, work, &lwork, &info);
    return info;
}

static inline integer c_dhsein( char job, char eigsrc, char initv,
                    logical* select, integer n, const double* h,
                    integer ldh, double* wr, const double* wi, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer* m, double* work,
                    integer* ifaill, integer* ifailr )
{
    integer info;
    FFUNC(dhsein,DHSEIN)(&job, &eigsrc, &initv, select, &n, h, &ldh, wr, wi, vl, &ldvl, vr, &ldvr, &mm, m, work,
                            ifaill, ifailr, &info);
    return info;
}

static inline integer c_dtrevc( char side, char howmny, logical* select,
                    integer n, const double* t, integer ldt, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer* m, double* work )
{
    integer info;
    FFUNC(dtrevc,DTREVC)(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, m, work, &info);
    return info;
}

static inline integer c_dtrsna( char job, char howmny, const logical* select,
                    integer n, const double* t, integer ldt,
                    const double* vl, integer ldvl, const double* vr,
                    integer ldvr, double* s, double* sep, integer mm,
                    integer* m, double* work, integer ldwork,
                    integer* iwork )
{
    integer info;
    FFUNC(dtrsna,DTRSNA)(&job, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, m, work, &ldwork, iwork, &info);
    return info;
}

static inline integer c_dtrexc( char compq, integer n, double* t, integer ldt,
                    double* q, integer ldq, integer* ifst,
                    integer* ilst, double* work )
{
    integer info;
    FFUNC(dtrexc,DTREXC)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, work, &info);
    return info;
}

static inline integer c_dtrsen( char job, char compq, const logical* select,
                    integer n, double* t, integer ldt, double* q,
                    integer ldq, double* wr, double* wi, integer* m,
                    double* s, double* sep, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dtrsen,DTRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, m, s, sep, work, &lwork, iwork, &liwork, &info);
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

static inline integer c_dgghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, double* a, integer lda, double* b,
                    integer ldb, double* q, integer ldq, double* z,
                    integer ldz )
{
    integer info;
    FFUNC(dgghrd,DGGHRD)(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info);
    return info;
}

static inline integer c_dggbal( char job, integer n, double* a, integer lda,
                    double* b, integer ldb, integer* ilo,
                    integer* ihi, double* lscale, double* rscale,
                    double* work )
{
    integer info;
    FFUNC(dggbal,DGGBAL)(&job, &n, a, &lda, b, &ldb, ilo, ihi, lscale, rscale, work, &info);
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

static inline integer c_dhgeqz( char job, char compq, char compz, integer n,
                    integer ilo, integer ihi, double* h,
                    integer ldh, double* t, integer ldt, double* alphar,
                    double* alphai, double* beta, double* q, integer ldq,
                    double* z, integer ldz, double* work, integer lwork )
{
    integer info;
    FFUNC(dhgeqz,DHGEQZ)(&job, &compq, &compz, &n, &ilo, &ihi, h, &ldh, t, &ldt, alphar, alphai, beta, q, &ldq, z, &ldz,
                            work, &lwork, &info);
    return info;
}

static inline integer c_dtgevc( char side, char howmny, const logical* select,
                    integer n, const double* s, integer lds,
                    const double* p, integer ldp, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer* m, double* work )
{
    integer info;
    FFUNC(dtgevc,DTGEVC)(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, m, work, &info);
    return info;
}

static inline integer c_dtgexc( logical wantq, logical wantz, integer n,
                    double* a, integer lda, double* b, integer ldb,
                    double* q, integer ldq, double* z, integer ldz,
                    integer* ifst, integer* ilst, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dtgexc,DTGEXC)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, work, &lwork, &info);
    return info;
}

static inline integer c_dtgsen( integer ijob, logical wantq,
                    logical wantz, const logical* select,
                    integer n, double* a, integer lda, double* b,
                    integer ldb, double* alphar, double* alphai,
                    double* beta, double* q, integer ldq, double* z,
                    integer ldz, integer* m, double* pl, double* pr,
                    double* dif, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dtgsen,DTGSEN)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alphar, alphai, beta, q, &ldq, z, &ldz,
                            m, pl, pr, dif, work, &lwork, iwork, &liwork, &info);
    return info;
}

static inline integer c_dtgsyl( char trans, integer ijob, integer m, integer n,
                    const double* a, integer lda, const double* b,
                    integer ldb, double* c, integer ldc,
                    const double* d, integer ldd, const double* e,
                    integer lde, double* f, integer ldf, double* scale,
                    double* dif, double* work, integer lwork,
                    integer* iwork )
{
    integer info;
    FFUNC(dtgsyl,DTGSYL)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work,
                            &lwork, iwork, &info);
    return info;
}

static inline integer c_dtgsna( char job, char howmny, const logical* select,
                    integer n, const double* a, integer lda,
                    const double* b, integer ldb, const double* vl,
                    integer ldvl, const double* vr, integer ldvr,
                    double* s, double* dif, integer mm, integer* m,
                    double* work, integer lwork, integer* iwork )
{
    integer info;
    FFUNC(dtgsna,DTGSNA)(&job, &howmny, select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, m, work,
                            &lwork, iwork, &info);
    return info;
}

static inline integer c_dggsvp( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, double* a, integer lda,
                    double* b, integer ldb, double tola, double tolb,
                    integer* k, integer* l, double* u, integer ldu,
                    double* v, integer ldv, double* q, integer ldq,
                    integer* iwork, double* tau, double* work )
{
    integer info;
    FFUNC(dggsvp,DGGSVP)(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, &tola, &tolb, k, l, u, &ldu, v, &ldv, q, &ldq,
                            iwork, tau, work, &info);
    return info;
}

static inline integer c_dtgsja( char jobu, char jobv, char jobq, integer m,
                    integer p, integer n, integer k, integer l,
                    double* a, integer lda, double* b, integer ldb,
                    double tola, double tolb, double* alpha, double* beta,
                    double* u, integer ldu, double* v, integer ldv,
                    double* q, integer ldq, double* work,
                    integer* ncycle )
{
    integer info;
    FFUNC(dtgsja,DTGSJA)(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, &tola, &tolb, alpha, beta, u, &ldu,
                            v, &ldv, q, &ldq, work, ncycle, &info);
    return info;
}

static inline integer c_dgels( char trans, integer m, integer n, integer nrhs,
                   double* a, integer lda, double* b, integer ldb,
                   double* work, integer lwork )
{
    integer info;
    FFUNC(dgels,DGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

static inline integer c_dgelsy( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb,
                    integer* jpvt, double rcond, integer* rank,
                    double* work, integer lwork )
{
    integer info;
    FFUNC(dgelsy,DGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work, &lwork, &info);
    return info;
}

static inline integer c_dgelss( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* s,
                    double rcond, integer* rank, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgelss,DGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, &info);
    return info;
}

static inline integer c_dgelsd( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* s,
                    double rcond, integer* rank, double* work,
                    integer lwork, integer* iwork )
{
    integer info;
    FFUNC(dgelsd,DGELSD)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, iwork, &info);
    return info;
}

static inline integer c_dgglse( integer m, integer n, integer p, double* a,
                    integer lda, double* b, integer ldb, double* c,
                    double* d, double* x, double* work, integer lwork )
{
    integer info;
    FFUNC(dgglse,DGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
    return info;
}

static inline integer c_dggglm( integer n, integer m, integer p, double* a,
                    integer lda, double* b, integer ldb, double* d,
                    double* x, double* y, double* work, integer lwork )
{
    integer info;
    FFUNC(dggglm,DGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
    return info;
}

static inline integer c_dsyev( char jobz, char uplo, integer n, double* a,
                   integer lda, double* w, double* work, integer lwork )
{
    integer info;
    FFUNC(dsyev,DSYEV)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    return info;
}

static inline integer c_dsyevd( char jobz, char uplo, integer n, double* a,
                    integer lda, double* w, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dsyevd,DSYEVD)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
    return info;
}

static inline integer c_dsyevx( char jobz, char range, char uplo, integer n,
                    double* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer* m, double* w, double* z, integer ldz,
                    double* work, integer lwork, integer* iwork,
                    integer* ifail )
{
    integer info;
    FFUNC(dsyevx,DSYEVX)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork,
                            iwork, ifail, &info);
    return info;
}

static inline integer c_dsyevr( char jobz, char range, char uplo, integer n,
                    double* a, integer lda, double vl, double vu,
                    integer il, integer iu, double abstol,
                    integer* m, double* w, double* z, integer ldz,
                    integer* isuppz, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dsyevr,DSYEVR)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    return info;
}

static inline integer c_dspev( char jobz, char uplo, integer n, double* ap, double* w,
                   double* z, integer ldz, double* work )
{
    integer info;
    FFUNC(dspev,DSPEV)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &info);
    return info;
}

static inline integer c_dspevd( char jobz, char uplo, integer n, double* ap,
                    double* w, double* z, integer ldz, double* work,
                    integer lwork, integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dspevd,DSPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

static inline integer c_dspevx( char jobz, char range, char uplo, integer n,
                    double* ap, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, double* work, integer* iwork,
                    integer* ifail )
{
    integer info;
    FFUNC(dspevx,DSPEVX)(&jobz, &range, &uplo, &n, ap, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, &info);
    return info;
}

static inline integer c_dsbev( char jobz, char uplo, integer n, integer kd,
                   double* ab, integer ldab, double* w, double* z,
                   integer ldz, double* work )
{
    integer info;
    FFUNC(dsbev,DSBEV)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info);
    return info;
}

static inline integer c_dsbevd( char jobz, char uplo, integer n, integer kd,
                    double* ab, integer ldab, double* w, double* z,
                    integer ldz, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dsbevd,DSBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

static inline integer c_dsbevx( char jobz, char range, char uplo, integer n,
                    integer kd, double* ab, integer ldab, double* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, double* work, integer* iwork,
                    integer* ifail )
{
    integer info;
    FFUNC(dsbevx,DSBEVX)(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                           work, iwork, ifail, &info);
    return info;
}

static inline integer c_dstev( char jobz, integer n, double* d, double* e, double* z,
                   integer ldz, double* work )
{
    integer info;
    FFUNC(dstev,DSTEV)(&jobz, &n, d, e, z, &ldz, work, &info);
    return info;
}

static inline integer c_dstevd( char jobz, integer n, double* d, double* e, double* z,
                    integer ldz, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dstevd,DSTEVD)(&jobz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

static inline integer c_dstevx( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, double* work, integer* iwork,
                    integer* ifail )
{
    integer info;
    FFUNC(dstevx,DSTEVX)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, &info);
    return info;
}

static inline integer c_dstevr( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, integer* isuppz,
                    double* work, integer lwork, integer* iwork,
                    integer liwork )
{
    integer info;
    FFUNC(dstevr,DSTEVR)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork,
                            iwork, &liwork, &info);
    return info;
}

static inline integer c_dgeev( char jobvl, char jobvr, integer n, double* a,
                   integer lda, double* wr, double* wi, double* vl,
                   integer ldvl, double* vr, integer ldvr, double* work,
                   integer lwork )
{
    integer info;
    FFUNC(dgeev,DGEEV)(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    return info;
}

static inline integer c_dgeevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, double* a, integer lda, double* wr,
                    double* wi, double* vl, integer ldvl, double* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    double* scale, double* abnrm, double* rconde,
                    double* rcondv, double* work, integer lwork,
                    integer* iwork )
{
    integer info;
    FFUNC(dgeevx,DGEEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, ilo, ihi, scale, abnrm,
                            rconde, rcondv, work, &lwork, iwork, &info);
    return info;
}

static inline integer c_dgesvd( char jobu, char jobvt, integer m, integer n,
                    double* a, integer lda, double* s, double* u,
                    integer ldu, double* vt, integer ldvt, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgesvd,DGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
    return info;
}

static inline integer c_dgesdd( char jobz, integer m, integer n, double* a,
                    integer lda, double* s, double* u, integer ldu,
                    double* vt, integer ldvt, double* work,
                    integer lwork, integer* iwork )
{
    integer info;
    FFUNC(dgesdd,DGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
    return info;
}

static inline integer c_dgejsv( char joba, char jobu, char jobv, char jobr, char jobt,
                    char jobp, integer m, integer n, double* a,
                    integer lda, double* sva, double* u, integer ldu,
                    double* v, integer ldv, double* work, integer lwork,
                    integer* iwork )
{
    integer info;
    FFUNC(dgejsv,DGEJSV)(&joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, a, &lda, sva, u, &ldu, v, &ldv, work, &lwork,
                            iwork, &info);
    return info;
}

static inline integer c_dgesvj( char joba, char jobu, char jobv, integer m,
                    integer n, double* a, integer lda, double* sva,
                    integer mv, double* v, integer ldv, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgesvj,DGESVJ)(&joba, &jobu, &jobv, &m, &n, a, &lda, sva, &mv, v, &ldv, work, &lwork, &info);
    return info;
}

static inline integer c_dggsvd( char jobu, char jobv, char jobq, integer m,
                    integer n, integer p, integer* k, integer* l,
                    double* a, integer lda, double* b, integer ldb,
                    double* alpha, double* beta, double* u, integer ldu,
                    double* v, integer ldv, double* q, integer ldq,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dggsvd,DGGSVD)(&jobu, &jobv, &jobq, &m, &n, &p, k, l, a, &lda, b, &ldb, alpha, beta, u, &ldu, v, &ldv, q, &ldq,
                            work, iwork, &info);
    return info;
}

static inline integer c_dsygv( integer itype, char jobz, char uplo, integer n,
                   double* a, integer lda, double* b, integer ldb,
                   double* w, double* work, integer lwork )
{
    integer info;
    FFUNC(dsygv,DSYGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info);
    return info;
}

static inline integer c_dsygvd( integer itype, char jobz, char uplo, integer n,
                    double* a, integer lda, double* b, integer ldb,
                    double* w, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dsygvd,DSYGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, iwork, &liwork, &info);
    return info;
}

static inline integer c_dsygvx( integer itype, char jobz, char range, char uplo,
                    integer n, double* a, integer lda, double* b,
                    integer ldb, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, double* work, integer lwork,
                    integer* iwork, integer* ifail )
{
    integer info;
    FFUNC(dsygvx,DSYGVX)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz,
                            work, &lwork, iwork, ifail, &info);
    return info;
}

static inline integer c_dspgv( integer itype, char jobz, char uplo, integer n,
                   double* ap, double* bp, double* w, double* z,
                   integer ldz, double* work )
{
    integer info;
    FFUNC(dspgv,DSPGV)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &info);
    return info;
}

static inline integer c_dspgvd( integer itype, char jobz, char uplo, integer n,
                    double* ap, double* bp, double* w, double* z,
                    integer ldz, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dspgvd,DSPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

static inline integer c_dspgvx( integer itype, char jobz, char range, char uplo,
                    integer n, double* ap, double* bp, double vl,
                    double vu, integer il, integer iu, double abstol,
                    integer* m, double* w, double* z, integer ldz,
                    double* work, integer* iwork, integer* ifail )
{
    integer info;
    FFUNC(dspgvx,DSPGVX)(&itype, &jobz, &range, &uplo, &n, ap, bp, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork,
                            ifail, &info);
    return info;
}

static inline integer c_dsbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, double* ab, integer ldab, double* bb,
                   integer ldbb, double* w, double* z, integer ldz,
                   double* work )
{
    integer info;
    FFUNC(dsbgv,DSBGV)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &info);
    return info;
}

static inline integer c_dsbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, double* ab, integer ldab, double* bb,
                    integer ldbb, double* w, double* z, integer ldz,
                    double* work, integer lwork, integer* iwork,
                    integer liwork )
{
    integer info;
    FFUNC(dsbgvd,DSBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

static inline integer c_dsbgvx( char jobz, char range, char uplo, integer n,
                    integer ka, integer kb, double* ab,
                    integer ldab, double* bb, integer ldbb, double* q,
                    integer ldq, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, double* work, integer* iwork,
                    integer* ifail )
{
    integer info;
    FFUNC(dsbgvx,DSBGVX)(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, &vl, &vu, &il, &iu, &abstol,
                            m, w, z, &ldz, work, iwork, ifail, &info);
    return info;
}

static inline integer c_dggev( char jobvl, char jobvr, integer n, double* a,
                   integer lda, double* b, integer ldb, double* alphar,
                   double* alphai, double* beta, double* vl, integer ldvl,
                   double* vr, integer ldvr, double* work,
                   integer lwork )
{
    integer info;
    FFUNC(dggev,DGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    return info;
}

static inline integer c_dggevx( char balanc, char jobvl, char jobvr, char sense,
                    integer n, double* a, integer lda, double* b,
                    integer ldb, double* alphar, double* alphai,
                    double* beta, double* vl, integer ldvl, double* vr,
                    integer ldvr, integer* ilo, integer* ihi,
                    double* lscale, double* rscale, double* abnrm,
                    double* bbnrm, double* rconde, double* rcondv, double* work,
                    integer lwork, integer* iwork, logical* bwork )
{
    integer info;
    FFUNC(dggevx,DGGEVX)(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr,
                            ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, iwork, bwork, &info);
    return info;
}

static inline void c_dsfrk( char transr, char uplo, char trans, integer n,
                   integer k, double alpha, const double* a,
                   integer lda, double beta, double* c )
{
    FFUNC(dsfrk,DSFRK)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c);
}

static inline void c_dtfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, double alpha,
                   const double* a, double* b, integer ldb )
{
    FFUNC(dtfsm,DTFSM)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb);
}

static inline integer c_dtfttp( char transr, char uplo, integer n, const double* arf,
                    double* ap )
{
    integer info;
    FFUNC(dtfttp,DTFTTP)(&transr, &uplo, &n, arf, ap, &info);
    return info;
}

static inline integer c_dtfttr( char transr, char uplo, integer n, const double* arf,
                    double* a, integer lda )
{
    integer info;
    FFUNC(dtfttr,DTFTTR)(&transr, &uplo, &n, arf, a, &lda, &info);
    return info;
}

static inline integer c_dtpttf( char transr, char uplo, integer n, const double* ap,
                    double* arf )
{
    integer info;
    FFUNC(dtpttf,DTPTTF)(&transr, &uplo, &n, ap, arf, &info);
    return info;
}

static inline integer c_dtpttr( char uplo, integer n, const double* ap, double* a,
                    integer lda )
{
    integer info;
    FFUNC(dtpttr,DTPTTR)(&uplo, &n, ap, a, &lda, &info);
    return info;
}

static inline integer c_dtrttf( char transr, char uplo, integer n, const double* a,
                    integer lda, double* arf )
{
    integer info;
    FFUNC(dtrttf,DTRTTF)(&transr, &uplo, &n, a, &lda, arf, &info);
    return info;
}

static inline integer c_dtrttp( char uplo, integer n, const double* a, integer lda,
                    double* ap )
{
    integer info;
    FFUNC(dtrttp,DTRTTP)(&uplo, &n, a, &lda, ap, &info);
    return info;
}

static inline integer c_dgeqrfp( integer m, integer n, double* a, integer lda,
                     double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dgeqrfp,DGEQRFP)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

#ifdef __cplusplus
}
#endif

/*
 * #define more familiar names for the C versions
 */
#define dgetrf c_dgetrf
#define dgbtrf c_dgbtrf
#define dgttrf c_dgttrf
#define dpotrf c_dpotrf
#define dpstrf c_dpstrf
#define dpftrf c_dpftrf
#define dpptrf c_dpptrf
#define dpbtrf c_dpbtrf
#define dpttrf c_dpttrf
#define dsytrf c_dsytrf
#define dsptrf c_dsptrf
#define dgetrs c_dgetrs
#define dgbtrs c_dgbtrs
#define dgttrs c_dgttrs
#define dpotrs c_dpotrs
#define dpftrs c_dpftrs
#define dpptrs c_dpptrs
#define dpbtrs c_dpbtrs
#define dpttrs c_dpttrs
#define dsytrs c_dsytrs
#define dsptrs c_dsptrs
#define dtrtrs c_dtrtrs
#define dtptrs c_dtptrs
#define dtbtrs c_dtbtrs
#define dgecon c_dgecon
#define dgbcon c_dgbcon
#define dgtcon c_dgtcon
#define dpocon c_dpocon
#define dppcon c_dppcon
#define dpbcon c_dpbcon
#define dptcon c_dptcon
#define dsycon c_dsycon
#define dspcon c_dspcon
#define dtrcon c_dtrcon
#define dtpcon c_dtpcon
#define dtbcon c_dtbcon
#define dgerfs c_dgerfs
#define dgerfsx c_dgerfsx
#define dgbrfs c_dgbrfs
#define dgbrfsx c_dgbrfsx
#define dgtrfs c_dgtrfs
#define dporfs c_dporfs
#define dporfsx c_dporfsx
#define dpprfs c_dpprfs
#define dpbrfs c_dpbrfs
#define dptrfs c_dptrfs
#define dsyrfs c_dsyrfs
#define dsyrfsx c_dsyrfsx
#define dsprfs c_dsprfs
#define dtrrfs c_dtrrfs
#define dtprfs c_dtprfs
#define dtbrfs c_dtbrfs
#define dgetri c_dgetri
#define dpotri c_dpotri
#define dpftri c_dpftri
#define dpptri c_dpptri
#define dsytri c_dsytri
#define dsptri c_dsptri
#define dtrtri c_dtrtri
#define dtftri c_dtftri
#define dtptri c_dtptri
#define dgeequ c_dgeequ
#define dgeequb c_dgeequb
#define dgbequ c_dgbequ
#define dgbequb c_dgbequb
#define dpoequ c_dpoequ
#define dpoequb c_dpoequb
#define dppequ c_dppequ
#define dpbequ c_dpbequ
#define dsyequb c_dsyequb
#define dgesv c_dgesv
#define dsgesv c_dsgesv
#define dgesvx c_dgesvx
#define dgesvxx c_dgesvxx
#define dgbsv c_dgbsv
#define dgbsvx c_dgbsvx
#define dgbsvxx c_dgbsvxx
#define dgtsv c_dgtsv
#define dgtsvx c_dgtsvx
#define dposv c_dposv
#define dsposv c_dsposv
#define dposvx c_dposvx
#define dposvxx c_dposvxx
#define dppsv c_dppsv
#define dppsvx c_dppsvx
#define dpbsv c_dpbsv
#define dpbsvx c_dpbsvx
#define dptsv c_dptsv
#define dptsvx c_dptsvx
#define dsysv c_dsysv
#define dsysvx c_dsysvx
#define dsysvxx c_dsysvxx
#define dspsv c_dspsv
#define dspsvx c_dspsvx
#define dgeqrf c_dgeqrf
#define dgeqpf c_dgeqpf
#define dgeqp3 c_dgeqp3
#define dorgqr c_dorgqr
#define dormqr c_dormqr
#define dgelqf c_dgelqf
#define dorglq c_dorglq
#define dormlq c_dormlq
#define dgeqlf c_dgeqlf
#define dorgql c_dorgql
#define dormql c_dormql
#define dgerqf c_dgerqf
#define dorgrq c_dorgrq
#define dormrq c_dormrq
#define dtzrzf c_dtzrzf
#define dormrz c_dormrz
#define dggqrf c_dggqrf
#define dggrqf c_dggrqf
#define dgebrd c_dgebrd
#define dgbbrd c_dgbbrd
#define dorgbr c_dorgbr
#define dormbr c_dormbr
#define dbdsqr c_dbdsqr
#define dbdsdc c_dbdsdc
#define dsytrd c_dsytrd
#define dorgtr c_dorgtr
#define dormtr c_dormtr
#define dsptrd c_dsptrd
#define dopgtr c_dopgtr
#define dopmtr c_dopmtr
#define dsbtrd c_dsbtrd
#define dsterf c_dsterf
#define dsteqr c_dsteqr
#define dstemr c_dstemr
#define dstedc c_dstedc
#define dstegr c_dstegr
#define dpteqr c_dpteqr
#define dstebz c_dstebz
#define dstein c_dstein
#define ddisna c_ddisna
#define dsygst c_dsygst
#define dspgst c_dspgst
#define dsbgst c_dsbgst
#define dpbstf c_dpbstf
#define dgehrd c_dgehrd
#define dorghr c_dorghr
#define dormhr c_dormhr
#define dgebal c_dgebal
#define dgebak c_dgebak
#define dhseqr c_dhseqr
#define dhsein c_dhsein
#define dtrevc c_dtrevc
#define dtrsna c_dtrsna
#define dtrexc c_dtrexc
#define dtrsen c_dtrsen
#define dtrsyl c_dtrsyl
#define dgghrd c_dgghrd
#define dggbal c_dggbal
#define dggbak c_dggbak
#define dhgeqz c_dhgeqz
#define dtgevc c_dtgevc
#define dtgexc c_dtgexc
#define dtgsen c_dtgsen
#define dtgsyl c_dtgsyl
#define dtgsna c_dtgsna
#define dggsvp c_dggsvp
#define dtgsja c_dtgsja
#define dgels c_dgels
#define dgelsy c_dgelsy
#define dgelss c_dgelss
#define dgelsd c_dgelsd
#define dgglse c_dgglse
#define dggglm c_dggglm
#define dsyev c_dsyev
#define dsyevd c_dsyevd
#define dsyevx c_dsyevx
#define dsyevr c_dsyevr
#define dspev c_dspev
#define dspevd c_dspevd
#define dspevx c_dspevx
#define dsbev c_dsbev
#define dsbevd c_dsbevd
#define dsbevx c_dsbevx
#define dstev c_dstev
#define dstevd c_dstevd
#define dstevx c_dstevx
#define dstevr c_dstevr
#define dgeev c_dgeev
#define dgeevx c_dgeevx
#define dgesvd c_dgesvd
#define dgesdd c_dgesdd
#define dgejsv c_dgejsv
#define dgesvj c_dgesvj
#define dggsvd c_dggsvd
#define dsygv c_dsygv
#define dsygvd c_dsygvd
#define dsygvx c_dsygvx
#define dspgv c_dspgv
#define dspgvd c_dspgvd
#define dspgvx c_dspgvx
#define dsbgv c_dsbgv
#define dsbgvd c_dsbgvd
#define dsbgvx c_dsbgvx
#define dggev c_dggev
#define dggevx c_dggevx
#define dsfrk c_dsfrk
#define dtfsm c_dtfsm
#define dtfttp c_dtfttp
#define dtfttr c_dtfttr
#define dtpttf c_dtpttf
#define dtpttr c_dtpttr
#define dtrttf c_dtrttf
#define dtrttp c_dtrttp
#define dgeqrfp c_dgeqrfp

#endif
