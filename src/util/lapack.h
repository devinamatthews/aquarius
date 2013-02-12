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
void FFUNC(dgetrf,DGETRF)( integer* m, integer* n, double* a, integer* lda,
                    integer* ipiv, integer *info );

void FFUNC(dgbtrf,DGBTRF)( integer* m, integer* n, integer* kl,
                    integer* ku, double* ab, integer* ldab,
                    integer* ipiv, integer *info );

void FFUNC(dgttrf,DGTTRF)( integer* n, double* dl, double* d, double* du,
                    double* du2, integer* ipiv, integer *info );

void FFUNC(dpotrf,DPOTRF)( char* uplo, integer* n, double* a, integer* lda,
                    integer *info );

void FFUNC(dpstrf,DPSTRF)( char* uplo, integer* n, double* a, integer* lda,
                    integer* piv, integer* rank, double* tol,
                    double* work, integer *info );

void FFUNC(dpftrf,DPFTRF)( char* transr, char* uplo, integer* n, double* a,
                    integer *info );

void FFUNC(dpptrf,DPPTRF)( char* uplo, integer* n, double* ap, integer *info );

void FFUNC(dpbtrf,DPBTRF)( char* uplo, integer* n, integer* kd, double* ab,
                    integer* ldab, integer *info );

void FFUNC(dpttrf,DPTTRF)( integer* n, double* d, double* e, integer *info );

void FFUNC(dsytrf,DSYTRF)( char* uplo, integer* n, double* a, integer* lda,
                    integer* ipiv, double* work, integer* lwork,
                    integer *info );

void FFUNC(dsptrf,DSPTRF)( char* uplo, integer* n, double* ap, integer* ipiv,
                    integer *info );

void FFUNC(dgetrs,DGETRS)( char* trans, integer* n, integer* nrhs,
                    const double* a, integer* lda, const integer* ipiv,
                    double* b, integer* ldb, integer *info );

void FFUNC(dgbtrs,DGBTRS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const double* ab, integer* ldab,
                    const integer* ipiv, double* b, integer* ldb,
                    integer *info );

void FFUNC(dgttrs,DGTTRS)( char* trans, integer* n, integer* nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* du2, const integer* ipiv, double* b,
                    integer* ldb, integer *info );

void FFUNC(dpotrs,DPOTRS)( char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, double* b,
                    integer* ldb, integer *info );

void FFUNC(dpftrs,DPFTRS)( char* transr, char* uplo, integer* n, integer* nrhs,
                    const double* a, double* b, integer* ldb,
                    integer *info );

void FFUNC(dpptrs,DPPTRS)( char* uplo, integer* n, integer* nrhs,
                    const double* ap, double* b, integer* ldb,
                    integer *info );

void FFUNC(dpbtrs,DPBTRS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const double* ab, integer* ldab, double* b,
                    integer* ldb, integer *info );

void FFUNC(dpttrs,DPTTRS)( integer* n, integer* nrhs, const double* d,
                    const double* e, double* b, integer* ldb,
                    integer *info );

void FFUNC(dsytrs,DSYTRS)( char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, const integer* ipiv,
                    double* b, integer* ldb, integer *info );

void FFUNC(dsptrs,DSPTRS)( char* uplo, integer* n, integer* nrhs,
                    const double* ap, const integer* ipiv, double* b,
                    integer* ldb, integer *info );

void FFUNC(dtrtrs,DTRTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const double* a, integer* lda,
                    double* b, integer* ldb, integer *info );

void FFUNC(dtptrs,DTPTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const double* ap, double* b,
                    integer* ldb, integer *info );

void FFUNC(dtbtrs,DTBTRS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const double* ab,
                    integer* ldab, double* b, integer* ldb,
                    integer *info );

void FFUNC(dgecon,DGECON)( char* norm, integer* n, const double* a, integer* lda,
                    double* anorm, double* rcond, double* work,
                    integer* iwork, integer *info );

void FFUNC(dgbcon,DGBCON)( char* norm, integer* n, integer* kl, integer* ku,
                    const double* ab, integer* ldab, const integer* ipiv,
                    double* anorm, double* rcond, double* work,
                    integer* iwork, integer *info );

void FFUNC(dgtcon,DGTCON)( char* norm, integer* n, const double* dl,
                    const double* d, const double* du, const double* du2,
                    const integer* ipiv, double* anorm, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(dpocon,DPOCON)( char* uplo, integer* n, const double* a, integer* lda,
                    double* anorm, double* rcond, double* work,
                    integer* iwork, integer *info );

void FFUNC(dppcon,DPPCON)( char* uplo, integer* n, const double* ap, double* anorm,
                    double* rcond, double* work, integer* iwork,
                    integer *info );

void FFUNC(dpbcon,DPBCON)( char* uplo, integer* n, integer* kd, const double* ab,
                    integer* ldab, double* anorm, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(dptcon,DPTCON)( integer* n, const double* d, const double* e,
                    double* anorm, double* rcond, double* work,
                    integer *info );

void FFUNC(dsycon,DSYCON)( char* uplo, integer* n, const double* a, integer* lda,
                    const integer* ipiv, double* anorm, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(dspcon,DSPCON)( char* uplo, integer* n, const double* ap,
                    const integer* ipiv, double* anorm, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(dtrcon,DTRCON)( char* norm, char* uplo, char* diag, integer* n,
                    const double* a, integer* lda, double* rcond,
                    double* work, integer* iwork, integer *info );

void FFUNC(dtpcon,DTPCON)( char* norm, char* uplo, char* diag, integer* n,
                    const double* ap, double* rcond, double* work,
                    integer* iwork, integer *info );

void FFUNC(dtbcon,DTBCON)( char* norm, char* uplo, char* diag, integer* n,
                    integer* kd, const double* ab, integer* ldab,
                    double* rcond, double* work, integer* iwork,
                    integer *info );

void FFUNC(dgerfs,DGERFS)( char* trans, integer* n, integer* nrhs,
                    const double* a, integer* lda, const double* af,
                    integer* ldaf, const integer* ipiv, const double* b,
                    integer* ldb, double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dgerfsx,DGERFSX)( char* trans, char* equed, integer* n, integer* nrhs,
                     const double* a, integer* lda, const double* af,
                     integer* ldaf, const integer* ipiv, const double* r,
                     const double* c, const double* b, integer* ldb,
                     double* x, integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(dgbrfs,DGBRFS)( char* trans, integer* n, integer* kl, integer* ku,
                    integer* nrhs, const double* ab, integer* ldab,
                    const double* afb, integer* ldafb,
                    const integer* ipiv, const double* b, integer* ldb,
                    double* x, integer* ldx, double* ferr, double* berr,
                    double* work, integer* iwork, integer *info );

void FFUNC(dgbrfsx,DGBRFSX)( char* trans, char* equed, integer* n, integer* kl,
                     integer* ku, integer* nrhs, const double* ab,
                     integer* ldab, const double* afb, integer* ldafb,
                     const integer* ipiv, const double* r, const double* c,
                     const double* b, integer* ldb, double* x,
                     integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(dgtrfs,DGTRFS)( char* trans, integer* n, integer* nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* dlf, const double* df, const double* duf,
                    const double* du2, const integer* ipiv, const double* b,
                    integer* ldb, double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dporfs,DPORFS)( char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, const double* af,
                    integer* ldaf, const double* b, integer* ldb,
                    double* x, integer* ldx, double* ferr, double* berr,
                    double* work, integer* iwork, integer *info );

void FFUNC(dporfsx,DPORFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const double* a, integer* lda, const double* af,
                     integer* ldaf, const double* s, const double* b,
                     integer* ldb, double* x, integer* ldx, double* rcond,
                     double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, double* work,
                     integer* iwork, integer *info );

void FFUNC(dpprfs,DPPRFS)( char* uplo, integer* n, integer* nrhs,
                    const double* ap, const double* afp, const double* b,
                    integer* ldb, double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dpbrfs,DPBRFS)( char* uplo, integer* n, integer* kd, integer* nrhs,
                    const double* ab, integer* ldab, const double* afb,
                    integer* ldafb, const double* b, integer* ldb,
                    double* x, integer* ldx, double* ferr, double* berr,
                    double* work, integer* iwork, integer *info );

void FFUNC(dptrfs,DPTRFS)( integer* n, integer* nrhs, const double* d,
                    const double* e, const double* df, const double* ef,
                    const double* b, integer* ldb, double* x,
                    integer* ldx, double* ferr, double* berr, double* work,
                    integer *info );

void FFUNC(dsyrfs,DSYRFS)( char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, const double* af,
                    integer* ldaf, const integer* ipiv, const double* b,
                    integer* ldb, double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dsyrfsx,DSYRFSX)( char* uplo, char* equed, integer* n, integer* nrhs,
                     const double* a, integer* lda, const double* af,
                     integer* ldaf, const integer* ipiv, const double* s,
                     const double* b, integer* ldb, double* x,
                     integer* ldx, double* rcond, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(dsprfs,DSPRFS)( char* uplo, integer* n, integer* nrhs,
                    const double* ap, const double* afp, const integer* ipiv,
                    const double* b, integer* ldb, double* x,
                    integer* ldx, double* ferr, double* berr, double* work,
                    integer* iwork, integer *info );

void FFUNC(dtrrfs,DTRRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const double* a, integer* lda,
                    const double* b, integer* ldb, const double* x,
                    integer* ldx, double* ferr, double* berr, double* work,
                    integer* iwork, integer *info );

void FFUNC(dtprfs,DTPRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* nrhs, const double* ap, const double* b,
                    integer* ldb, const double* x, integer* ldx,
                    double* ferr, double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dtbrfs,DTBRFS)( char* uplo, char* trans, char* diag, integer* n,
                    integer* kd, integer* nrhs, const double* ab,
                    integer* ldab, const double* b, integer* ldb,
                    const double* x, integer* ldx, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dgetri,DGETRI)( integer* n, double* a, integer* lda,
                    const integer* ipiv, double* work, integer* lwork,
                    integer *info );

void FFUNC(dpotri,DPOTRI)( char* uplo, integer* n, double* a, integer* lda,
                    integer *info );

void FFUNC(dpftri,DPFTRI)( char* transr, char* uplo, integer* n, double* a,
                    integer *info );

void FFUNC(dpptri,DPPTRI)( char* uplo, integer* n, double* ap, integer *info );

void FFUNC(dsytri,DSYTRI)( char* uplo, integer* n, double* a, integer* lda,
                    const integer* ipiv, double* work, integer *info );

void FFUNC(dsptri,DSPTRI)( char* uplo, integer* n, double* ap,
                    const integer* ipiv, double* work, integer *info );

void FFUNC(dtrtri,DTRTRI)( char* uplo, char* diag, integer* n, double* a,
                    integer* lda, integer *info );

void FFUNC(dtftri,DTFTRI)( char* transr, char* uplo, char* diag, integer* n,
                    double* a, integer *info );

void FFUNC(dtptri,DTPTRI)( char* uplo, char* diag, integer* n, double* ap,
                    integer *info );

void FFUNC(dgeequ,DGEEQU)( integer* m, integer* n, const double* a,
                    integer* lda, double* r, double* c, double* rowcnd,
                    double* colcnd, double* amax, integer *info );

void FFUNC(dgeequb,DGEEQUB)( integer* m, integer* n, const double* a,
                     integer* lda, double* r, double* c, double* rowcnd,
                     double* colcnd, double* amax, integer *info );

void FFUNC(dgbequ,DGBEQU)( integer* m, integer* n, integer* kl,
                    integer* ku, const double* ab, integer* ldab,
                    double* r, double* c, double* rowcnd, double* colcnd,
                    double* amax, integer *info );

void FFUNC(dgbequb,DGBEQUB)( integer* m, integer* n, integer* kl,
                     integer* ku, const double* ab, integer* ldab,
                     double* r, double* c, double* rowcnd, double* colcnd,
                     double* amax, integer *info );

void FFUNC(dpoequ,DPOEQU)( integer* n, const double* a, integer* lda, double* s,
                    double* scond, double* amax, integer *info );

void FFUNC(dpoequb,DPOEQUB)( integer* n, const double* a, integer* lda, double* s,
                     double* scond, double* amax, integer *info );

void FFUNC(dppequ,DPPEQU)( char* uplo, integer* n, const double* ap, double* s,
                    double* scond, double* amax, integer *info );

void FFUNC(dpbequ,DPBEQU)( char* uplo, integer* n, integer* kd, const double* ab,
                    integer* ldab, double* s, double* scond, double* amax,
                    integer *info );

void FFUNC(dsyequb,DSYEQUB)( char* uplo, integer* n, const double* a,
                     integer* lda, double* s, double* scond, double* amax,
                     double* work, integer *info );

void FFUNC(dgesv,DGESV)( integer* n, integer* nrhs, double* a, integer* lda,
                   integer* ipiv, double* b, integer* ldb,
                   integer *info );

void FFUNC(dsgesv,DSGESV)( integer* n, integer* nrhs, double* a, integer* lda,
                    integer* ipiv, double* b, integer* ldb, double* x,
                    integer* ldx, double* work, float* swork,
                    integer* iter, integer *info );

void FFUNC(dgesvx,DGESVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    double* a, integer* lda, double* af, integer* ldaf,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer* ldb, double* x, integer* ldx,
                    double* rcond, double* ferr, double* berr, double* work,
                    integer* iwork, integer *info );

void FFUNC(dgesvxx,DGESVXX)( char* fact, char* trans, integer* n, integer* nrhs,
                     double* a, integer* lda, double* af, integer* ldaf,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer* ldb, double* x, integer* ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(dgbsv,DGBSV)( integer* n, integer* kl, integer* ku,
                   integer* nrhs, double* ab, integer* ldab,
                   integer* ipiv, double* b, integer* ldb,
                   integer *info );

void FFUNC(dgbsvx,DGBSVX)( char* fact, char* trans, integer* n, integer* kl,
                    integer* ku, integer* nrhs, double* ab,
                    integer* ldab, double* afb, integer* ldafb,
                    integer* ipiv, char* equed, double* r, double* c,
                    double* b, integer* ldb, double* x, integer* ldx,
                    double* rcond, double* ferr, double* berr, double* work,
                    integer* iwork, integer *info );

void FFUNC(dgbsvxx,DGBSVXX)( char* fact, char* trans, integer* n, integer* kl,
                     integer* ku, integer* nrhs, double* ab,
                     integer* ldab, double* afb, integer* ldafb,
                     integer* ipiv, char* equed, double* r, double* c,
                     double* b, integer* ldb, double* x, integer* ldx,
                     double* rcond, double* rpvgrw, double* berr,
                     integer* n_err_bnds, double* err_bnds_norm,
                     double* err_bnds_comp, integer* nparams, double* params,
                     double* work, integer* iwork, integer *info );

void FFUNC(dgtsv,DGTSV)( integer* n, integer* nrhs, double* dl, double* d,
                   double* du, double* b, integer* ldb, integer *info );

void FFUNC(dgtsvx,DGTSVX)( char* fact, char* trans, integer* n, integer* nrhs,
                    const double* dl, const double* d, const double* du,
                    double* dlf, double* df, double* duf, double* du2,
                    integer* ipiv, const double* b, integer* ldb,
                    double* x, integer* ldx, double* rcond, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dposv,DPOSV)( char* uplo, integer* n, integer* nrhs, double* a,
                   integer* lda, double* b, integer* ldb,
                   integer *info );

void FFUNC(dsposv,DSPOSV)( char* uplo, integer* n, integer* nrhs, double* a,
                    integer* lda, double* b, integer* ldb, double* x,
                    integer* ldx, double* work, float* swork,
                    integer* iter, integer *info );

void FFUNC(dposvx,DPOSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    double* a, integer* lda, double* af, integer* ldaf,
                    char* equed, double* s, double* b, integer* ldb,
                    double* x, integer* ldx, double* rcond, double* ferr,
                    double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dposvxx,DPOSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     double* a, integer* lda, double* af, integer* ldaf,
                     char* equed, double* s, double* b, integer* ldb,
                     double* x, integer* ldx, double* rcond, double* rpvgrw,
                     double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, double* work,
                     integer* iwork, integer *info );

void FFUNC(dppsv,DPPSV)( char* uplo, integer* n, integer* nrhs, double* ap,
                   double* b, integer* ldb, integer *info );

void FFUNC(dppsvx,DPPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    double* ap, double* afp, char* equed, double* s, double* b,
                    integer* ldb, double* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dpbsv,DPBSV)( char* uplo, integer* n, integer* kd, integer* nrhs,
                   double* ab, integer* ldab, double* b, integer* ldb,
                   integer *info );

void FFUNC(dpbsvx,DPBSVX)( char* fact, char* uplo, integer* n, integer* kd,
                    integer* nrhs, double* ab, integer* ldab, double* afb,
                    integer* ldafb, char* equed, double* s, double* b,
                    integer* ldb, double* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer* iwork,
                    integer *info );

void FFUNC(dptsv,DPTSV)( integer* n, integer* nrhs, double* d, double* e,
                   double* b, integer* ldb, integer *info );

void FFUNC(dptsvx,DPTSVX)( char* fact, integer* n, integer* nrhs,
                    const double* d, const double* e, double* df, double* ef,
                    const double* b, integer* ldb, double* x,
                    integer* ldx, double* rcond, double* ferr, double* berr,
                    double* work, integer *info );

void FFUNC(dsysv,DSYSV)( char* uplo, integer* n, integer* nrhs, double* a,
                   integer* lda, integer* ipiv, double* b,
                   integer* ldb, double* work, integer* lwork,
                   integer *info );

void FFUNC(dsysvx,DSYSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const double* a, integer* lda, double* af,
                    integer* ldaf, integer* ipiv, const double* b,
                    integer* ldb, double* x, integer* ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(dsysvxx,DSYSVXX)( char* fact, char* uplo, integer* n, integer* nrhs,
                     double* a, integer* lda, double* af, integer* ldaf,
                     integer* ipiv, char* equed, double* s, double* b,
                     integer* ldb, double* x, integer* ldx, double* rcond,
                     double* rpvgrw, double* berr, integer* n_err_bnds,
                     double* err_bnds_norm, double* err_bnds_comp,
                     integer* nparams, double* params, double* work,
                     integer* iwork, integer *info );

void FFUNC(dspsv,DSPSV)( char* uplo, integer* n, integer* nrhs, double* ap,
                   integer* ipiv, double* b, integer* ldb,
                   integer *info );

void FFUNC(dspsvx,DSPSVX)( char* fact, char* uplo, integer* n, integer* nrhs,
                    const double* ap, double* afp, integer* ipiv,
                    const double* b, integer* ldb, double* x,
                    integer* ldx, double* rcond, double* ferr, double* berr,
                    double* work, integer* iwork, integer *info );

void FFUNC(dgeqrf,DGEQRF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(dgeqpf,DGEQPF)( integer* m, integer* n, double* a, integer* lda,
                    integer* jpvt, double* tau, double* work,
                    integer *info );

void FFUNC(dgeqp3,DGEQP3)( integer* m, integer* n, double* a, integer* lda,
                    integer* jpvt, double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(dorgqr,DORGQR)( integer* m, integer* n, integer* k, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(dormqr,DORMQR)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(dgelqf,DGELQF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(dorglq,DORGLQ)( integer* m, integer* n, integer* k, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(dormlq,DORMLQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(dgeqlf,DGEQLF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(dorgql,DORGQL)( integer* m, integer* n, integer* k, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(dormql,DORMQL)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(dgerqf,DGERQF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(dorgrq,DORGRQ)( integer* m, integer* n, integer* k, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(dormrq,DORMRQ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(dtzrzf,DTZRZF)( integer* m, integer* n, double* a, integer* lda,
                    double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(dormrz,DORMRZ)( char* side, char* trans, integer* m, integer* n,
                    integer* k, integer* l, const double* a,
                    integer* lda, const double* tau, double* c,
                    integer* ldc, double* work, integer* lwork,
                    integer *info );

void FFUNC(dggqrf,DGGQRF)( integer* n, integer* m, integer* p, double* a,
                    integer* lda, double* taua, double* b, integer* ldb,
                    double* taub, double* work, integer* lwork,
                    integer *info );

void FFUNC(dggrqf,DGGRQF)( integer* m, integer* p, integer* n, double* a,
                    integer* lda, double* taua, double* b, integer* ldb,
                    double* taub, double* work, integer* lwork,
                    integer *info );

void FFUNC(dgebrd,DGEBRD)( integer* m, integer* n, double* a, integer* lda,
                    double* d, double* e, double* tauq, double* taup,
                    double* work, integer* lwork, integer *info );

void FFUNC(dgbbrd,DGBBRD)( char* vect, integer* m, integer* n, integer* ncc,
                    integer* kl, integer* ku, double* ab,
                    integer* ldab, double* d, double* e, double* q,
                    integer* ldq, double* pt, integer* ldpt, double* c,
                    integer* ldc, double* work, integer *info );

void FFUNC(dorgbr,DORGBR)( char* vect, integer* m, integer* n, integer* k,
                    double* a, integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(dormbr,DORMBR)( char* vect, char* side, char* trans, integer* m,
                    integer* n, integer* k, const double* a,
                    integer* lda, const double* tau, double* c,
                    integer* ldc, double* work, integer* lwork,
                    integer *info );

void FFUNC(dbdsqr,DBDSQR)( char* uplo, integer* n, integer* ncvt,
                    integer* nru, integer* ncc, double* d, double* e,
                    double* vt, integer* ldvt, double* u, integer* ldu,
                    double* c, integer* ldc, double* work,
                    integer *info );

void FFUNC(dbdsdc,DBDSDC)( char* uplo, char* compq, integer* n, double* d,
                    double* e, double* u, integer* ldu, double* vt,
                    integer* ldvt, double* q, integer* iq, double* work,
                    integer* iwork, integer *info );

void FFUNC(dsytrd,DSYTRD)( char* uplo, integer* n, double* a, integer* lda,
                    double* d, double* e, double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(dorgtr,DORGTR)( char* uplo, integer* n, double* a, integer* lda,
                    const double* tau, double* work, integer* lwork,
                    integer *info );

void FFUNC(dormtr,DORMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const double* a, integer* lda,
                    const double* tau, double* c, integer* ldc, double* work,
                    integer* lwork, integer *info );

void FFUNC(dsptrd,DSPTRD)( char* uplo, integer* n, double* ap, double* d, double* e,
                    double* tau, integer *info );

void FFUNC(dopgtr,DOPGTR)( char* uplo, integer* n, const double* ap,
                    const double* tau, double* q, integer* ldq, double* work,
                    integer *info );

void FFUNC(dopmtr,DOPMTR)( char* side, char* uplo, char* trans, integer* m,
                    integer* n, const double* ap, const double* tau,
                    double* c, integer* ldc, double* work,
                    integer *info );

void FFUNC(dsbtrd,DSBTRD)( char* vect, char* uplo, integer* n, integer* kd,
                    double* ab, integer* ldab, double* d, double* e,
                    double* q, integer* ldq, double* work,
                    integer *info );

void FFUNC(dsterf,DSTERF)( integer* n, double* d, double* e, integer *info );

void FFUNC(dsteqr,DSTEQR)( char* compz, integer* n, double* d, double* e, double* z,
                    integer* ldz, double* work, integer *info );

void FFUNC(dstemr,DSTEMR)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, integer* m, double* w, double* z,
                    integer* ldz, integer* nzc, integer* isuppz,
                    logical* tryrac, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dstedc,DSTEDC)( char* compz, integer* n, double* d, double* e, double* z,
                    integer* ldz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dstegr,DSTEGR)( char* jobz, char* range, integer* n, double* d,
                    double* e, double* vl, double* vu, integer* il,
                    integer* iu, double* abstol, integer* m, double* w,
                    double* z, integer* ldz, integer* isuppz,
                    double* work, integer* lwork, integer* iwork,
                    integer* liwork, integer *info );

void FFUNC(dpteqr,DPTEQR)( char* compz, integer* n, double* d, double* e, double* z,
                    integer* ldz, double* work, integer *info );

void FFUNC(dstebz,DSTEBZ)( char* range, char* order, integer* n, double* vl,
                    double* vu, integer* il, integer* iu, double* abstol,
                    const double* d, const double* e, integer* m,
                    integer* nsplit, double* w, integer* iblock,
                    integer* isplit, double* work, integer* iwork,
                    integer *info );

void FFUNC(dstein,DSTEIN)( integer* n, const double* d, const double* e,
                    integer* m, const double* w, const integer* iblock,
                    const integer* isplit, double* z, integer* ldz,
                    double* work, integer* iwork, integer* ifailv,
                    integer *info );

void FFUNC(ddisna,DDISNA)( char* job, integer* m, integer* n, const double* d,
                    double* sep, integer *info );

void FFUNC(dsygst,DSYGST)( integer* itype, char* uplo, integer* n, double* a,
                    integer* lda, const double* b, integer* ldb,
                    integer *info );

void FFUNC(dspgst,DSPGST)( integer* itype, char* uplo, integer* n, double* ap,
                    const double* bp, integer *info );

void FFUNC(dsbgst,DSBGST)( char* vect, char* uplo, integer* n, integer* ka,
                    integer* kb, double* ab, integer* ldab,
                    const double* bb, integer* ldbb, double* x,
                    integer* ldx, double* work, integer *info );

void FFUNC(dpbstf,DPBSTF)( char* uplo, integer* n, integer* kb, double* bb,
                    integer* ldbb, integer *info );

void FFUNC(dgehrd,DGEHRD)( integer* n, integer* ilo, integer* ihi, double* a,
                    integer* lda, double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(dorghr,DORGHR)( integer* n, integer* ilo, integer* ihi, double* a,
                    integer* lda, const double* tau, double* work,
                    integer* lwork, integer *info );

void FFUNC(dormhr,DORMHR)( char* side, char* trans, integer* m, integer* n,
                    integer* ilo, integer* ihi, const double* a,
                    integer* lda, const double* tau, double* c,
                    integer* ldc, double* work, integer* lwork,
                    integer *info );

void FFUNC(dgebal,DGEBAL)( char* job, integer* n, double* a, integer* lda,
                    integer* ilo, integer* ihi, double* scale,
                    integer *info );

void FFUNC(dgebak,DGEBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const double* scale, integer* m,
                    double* v, integer* ldv, integer *info );

void FFUNC(dhseqr,DHSEQR)( char* job, char* compz, integer* n, integer* ilo,
                    integer* ihi, double* h, integer* ldh, double* wr,
                    double* wi, double* z, integer* ldz, double* work,
                    integer* lwork, integer *info );

void FFUNC(dhsein,DHSEIN)( char* job, char* eigsrc, char* initv,
                    logical* select, integer* n, const double* h,
                    integer* ldh, double* wr, const double* wi, double* vl,
                    integer* ldvl, double* vr, integer* ldvr,
                    integer* mm, integer* m, double* work,
                    integer* ifaill, integer* ifailr, integer *info );

void FFUNC(dtrevc,DTREVC)( char* side, char* howmny, logical* select,
                    integer* n, const double* t, integer* ldt, double* vl,
                    integer* ldvl, double* vr, integer* ldvr,
                    integer* mm, integer* m, double* work,
                    integer *info );

void FFUNC(dtrsna,DTRSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const double* t, integer* ldt,
                    const double* vl, integer* ldvl, const double* vr,
                    integer* ldvr, double* s, double* sep, integer* mm,
                    integer* m, double* work, integer* ldwork,
                    integer* iwork, integer *info );

void FFUNC(dtrexc,DTREXC)( char* compq, integer* n, double* t, integer* ldt,
                    double* q, integer* ldq, integer* ifst,
                    integer* ilst, double* work, integer *info );

void FFUNC(dtrsen,DTRSEN)( char* job, char* compq, const logical* select,
                    integer* n, double* t, integer* ldt, double* q,
                    integer* ldq, double* wr, double* wi, integer* m,
                    double* s, double* sep, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dtrsyl,DTRSYL)( char* trana, char* tranb, integer* isgn, integer* m,
                    integer* n, const double* a, integer* lda,
                    const double* b, integer* ldb, double* c,
                    integer* ldc, double* scale, integer *info );

void FFUNC(dgghrd,DGGHRD)( char* compq, char* compz, integer* n, integer* ilo,
                    integer* ihi, double* a, integer* lda, double* b,
                    integer* ldb, double* q, integer* ldq, double* z,
                    integer* ldz, integer *info );

void FFUNC(dggbal,DGGBAL)( char* job, integer* n, double* a, integer* lda,
                    double* b, integer* ldb, integer* ilo,
                    integer* ihi, double* lscale, double* rscale,
                    double* work, integer *info );

void FFUNC(dggbak,DGGBAK)( char* job, char* side, integer* n, integer* ilo,
                    integer* ihi, const double* lscale, const double* rscale,
                    integer* m, double* v, integer* ldv,
                    integer *info );

void FFUNC(dhgeqz,DHGEQZ)( char* job, char* compq, char* compz, integer* n,
                    integer* ilo, integer* ihi, double* h,
                    integer* ldh, double* t, integer* ldt, double* alphar,
                    double* alphai, double* beta, double* q, integer* ldq,
                    double* z, integer* ldz, double* work, integer* lwork,
                    integer *info );

void FFUNC(dtgevc,DTGEVC)( char* side, char* howmny, const logical* select,
                    integer* n, const double* s, integer* lds,
                    const double* p, integer* ldp, double* vl,
                    integer* ldvl, double* vr, integer* ldvr,
                    integer* mm, integer* m, double* work,
                    integer *info );

void FFUNC(dtgexc,DTGEXC)( logical* wantq, logical* wantz, integer* n,
                    double* a, integer* lda, double* b, integer* ldb,
                    double* q, integer* ldq, double* z, integer* ldz,
                    integer* ifst, integer* ilst, double* work,
                    integer* lwork, integer *info );

void FFUNC(dtgsen,DTGSEN)( integer* ijob, logical* wantq,
                    logical* wantz, const logical* select,
                    integer* n, double* a, integer* lda, double* b,
                    integer* ldb, double* alphar, double* alphai,
                    double* beta, double* q, integer* ldq, double* z,
                    integer* ldz, integer* m, double* pl, double* pr,
                    double* dif, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dtgsyl,DTGSYL)( char* trans, integer* ijob, integer* m, integer* n,
                    const double* a, integer* lda, const double* b,
                    integer* ldb, double* c, integer* ldc,
                    const double* d, integer* ldd, const double* e,
                    integer* lde, double* f, integer* ldf, double* scale,
                    double* dif, double* work, integer* lwork,
                    integer* iwork, integer *info );

void FFUNC(dtgsna,DTGSNA)( char* job, char* howmny, const logical* select,
                    integer* n, const double* a, integer* lda,
                    const double* b, integer* ldb, const double* vl,
                    integer* ldvl, const double* vr, integer* ldvr,
                    double* s, double* dif, integer* mm, integer* m,
                    double* work, integer* lwork, integer* iwork,
                    integer *info );

void FFUNC(dggsvp,DGGSVP)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, double* a, integer* lda,
                    double* b, integer* ldb, double* tola, double* tolb,
                    integer* k, integer* l, double* u, integer* ldu,
                    double* v, integer* ldv, double* q, integer* ldq,
                    integer* iwork, double* tau, double* work,
                    integer *info );

void FFUNC(dtgsja,DTGSJA)( char* jobu, char* jobv, char* jobq, integer* m,
                    integer* p, integer* n, integer* k, integer* l,
                    double* a, integer* lda, double* b, integer* ldb,
                    double* tola, double* tolb, double* alpha, double* beta,
                    double* u, integer* ldu, double* v, integer* ldv,
                    double* q, integer* ldq, double* work,
                    integer* ncycle, integer *info );

void FFUNC(dgels,DGELS)( char* trans, integer* m, integer* n, integer* nrhs,
                   double* a, integer* lda, double* b, integer* ldb,
                   double* work, integer* lwork, integer *info );

void FFUNC(dgelsy,DGELSY)( integer* m, integer* n, integer* nrhs, double* a,
                    integer* lda, double* b, integer* ldb,
                    integer* jpvt, double* rcond, integer* rank,
                    double* work, integer* lwork, integer *info );

void FFUNC(dgelss,DGELSS)( integer* m, integer* n, integer* nrhs, double* a,
                    integer* lda, double* b, integer* ldb, double* s,
                    double* rcond, integer* rank, double* work,
                    integer* lwork, integer *info );

void FFUNC(dgelsd,DGELSD)( integer* m, integer* n, integer* nrhs, double* a,
                    integer* lda, double* b, integer* ldb, double* s,
                    double* rcond, integer* rank, double* work,
                    integer* lwork, integer* iwork, integer *info );

void FFUNC(dgglse,DGGLSE)( integer* m, integer* n, integer* p, double* a,
                    integer* lda, double* b, integer* ldb, double* c,
                    double* d, double* x, double* work, integer* lwork,
                    integer *info );

void FFUNC(dggglm,DGGGLM)( integer* n, integer* m, integer* p, double* a,
                    integer* lda, double* b, integer* ldb, double* d,
                    double* x, double* y, double* work, integer* lwork,
                    integer *info );

void FFUNC(dsyev,DSYEV)( char* jobz, char* uplo, integer* n, double* a,
                   integer* lda, double* w, double* work, integer* lwork,
                   integer *info );

void FFUNC(dsyevd,DSYEVD)( char* jobz, char* uplo, integer* n, double* a,
                    integer* lda, double* w, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dsyevx,DSYEVX)( char* jobz, char* range, char* uplo, integer* n,
                    double* a, integer* lda, double* vl, double* vu,
                    integer* il, integer* iu, double* abstol,
                    integer* m, double* w, double* z, integer* ldz,
                    double* work, integer* lwork, integer* iwork,
                    integer* ifail, integer *info );

void FFUNC(dsyevr,DSYEVR)( char* jobz, char* range, char* uplo, integer* n,
                    double* a, integer* lda, double* vl, double* vu,
                    integer* il, integer* iu, double* abstol,
                    integer* m, double* w, double* z, integer* ldz,
                    integer* isuppz, double* work, integer* lwork,
                    integer* iwork, integer* liwork, integer *info );

void FFUNC(dspev,DSPEV)( char* jobz, char* uplo, integer* n, double* ap, double* w,
                   double* z, integer* ldz, double* work, integer *info );

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
inline integer c_dgetrf( integer m, integer n, double* a, integer lda,
                    integer* ipiv )
{
    integer info;
    FFUNC(dgetrf,DGETRF)(&m, &n, a, &lda, ipiv, &info);
    return info;
}

inline integer c_dgbtrf( integer m, integer n, integer kl,
                    integer ku, double* ab, integer ldab,
                    integer* ipiv )
{
    integer info;
    FFUNC(dgbtrf,DGBTRF)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

inline integer c_dgttrf( integer n, double* dl, double* d, double* du,
                    double* du2, integer* ipiv )
{
    integer info;
    FFUNC(dgttrf,DGTTRF)(&n, dl, d, du, du2, ipiv, &info);
    return info;
}

inline integer c_dpotrf( char uplo, integer n, double* a, integer lda )
{
    integer info;
    FFUNC(dpotrf,DPOTRF)(&uplo, &n, a, &lda, &info);
    return info;
}

inline integer c_dpstrf( char uplo, integer n, double* a, integer lda,
                    integer* piv, integer* rank, double tol,
                    double* work )
{
    integer info;
    FFUNC(dpstrf,DPSTRF)(&uplo, &n, a, &lda, piv, rank, &tol, work, &info);
    return info;
}

inline integer c_dpftrf( char transr, char uplo, integer n, double* a )
{
    integer info;
    FFUNC(dpftrf,DPFTRF)(&transr, &uplo, &n, a, &info);
    return info;
}

inline integer c_dpptrf( char uplo, integer n, double* ap )
{
    integer info;
    FFUNC(dpptrf,DPPTRF)(&uplo, &n, ap, &info);
    return info;
}

inline integer c_dpbtrf( char uplo, integer n, integer kd, double* ab,
                    integer ldab )
{
    integer info;
    FFUNC(dpbtrf,DPBTRF)(&uplo, &n, &kd, ab, &ldab, &info);
    return info;
}

inline integer c_dpttrf( integer n, double* d, double* e )
{
    integer info;
    FFUNC(dpttrf,DPTTRF)(&n, d, e, &info);
    return info;
}

inline integer c_dsytrf( char uplo, integer n, double* a, integer lda,
                    integer* ipiv, double* work, integer* lwork )
{
    integer info;
    FFUNC(dsytrf,DSYTRF)(&uplo, &n, a, &lda, ipiv, work, lwork, &info);
    return info;
}

inline integer c_dsptrf( char uplo, integer n, double* ap, integer* ipiv )
{
    integer info;
    FFUNC(dsptrf,DSPTRF)(&uplo, &n, ap, ipiv, &info);
    return info;
}

inline integer c_dgetrs( char trans, integer n, integer nrhs,
                    const double* a, integer lda, const integer* ipiv,
                    double* b, integer ldb )
{
    integer info;
    FFUNC(dgetrs,DGETRS)(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

inline integer c_dgbtrs( char trans, integer n, integer kl, integer ku,
                    integer nrhs, const double* ab, integer ldab,
                    const integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dgbtrs,DGBTRS)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

inline integer c_dgttrs( char trans, integer n, integer nrhs,
                    const double* dl, const double* d, const double* du,
                    const double* du2, const integer* ipiv, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dgttrs,DGTTRS)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);
    return info;
}

inline integer c_dpotrs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dpotrs,DPOTRS)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

inline integer c_dpftrs( char transr, char uplo, integer n, integer nrhs,
                    const double* a, double* b, integer ldb )
{
    integer info;
    FFUNC(dpftrs,DPFTRS)(&transr, &uplo, &n, &nrhs, a, b, &ldb, &info);
    return info;
}

inline integer c_dpptrs( char uplo, integer n, integer nrhs,
                    const double* ap, double* b, integer ldb )
{
    integer info;
    FFUNC(dpptrs,DPPTRS)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

inline integer c_dpbtrs( char uplo, integer n, integer kd, integer nrhs,
                    const double* ab, integer ldab, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dpbtrs,DPBTRS)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

inline integer c_dpttrs( integer n, integer nrhs, const double* d,
                    const double* e, double* b, integer ldb )
{
    integer info;
    FFUNC(dpttrs,DPTTRS)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

inline integer c_dsytrs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const integer* ipiv,
                    double* b, integer ldb )
{
    integer info;
    FFUNC(dsytrs,DSYTRS)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

inline integer c_dsptrs( char uplo, integer n, integer nrhs,
                    const double* ap, const integer* ipiv, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dsptrs,DSPTRS)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

inline integer c_dtrtrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* a, integer lda,
                    double* b, integer ldb )
{
    integer info;
    FFUNC(dtrtrs,DTRTRS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

inline integer c_dtptrs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* ap, double* b,
                    integer ldb )
{
    integer info;
    FFUNC(dtptrs,DTPTRS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

inline integer c_dtbtrs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const double* ab,
                    integer ldab, double* b, integer ldb )
{
    integer info;
    FFUNC(dtbtrs,DTBTRS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

inline integer c_dgecon( char norm, integer n, const double* a, integer lda,
                    double anorm, double* rcond, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dgecon,DGECON)(&norm, &n, a, &lda, &anorm, rcond, work, iwork, &info);
    return info;
}

inline integer c_dgbcon( char norm, integer n, integer kl, integer ku,
                    const double* ab, integer ldab, const integer* ipiv,
                    double anorm, double* rcond, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dgbcon,DGBCON)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, iwork, &info);
    return info;
}

inline integer c_dgtcon( char norm, integer n, const double* dl,
                    const double* d, const double* du, const double* du2,
                    const integer* ipiv, double anorm, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dgtcon,DGTCON)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, iwork, &info);
    return info;
}

inline integer c_dpocon( char uplo, integer n, const double* a, integer lda,
                    double anorm, double* rcond, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dpocon,DPOCON)(&uplo, &n, a, &lda, &anorm, rcond, work, iwork, &info);
    return info;
}

inline integer c_dppcon( char uplo, integer n, const double* ap, double anorm,
                    double* rcond, double* work, integer* iwork )
{
    integer info;
    FFUNC(dppcon,DPPCON)(&uplo, &n, ap, &anorm, rcond, work, iwork, &info);
    return info;
}

inline integer c_dpbcon( char uplo, integer n, integer kd, const double* ab,
                    integer ldab, double anorm, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dpbcon,DPBCON)(&uplo, &n, &kd, ab, &ldab, &anorm, rcond, work, iwork, &info);
    return info;
}

inline integer c_dptcon( integer n, const double* d, const double* e,
                    double anorm, double* rcond, double* work )
{
    integer info;
    FFUNC(dptcon,DPTCON)(&n, d, e, &anorm, rcond, work, &info);
    return info;
}

inline integer c_dsycon( char uplo, integer n, const double* a, integer lda,
                    const integer* ipiv, double anorm, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dsycon,DSYCON)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, iwork, &info);
    return info;
}

inline integer c_dspcon( char uplo, integer n, const double* ap,
                    const integer* ipiv, double anorm, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dspcon,DSPCON)(&uplo, &n, ap, ipiv, &anorm, rcond, work, iwork, &info);
    return info;
}

inline integer c_dtrcon( char norm, char uplo, char diag, integer n,
                    const double* a, integer lda, double* rcond,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dtrcon,DTRCON)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, iwork, &info);
    return info;
}

inline integer c_dtpcon( char norm, char uplo, char diag, integer n,
                    const double* ap, double* rcond, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dtpcon,DTPCON)(&norm, &uplo, &diag, &n, ap, rcond, work, iwork, &info);
    return info;
}

inline integer c_dtbcon( char norm, char uplo, char diag, integer n,
                    integer kd, const double* ab, integer ldab,
                    double* rcond, double* work, integer* iwork )
{
    integer info;
    FFUNC(dtbcon,DTBCON)(&norm, &uplo, &diag, &n, &kd, ab, &ldab, rcond, work, iwork, &info);
    return info;
}

inline integer c_dgerfs( char trans, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dgerfs,DGERFS)(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dgerfsx( char trans, char equed, integer n, integer nrhs,
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

inline integer c_dgbrfs( char trans, integer n, integer kl, integer ku,
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

inline integer c_dgbrfsx( char trans, char equed, integer n, integer kl,
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

inline integer c_dgtrfs( char trans, integer n, integer nrhs,
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

inline integer c_dporfs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dporfs,DPORFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dporfsx( char uplo, char equed, integer n, integer nrhs,
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

inline integer c_dpprfs( char uplo, integer n, integer nrhs,
                    const double* ap, const double* afp, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dpprfs,DPPRFS)(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dpbrfs( char uplo, integer n, integer kd, integer nrhs,
                    const double* ab, integer ldab, const double* afb,
                    integer ldafb, const double* b, integer ldb,
                    double* x, integer ldx, double* ferr, double* berr,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dpbrfs,DPBRFS)(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dptrfs( integer n, integer nrhs, const double* d,
                    const double* e, const double* df, const double* ef,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* ferr, double* berr, double* work )
{
    integer info;
    FFUNC(dptrfs,DPTRFS)(&n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, &info);
    return info;
}

inline integer c_dsyrfs( char uplo, integer n, integer nrhs,
                    const double* a, integer lda, const double* af,
                    integer ldaf, const integer* ipiv, const double* b,
                    integer ldb, double* x, integer ldx, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dsyrfs,DSYRFS)(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dsyrfsx( char uplo, char equed, integer n, integer nrhs,
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

inline integer c_dsprfs( char uplo, integer n, integer nrhs,
                    const double* ap, const double* afp, const integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* ferr, double* berr, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dsprfs,DSPRFS)(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dtrrfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* a, integer lda,
                    const double* b, integer ldb, const double* x,
                    integer ldx, double* ferr, double* berr, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dtrrfs,DTRRFS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dtprfs( char uplo, char trans, char diag, integer n,
                    integer nrhs, const double* ap, const double* b,
                    integer ldb, const double* x, integer ldx,
                    double* ferr, double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dtprfs,DTPRFS)(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dtbrfs( char uplo, char trans, char diag, integer n,
                    integer kd, integer nrhs, const double* ab,
                    integer ldab, const double* b, integer ldb,
                    const double* x, integer ldx, double* ferr,
                    double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dtbrfs,DTBRFS)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dgetri( integer n, double* a, integer lda,
                    const integer* ipiv, double* work, integer lwork )
{
    integer info;
    FFUNC(dgetri,DGETRI)(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

inline integer c_dpotri( char uplo, integer n, double* a, integer lda )
{
    integer info;
    FFUNC(dpotri,DPOTRI)(&uplo, &n, a, &lda, &info);
    return info;
}

inline integer c_dpftri( char transr, char uplo, integer n, double* a )
{
    integer info;
    FFUNC(dpftri,DPFTRI)(&transr, &uplo, &n, a, &info);
    return info;
}

inline integer c_dpptri( char uplo, integer n, double* ap )
{
    integer info;
    FFUNC(dpptri,DPPTRI)(&uplo, &n, ap, &info);
    return info;
}

inline integer c_dsytri( char uplo, integer n, double* a, integer lda,
                    const integer* ipiv, double* work )
{
    integer info;
    FFUNC(dsytri,DSYTRI)(&uplo, &n, a, &lda, ipiv, work, &info);
    return info;
}

inline integer c_dsptri( char uplo, integer n, double* ap,
                    const integer* ipiv, double* work )
{
    integer info;
    FFUNC(dsptri,DSPTRI)(&uplo, &n, ap, ipiv, work, &info);
    return info;
}

inline integer c_dtrtri( char uplo, char diag, integer n, double* a,
                    integer lda )
{
    integer info;
    FFUNC(dtrtri,DTRTRI)(&uplo, &diag, &n, a, &lda, &info);
    return info;
}

inline integer c_dtftri( char transr, char uplo, char diag, integer n,
                    double* a )
{
    integer info;
    FFUNC(dtftri,DTFTRI)(&transr, &uplo, &diag, &n, a, &info);
    return info;
}

inline integer c_dtptri( char uplo, char diag, integer n, double* ap )
{
    integer info;
    FFUNC(dtptri,DTPTRI)(&uplo, &diag, &n, ap, &info);
    return info;
}

inline integer c_dgeequ( integer m, integer n, const double* a,
                    integer lda, double* r, double* c, double* rowcnd,
                    double* colcnd, double* amax )
{
    integer info;
    FFUNC(dgeequ,DGEEQU)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

inline integer c_dgeequb( integer m, integer n, const double* a,
                     integer lda, double* r, double* c, double* rowcnd,
                     double* colcnd, double* amax )
{
    integer info;
    FFUNC(dgeequb,DGEEQUB)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

inline integer c_dgbequ( integer m, integer n, integer kl,
                    integer ku, const double* ab, integer ldab,
                    double* r, double* c, double* rowcnd, double* colcnd,
                    double* amax )
{
    integer info;
    FFUNC(dgbequ,DGBEQU)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

inline integer c_dgbequb( integer m, integer n, integer kl,
                     integer ku, const double* ab, integer ldab,
                     double* r, double* c, double* rowcnd, double* colcnd,
                     double* amax )
{
    integer info;
    FFUNC(dgbequb,DGBEQUB)(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
    return info;
}

inline integer c_dpoequ( integer n, const double* a, integer lda, double* s,
                    double* scond, double* amax )
{
    integer info;
    FFUNC(dpoequ,DPOEQU)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

inline integer c_dpoequb( integer n, const double* a, integer lda, double* s,
                     double* scond, double* amax )
{
    integer info;
    FFUNC(dpoequb,DPOEQUB)(&n, a, &lda, s, scond, amax, &info);
    return info;
}

inline integer c_dppequ( char uplo, integer n, const double* ap, double* s,
                    double* scond, double* amax )
{
    integer info;
    FFUNC(dppequ,DPPEQU)(&uplo, &n, ap, s, scond, amax, &info);
    return info;
}

inline integer c_dpbequ( char uplo, integer n, integer kd, const double* ab,
                    integer ldab, double* s, double* scond, double* amax )
{
    integer info;
    FFUNC(dpbequ,DPBEQU)(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info);
    return info;
}

inline integer c_dsyequb( char uplo, integer n, const double* a,
                     integer lda, double* s, double* scond, double* amax,
                     double* work )
{
    integer info;
    FFUNC(dsyequb,DSYEQUB)(&uplo, &n, a, &lda, s, scond, amax, work, &info);
    return info;
}

inline integer c_dgesv( integer n, integer nrhs, double* a, integer lda,
                   integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dgesv,DGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

inline integer c_dsgesv( integer n, integer nrhs, double* a, integer lda,
                    integer* ipiv, double* b, integer ldb, double* x,
                    integer ldx, double* work, float* swork,
                    integer* iter )
{
    integer info;
    FFUNC(dsgesv,DSGESV)(&n, &nrhs, a, &lda, ipiv, b, &ldb, x, &ldx, work, swork, iter, &info);
    return info;
}

inline integer c_dgesvx( char fact, char trans, integer n, integer nrhs,
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

inline integer c_dgesvxx( char fact, char trans, integer n, integer nrhs,
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

inline integer c_dgbsv( integer n, integer kl, integer ku,
                   integer nrhs, double* ab, integer ldab,
                   integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dgbsv,DGBSV)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

inline integer c_dgbsvx( char fact, char trans, integer n, integer kl,
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

inline integer c_dgbsvxx( char fact, char trans, integer n, integer kl,
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

inline integer c_dgtsv( integer n, integer nrhs, double* dl, double* d,
                   double* du, double* b, integer ldb )
{
    integer info;
    FFUNC(dgtsv,DGTSV)(&n, &nrhs, dl, d, du, b, &ldb, &info);
    return info;
}

inline integer c_dgtsvx( char fact, char trans, integer n, integer nrhs,
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

inline integer c_dposv( char uplo, integer n, integer nrhs, double* a,
                   integer lda, double* b, integer ldb )
{
    integer info;
    FFUNC(dposv,DPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

inline integer c_dsposv( char uplo, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* x,
                    integer ldx, double* work, float* swork,
                    integer* iter )
{
    integer info;
    FFUNC(dsposv,DSPOSV)(&uplo, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, work, swork, iter, &info);
    return info;
}

inline integer c_dposvx( char fact, char uplo, integer n, integer nrhs,
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

inline integer c_dposvxx( char fact, char uplo, integer n, integer nrhs,
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

inline integer c_dppsv( char uplo, integer n, integer nrhs, double* ap,
                   double* b, integer ldb )
{
    integer info;
    FFUNC(dppsv,DPPSV)(&uplo, &n, &nrhs, ap, b, &ldb, &info);
    return info;
}

inline integer c_dppsvx( char fact, char uplo, integer n, integer nrhs,
                    double* ap, double* afp, char* equed, double* s, double* b,
                    integer ldb, double* x, integer ldx, double* rcond,
                    double* ferr, double* berr, double* work, integer* iwork )
{
    integer info;
    FFUNC(dppsvx,DPPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dpbsv( char uplo, integer n, integer kd, integer nrhs,
                   double* ab, integer ldab, double* b, integer ldb )
{
    integer info;
    FFUNC(dpbsv,DPBSV)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info);
    return info;
}

inline integer c_dpbsvx( char fact, char uplo, integer n, integer kd,
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

inline integer c_dptsv( integer n, integer nrhs, double* d, double* e,
                   double* b, integer ldb )
{
    integer info;
    FFUNC(dptsv,DPTSV)(&n, &nrhs, d, e, b, &ldb, &info);
    return info;
}

inline integer c_dptsvx( char fact, integer n, integer nrhs,
                    const double* d, const double* e, double* df, double* ef,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* rcond, double* ferr, double* berr,
                    double* work )
{
    integer info;
    FFUNC(dptsvx,DPTSVX)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, &info);
    return info;
}

inline integer c_dsysv( char uplo, integer n, integer nrhs, double* a,
                   integer lda, integer* ipiv, double* b,
                   integer ldb, double* work, integer lwork )
{
    integer info;
    FFUNC(dsysv,DSYSV)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info);
    return info;
}

inline integer c_dsysvx( char fact, char uplo, integer n, integer nrhs,
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

inline integer c_dsysvxx( char fact, char uplo, integer n, integer nrhs,
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

inline integer c_dspsv( char uplo, integer n, integer nrhs, double* ap,
                   integer* ipiv, double* b, integer ldb )
{
    integer info;
    FFUNC(dspsv,DSPSV)(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info);
    return info;
}

inline integer c_dspsvx( char fact, char uplo, integer n, integer nrhs,
                    const double* ap, double* afp, integer* ipiv,
                    const double* b, integer ldb, double* x,
                    integer ldx, double* rcond, double* ferr, double* berr,
                    double* work, integer* iwork )
{
    integer info;
    FFUNC(dspsvx,DSPSVX)(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, &info);
    return info;
}

inline integer c_dgeqrf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dgeqrf,DGEQRF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dgeqpf( integer m, integer n, double* a, integer lda,
                    integer* jpvt, double* tau, double* work )
{
    integer info;
    FFUNC(dgeqpf,DGEQPF)(&m, &n, a, &lda, jpvt, tau, work, &info);
    return info;
}

inline integer c_dgeqp3( integer m, integer n, double* a, integer lda,
                    integer* jpvt, double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgeqp3,DGEQP3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, &info);
    return info;
}

inline integer c_dorgqr( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorgqr,DORGQR)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dormqr( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormqr,DORMQR)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

inline integer c_dgelqf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dgelqf,DGELQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dorglq( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorglq,DORGLQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dormlq( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormlq,DORMLQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

inline integer c_dgeqlf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dgeqlf,DGEQLF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dorgql( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorgql,DORGQL)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dormql( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormql,DORMQL)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

inline integer c_dgerqf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dgerqf,DGERQF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dorgrq( integer m, integer n, integer k, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorgrq,DORGRQ)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dormrq( char side, char trans, integer m, integer n,
                    integer k, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormrq,DORMRQ)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

inline integer c_dtzrzf( integer m, integer n, double* a, integer lda,
                    double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dtzrzf,DTZRZF)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dormrz( char side, char trans, integer m, integer n,
                    integer k, integer l, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc, double* work, integer lwork )
{
    integer info;
    FFUNC(dormrz,DORMRZ)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

inline integer c_dggqrf( integer n, integer m, integer p, double* a,
                    integer lda, double* taua, double* b, integer ldb,
                    double* taub, double* work, integer lwork )
{
    integer info;
    FFUNC(dggqrf,DGGQRF)(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    return info;
}

inline integer c_dggrqf( integer m, integer p, integer n, double* a,
                    integer lda, double* taua, double* b, integer ldb,
                    double* taub, double* work, integer lwork )
{
    integer info;
    FFUNC(dggrqf,DGGRQF)(&m, &p, &n, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
    return info;
}

inline integer c_dgebrd( integer m, integer n, double* a, integer lda,
                    double* d, double* e, double* tauq, double* taup,
                    double* work, integer lwork )
{
    integer info;
    FFUNC(dgebrd,DGEBRD)(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
    return info;
}

inline integer c_dgbbrd( char vect, integer m, integer n, integer ncc,
                    integer kl, integer ku, double* ab,
                    integer ldab, double* d, double* e, double* q,
                    integer ldq, double* pt, integer ldpt, double* c,
                    integer ldc, double* work )
{
    integer info;
    FFUNC(dgbbrd,DGBBRD)(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info);
    return info;
}

inline integer c_dorgbr( char vect, integer m, integer n, integer k,
                    double* a, integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorgbr,DORGBR)(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dormbr( char vect, char side, char trans, integer m,
                    integer n, integer k, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc, double* work, integer lwork )
{
    integer info;
    FFUNC(dormbr,DORMBR)(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

inline integer c_dbdsqr( char uplo, integer n, integer ncvt,
                    integer nru, integer ncc, double* d, double* e,
                    double* vt, integer ldvt, double* u, integer ldu,
                    double* c, integer ldc, double* work )
{
    integer info;
    FFUNC(dbdsqr,DBDSQR)(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info);
    return info;
}

inline integer c_dbdsdc( char uplo, char compq, integer n, double* d,
                    double* e, double* u, integer ldu, double* vt,
                    integer ldvt, double* q, integer* iq, double* work,
                    integer* iwork )
{
    integer info;
    FFUNC(dbdsdc,DBDSDC)(&uplo, &compq, &n, d, e, u, &ldu, vt, &ldvt, q, iq, work, iwork, &info);
    return info;
}

inline integer c_dsytrd( char uplo, integer n, double* a, integer lda,
                    double* d, double* e, double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dsytrd,DSYTRD)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info);
    return info;
}

inline integer c_dorgtr( char uplo, integer n, double* a, integer lda,
                    const double* tau, double* work, integer lwork )
{
    integer info;
    FFUNC(dorgtr,DORGTR)(&uplo, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dormtr( char side, char uplo, char trans, integer m,
                    integer n, const double* a, integer lda,
                    const double* tau, double* c, integer ldc, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dormtr,DORMTR)(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

inline integer c_dsptrd( char uplo, integer n, double* ap, double* d, double* e,
                    double* tau )
{
    integer info;
    FFUNC(dsptrd,DSPTRD)(&uplo, &n, ap, d, e, tau, &info);
    return info;
}

inline integer c_dopgtr( char uplo, integer n, const double* ap,
                    const double* tau, double* q, integer ldq, double* work )
{
    integer info;
    FFUNC(dopgtr,DOPGTR)(&uplo, &n, ap, tau, q, &ldq, work, &info);
    return info;
}

inline integer c_dopmtr( char side, char uplo, char trans, integer m,
                    integer n, const double* ap, const double* tau,
                    double* c, integer ldc, double* work )
{
    integer info;
    FFUNC(dopmtr,DOPMTR)(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info);
    return info;
}

inline integer c_dsbtrd( char vect, char uplo, integer n, integer kd,
                    double* ab, integer ldab, double* d, double* e,
                    double* q, integer ldq, double* work )
{
    integer info;
    FFUNC(dsbtrd,DSBTRD)(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info);
    return info;
}

inline integer c_dsterf( integer n, double* d, double* e )
{
    integer info;
    FFUNC(dsterf,DSTERF)(&n, d, e, &info);
    return info;
}

inline integer c_dsteqr( char compz, integer n, double* d, double* e, double* z,
                    integer ldz, double* work )
{
    integer info;
    FFUNC(dsteqr,DSTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    return info;
}

inline integer c_dstemr( char jobz, char range, integer n, double* d,
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

inline integer c_dstedc( char compz, integer n, double* d, double* e, double* z,
                    integer ldz, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dstedc,DSTEDC)(&compz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline integer c_dstegr( char jobz, char range, integer n, double* d,
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

inline integer c_dpteqr( char compz, integer n, double* d, double* e, double* z,
                    integer ldz, double* work )
{
    integer info;
    FFUNC(dpteqr,DPTEQR)(&compz, &n, d, e, z, &ldz, work, &info);
    return info;
}

inline integer c_dstebz( char range, char order, integer n, double vl,
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

inline integer c_dstein( integer n, const double* d, const double* e,
                    integer m, const double* w, const integer* iblock,
                    const integer* isplit, double* z, integer ldz,
                    double* work, integer* iwork, integer* ifailv )
{
    integer info;
    FFUNC(dstein,DSTEIN)(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifailv, &info);
    return info;
}

inline integer c_ddisna( char job, integer m, integer n, const double* d,
                    double* sep )
{
    integer info;
    FFUNC(ddisna,DDISNA)(&job, &m, &n, d, sep, &info);
    return info;
}

inline integer c_dsygst( integer itype, char uplo, integer n, double* a,
                    integer lda, const double* b, integer ldb )
{
    integer info;
    FFUNC(dsygst,DSYGST)(&itype, &uplo, &n, a, &lda, b, &ldb, &info);
    return info;
}

inline integer c_dspgst( integer itype, char uplo, integer n, double* ap,
                    const double* bp )
{
    integer info;
    FFUNC(dspgst,DSPGST)(&itype, &uplo, &n, ap, bp, &info);
    return info;
}

inline integer c_dsbgst( char vect, char uplo, integer n, integer ka,
                    integer kb, double* ab, integer ldab,
                    const double* bb, integer ldbb, double* x,
                    integer ldx, double* work )
{
    integer info;
    FFUNC(dsbgst,DSBGST)(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, &info);
    return info;
}

inline integer c_dpbstf( char uplo, integer n, integer kb, double* bb,
                    integer ldbb )
{
    integer info;
    FFUNC(dpbstf,DPBSTF)(&uplo, &n, &kb, bb, &ldbb, &info);
    return info;
}

inline integer c_dgehrd( integer n, integer ilo, integer ihi, double* a,
                    integer lda, double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgehrd,DGEHRD)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dorghr( integer n, integer ilo, integer ihi, double* a,
                    integer lda, const double* tau, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dorghr,DORGHR)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
    return info;
}

inline integer c_dormhr( char side, char trans, integer m, integer n,
                    integer ilo, integer ihi, const double* a,
                    integer lda, const double* tau, double* c,
                    integer ldc, double* work, integer lwork )
{
    integer info;
    FFUNC(dormhr,DORMHR)(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info);
    return info;
}

inline integer c_dgebal( char job, integer n, double* a, integer lda,
                    integer* ilo, integer* ihi, double* scale )
{
    integer info;
    FFUNC(dgebal,DGEBAL)(&job, &n, a, &lda, ilo, ihi, scale, &info);
    return info;
}

inline integer c_dgebak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* scale, integer m,
                    double* v, integer ldv )
{
    integer info;
    FFUNC(dgebak,DGEBAK)(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info);
    return info;
}

inline integer c_dhseqr( char job, char compz, integer n, integer ilo,
                    integer ihi, double* h, integer ldh, double* wr,
                    double* wi, double* z, integer ldz, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dhseqr,DHSEQR)(&job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, work, &lwork, &info);
    return info;
}

inline integer c_dhsein( char job, char eigsrc, char initv,
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

inline integer c_dtrevc( char side, char howmny, logical* select,
                    integer n, const double* t, integer ldt, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer* m, double* work )
{
    integer info;
    FFUNC(dtrevc,DTREVC)(&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, m, work, &info);
    return info;
}

inline integer c_dtrsna( char job, char howmny, const logical* select,
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

inline integer c_dtrexc( char compq, integer n, double* t, integer ldt,
                    double* q, integer ldq, integer* ifst,
                    integer* ilst, double* work )
{
    integer info;
    FFUNC(dtrexc,DTREXC)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, work, &info);
    return info;
}

inline integer c_dtrsen( char job, char compq, const logical* select,
                    integer n, double* t, integer ldt, double* q,
                    integer ldq, double* wr, double* wi, integer* m,
                    double* s, double* sep, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dtrsen,DTRSEN)(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, m, s, sep, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline integer c_dtrsyl( char trana, char tranb, integer isgn, integer m,
                    integer n, const double* a, integer lda,
                    const double* b, integer ldb, double* c,
                    integer ldc, double* scale )
{
    integer info;
    FFUNC(dtrsyl,DTRSYL)(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, &info);
    return info;
}

inline integer c_dgghrd( char compq, char compz, integer n, integer ilo,
                    integer ihi, double* a, integer lda, double* b,
                    integer ldb, double* q, integer ldq, double* z,
                    integer ldz )
{
    integer info;
    FFUNC(dgghrd,DGGHRD)(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info);
    return info;
}

inline integer c_dggbal( char job, integer n, double* a, integer lda,
                    double* b, integer ldb, integer* ilo,
                    integer* ihi, double* lscale, double* rscale,
                    double* work )
{
    integer info;
    FFUNC(dggbal,DGGBAL)(&job, &n, a, &lda, b, &ldb, ilo, ihi, lscale, rscale, work, &info);
    return info;
}

inline integer c_dggbak( char job, char side, integer n, integer ilo,
                    integer ihi, const double* lscale, const double* rscale,
                    integer m, double* v, integer ldv )
{
    integer info;
    FFUNC(dggbak,DGGBAK)(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info);
    return info;
}

inline integer c_dhgeqz( char job, char compq, char compz, integer n,
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

inline integer c_dtgevc( char side, char howmny, const logical* select,
                    integer n, const double* s, integer lds,
                    const double* p, integer ldp, double* vl,
                    integer ldvl, double* vr, integer ldvr,
                    integer mm, integer* m, double* work )
{
    integer info;
    FFUNC(dtgevc,DTGEVC)(&side, &howmny, select, &n, s, &lds, p, &ldp, vl, &ldvl, vr, &ldvr, &mm, m, work, &info);
    return info;
}

inline integer c_dtgexc( logical wantq, logical wantz, integer n,
                    double* a, integer lda, double* b, integer ldb,
                    double* q, integer ldq, double* z, integer ldz,
                    integer* ifst, integer* ilst, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dtgexc,DTGEXC)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, work, &lwork, &info);
    return info;
}

inline integer c_dtgsen( integer ijob, logical wantq,
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

inline integer c_dtgsyl( char trans, integer ijob, integer m, integer n,
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

inline integer c_dtgsna( char job, char howmny, const logical* select,
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

inline integer c_dggsvp( char jobu, char jobv, char jobq, integer m,
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

inline integer c_dtgsja( char jobu, char jobv, char jobq, integer m,
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

inline integer c_dgels( char trans, integer m, integer n, integer nrhs,
                   double* a, integer lda, double* b, integer ldb,
                   double* work, integer lwork )
{
    integer info;
    FFUNC(dgels,DGELS)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

inline integer c_dgelsy( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb,
                    integer* jpvt, double rcond, integer* rank,
                    double* work, integer lwork )
{
    integer info;
    FFUNC(dgelsy,DGELSY)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, rank, work, &lwork, &info);
    return info;
}

inline integer c_dgelss( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* s,
                    double rcond, integer* rank, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgelss,DGELSS)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, &info);
    return info;
}

inline integer c_dgelsd( integer m, integer n, integer nrhs, double* a,
                    integer lda, double* b, integer ldb, double* s,
                    double rcond, integer* rank, double* work,
                    integer lwork, integer* iwork )
{
    integer info;
    FFUNC(dgelsd,DGELSD)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, iwork, &info);
    return info;
}

inline integer c_dgglse( integer m, integer n, integer p, double* a,
                    integer lda, double* b, integer ldb, double* c,
                    double* d, double* x, double* work, integer lwork )
{
    integer info;
    FFUNC(dgglse,DGGLSE)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
    return info;
}

inline integer c_dggglm( integer n, integer m, integer p, double* a,
                    integer lda, double* b, integer ldb, double* d,
                    double* x, double* y, double* work, integer lwork )
{
    integer info;
    FFUNC(dggglm,DGGGLM)(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
    return info;
}

inline integer c_dsyev( char jobz, char uplo, integer n, double* a,
                   integer lda, double* w, double* work, integer lwork )
{
    integer info;
    FFUNC(dsyev,DSYEV)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    return info;
}

inline integer c_dsyevd( char jobz, char uplo, integer n, double* a,
                    integer lda, double* w, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dsyevd,DSYEVD)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline integer c_dsyevx( char jobz, char range, char uplo, integer n,
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

inline integer c_dsyevr( char jobz, char range, char uplo, integer n,
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

inline integer c_dspev( char jobz, char uplo, integer n, double* ap, double* w,
                   double* z, integer ldz, double* work )
{
    integer info;
    FFUNC(dspev,DSPEV)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &info);
    return info;
}

inline integer c_dspevd( char jobz, char uplo, integer n, double* ap,
                    double* w, double* z, integer ldz, double* work,
                    integer lwork, integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dspevd,DSPEVD)(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline integer c_dspevx( char jobz, char range, char uplo, integer n,
                    double* ap, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, double* work, integer* iwork,
                    integer* ifail )
{
    integer info;
    FFUNC(dspevx,DSPEVX)(&jobz, &range, &uplo, &n, ap, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, &info);
    return info;
}

inline integer c_dsbev( char jobz, char uplo, integer n, integer kd,
                   double* ab, integer ldab, double* w, double* z,
                   integer ldz, double* work )
{
    integer info;
    FFUNC(dsbev,DSBEV)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info);
    return info;
}

inline integer c_dsbevd( char jobz, char uplo, integer n, integer kd,
                    double* ab, integer ldab, double* w, double* z,
                    integer ldz, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dsbevd,DSBEVD)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline integer c_dsbevx( char jobz, char range, char uplo, integer n,
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

inline integer c_dstev( char jobz, integer n, double* d, double* e, double* z,
                   integer ldz, double* work )
{
    integer info;
    FFUNC(dstev,DSTEV)(&jobz, &n, d, e, z, &ldz, work, &info);
    return info;
}

inline integer c_dstevd( char jobz, integer n, double* d, double* e, double* z,
                    integer ldz, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dstevd,DSTEVD)(&jobz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline integer c_dstevx( char jobz, char range, integer n, double* d,
                    double* e, double vl, double vu, integer il,
                    integer iu, double abstol, integer* m, double* w,
                    double* z, integer ldz, double* work, integer* iwork,
                    integer* ifail )
{
    integer info;
    FFUNC(dstevx,DSTEVX)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, &info);
    return info;
}

inline integer c_dstevr( char jobz, char range, integer n, double* d,
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

inline integer c_dgeev( char jobvl, char jobvr, integer n, double* a,
                   integer lda, double* wr, double* wi, double* vl,
                   integer ldvl, double* vr, integer ldvr, double* work,
                   integer lwork )
{
    integer info;
    FFUNC(dgeev,DGEEV)(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    return info;
}

inline integer c_dgeevx( char balanc, char jobvl, char jobvr, char sense,
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

inline integer c_dgesvd( char jobu, char jobvt, integer m, integer n,
                    double* a, integer lda, double* s, double* u,
                    integer ldu, double* vt, integer ldvt, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgesvd,DGESVD)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
    return info;
}

inline integer c_dgesdd( char jobz, integer m, integer n, double* a,
                    integer lda, double* s, double* u, integer ldu,
                    double* vt, integer ldvt, double* work,
                    integer lwork, integer* iwork )
{
    integer info;
    FFUNC(dgesdd,DGESDD)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
    return info;
}

inline integer c_dgejsv( char joba, char jobu, char jobv, char jobr, char jobt,
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

inline integer c_dgesvj( char joba, char jobu, char jobv, integer m,
                    integer n, double* a, integer lda, double* sva,
                    integer mv, double* v, integer ldv, double* work,
                    integer lwork )
{
    integer info;
    FFUNC(dgesvj,DGESVJ)(&joba, &jobu, &jobv, &m, &n, a, &lda, sva, &mv, v, &ldv, work, &lwork, &info);
    return info;
}

inline integer c_dggsvd( char jobu, char jobv, char jobq, integer m,
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

inline integer c_dsygv( integer itype, char jobz, char uplo, integer n,
                   double* a, integer lda, double* b, integer ldb,
                   double* w, double* work, integer lwork )
{
    integer info;
    FFUNC(dsygv,DSYGV)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info);
    return info;
}

inline integer c_dsygvd( integer itype, char jobz, char uplo, integer n,
                    double* a, integer lda, double* b, integer ldb,
                    double* w, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dsygvd,DSYGVD)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline integer c_dsygvx( integer itype, char jobz, char range, char uplo,
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

inline integer c_dspgv( integer itype, char jobz, char uplo, integer n,
                   double* ap, double* bp, double* w, double* z,
                   integer ldz, double* work )
{
    integer info;
    FFUNC(dspgv,DSPGV)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &info);
    return info;
}

inline integer c_dspgvd( integer itype, char jobz, char uplo, integer n,
                    double* ap, double* bp, double* w, double* z,
                    integer ldz, double* work, integer lwork,
                    integer* iwork, integer liwork )
{
    integer info;
    FFUNC(dspgvd,DSPGVD)(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline integer c_dspgvx( integer itype, char jobz, char range, char uplo,
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

inline integer c_dsbgv( char jobz, char uplo, integer n, integer ka,
                   integer kb, double* ab, integer ldab, double* bb,
                   integer ldbb, double* w, double* z, integer ldz,
                   double* work )
{
    integer info;
    FFUNC(dsbgv,DSBGV)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &info);
    return info;
}

inline integer c_dsbgvd( char jobz, char uplo, integer n, integer ka,
                    integer kb, double* ab, integer ldab, double* bb,
                    integer ldbb, double* w, double* z, integer ldz,
                    double* work, integer lwork, integer* iwork,
                    integer liwork )
{
    integer info;
    FFUNC(dsbgvd,DSBGVD)(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, iwork, &liwork, &info);
    return info;
}

inline integer c_dsbgvx( char jobz, char range, char uplo, integer n,
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

inline integer c_dggev( char jobvl, char jobvr, integer n, double* a,
                   integer lda, double* b, integer ldb, double* alphar,
                   double* alphai, double* beta, double* vl, integer ldvl,
                   double* vr, integer ldvr, double* work,
                   integer lwork )
{
    integer info;
    FFUNC(dggev,DGGEV)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    return info;
}

inline integer c_dggevx( char balanc, char jobvl, char jobvr, char sense,
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

inline void c_dsfrk( char transr, char uplo, char trans, integer n,
                   integer k, double alpha, const double* a,
                   integer lda, double beta, double* c )
{
    FFUNC(dsfrk,DSFRK)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c);
}

inline void c_dtfsm( char transr, char side, char uplo, char trans,
                   char diag, integer m, integer n, double alpha,
                   const double* a, double* b, integer ldb )
{
    FFUNC(dtfsm,DTFSM)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb);
}

inline integer c_dtfttp( char transr, char uplo, integer n, const double* arf,
                    double* ap )
{
    integer info;
    FFUNC(dtfttp,DTFTTP)(&transr, &uplo, &n, arf, ap, &info);
    return info;
}

inline integer c_dtfttr( char transr, char uplo, integer n, const double* arf,
                    double* a, integer lda )
{
    integer info;
    FFUNC(dtfttr,DTFTTR)(&transr, &uplo, &n, arf, a, &lda, &info);
    return info;
}

inline integer c_dtpttf( char transr, char uplo, integer n, const double* ap,
                    double* arf )
{
    integer info;
    FFUNC(dtpttf,DTPTTF)(&transr, &uplo, &n, ap, arf, &info);
    return info;
}

inline integer c_dtpttr( char uplo, integer n, const double* ap, double* a,
                    integer lda )
{
    integer info;
    FFUNC(dtpttr,DTPTTR)(&uplo, &n, ap, a, &lda, &info);
    return info;
}

inline integer c_dtrttf( char transr, char uplo, integer n, const double* a,
                    integer lda, double* arf )
{
    integer info;
    FFUNC(dtrttf,DTRTTF)(&transr, &uplo, &n, a, &lda, arf, &info);
    return info;
}

inline integer c_dtrttp( char uplo, integer n, const double* a, integer lda,
                    double* ap )
{
    integer info;
    FFUNC(dtrttp,DTRTTP)(&uplo, &n, a, &lda, ap, &info);
    return info;
}

inline integer c_dgeqrfp( integer m, integer n, double* a, integer lda,
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
