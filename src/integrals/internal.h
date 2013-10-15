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

#ifndef _AQUARIUS_INTEGRALS_INTERNAL_H_
#define _AQUARIUS_INTEGRALS_INTERNAL_H_

#ifdef _OPENMP
#include "omp.h"
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

#include <assert.h>

#include "util/util.h"
#include "util/math_ext.h"
#include "util/lapack.h"
#include "memory/memory.h"
#include "time/c_time.h"

#ifdef __cplusplus
extern "C"
{
#endif

// elementdata.c

typedef struct
{
    double charge;      /**< Charge of this element in au. */
    const int Z;        /**< Atomic number. */
    const int isotope;  /**< Isotope number of this element. This is not e.g. the nucleon number, but an index
                             into the spin, mass, and nucleon arrays. */
    const char symbol[3]; /**< Symbol of this element, e.g. H or Li. */
    const char name[14];   /**< Full name of this element. */
    const int spin[25];    /**< Spin of each isotope for this element, given in units of \f$\frac{\hbar}{2}\f$. */
    const double mass[25]; /**< Mass of each isotope for this element, given in units of amu. */
    const int nucleon[25]; /**< Nucleon number (P+N) of each isotope for this element. */
} element_t;

extern const element_t elements[];

// ishprim.c

void ishprim(int la, int lb, int lc, int ld,
             const double* posa, const double* posb, const double* posc, const double* posd,
             double za, double zb, double zc, double zd, double* integrals,
             double* xtable);

// momprim.c

void momprim(int la, int lb, const double* posa, const double* posb, double za, double zb,
             int lc, const double* posc, double* integrals, double* table);

// osprim.c

void osprim(int la, int lb, int lc, int ld,
            const double* posa, const double* posb, const double* posc, const double* posd,
            double za, double zb, double zc, double zd, double* integrals,
            double* xtable);

// osinv.c

void osinv(int la, int lb, int lc, int ld,
           const double* posa, const double* posb, const double* posc, const double* posd,
           int na, int nb, int nc, int nd,
           const double* za, const double* zb, const double* zc, const double* zd, double* integrals);

// vrr.c

void vrr(int le0, int le1, int lf0, int lf1,
         const double* posa, const double* posb, const double* posc, const double* posd,
         int na, int nb, int nc, int nd,
         const double* za, const double* zb, const double* zc, const double* zd, double* integrals);

// hrr.c

void hrr(int la, int lb,
         const double* posa, const double* posb,
         double* integrals1, double* integrals2);

void hrr1(int la0, int la1, int lb,
          const double* posa, const double* posb,
          double* integrals1, double* integrals2);

// fmgamma.c

double fm(double T, int m);
double fmtaylor(double T, int m);
double fmasymptotic(double T, int m);
void fmrecursive(double T, int n, double* array);
void calcfmtable();

// keiprim.c

void keiprim(int la, int lb, const double* posa, const double* posb,
             double za, double zb, double* integrals,
             double* stable, double* ttable);

// naiprim.c

void naiprim(int la, int lb, const double* posa, const double* posb, const double* posc,
             double za, double zb, double charge, double* integrals,
             double* gtable);

// oviprim.c

void oviprim(int la, int lb, const double* posa, const double* posb,
             double za, double zb, double* integrals,
             double* stable);

// rys.c

void rysquad(double T, int n, double* rt, double* wt);
void quadrature(int n, double* a, double* b, double mu0, double* rt, double* wt);

// blasx.c

void zero(size_t n, double* a, size_t inca);
void dxypz(size_t n, double alpha, const double* a, size_t inca,
                                   const double* b, size_t incb,
                     double  beta,       double* c, size_t incc);

#ifdef __cplusplus
}
#endif

#endif
