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

#ifndef _AQUARIUS_SLIDE_INTEGRALS_H_
#define _AQUARIUS_SLIDE_INTEGRALS_H_

// ish.c

void ishprim(const int la, const int lb, const int lc, const int ld,
             const double* posa, const double* posb, const double* posc, const double* posd,
             const double za, const double zb, const double zc, const double zd, double* integrals,
             double* xtable);

// moment.c

void momprim(const int la, const int lb, const double* posa, const double* posb, const double za, const double zb,
             const int lc, const double* posc, double* integrals, double* table);

// os.c

int oseri(const int la, const int lb, const int lc, const int ld,
          const double* posa, const double* posb, const double* posc, const double* posd,
          const int na, const int nb, const int nc, const int nd,
          const double* za, const double* zb, const double* zc, const double* zd, double* integrals);

void osprim(const int la, const int lb, const int lc, const int ld,
            const double* posa, const double* posb, const double* posc, const double* posd,
            const double za, const double zb, const double zc, const double zd, double* integrals,
            double* xtable);

// osinv.c

int osinv(const int la, const int lb, const int lc, const int ld,
          const double* posa, const double* posb, const double* posc, const double* posd,
          const int na, const int nb, const int nc, const int nd,
          const double* za, const double* zb, const double* zc, const double* zd, double* integrals);

// vrr.c

int vrr(const int le0, const int le1, const int lf0, const int lf1,
        const double* posa, const double* posb, const double* posc, const double* posd,
        const int na, const int nb, const int nc, const int nd,
        const double* za, const double* zb, const double* zc, const double* zd, double* integrals);

// hrr.c

void hrr(const int la, const int lb,
         const double* posa, const double* posb,
         double* integrals1, double* integrals2);

void hrr1(const int la0, const int la1, const int lb,
          const double* posa, const double* posb,
          double* integrals1, double* integrals2);

// fmgamma.c

double fm(const double T, const int m);
double fmtaylor(const double T, const int m);
double fmasymptotic(const double T, const int m);
void fmrecursive(const double T, const int n, double* array);
void calcfmtable();

// kei.c

void keiprim(const int la, const int lb, const double* posa, const double* posb,
             const double za, const double zb, double* integrals,
             double* stable, double* ttable);

// nai.c

void naiprim(const int la, const int lb, const double* posa, const double* posb, const double* posc,
             const double za, const double zb, const int charge, double* integrals,
             double* gtable);

// ovi.c

void oviprim(const int la, const int lb, const double* posa, const double* posb,
             const double za, const double zb, double* integrals,
             double* stable);

// rys.c

void rysquad(const double T, const int n, double* rt, double* wt);
void quadrature(const int n, double* a, double* b, const double mu0, double* rt, double* wt);

#endif
