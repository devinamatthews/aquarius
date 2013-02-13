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

#ifndef _AQUARIUS_SLIDE_INTERNAL_H_
#define _AQUARIUS_SLIDE_INTERNAL_H_

#include "slide.h"
#include "memory/memory.h"
#include "util/util.h"
#include "util/blas.h"
#include "util/lapack.h"
#include "integrals.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#define MODULO(a,b) ((a)%(b) < 0 ? (a)%(b)+(b) : (a)%(b))

#define MAKEIDX(irrep,idx) (uint16_t)((((irrep)&0x7)<<13)|((idx)&0x1FFF))

#ifdef __cplusplus
extern "C"
{
#endif

// symdata.c

extern const int dirprd[8][8];
extern const int pairstabs[8][8];
extern const int allorders[8];
extern const int allops[8][8];
extern const int allstaborders[8][8];
extern const char* allirreplabels[8][8];
extern const int allchars[8][8][8];
extern const int allstabs[8][8][8];
extern const int alldcrindexs[8][8][8];
extern const int alldcrdegs[8][8][8];

// elementdata.c

extern const element_t elements[];

// symmetry.c

/**
 * contains parity, labeling, and transformation information for cartesian and spherical harmonic
 * functions of a given total angular momentum
 */
typedef struct
{
    int nspher, ncart;
    int* spherparity[8];
    int* cartparity[8];
    double* cart2spher;
    char** cartlabels;
    char** spherlabels;
} angular_data_t;

extern const char* angmomlabels[9];
extern const char* oplabels[8];
extern const char* grouplabels[8];
extern char* irreplabels[8];
extern int group, order;
extern int ops[8], staborder[8];
extern int chars[8][8], stabs[8][8], dcrindex[8][8], dcrdeg[8][8];
extern angular_data_t* angmom[SLIDE_MAX_L + 1];
extern double geomtol;
extern int ordering;
extern int cartesian_order[9][55];
extern int spherical_order[9][55];
extern int cartesian_order_inverse[9][55];
extern int spherical_order_inverse[9][55];

angular_data_t* newangulardata(const int L);
void freeangulardata(angular_data_t* ang);

/**
 * Binomial coefficient \f$\left(\begin{array}{c}a\\b\end{array}\right) = \frac{a!}{b!(a-b)!},\;0 \le b \le a\f$.
 */
int binom(const int a, const int b);

/**
 * \f$n!\f$
 */
int fact(const int n);

/**
 * \f$n!!\f$
 */
int dfact(const int n);

/**
 * Normalized cartesian-to-spherical transformation matrix elements.
 *
 * See e.g. Schlegel, H.B.; Frisch, M.J., Int. J. Quant. Chem. 54, 83 (1995);
 * the final expression is:
 *
 * \f[
 * \phi_s(n,l,m) = \sum\limits_{l_x+l_y+l_z=n} C(l,m,l_x,l_y,l_z) \phi_c(l_x,l_y,l_z)
 * \f]
 *
 * \f[
 * C(l,m,l_x,l_y,l_z) = \left(\frac{\left(2-\delta_{m0}\right) \left(\begin{array}{c}2 l_x\\l_x\end{array}\right)
 *                                  \left(\begin{array}{c}2 l_y\\l_y\end{array}\right) \left(\begin{array}{c}2 l_z\\l_z\end{array}\right)
 *                                  \left(\begin{array}{c}l+|m|\\|m|\end{array}\right)}
 *                                 {\left(\begin{array}{c}2 l\\l\end{array}\right) \left(\begin{array}{c}l\\|m|\end{array}\right)
 *                                  \left(\begin{array}{c}l\\l_x\end{array}\right) \left(\begin{array}{c}l-l_x\\l_y\end{array}\right)}\right)^{1/2}
 *
 *                      \begin{array}{c}\cos\\or\\\sin\end{array} \left(\frac{\left(|m|-l_x\right) \pi}{2}\right)
 *
 *                      \sum\limits_{i=0}^{(l-|m|)/2} \sum\limits_{k=0}^{j} \frac{(-1)^{i+k}}{2^l}
 *
 *                      \left(\begin{array}{c}2l-2i\\l+|m|\end{array}\right) \left(\begin{array}{c}l\\i\end{array}\right)
 *                      \left(\begin{array}{c}i\\j\end{array}\right) \left(\begin{array}{c}j\\k\end{array}\right)
 *                      \left(\begin{array}{c}|m|\\l_x-2k\end{array}\right)
 * \f]
 *
 * where \f$\left(\begin{array}{c}a\\b\end{array}\right)\f$ is a binomial coefficient, \f$\delta_{ab}\f$ is the Kronecker delta,
 * \f$j = (l_x+l_y-|m|)/2\f$ integer, and cos is taken for \f$m >= 0\f$, sin for \f$m < 0\f$.
 */
double cartcoef(const int l, int m, int lx, int ly, int lz);

/**
 * Apply one of the 8 symmetry operations to all of the points pos.
 * Only operations of D2h are considered, so operations only flip the sign of x, y, and/or z.
 */
void applysymop(const int n, double (*pos)[3], const int op);

/**
 * Return true if for each point in pos, the operation op maps this point to some other point (or itself) in pos.
 */
bool checksymop(const int n, const double (*pos)[3], const int op);

/**
 * Find the stabilizer for the given center.
 */
int getstab(const double pos[3]);

// memory.c

typedef struct
{
	int n;
	double** w;
} workspace_t;

int reserve_workspaces(workspace_t* w, size_t size);
double* active_workspace(workspace_t w);
void release_workspaces(workspace_t w);

// transform.c

void transpose(const int m, const int n, const double alpha, const double* A, const int lda,
		                                 const double  beta,       double *B, const int ldb);

/**
 * transform a block of integrals (integrals of all angular momenta types, all primitives) to the symmetry-adapted basis
 * integrals must be ordered e.g. [abcdijkl] with the fast-running abcd indices over contraction and ijkl over angular momentum type
 *
 * this is just the icky part of the transformation which accumulates onto the sointegrals buffer using characters and parity coefficients
 * the driver code must loop over the DCRs
 *
 * the sointegrals buffer is larger than aointegrals sinces of multiple irreps can be generated from one ao
 */
void ao2so2(const shell_t* a, const shell_t* b, const int r, double* aointegrals, double* sointegrals);
void ao2so4(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int r, const int t, const int st, double* aointegrals, double* sointegrals);
void cart2spher2r(const shell_t* a, const shell_t* b, const int nprim, double* buf1, double* buf2);
void cart2spher2l(const shell_t* a, const shell_t* b, const int nprim, double* buf1, double* buf2);
void cart2spher4r(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int nprim, double* buf1, double* buf2);
void cart2spher4l(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int nprim, double* buf1, double* buf2);
void prim2contr2r(const shell_t* a, const shell_t* b, const int nspher, double* buf1, double* buf2);
void prim2contr2l(const shell_t* a, const shell_t* b, const int nspher, double* buf1, double* buf2);
void prim2contr4r(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int nspher, double* buf1, double* buf2);
void prim2contr4l(const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d, const int nspher, double* buf1, double* buf2);

// blasx.c

void zero(const size_t n, double* a, const size_t inca);

void dxypz(const size_t n, const double alpha, const double* a, const size_t inca,
		                                       const double* b, const size_t incb,
		                   const double  beta,       double* c, const size_t incc);

#ifdef __cplusplus
}
#endif

#endif
