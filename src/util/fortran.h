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

#ifndef _AQUARIUS_UTIL_FORTRAN_H_
#define _AQUARIUS_UTIL_FORTRAN_H_

#define LOWER_NO_UNDERSCORE 1
#define LOWER_UNDERSCORE 2
#define UPPER_NO_UNDERSCORE 3
#define UPPER_UNDERSCORE 4

#if F77_NAME == UPPER_NO_UNDERSCORE
#define FFUNC(name,NAME) NAME
#elif F77_NAME == UPPER_UNDERSCORE
#define FFUNC(name,NAME) NAME##_
#elif F77_NAME == LOWER_NO_UNDERSCORE
#define FFUNC(name,NAME) name
#elif F77_NAME == LOWER_UNDERSCORE
#define FFUNC(name,NAME) name##_
#else
#error "Unsupported Fortran naming convention"
#endif

#if FORTRAN_INTEGER_SIZE == 4
typedef int32_t integer;
#elif FORTRAN_INTEGER_SIZE == 8
typedef int64_t integer;
#else
#error "Unsupported Fortran integer size"
#endif

typedef integer logical;
typedef float real;

#if defined(__cplusplus)
#include <complex>
typedef std::complex<float> scomplex;
typedef std::complex<double> dcomplex;
#elif __STDC_VERSION__ == 199901L
#include <complex.h>
typedef float complex scomplex;
typedef double complex dcomplex;
#else
typedef struct scomplex_ { float real,imag; } scomplex;
typedef struct dcomplex_ { double real,imag; } dcomplex;
#endif

#endif
