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

#include "internal.h"

void zero(const size_t n, double* a, const size_t inca)
{
	if (inca == 1)
	{
		memset(a, 0, sizeof(double)*n);
	}
	else
	{
		size_t ia = 0;
		for (size_t i = 0;i < n;i++)
		{
			a[ia] = 0.0;
			ia += inca;
		}
	}
}

void dxypz(const size_t n, const double alpha, const double* restrict a, const size_t inca,
		                                        const double* restrict b, const size_t incb,
		                    const double  beta,       double* restrict c, const size_t incc)
{
	if (inca == 1 && incb == 1 && incc == 1)
	{
		if (alpha == 1.0)
		{
			if (beta == 0.0)
			{
				for (size_t i = 0;i < n;i++)
				{
					c[i] = a[i]*b[i];
				}
			}
			else
			{
				for (size_t i = 0;i < n;i++)
				{
					c[i] = beta*c[i] + a[i]*b[i];
				}
			}
		}
		else
		{
			if (beta == 0.0)
			{
				for (size_t i = 0;i < n;i++)
				{
					c[i] = alpha*a[i]*b[i];
				}
			}
			else
			{
				for (size_t i = 0;i < n;i++)
				{
					c[i] = beta*c[i] + alpha*a[i]*b[i];
				}
			}
		}
	}
	else
	{
		size_t ia = 0;
		size_t ib = 0;
		size_t ic = 0;
		for (size_t i = 0;i < n;i++)
		{
			c[ic] = beta*c[ic] + alpha*a[ia]*b[ib];
			ia += inca;
			ib += incb;
			ic += incc;
		}
	}
}
