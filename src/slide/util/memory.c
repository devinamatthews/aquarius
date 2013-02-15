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

#include "internal.h"

int reserve_workspaces(workspace_t* w, size_t size)
{
	if (w == NULL)
	{
		assert(0);
		return -1;
	}

	#ifdef ENABLE_OPENMP
	w->n = omp_get_max_threads();
	#else //ENABLE_OPENMP
	w->n = 1;
	#endif //ENABLE_OPENMP

	w->w = MALLOC(double*, w->n);
	if (w->w == NULL)
	{
		assert(0);
		return -1;
	}

	for (int i = 0;i < w->n;i++)
	{
		w->w[i] = MALLOC(double, size);
		if (w->w[i] == NULL)
		{
			assert(0);
			for (int j = i-1;j >= 0;j--) FREE(w->w[j]);
			FREE(w->w);
			return -1;
		}
	}

	return 0;
}

double* active_workspace(workspace_t w)
{
	#ifdef ENABLE_OPENMP
	return w.w[omp_get_thread_num()];
	#else //ENABLE_OPENMP
	return w.w[0];
	#endif //ENABLE_OPENMP
}

void release_workspaces(workspace_t w)
{
	for (int i = 0;i < w.n;i++)
	{
		FREE(w.w[i]);
	}
	FREE(w.w);
}
