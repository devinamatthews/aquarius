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

#ifndef _AQUARIUS_MEMORY_H_
#define _AQUARIUS_MEMORY_H_

#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "util/util.h"

#define ALIGNMENT 16

#define ALIGN(x) ((x)+((ALIGNMENT-((x)&(ALIGNMENT-1)))&(ALIGNMENT-1)))

#ifdef __INTEL_COMPILER

#define ALIGNED_LOOP(x) \
_Pragma("ivdep") \
_Pragma("vector aligned") \
for (x)

#else

#define ALIGNED_LOOP(x) \
for (x)

#endif

#define MALLOC(type, size) (type*)aq_malloc(sizeof(type)*(size), __FILE__, __LINE__, 0)

#define REALLOC(type, ptr, size) (type*)aq_realloc(ptr, sizeof(type)*(size), __FILE__, __LINE__, 0)

#define SAFE_MALLOC(type, size) (type*)aq_malloc(sizeof(type)*(size), __FILE__, __LINE__, 1)

#define SAFE_REALLOC(type, ptr, size) (type*)aq_realloc(ptr, sizeof(type)*(size), __FILE__, __LINE__, 1)

#define FREE(ptr) aq_free(ptr, __FILE__, __LINE__)

#ifdef __cplusplus
extern "C"
{
#endif

extern size_t mem_limit, mem_used;

void set_memory_limit(const size_t bytes);

size_t get_memory_limit();

size_t get_memory_used();

void* aq_malloc(const size_t size, const char* who, const int where, const int bailout);

void* aq_realloc(void* ptr, const size_t size, const char* who, const int where, const int bailout);

void aq_free(void* ptr, const char* who, const int where);

#ifdef __cplusplus
}
#endif

#endif
