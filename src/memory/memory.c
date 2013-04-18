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

#include "memory.h"

#include <stdint.h>

size_t mem_limit = 0, mem_used = 0;

void set_memory_limit(const size_t bytes)
{
    mem_limit = bytes;
}

size_t get_memory_limit()
{
    return mem_limit;
}

size_t get_memory_used()
{
    return mem_used;
}

void* aq_malloc(const size_t size_, const char *file, const int line, const int bailout)
{
	assert(ALIGNMENT >= sizeof(size_t));

	size_t size = MAX(size_,1);

	if (mem_limit > 0 && mem_used+size+ALIGNMENT > mem_limit)
    {
        if (bailout) ERROR_AT(file, line, "Could not allocate memory");
	    return NULL;
    }

    void* mem;
    if (posix_memalign(&mem, ALIGNMENT, size+ALIGNMENT) != 0)
    {
        if (bailout) ERROR_AT(file, line, "Could not allocate memory");
        return NULL;
    }

    mem_used += size+ALIGNMENT;
    *(size_t*)mem = size;

    return (void*)((intptr_t)mem+ALIGNMENT);
}

void* aq_realloc(void* ptr, const size_t size_, const char *file, const int line, const int bailout)
{
    size_t size = MAX(size_,1);

	if (ptr == NULL) return aq_malloc(size, file, line, bailout);
	if (size == 0)
	{
		aq_free(ptr, file, line);
		return NULL;
	}

	ptr -= ALIGNMENT;
	size_t old_size = *(size_t*)ptr;

	if (mem_limit > 0 && mem_used+size-old_size > mem_limit)
    {
        if (bailout) ERROR_AT(file, line, "Could not allocate memory");
	    return NULL;
    }

    void* mem = aq_malloc(size, file, line, bailout);
    if (mem == NULL)
    {
        if (bailout) ERROR_AT(file, line, "Could not allocate memory");
        return NULL;
    }

    memcpy(mem+ALIGNMENT, ptr+ALIGNMENT, MIN(size,old_size));

    aq_free(ptr, file, line);

    mem_used += size-old_size;
    *(size_t*)mem = size;

    return mem+ALIGNMENT;
}

void aq_free(void* ptr, const char *file, const int line)
{
	if (ptr == NULL) return;

	ptr -= ALIGNMENT;
	size_t old_size = *(size_t*)ptr;

    free(ptr);

    mem_used -= old_size+ALIGNMENT;
}
