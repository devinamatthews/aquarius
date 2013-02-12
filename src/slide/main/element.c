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

const element_t* SLIDE_get_element(const char* symbol)
{
    char rs[3] = "  ";
    int i;

    if (symbol[1] == ' ' || symbol[1] == '\0')
    {
        rs[1] = (char)toupper(symbol[0]);
    }
    else if (symbol[0] == ' ')
    {
        rs[1] = (char)toupper(symbol[1]);
    }
    else
    {
        rs[0] = (char)toupper(symbol[0]);
        rs[1] = (char)tolower(symbol[1]);
    }

    for (i = 0;i < 102;i++)
    {
        if (strncmp(elements[i].symbol, rs, 2) == 0) return &elements[i];
    }

    return NULL;
}

element_t* SLIDE_copy_element(const element_t* element)
{
    element_t* e;

    e = MALLOC(element_t, 1);
    if (e == NULL) return NULL;

    memcpy(e, element, sizeof(element_t));

    return e;
}

void SLIDE_free_element(element_t* element)
{
	FREE(element);
}
