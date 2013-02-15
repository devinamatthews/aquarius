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

center_t* SLIDE_new_center(const double pos[3], const element_t* element)
{
    int i, j;
    center_t* c = MALLOC(center_t, 1);
    int dcr;
    double xyz[3];

    if (c == NULL) return NULL;

    c->stabilizer = getstab(pos);
    c->degeneracy = order/staborder[c->stabilizer];
    c->element = (element_t*)element;

    dcr = dcrindex[c->stabilizer][c->stabilizer];

    for (i = 0;i < c->degeneracy;i++)
    {
        memcpy(c->centers[i], pos, sizeof(double)*3);
        applysymop(1, &c->centers[i], stabs[dcr][i]);
    }

    for (j = 0;j < order;j++)
    {
        for (i = 0;i < c->degeneracy;i++)
        {
            memcpy(xyz, pos, sizeof(double)*3);
            applysymop(1, (double(*)[3])xyz, ops[j]);
            if (dist(xyz, c->centers[i]) < geomtol)
            {
                c->centermap[j] = i;
            }
        }
    }

    return c;
}

void SLIDE_free_center(center_t* center)
{
    FREE(center);
}
