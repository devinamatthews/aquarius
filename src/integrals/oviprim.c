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

/*
 * Calculate OVIs with the algorithm of Ishida
 *  K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 */
void oviprim(int la, int lb, const double* posa, const double* posb, double za, double zb,
             double* integrals, double* stable_)
{
    const double PI_32 = 5.5683279968317078452848179821188;

    double posp[3];
    double zp;
    double afac, bfac;
    double A0;
    int a, b, i, xyz;
    int ax, ay, az, bx, by, bz;
    double (*stable)[la+1][lb+1] = (double(*)[la+1][lb+1])stable_;

    zp = za + zb;
    A0 = PI_32 * exp(-za * zb * dist2(posa, posb) / zp) / pow(zp, 1.5);

    for (xyz = 0;xyz < 3;xyz++)
    {
        posp[xyz] = (posa[xyz] * za + posb[xyz] * zb) / zp;

        afac = posp[xyz] - posa[xyz];
        bfac = posp[xyz] - posb[xyz];

        stable[xyz][0][0] = 1.0;

        for (a = 0;a <= la;a++)
        {
            if (a < la)
            {
                stable[xyz][a + 1][0] = afac * stable[xyz][a][0];
                if (a > 0) stable[xyz][a + 1][0] += a * stable[xyz][a - 1][0] / (2 * zp);
            }

            for (b = 0;b < lb;b++)
            {
                stable[xyz][a][b + 1] = bfac * stable[xyz][a][b];
                if (a > 0) stable[xyz][a][b + 1] += a * stable[xyz][a - 1][b] / (2 * zp);
                if (b > 0) stable[xyz][a][b + 1] += b * stable[xyz][a][b - 1] / (2 * zp);
            }
        }
    }

    for (a = 0;a <= la;a++)
    {
        for (b = 0;b <= lb;b++)
        {
            stable[0][a][b] *= A0;
        }
    }

    i = 0;
    for (bx = 0;bx <= lb;bx++)
    {
        for (by = 0;by <= lb - bx;by++)
        {
            bz = lb - bx - by;
            for (ax = 0;ax <= la;ax++)
            {
                for (ay = 0;ay <= la - ax;ay++)
                {
                    az = la - ax - ay;
                    integrals[i] = stable[0][ax][bx] * stable[1][ay][by] * stable[2][az][bz];
                    i++;
                }
            }
        }
    }
}
