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

static void filltable(double* table, const int la, const int lb, const int lc,
                      const double afac, const double bfac, const double cfac, const double sfac,
                      const int ainc, const int binc, const int cinc);

/*
 * Calculate moment integrals with the algorithm of Obara and Saika
 *  S. Obara; A. Saika, J. Chem. Phys. 84, 3963 (1986)
 */
void momprim(const int la, const int lb, const double* posa, const double* posb, const double za, const double zb,
             const int lc, const double* posc, double* integrals, double* table)
{
    const double PI_32 = 5.5683279968317078452848179821188;

    double zp = za + zb;
    double A0 = PI_32 * exp(-za * zb * dist2(posa, posb) / zp) / pow(zp, 1.5);

    double posp[3];
    posp[0] = (posa[0]*za + posb[0]*zb)/zp;
    posp[1] = (posa[1]*za + posb[1]*zb)/zp;
    posp[2] = (posa[2]*za + posb[2]*zb)/zp;

    double afac[3], bfac[3], cfac[3];
    afac[0] = posp[0] - posa[0];
    afac[1] = posp[1] - posa[1];
    afac[2] = posp[2] - posa[2];
    bfac[0] = posp[0] - posb[0];
    bfac[1] = posp[1] - posb[1];
    bfac[2] = posp[2] - posb[2];
    cfac[0] = posp[0] - posc[0];
    cfac[1] = posp[1] - posc[1];
    cfac[2] = posp[2] - posc[2];

    double sfac = 1.0/(2*zp);

    table[0] = A0;

    /*
     * table is dimensioned as table[lc+1][lb+1][la+1]
     * (e.g. table[k][j][i] ~ (i|k|j))
     *
     * after table is initially filled, table[lc][lb][la] is always
     * a valid integral for some xa+ya+za=la, xb+yb+zb=lb, etc.
     */
    int ainc = 1;
    int binc = ainc*(la+1);
    int cinc = binc*(lb+1);

    double* target = &table[la*ainc+lb*binc+lc*cinc];

    /*
     * integrals is dimensioned as integrals[N(lc)][N(lb)][N(la)]
     * where N(l) is the number of ways to partition l units of angular momentum
     * amongst x,y,z (N(l) = (l+1)(l+2)/2)
     */
    int iinc = 1;
    int jinc = iinc*(la+1)*(la+2)/2;
    int kinc = jinc*(lb+1)*(lb+2)/2;
    int nint = kinc*(lc+1)*(lc+2)/2;

    // fill table with x
    filltable(table, la, lb, lc,
              afac[0], bfac[0], cfac[0], sfac,
              ainc, binc, cinc);

    // loop over all possible distributions of x momenta
    double* table1 = target;
    double* integral = &integrals[nint-1];
    for (int cx = lc;cx >= 0;cx--)
    {
        for (int bx = lb;bx >= 0;bx--)
        {
            for (int ax = la;ax >= 0;ax--)
            {
                // and fill remainder with y from that point
                filltable(table1, la-ax, lb-bx, lc-cx,
                          afac[1], bfac[1], cfac[1], sfac,
                          ainc, binc, cinc);

                // loop over all possible distirubtions of y momenta given x
                double* table2 = target;
                for (int cy = lc-cx;cy >= 0;cy--)
                {
                    for (int by = lb-bx;by >= 0;by--)
                    {
                        for (int ay = la-ax;ay >= 0;ay--)
                        {
                            // and fill remainder with z from that point
                            filltable(table2, la-ax-ay, lb-bx-by, lc-cx-cy,
                                      afac[2], bfac[2], cfac[2], sfac,
                                      ainc, binc, cinc);

                            *integral = *target;

                            table2 -= ainc;
                            integral -= iinc;
                        }
                        table2 -= binc-ainc*(la-ax+1);
                        integral -= jinc-iinc*(la-ax+1);
                    }
                    table2 -= cinc-binc*(lb-bx+1);
                    integral -= kinc-jinc*(lb-bx+1);
                }
                table1 -= ainc;
                integral -= iinc*(la-ax+1)-kinc*(lc-cx+1);
            }
            table1 -= binc-ainc*(la+1);
            integral -= jinc*(lb-bx+1)-iinc*(la+1)*(la+2)/2;
        }
        table1 -= cinc-binc*(lb+1);
        integral -= kinc*(lc-cx+1)-jinc*(lb+1)*(lb+2)/2;
    }
}

static void filltable(double* table, const int la, const int lb, const int lc,
                      const double afac, const double bfac, const double cfac, const double sfac,
                      const int ainc, const int binc, const int cinc)
{
    for (int c = 0;c <= lc;c++)
    {
        if (c < lc)
        {
            table[cinc] = cfac*table[0];

            if (c > 0)
            {
                table[cinc] += c*sfac*table[-cinc];
            }
        }

        for (int b = 0;b <= lb;b++)
        {
            if (b < lb)
            {
                table[binc] = bfac*table[0];

                if (b > 0)
                {
                    table[binc] += b*sfac*table[-binc];
                }

                if (c > 0)
                {
                    table[binc] += c*sfac*table[-cinc];
                }
            }

            for (int a = 0;a <= la;a++)
            {
                if (a < la)
                {
                        table[ainc] = afac*table[0];

                        if (a > 0)
                        {
                            table[ainc] += a*sfac*table[-ainc];
                        }

                        if (b > 0)
                        {
                            table[ainc] += b*sfac*table[-binc];
                        }

                        if (c > 0)
                        {
                            table[ainc] += c*sfac*table[-cinc];
                        }
                }

                table += ainc;
            }
            table += binc-ainc*(la+1);
        }
        table += cinc-binc*(lb+1);
    }
}
