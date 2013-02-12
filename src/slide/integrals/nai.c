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

void fillgtable(double* table, const int la, const int lb,
                const double afac, const double bfac, const double cfac,
                const double sfac, const int ainc, const int binc);

/*
 * Calculate NAIs with the Rys Polynomial algorithm of Ishida
 *  K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
 *  Ishida, K., J. Chem. Phys., 98, 2176 (1993)
 */
void naiprim(const int la, const int lb, const double* posa, const double* posb, const double* posc,
             const double za, const double zb, const int charge, double* integrals,
             double* gtable)
{
    const double PI = 3.1415926535897932384626433832795;
    const int vmax = la + lb;

    double posp[3];
    double zp;
    double afac[3], bfac[3], cfac[3];
    double sfac;
    double A0, Z;
    int i, nint;
    int ax, ay, bx, by;
    int ainc, binc, iinc, jinc;
    double* table1;
    double* table2;
    double* target;
    double* integral;

    ainc = (vmax + 1);
    binc = ainc * (la + 1);

    iinc = 1;
    jinc = iinc * (la + 1) * (la + 2) / 2;
    nint = jinc * (lb + 1) * (lb + 2) / 2;

    zp = za + zb;

    posp[0] = (posa[0] * za + posb[0] * zb) / zp;
    posp[1] = (posa[1] * za + posb[1] * zb) / zp;
    posp[2] = (posa[2] * za + posb[2] * zb) / zp;

    afac[0] = posp[0] - posa[0];
    afac[1] = posp[1] - posa[1];
    afac[2] = posp[2] - posa[2];
    bfac[0] = posp[0] - posb[0];
    bfac[1] = posp[1] - posb[1];
    bfac[2] = posp[2] - posb[2];
    cfac[0] = posp[0] - posc[0];
    cfac[1] = posp[1] - posc[1];
    cfac[2] = posp[2] - posc[2];
    sfac = 0.5 / zp;

    A0 = -charge * 2 * PI * exp(-za * zb * dist2(posa, posb) / zp) / zp;
    Z = dist2(posp, posc) * zp;

    //fmrecursive(Z, vmax, (double*)gtable);
    for (i = 0;i <= vmax;i++)
        gtable[i] = fm(Z, i) * A0;

    target = &gtable[la * ainc + lb * binc];

    // fill table with x
    fillgtable(&gtable[0], la, lb, afac[0], bfac[0], cfac[0], sfac, ainc, binc);

    // loop over all possible distributions of x momenta
    table1 = target;
    integral = &integrals[nint - 1];
    for (bx = lb;bx >= 0;bx--)
    {
        for (ax = la;ax >= 0;ax--)
        {
            // and fill remainder with y from that point
            fillgtable(table1, la - ax, lb - bx, afac[1], bfac[1], cfac[1], sfac, ainc, binc);

            // loop over all possible distributions of y momenta given x
            table2 = target;
            for (by = lb - bx;by >= 0;by--)
            {
                for (ay = la - ax;ay >= 0;ay--)
                {
                    // and fill remainder with z from that point
                    fillgtable(table2, la - ax - ay, lb - bx - by, afac[2], bfac[2], cfac[2], sfac, ainc, binc);

                    *integral += *target;

                    table2 -= ainc;
                    integral -= iinc;
                }
                table2 -= binc - ainc * (la - ax + 1);
                integral -= jinc - iinc * (la - ax + 1);
            }
            table1 -= ainc;
            integral -= iinc * (la - ax + 1) - jinc * (lb - bx + 1);
        }
        table1 -= binc - ainc * (la + 1);
        integral -= jinc * (lb - bx + 1) - iinc * (la + 1) * (la + 2) / 2;
    }
}

void fillgtable(double* table, const int la, const int lb, const double afac, const double bfac, const double cfac, const double sfac, const int ainc, const int binc)
{
    int a, b;
    int v;
    int vmax = la + lb;

    for (b = 0;b <= lb;b++)
    {
        if (b < lb)
        {
            for (v = 0;v < vmax - b;v++)
            {
                *(table + v + binc) = bfac * (*(table + v)) - cfac * (*(table + v + 1));

                if (b > 0)
                {
                    *(table + v + binc) += b * sfac * (*(table + v - binc) - *(table + v - binc + 1));
                }
            }
        }

        for (a = 0;a <= la;a++)
        {
            if (a < la)
            {
                for (v = 0;v < vmax - a - b;v++)
                {
                    *(table + v + ainc) = afac * (*(table + v)) - cfac * (*(table + v + 1));

                    if (a > 0)
                    {
                        *(table + v + ainc) += a * sfac * (*(table + v - ainc) - *(table + v - ainc + 1));
                    }

                    if (b > 0)
                    {
                        *(table + v + ainc) += b * sfac * (*(table + v - binc) - *(table + v - binc + 1));
                    }
                }
            }
            table += ainc;
        }
        table += binc - ainc * (la + 1);
    }
}
