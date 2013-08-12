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
 * Calculate ERIs with the Rys Polynomial algorithm of Ishida
 *  Ishida, K. J. Chem. Phys. 95, 5198-205 (1991)
 *  Ishida, K. J. Chem. Phys. 98, 2176 (1993)
 */
void ishprim(int la, int lb, int lc, int ld,
             const double* posa, const double* posb, const double* posc, const double* posd,
             double za, double zb, double zc, double zd, double* restrict integrals,
             double* restrict xtable_)
{
    const double PI_52 = 17.493418327624862846262821679872;
    int nrys = (la + lb + lc + ld) / 2+1;

    double rts[nrys], wts[nrys];
    double posp[3], posq[3];
    double zp, zq;
    double pfac, qfac, afac, bfac, cfac, dfac;
    double aafac, bbfac, ccfac, ddfac, s1fac, s2fac, gfac;
    double A0, Z;
    int a, b, c, d, i, v, xyz;
    int ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
    double (*xtable)[lb+1][lc+1][ld+1][(la+lb+lc+ld)/2+1][3] = (double(*)[lb+1][lc+1][ld+1][(la+lb+lc+ld)/2+1][3])xtable_;

    zp = za + zb;
    zq = zc + zd;

    posp[0] = (posa[0] * za + posb[0] * zb) / zp;
    posp[1] = (posa[1] * za + posb[1] * zb) / zp;
    posp[2] = (posa[2] * za + posb[2] * zb) / zp;
    posq[0] = (posc[0] * zc + posd[0] * zd) / zq;
    posq[1] = (posc[1] * zc + posd[1] * zd) / zq;
    posq[2] = (posc[2] * zc + posd[2] * zd) / zq;

    A0 = 2 * PI_52 * exp(-za * zb * dist2(posa, posb) / zp - zc * zd * dist2(posc, posd) / zq) / (sqrt(zp + zq) * zp * zq);
    Z = dist2(posp, posq) * zp * zq / (zp + zq);
    rysquad(Z, nrys, rts, wts);

    for (xyz = 0;xyz < 3;xyz++)
    {
        pfac = zq * (posq[xyz] - posp[xyz]) / (zp + zq);
        qfac = zp * (posp[xyz] - posq[xyz]) / (zp + zq);

        afac = posp[xyz] - posa[xyz];
        bfac = posp[xyz] - posb[xyz];
        cfac = posq[xyz] - posc[xyz];
        dfac = posq[xyz] - posd[xyz];

        for (v = 0;v < nrys;v++)
        {
            aafac = afac + pfac * rts[v];
            bbfac = bfac + pfac * rts[v];
            ccfac = cfac + qfac * rts[v];
            ddfac = dfac + qfac * rts[v];

            /*
             * there is a typo in Ishida (JCP v98), the definition of G after Eq. 5 should read G = \xi s_i^2
             */
            gfac = rts[v] / (2 * (zp + zq));
            s1fac = (0.5 - zq * gfac) / zp;
            s2fac = (0.5 - zp * gfac) / zq;

            xtable[0][0][0][0][v][xyz] = 1.0;

            for (a = 0;a <= la;a++)
            {
                if (a < la)
                {
                               xtable[0][0][0][a+1][v][xyz] =      aafac * xtable[0][0][0][a]  [v][xyz];
                    if (a > 0) xtable[0][0][0][a+1][v][xyz] += a * s1fac * xtable[0][0][0][a-1][v][xyz];
                }

                for (b = 0;b <= lb;b++)
                {
                    if (b < lb)
                    {
                                   xtable[0][0][b+1][a][v][xyz] =      bbfac * xtable[0][0][b][a]  [v][xyz];
                        if (a > 0) xtable[0][0][b+1][a][v][xyz] += a * s1fac * xtable[0][0][b][a-1][v][xyz];
                        if (b > 0) xtable[0][0][b+1][a][v][xyz] += b * s1fac * xtable[0][0][b-1][a][v][xyz];
                    }

                    for (c = 0;c <= lc;c++)
                    {
                        if (c < lc)
                        {
                                       xtable[0][c+1][b][a][v][xyz] =      ccfac * xtable[0][c][b][a]  [v][xyz];
                            if (a > 0) xtable[0][c+1][b][a][v][xyz] += a *  gfac * xtable[0][c][b][a-1][v][xyz];
                            if (b > 0) xtable[0][c+1][b][a][v][xyz] += b *  gfac * xtable[0][c][b-1][a][v][xyz];
                            if (c > 0) xtable[0][c+1][b][a][v][xyz] += c * s2fac * xtable[0][c-1][b][a][v][xyz];
                        }

                        for (d = 0;d < ld;d++)
                        {
                                       xtable[d+1][c][b][a][v][xyz] =      ddfac * xtable[d][c][b][a]  [v][xyz];
                            if (a > 0) xtable[d+1][c][b][a][v][xyz] += a *  gfac * xtable[d][c][b][a-1][v][xyz];
                            if (b > 0) xtable[d+1][c][b][a][v][xyz] += b *  gfac * xtable[d][c][b-1][a][v][xyz];
                            if (c > 0) xtable[d+1][c][b][a][v][xyz] += c * s2fac * xtable[d][c-1][b][a][v][xyz];
                            if (d > 0) xtable[d+1][c][b][a][v][xyz] += d * s2fac * xtable[d-1][c][b][a][v][xyz];
                        }
                    }
                }
            }
        }
    }

    for (d = 0;d <= ld;d++)
    {
        for (c = 0;c <= lc;c++)
        {
            for (b = 0;b <= lb;b++)
            {
                for (a = 0;a <= la;a++)
                {
                    for (v = 0;v < nrys;v++)
                    {
                        xtable[d][c][b][a][v][0] *= wts[v] * A0;
                    }
                }
            }
        }
    }

    i = 0;
    for (dx = 0;dx <= ld;dx++)
    {
        for (dy = 0;dy <= ld - dx;dy++)
        {
            dz = ld - dx - dy;
            for (cx = 0;cx <= lc;cx++)
            {
                for (cy = 0;cy <= lc - cx;cy++)
                {
                    cz = lc - cx - cy;
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
                                    integrals[i] = 0.0;
                                    for (v = 0;v < nrys;v++)
                                    {
                                        integrals[i] += xtable[dx][cx][bx][ax][v][0] *
                                                        xtable[dy][cy][by][ay][v][1] *
                                                        xtable[dz][cz][bz][az][v][2];
                                    }
                                    i++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
