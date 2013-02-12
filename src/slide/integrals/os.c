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

static void filltable(double* table, const int la, const int lb, const int lc, const int ld,
        const double afac, const double bfac, const double cfac, const double dfac, const double pfac, const double qfac,
        const double s1fac, const double t1fac, const double s2fac, const double t2fac, const double gfac,
        const int ainc, const int binc, const int cinc, const int dinc);

int oseri(const int la, const int lb, const int lc, const int ld,
          const double* posa, const double* posb, const double* posc, const double* posd,
          const int na, const int nb, const int nc, const int nd,
          const double* za, const double* zb, const double* zc, const double* zd, double* restrict integrals)
{
	size_t worksize = (la+1)*(lb+1)*(lc+1)*(ld+1)*(la+lb+lc+ld+1);
	workspace_t work;
    if (reserve_workspaces(&work, worksize) == -1) return -1;

    size_t aosize1 = ((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)*
    		         ((lc+1)*(lc+2)/2)*((ld+1)*(ld+2)/2);

	#ifdef ENABLE_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif //ENABLE_OPENMP
	for (int m = 0;m < na*nb*nc*nd;m++)
	{
		int h = m / (na*nb*nc);
		int n = m - h*na*nb*nc;
		int g = n / (na*nb);
		int o = n - g*na*nb;
		int f = o / na;
		int e = o - f*na;

		size_t idx = aosize1*m;

		osprim(la, lb, lc, ld,
			   posa, posb, posc, posd,
			   za[e], zb[f], zc[g], zd[h],
			   &integrals[idx], active_workspace(work));
	}

	release_workspaces(work);

	return 0;
}

/**
 * Calculate ERIs with the recursive algorithm of Obara and Saika
 *  S. Obara; A. Saika, J. Chem. Phys. 84, 3963 (1986)
 */
void osprim(const int la, const int lb, const int lc, const int ld,
            const double* posa, const double* posb, const double* posc, const double* posd,
            const double za, const double zb, const double zc, const double zd, double* restrict integrals,
            double* restrict xtable)
{
    const double TWO_PI_52 = 34.98683665524972497; // 2*pi^(5/2)
    const int vmax = la+lb+lc+ld;

    double zp = za+zb;
    double zq = zc+zd;

    double A0 = TWO_PI_52*exp(-za*zb*dist2(posa, posb)/zp
    		                  -zc*zd*dist2(posc, posd)/zq)/(zp*zq*sqrt(zp+zq));

    /*
    if (fabs(A0) < 1e-14)
    {
    	dzero(((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)*
    		  ((lc+1)*(lc+2)/2)*((ld+1)*(ld+2)/2), integrals, 1);
    	return;
    }
    */

    double posp[3], posq[3], posw[3];
    posp[0] = (posa[0]*za + posb[0]*zb)/zp;
    posp[1] = (posa[1]*za + posb[1]*zb)/zp;
    posp[2] = (posa[2]*za + posb[2]*zb)/zp;
    posq[0] = (posc[0]*zc + posd[0]*zd)/zq;
    posq[1] = (posc[1]*zc + posd[1]*zd)/zq;
    posq[2] = (posc[2]*zc + posd[2]*zd)/zq;
    posw[0] = (posp[0]*zp + posq[0]*zq)/(zp+zq);
    posw[1] = (posp[1]*zp + posq[1]*zq)/(zp+zq);
    posw[2] = (posp[2]*zp + posq[2]*zq)/(zp+zq);

    double afac[3], bfac[3], cfac[3], dfac[3], pfac[3], qfac[3];
    afac[0] = posp[0] - posa[0];
    afac[1] = posp[1] - posa[1];
    afac[2] = posp[2] - posa[2];
    bfac[0] = posp[0] - posb[0];
    bfac[1] = posp[1] - posb[1];
    bfac[2] = posp[2] - posb[2];
    cfac[0] = posq[0] - posc[0];
    cfac[1] = posq[1] - posc[1];
    cfac[2] = posq[2] - posc[2];
    dfac[0] = posq[0] - posd[0];
    dfac[1] = posq[1] - posd[1];
    dfac[2] = posq[2] - posd[2];
    pfac[0] = posw[0] - posp[0];
    pfac[1] = posw[1] - posp[1];
    pfac[2] = posw[2] - posp[2];
    qfac[0] = posw[0] - posq[0];
    qfac[1] = posw[1] - posq[1];
    qfac[2] = posw[2] - posq[2];

    double s1fac = 1.0/(2*zp);
    double s2fac = 1.0/(2*zq);
    double gfac = 1.0/(2*(zp+zq));
    double t1fac = -gfac*zq/zp;
    double t2fac = -gfac*zp/zq;

    double Z = dist2(posp,posq)*zp*zq/(zp+zq);

    /*
     * Get (00|00)[m] and scale by A0
     */
    //fmrecursive(Z, vmax, xtable);
    for (int v = 0;v <= vmax;v++)
    {
        //xtable[v] *= A0;
        xtable[v] = fm(Z, v)*A0;
    }

    /*
     * xtable is dimensioned as xtable[ld+1][lc+1][lb+1][la+1][la+lb+lc+ld+1]
     * (e.g. xtable[l][k][j][i][m] ~ (ij|kl)[m])
     *
     * after xtable is initially filled, xtable[ld][lc][lb][la][0] is always
     * a valid integral for some xa+ya+za=la, xb+yb+zb=lb, etc.
     */
    int ainc = (vmax+1);
    int binc = ainc*(la+1);
    int cinc = binc*(lb+1);
    int dinc = cinc*(lc+1);

    double* target = &xtable[la*ainc+lb*binc+lc*cinc+ld*dinc];

    /*
     * integrals is dimensioned as integrals[N(ld)][N(lc)][N(lb)][N(la)]
     * where N(l) is the number of ways to partition l units of angular momentum
     * amongst x,y,z (N(l) = (l+1)(l+2)/2)
     */
    int iinc = 1;
    int jinc = iinc*(la+1)*(la+2)/2;
    int kinc = jinc*(lb+1)*(lb+2)/2;
    int linc = kinc*(lc+1)*(lc+2)/2;
    int nint = linc*(ld+1)*(ld+2)/2;

    // fill table with x
    filltable(xtable, la, lb, lc, ld,
              afac[0], bfac[0], cfac[0], dfac[0], pfac[0], qfac[0],
              s1fac, t1fac, s2fac, t2fac, gfac,
              ainc, binc, cinc, dinc);

    // loop over all possible distributions of x momenta
    double* table1 = target;
    double* integral = &integrals[nint-1];
    for (int dx = ld;dx >= 0;dx--)
    {
        for (int cx = lc;cx >= 0;cx--)
        {
            for (int bx = lb;bx >= 0;bx--)
            {
                for (int ax = la;ax >= 0;ax--)
                {
                    // and fill remainder with y from that point
                    filltable(table1, la-ax, lb-bx, lc-cx, ld-dx,
                              afac[1], bfac[1], cfac[1], dfac[1], pfac[1], qfac[1],
                              s1fac, t1fac, s2fac, t2fac, gfac,
                              ainc, binc, cinc, dinc);

                    // loop over all possible distirubtions of y momenta given x
                    double* table2 = target;
                    for (int dy = ld-dx;dy >= 0;dy--)
                    {
                        for (int cy = lc-cx;cy >= 0;cy--)
                        {
                            for (int by = lb-bx;by >= 0;by--)
                            {
                                for (int ay = la-ax;ay >= 0;ay--)
                                {
                                    // and fill remainder with z from that point
                                    filltable(table2, la-ax-ay, lb-bx-by, lc-cx-cy, ld-dx-dy,
                                              afac[2], bfac[2], cfac[2], dfac[2], pfac[2], qfac[2],
                                              s1fac, t1fac, s2fac, t2fac, gfac,
                                              ainc, binc, cinc, dinc);

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
                        table2 -= dinc-cinc*(lc-cx+1);
                        integral -= linc-kinc*(lc-cx+1);
                    }
                    table1 -= ainc;
                    integral -= iinc*(la-ax+1)-linc*(ld-dx+1);
                }
                table1 -= binc-ainc*(la+1);
                integral -= jinc*(lb-bx+1)-iinc*(la+1)*(la+2)/2;
            }
            table1 -= cinc-binc*(lb+1);
            integral -= kinc*(lc-cx+1)-jinc*(lb+1)*(lb+2)/2;
        }
        table1 -= dinc-cinc*(lc+1);
        integral -= linc*(ld-dx+1)-kinc*(lc+1)*(lc+2)/2;
    }
}

static void filltable(double* table, const int la, const int lb, const int lc, const int ld,
        const double afac, const double bfac, const double cfac, const double dfac, const double pfac, const double qfac,
        const double s1fac, const double t1fac, const double s2fac, const double t2fac, const double gfac,
        const int ainc, const int binc, const int cinc, const int dinc)
{
    const int vmax = la+lb+lc+ld;

    for (int d = 0;d <= ld;d++)
    {
        if (d < ld)
        {
            for (int v = 0;v < vmax-d;v++)
            {
                table[v+dinc] = dfac*(table[v  ]) +
                                qfac*(table[v+1]);

                if (d > 0)
                {
                    table[v+dinc] += d*(s2fac*(table[v-dinc  ]) +
                                        t2fac*(table[v-dinc+1]));
                }
            }
        }

        for (int c = 0;c <= lc;c++)
        {
            if (c < lc)
            {
                for (int v = 0;v < vmax-c-d;v++)
                {
                    table[v+cinc] = cfac*(table[v  ]) +
                                    qfac*(table[v+1]);

                    if (c > 0)
                    {
                        table[v+cinc] += c*(s2fac*(table[v-cinc  ]) +
                                            t2fac*(table[v-cinc+1]));
                    }

                    if (d > 0)
                    {
                        table[v+cinc] += d*(s2fac*(table[v-dinc  ]) +
                                            t2fac*(table[v-dinc+1]));
                    }
                }
            }

            for (int b = 0;b <= lb;b++)
            {
                if (b < lb)
                {
                    for (int v = 0;v < vmax-b-c-d;v++)
                    {
                        table[v+binc] = bfac*(table[v  ]) +
                                        pfac*(table[v+1]);

                        if (c > 0)
                        {
                            table[v+binc] += c*gfac*(table[v-cinc+1]);
                        }

                        if (d > 0)
                        {
                            table[v+binc] += d*gfac*(table[v-dinc+1]);
                        }

                        if (b > 0)
                        {
                            table[v+binc] += b*(s1fac*(table[v-binc  ]) +
                                                t1fac*(table[v-binc+1]));
                        }
                    }
                }

                for (int a = 0;a <= la;a++)
                {
                    if (a < la)
                    {
                        for (int v = 0;v < vmax-a-b-c-d;v++)
                        {
                            table[v+ainc] = afac*(table[v  ]) +
                                            pfac*(table[v+1]);

                            if (c > 0)
                            {
                                table[v+ainc] += c*gfac*(table[v-cinc+1]);
                            }

                            if (d > 0)
                            {
                                table[v+ainc] += d*gfac*(table[v-dinc+1]);
                            }

                            if (a > 0)
                            {
                                table[v+ainc] += a*(s1fac*(table[v-ainc  ]) +
                                                    t1fac*(table[v-ainc+1]));
                            }

                            if (b > 0)
                            {
                                table[v+ainc] += b*(s1fac*(table[v-binc  ]) +
                                                    t1fac*(table[v-binc+1]));
                            }
                        }
                    }

                    table += ainc;
                }
                table += binc-ainc*(la+1);
            }
            table += cinc-binc*(lb+1);
        }
        table += dinc-cinc*(lc+1);
    }
}
