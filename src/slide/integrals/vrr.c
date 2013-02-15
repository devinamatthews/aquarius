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

static void filltable(const int np, double* table, const int le, const int lf,
                      const double*  efac, const double*  ffac,
                      const double*  pfac, const double*  qfac,
                      const double* s1fac, const double* t1fac,
                      const double* s2fac, const double* t2fac, const double* gfac,
                      const int vinc, const int einc, const int finc);

int vrr(const int le0, const int le1, const int lf0, const int lf1,
        const double* posa, const double* posb, const double* posc, const double* posd,
        const int na, const int nb, const int nc, const int nd,
        const double* za, const double* zb, const double* zc, const double* zd, double* restrict integrals)
{
    const double PI_52 = 17.493418327624862846262821679872; // \pi^(5/2)
    const int vmax = le1+lf1;

    //size_t aosize1 = (((le1+1)*(le1+2)*(le1+3)/6)-(le0*(le0+1)*(le0+2)/6))*
    //                 (((lf1+1)*(lf1+2)*(lf1+3)/6)-(lf0*(lf0+1)*(lf0+2)/6));
    size_t aosize3 = na*nb*nc*nd;

    /*
     * xtable is dimensioned as xtable[lf1+1][le1+1][le1+lf1+1][prim]
     * (e.g. xtable[k][i][m][...] = (i0|k0)[m][...])
     *
     * after xtable is initially filled, xtable[lf][le][0][...] is always
     * a valid integral for some xe+ye+ze=le, xf+yf+zf=lf, etc.
     */
    int vinc = ALIGN(aosize3);
    int einc = vinc*(vmax+1);
    int finc = einc*(le1+1);

	size_t worksize = finc*(lf1+1);
	double* xtable = MALLOC(double, worksize);
    if (xtable == NULL) return -1;

    double *efac, *ffac, *pfac, *qfac, *s1fac, *s2fac, *t1fac, *t2fac, *gfac;
     efac = MALLOC(double, 3*aosize3); if ( efac == NULL) return -1;
     ffac = MALLOC(double, 3*aosize3); if ( ffac == NULL) return -1;
     pfac = MALLOC(double, 3*aosize3); if ( pfac == NULL) return -1;
     qfac = MALLOC(double, 3*aosize3); if ( qfac == NULL) return -1;
    s1fac = MALLOC(double,   aosize3); if (s1fac == NULL) return -1;
    s2fac = MALLOC(double,   aosize3); if (s2fac == NULL) return -1;
    t1fac = MALLOC(double,   aosize3); if (t1fac == NULL) return -1;
    t2fac = MALLOC(double,   aosize3); if (t2fac == NULL) return -1;
     gfac = MALLOC(double,   aosize3); if ( gfac == NULL) return -1;

    int i = 0;
    for (int h = 0;h < nd;h++)
    {
        for (int g = 0;g < nc;g++)
        {
            for (int f = 0;f < nb;f++)
            {
                for (int e = 0;e < na;e++)
                {
                    double zp = za[e]+zb[f];
                    double zq = zc[g]+zd[h];

                    double posp[3], posq[3], posw[3];
                    posp[0] = (posa[0]*za[e] + posb[0]*zb[f])/zp;
                    posp[1] = (posa[1]*za[e] + posb[1]*zb[f])/zp;
                    posp[2] = (posa[2]*za[e] + posb[2]*zb[f])/zp;
                    posq[0] = (posc[0]*zc[g] + posd[0]*zd[h])/zq;
                    posq[1] = (posc[1]*zc[g] + posd[1]*zd[h])/zq;
                    posq[2] = (posc[2]*zc[g] + posd[2]*zd[h])/zq;
                    posw[0] = (posp[0]*zp + posq[0]*zq)/(zp+zq);
                    posw[1] = (posp[1]*zp + posq[1]*zq)/(zp+zq);
                    posw[2] = (posp[2]*zp + posq[2]*zq)/(zp+zq);

                    efac[i+0*aosize3] = posp[0] - posa[0];
                    efac[i+1*aosize3] = posp[1] - posa[1];
                    efac[i+2*aosize3] = posp[2] - posa[2];
                    ffac[i+0*aosize3] = posq[0] - posc[0];
                    ffac[i+1*aosize3] = posq[1] - posc[1];
                    ffac[i+2*aosize3] = posq[2] - posc[2];
                    pfac[i+0*aosize3] = posw[0] - posp[0];
                    pfac[i+1*aosize3] = posw[1] - posp[1];
                    pfac[i+2*aosize3] = posw[2] - posp[2];
                    qfac[i+0*aosize3] = posw[0] - posq[0];
                    qfac[i+1*aosize3] = posw[1] - posq[1];
                    qfac[i+2*aosize3] = posw[2] - posq[2];

                    s1fac[i] = 1.0/(2*zp);
                    s2fac[i] = 1.0/(2*zq);
                     gfac[i] = 1.0/(2*(zp+zq));
                    t1fac[i] = -gfac[i]*zq/zp;
                    t2fac[i] = -gfac[i]*zp/zq;

                    double A0 = 2*PI_52*exp(-za[e]*zb[f]*dist2(posa, posb)/zp
                    		                -zc[g]*zd[h]*dist2(posc, posd)/zq)/(sqrt(zp+zq)*zp*zq);
                    double Z = dist2(posp,posq)*zp*zq/(zp+zq);

                    /*
                     * Get (00|00)[m] and scale by A0
                     */
                    double ssssm[vmax+1];
                    fmrecursive(Z, vmax, ssssm);
                    for (int v = 0;v <= vmax;v++) xtable[i+v*vinc] = A0*ssssm[v];

                    i++;
                }
            }
		}
    }

    /*
     * integrals is dimensioned as integrals[N(ld)][N(lc)][N(lb)][N(la)]
     * where N(l) is the number of ways to partition l units of angular momentum
     * amongst x,y,z (N(l) = (l+1)(l+2)/2)
     */
    int iinc = aosize3;
    int kinc = iinc*(((le1+1)*(le1+2)*(le1+3)/6)-(le0*(le0+1)*(le0+2)/6));
    int nint = kinc*(((lf1+1)*(lf1+2)*(lf1+3)/6)-(lf0*(lf0+1)*(lf0+2)/6));

    double* integral[le1+1-le0][lf1+1-lf0];
    double* target[le1+1-le0][lf1+1-lf0];

    integral[0][0] = integrals+nint-1;
    target[0][0] = xtable+le1*einc+lf1*finc;
    for (int k = 0;k < lf1+1-lf0;k++)
    {
    	int f = lf1-k;

    	if (k > 0)
    	{
    		integral[0][k] = integral[0][k-1]-kinc*(f+2)*(f+3)/2;
    		target[0][k] = target[0][k-1]-finc;
    	}

    	for (int i = 1;i < le1+1-le0;i++)
    	{
    		int e = le1-i;
    		integral[i][k] = integral[i-1][k]-iinc*(e+2)*(e+3)/2;
    		target[i][k] = target[i-1][k]-einc;
    	}
    }

    // fill table with x
    filltable(aosize3, xtable, le1, lf1,
              &efac[0*aosize3], &ffac[0*aosize3],
              &pfac[0*aosize3], &qfac[0*aosize3],
              s1fac, t1fac, s2fac, t2fac, gfac,
              vinc, einc, finc);

    // loop over all possible distributions of x momenta
    double* table1 = target[0][0];
    for (int fx = lf1;fx >= 0;fx--)
    {
        for (int ex = le1;ex >= 0;ex--)
        {
			// and fill remainder with y from that point
			filltable(aosize3, table1, le1-ex, lf1-fx,
					  &efac[1*aosize3], &ffac[1*aosize3],
					  &pfac[1*aosize3], &qfac[1*aosize3],
					  s1fac, t1fac, s2fac, t2fac, gfac,
					  vinc, einc, finc);

			// loop over all possible distirubtions of y momenta given x
			double* table2 = target[0][0];
			for (int fy = lf1-fx;fy >= 0;fy--)
			{
				for (int ey = le1-ex;ey >= 0;ey--)
				{
					// and fill remainder with z from that point
					filltable(aosize3, table2, le1-ex-ey, lf1-fx-fy,
							  &efac[2*aosize3], &ffac[2*aosize3],
							  &pfac[2*aosize3], &qfac[2*aosize3],
							  s1fac, t1fac, s2fac, t2fac, gfac,
							  vinc, einc, finc);

					for (int k = 0;k <= lf1-MAX(fx+fy,lf0);k++)
					{
						for (int i = 0;i <= le1-MAX(ex+ey,le0);i++)
						{
							dcopy(aosize3, target[i][k], 1, integral[i][k], 1);
							integral[i][k] -= iinc;
						}
					}

					table2 -= einc;
				}

                for (int k = 0;k <= lf1-MAX(fx+fy,lf0);k++)
                {
                    for (int i = 0;i <= le1-MAX(ex,le0);i++)
                    {
                        int e = le1-i;
                        integral[i][k] -= kinc-iinc*(e-ex+1);
                    }
                }

				table2 -= finc-einc*(le1-ex+1);
			}

            for (int k = 0;k <= lf1-MAX(fx,lf0);k++)
            {
                int f = lf1-k;
                for (int i = 0;i <= le1-MAX(ex,le0);i++)
                {
                    int e = le1-i;
                    integral[i][k] -= iinc*(e-ex+1)-kinc*(f-fx+1);
                }
            }

			table1 -= einc;
        }

        for (int k = 0;k <= lf1-MAX(fx,lf0);k++)
        {
            int f = lf1-k;
            for (int i = 0;i <= le1-le0;i++)
            {
                int e = le1-i;
                integral[i][k] -= kinc*(f-fx+1)-iinc*(e+1)*(e+2)/2;
            }
        }

        table1 -= finc-einc*(le1+1);
    }

    FREE(xtable);
    FREE( efac);
    FREE( ffac);
    FREE( pfac);
    FREE( qfac);
    FREE(s1fac);
    FREE(s2fac);
    FREE(t1fac);
    FREE(t2fac);
    FREE( gfac);

    return 0;
}

static void filltable(const int np, double* restrict table, const int le, const int lf,
					  const double* restrict  efac, const double* restrict  ffac,
					  const double* restrict  pfac, const double* restrict  qfac,
					  const double* restrict s1fac, const double* restrict t1fac,
					  const double* restrict s2fac, const double* restrict t2fac, const double* restrict gfac,
					  const int vinc, const int einc, const int finc)
{
    const int vmax = le+lf;

    int e, f;

	f = 0;
	if (f < lf)
	{
		for (int v = 0;v < vmax-f;v++)
		{
			ALIGNED_LOOP(int p = 0;p < np;p++)
			{
				table[p+finc] = ffac[p]*table[p] +
								qfac[p]*table[p+vinc];
			}
			table += vinc;
		}
		table -= vinc*(vmax-f);
	}

	e = 0;
	if (e < le)
	{
		for (int v = 0;v < vmax-e-f;v++)
		{
			ALIGNED_LOOP(int p = 0;p < np;p++)
			{
				table[p+einc] = efac[p]*table[p] +
								pfac[p]*table[p+vinc];
			}
			table += vinc;
		}
		table -= vinc*(vmax-e-f);
	}
	table += einc;

	for (e = 1;e <= le;e++)
	{
		if (e < le)
		{
			for (int v = 0;v < vmax-e-f;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+einc] = efac[p]*table[p] +
									pfac[p]*table[p+vinc] +
					                e*(s1fac[p]*table[p-einc] +
									   t1fac[p]*table[p-einc+vinc]);
				}
				table += vinc;
			}
			table -= vinc*(vmax-e-f);
		}
		table += einc;
	}
	table += finc-einc*(le+1);

	for (f = 1;f <= lf;f++)
	{
		if (f < lf)
		{
			for (int v = 0;v < vmax-f;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+finc] = ffac[p]*table[p] +
									qfac[p]*table[p+vinc] +
						            f*(s2fac[p]*table[p-finc] +
									   t2fac[p]*table[p-finc+vinc]);
				}
				table += vinc;
			}
			table -= vinc*(vmax-f);
		}

		e = 0;
		if (e < le)
		{
			for (int v = 0;v < vmax-e-f;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+einc] = efac[p]*table[p] +
									pfac[p]*table[p+vinc] +
						            f*gfac[p]*table[p-finc];
				}
				table += vinc;
			}
			table -= vinc*(vmax-e-f);
		}
		table += einc;

		for (e = 1;e <= le;e++)
		{
			if (e < le)
			{
				for (int v = 0;v < vmax-e-f;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+einc] = efac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            f*gfac[p]*table[p-finc] +
							            e*(s1fac[p]*table[p-einc] +
										   t1fac[p]*table[p-einc+vinc]);
					}
					table += vinc;
				}
				table -= vinc*(vmax-e-f);
			}
			table += einc;
		}
		table += finc-einc*(le+1);
	}
}
