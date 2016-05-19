#include "internal.h"

static void filltable(int np, double* table, int la, int lb, int lc, int ld,
					  const double* afac, const double* bfac,
					  const double* cfac, const double* dfac,
					  const double* pfac, const double* qfac,
					  const double* s1fac, const double* t1fac,
					  const double* s2fac, const double* t2fac, const double* gfac,
					  int vinc, int ainc, int binc, int cinc, int dinc);

void osinv(int la, int lb, int lc, int ld,
           const double* posa, const double* posb, const double* posc, const double* posd,
           int na, int nb, int nc, int nd,
           const double* za, const double* zb, const double* zc, const double* zd, double* restrict integrals)
{
    const double PI_52 = 17.493418327624862846262821679872; // \pi^(5/2)
    int vmax = la+lb+lc+ld;

    size_t aosize1 = ((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)*
    		         ((lc+1)*(lc+2)/2)*((ld+1)*(ld+2)/2);
    size_t aosize3 = na*nb*nc*nd;

    /*
     * xtable is dimensioned as xtable[ld+1][lc+1][lb+1][la+1][la+lb+lc+ld+1][prim]
     * (e.g. xtable[l][k][j][i][m] ~ (ij|kl)[m])
     *
     * after xtable is initially filled, xtable[ld][lc][lb][la][0] is always
     * a valid integral for some xa+ya+za=la, xb+yb+zb=lb, etc.
     */
    int vinc = ALIGN(aosize3);
    int ainc = vinc*(vmax+1);
    int binc = ainc*(la+1);
    int cinc = binc*(lb+1);
    int dinc = cinc*(lc+1);

	size_t worksize = dinc*(ld+1);
	double* xtable = SAFE_MALLOC(double, worksize);

    double *afac, *bfac, *cfac, *dfac, *pfac, *qfac, *s1fac, *s2fac, *t1fac, *t2fac, *gfac;
     afac = SAFE_MALLOC(double, 3*aosize3);
     bfac = SAFE_MALLOC(double, 3*aosize3);
     cfac = SAFE_MALLOC(double, 3*aosize3);
     dfac = SAFE_MALLOC(double, 3*aosize3);
     pfac = SAFE_MALLOC(double, 3*aosize3);
     qfac = SAFE_MALLOC(double, 3*aosize3);
    s1fac = SAFE_MALLOC(double,   aosize3);
    s2fac = SAFE_MALLOC(double,   aosize3);
    t1fac = SAFE_MALLOC(double,   aosize3);
    t2fac = SAFE_MALLOC(double,   aosize3);
     gfac = SAFE_MALLOC(double,   aosize3);

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

                    afac[i+0*aosize3] = posp[0] - posa[0];
                    afac[i+1*aosize3] = posp[1] - posa[1];
                    afac[i+2*aosize3] = posp[2] - posa[2];
                    bfac[i+0*aosize3] = posp[0] - posb[0];
                    bfac[i+1*aosize3] = posp[1] - posb[1];
                    bfac[i+2*aosize3] = posp[2] - posb[2];
                    cfac[i+0*aosize3] = posq[0] - posc[0];
                    cfac[i+1*aosize3] = posq[1] - posc[1];
                    cfac[i+2*aosize3] = posq[2] - posc[2];
                    dfac[i+0*aosize3] = posq[0] - posd[0];
                    dfac[i+1*aosize3] = posq[1] - posd[1];
                    dfac[i+2*aosize3] = posq[2] - posd[2];
                    pfac[i+0*aosize3] = posw[0] - posp[0];
                    pfac[i+1*aosize3] = posw[1] - posp[1];
                    pfac[i+2*aosize3] = posw[2] - posp[2];
                    qfac[i+0*aosize3] = posw[0] - posq[0];
                    qfac[i+1*aosize3] = posw[1] - posq[1];
                    qfac[i+2*aosize3] = posw[2] - posq[2];

                    s1fac[i] = 0.5/zp;
                    s2fac[i] = 0.5/zq;
                     gfac[i] = 0.5/(zp+zq);
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
    PROFILE_FLOPS((int64_t)na*nb*nc*nd*(vmax+56+17*DIV_FLOPS+SQRT_FLOPS+EXP_FLOPS));

    // fill table with x
    filltable(aosize3, xtable, la, lb, lc, ld,
              &afac[0*aosize3], &bfac[0*aosize3],
              &cfac[0*aosize3], &dfac[0*aosize3],
              &pfac[0*aosize3], &qfac[0*aosize3],
              s1fac, t1fac, s2fac, t2fac, gfac,
              vinc, ainc, binc, cinc, dinc);

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
                    filltable(aosize3, table1, la-ax, lb-bx, lc-cx, ld-dx,
                              &afac[1*aosize3], &bfac[1*aosize3],
                              &cfac[1*aosize3], &dfac[1*aosize3],
                              &pfac[1*aosize3], &qfac[1*aosize3],
                              s1fac, t1fac, s2fac, t2fac, gfac,
                              vinc, ainc, binc, cinc, dinc);

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
                                    filltable(aosize3, table2, la-ax-ay, lb-bx-by, lc-cx-cy, ld-dx-dy,
                                              &afac[2*aosize3], &bfac[2*aosize3],
                                              &cfac[2*aosize3], &dfac[2*aosize3],
                                              &pfac[2*aosize3], &qfac[2*aosize3],
                                              s1fac, t1fac, s2fac, t2fac, gfac,
                                              vinc, ainc, binc, cinc, dinc);

                                    dcopy(aosize3, target, 1, integral, aosize1);

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

    FREE(xtable);
    FREE( afac);
    FREE( bfac);
    FREE( cfac);
    FREE( dfac);
    FREE( pfac);
    FREE( qfac);
    FREE(s1fac);
    FREE(s2fac);
    FREE(t1fac);
    FREE(t2fac);
    FREE( gfac);
}

static void filltable(int np, double* restrict table, int la, int lb, int lc, int ld,
					  const double* restrict  afac, const double* restrict  bfac,
					  const double* restrict  cfac, const double* restrict  dfac,
					  const double* restrict  pfac, const double* restrict  qfac,
					  const double* restrict s1fac, const double* restrict t1fac,
					  const double* restrict s2fac, const double* restrict t2fac, const double* restrict gfac,
					  int vinc, int ainc, int binc, int cinc, int dinc)
{
    int vmax = la+lb+lc+ld;

    int a, b, c, d;

    d = 0;
	if (d < ld)
	{
	    PROFILE_FLOPS(3*np*(vmax-d));
		for (int v = 0;v < vmax-d;v++)
		{
			ALIGNED_LOOP(int p = 0;p < np;p++)
			{
				table[p+dinc] = dfac[p]*table[p] +
								qfac[p]*table[p+vinc];
			}
			table += vinc;
		}
		table -= vinc*(vmax-d);
	}

	c = 0;
	if (c < lc)
	{
        PROFILE_FLOPS(3*np*(vmax-c-d));
		for (int v = 0;v < vmax-c-d;v++)
		{
			ALIGNED_LOOP(int p = 0;p < np;p++)
			{
				table[p+cinc] = cfac[p]*table[p] +
								qfac[p]*table[p+vinc];
			}
			table += vinc;
		}
		table -= vinc*(vmax-c-d);
	}

	b = 0;
	if (b < lb)
	{
        PROFILE_FLOPS(3*np*(vmax-b-c-d));
		for (int v = 0;v < vmax-b-c-d;v++)
		{
			ALIGNED_LOOP(int p = 0;p < np;p++)
			{
				table[p+binc] = bfac[p]*table[p] +
								pfac[p]*table[p+vinc];
			}
			table += vinc;
		}
		table -= vinc*(vmax-b-c-d);
	}

	a = 0;
	if (a < la)
	{
        PROFILE_FLOPS(3*np*(vmax-a-b-c-d));
		for (int v = 0;v < vmax-a-b-c-d;v++)
		{
			ALIGNED_LOOP(int p = 0;p < np;p++)
			{
				table[p+ainc] = afac[p]*table[p] +
								pfac[p]*table[p+vinc];
			}
			table += vinc;
		}
		table -= vinc*(vmax-a-b-c-d);
	}
	table += ainc;

	for (a = 1;a <= la;a++)
	{
		if (a < la)
		{
	        PROFILE_FLOPS(8*np*(vmax-a-b-c-d));
			for (int v = 0;v < vmax-a-b-c-d;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+ainc] = afac[p]*table[p] +
									pfac[p]*table[p+vinc] +
					                a*(s1fac[p]*table[p-ainc] +
									   t1fac[p]*table[p-ainc+vinc]);
				}
				table += vinc;
			}
			table -= vinc*(vmax-a-b-c-d);
		}
		table += ainc;
	}
	table += binc-ainc*(la+1);

	for (b = 1;b <= lb;b++)
	{
		if (b < lb)
		{
	        PROFILE_FLOPS(8*np*(vmax-b-c-d));
			for (int v = 0;v < vmax-b-c-d;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+binc] = bfac[p]*table[p] +
									pfac[p]*table[p+vinc] +
						            b*(s1fac[p]*table[p-binc] +
									   t1fac[p]*table[p-binc+vinc]);
				}
				table += vinc;
			}
			table -= vinc*(vmax-b-c-d);
		}

		a = 0;
		if (a < la)
		{
	        PROFILE_FLOPS(8*np*(vmax-a-b-c-d));
			for (int v = 0;v < vmax-a-b-c-d;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+ainc] = afac[p]*table[p] +
									pfac[p]*table[p+vinc] +
						            b*(s1fac[p]*table[p-binc] +
									   t1fac[p]*table[p-binc+vinc]);
				}
				table += vinc;
			}
			table -= vinc*(vmax-a-b-c-d);
		}
		table += ainc;

		for (a = 1;a <= la;a++)
		{
			if (a < la)
			{
		        PROFILE_FLOPS(13*np*(vmax-a-b-c-d));
				for (int v = 0;v < vmax-a-b-c-d;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+ainc] = afac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            a*(s1fac[p]*table[p-ainc] +
										   t1fac[p]*table[p-ainc+vinc]) +
						                b*(s1fac[p]*table[p-binc] +
									       t1fac[p]*table[p-binc+vinc]);
					}
					table += vinc;
				}
				table -= vinc*(vmax-a-b-c-d);
			}
			table += ainc;
		}
		table += binc-ainc*(la+1);
	}
	table += cinc-binc*(lb+1);

	for (c = 1;c <= lc;c++)
	{
		if (c < lc)
		{
	        PROFILE_FLOPS(8*np*(vmax-c-d));
			for (int v = 0;v < vmax-c-d;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+cinc] = cfac[p]*table[p] +
									qfac[p]*table[p+vinc] +
						            c*(s2fac[p]*table[p-cinc] +
									   t2fac[p]*table[p-cinc+vinc]);
				}
				table += vinc;
			}
			table -= vinc*(vmax-c-d);
		}

		b = 0;
		if (b < lb)
		{
	        PROFILE_FLOPS(6*np*(vmax-b-c-d));
			for (int v = 0;v < vmax-b-c-d;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+binc] = bfac[p]*table[p] +
									pfac[p]*table[p+vinc] +
						            c*gfac[p]*table[p-cinc];
				}
				table += vinc;
			}
			table -= vinc*(vmax-b-c-d);
		}

		a = 0;
		if (a < la)
		{
	        PROFILE_FLOPS(6*np*(vmax-a-b-c-d));
			for (int v = 0;v < vmax-a-b-c-d;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+ainc] = afac[p]*table[p] +
									pfac[p]*table[p+vinc] +
						            c*gfac[p]*table[p-cinc];
				}
				table += vinc;
			}
			table -= vinc*(vmax-a-b-c-d);
		}
		table += ainc;

		for (a = 1;a <= la;a++)
		{
			if (a < la)
			{
		        PROFILE_FLOPS(11*np*(vmax-a-b-c-d));
				for (int v = 0;v < vmax-a-b-c-d;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+ainc] = afac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            c*gfac[p]*table[p-cinc] +
							            a*(s1fac[p]*table[p-ainc] +
										   t1fac[p]*table[p-ainc+vinc]);
					}
					table += vinc;
				}
				table -= vinc*(vmax-a-b-c-d);
			}
			table += ainc;
		}
		table += binc-ainc*(la+1);

		for (b = 1;b <= lb;b++)
		{
			if (b < lb)
			{
		        PROFILE_FLOPS(11*np*(vmax-b-c-d));
				for (int v = 0;v < vmax-b-c-d;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+binc] = bfac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            c*gfac[p]*table[p-cinc] +
						                b*(s1fac[p]*table[p-binc] +
									       t1fac[p]*table[p-binc+vinc]);
					}
					table += vinc;
				}
				table -= vinc*(vmax-b-c-d);
			}

			a = 0;
			if (a < la)
			{
		        PROFILE_FLOPS(11*np*(vmax-a-b-c-d));
				for (int v = 0;v < vmax-a-b-c-d;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+ainc] = afac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            c*gfac[p]*table[p-cinc] +
							            b*(s1fac[p]*table[p-binc] +
										   t1fac[p]*table[p-binc+vinc]);
					}
					table += vinc;
				}
				table -= vinc*(vmax-a-b-c-d);
			}
			table += ainc;

			for (a = 1;a <= la;a++)
			{
				if (a < la)
				{
			        PROFILE_FLOPS(16*np*(vmax-a-b-c-d));
					for (int v = 0;v < vmax-a-b-c-d;v++)
					{
						ALIGNED_LOOP(int p = 0;p < np;p++)
						{
							table[p+ainc] = afac[p]*table[p] +
											pfac[p]*table[p+vinc] +
								            c*gfac[p]*table[p-cinc] +
								            a*(s1fac[p]*table[p-ainc] +
											   t1fac[p]*table[p-ainc+vinc]) +
								            b*(s1fac[p]*table[p-binc] +
											   t1fac[p]*table[p-binc+vinc]);
						}
						table += vinc;
					}
					table -= vinc*(vmax-a-b-c-d);
				}
				table += ainc;
			}
			table += binc-ainc*(la+1);
		}
		table += cinc-binc*(lb+1);
	}
	table += dinc-cinc*(lc+1);

    for (d = 1;d <= ld;d++)
    {
        if (d < ld)
        {
            PROFILE_FLOPS(8*np*(vmax-d));
            for (int v = 0;v < vmax-d;v++)
            {
				ALIGNED_LOOP(int p = 0;p < np;p++)
            	{
            		table[p+dinc] = dfac[p]*table[p] +
            		                qfac[p]*table[p+vinc] +
	            		            d*(s2fac[p]*table[p-dinc] +
	            				       t2fac[p]*table[p-dinc+vinc]);
            	}
                table += vinc;
            }
            table -= vinc*(vmax-d);
        }

        c = 0;
		if (c < lc)
		{
	        PROFILE_FLOPS(8*np*(vmax-c-d));
			for (int v = 0;v < vmax-c-d;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+cinc] = cfac[p]*table[p] +
									qfac[p]*table[p+vinc] +
						            d*(s2fac[p]*table[p-dinc] +
									   t2fac[p]*table[p-dinc+vinc]);
				}
				table += vinc;
			}
			table -= vinc*(vmax-c-d);
		}

		b = 0;
		if (b < lb)
		{
	        PROFILE_FLOPS(6*np*(vmax-b-c-d));
			for (int v = 0;v < vmax-b-c-d;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+binc] = bfac[p]*table[p] +
									pfac[p]*table[p+vinc] +
						            d*gfac[p]*table[p-dinc];
				}
				table += vinc;
			}
			table -= vinc*(vmax-b-c-d);
		}

		a = 0;
		if (a < la)
		{
	        PROFILE_FLOPS(6*np*(vmax-a-b-c-d));
			for (int v = 0;v < vmax-a-b-c-d;v++)
			{
				ALIGNED_LOOP(int p = 0;p < np;p++)
				{
					table[p+ainc] = afac[p]*table[p] +
									pfac[p]*table[p+vinc] +
						            d*gfac[p]*table[p-dinc];
				}
				table += vinc;
			}
			table -= vinc*(vmax-a-b-c-d);
		}
		table += ainc;

		for (a = 1;a <= la;a++)
		{
			if (a < la)
			{
		        PROFILE_FLOPS(11*np*(vmax-a-b-c-d));
				for (int v = 0;v < vmax-a-b-c-d;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+ainc] = afac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            d*gfac[p]*table[p-dinc] +
							            a*(s1fac[p]*table[p-ainc] +
										   t1fac[p]*table[p-ainc+vinc]);
					}
					table += vinc;
				}
				table -= vinc*(vmax-a-b-c-d);
			}
			table += ainc;
		}
		table += binc-ainc*(la+1);

		for (b = 1;b <= lb;b++)
		{
			if (b < lb)
			{
		        PROFILE_FLOPS(11*np*(vmax-b-c-d));
				for (int v = 0;v < vmax-b-c-d;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+binc] = bfac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            d*gfac[p]*table[p-dinc] +
							            b*(s1fac[p]*table[p-binc] +
										   t1fac[p]*table[p-binc+vinc]);
					}
					table += vinc;
				}
				table -= vinc*(vmax-b-c-d);
			}

			a = 0;
			if (a < la)
			{
		        PROFILE_FLOPS(11*np*(vmax-a-b-c-d));
				for (int v = 0;v < vmax-a-b-c-d;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+ainc] = afac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            d*gfac[p]*table[p-dinc] +
							            b*(s1fac[p]*table[p-binc] +
										   t1fac[p]*table[p-binc+vinc]);
					}
					table += vinc;
				}
				table -= vinc*(vmax-a-b-c-d);
			}
			table += ainc;

			for (a = 1;a <= la;a++)
			{
				if (a < la)
				{
			        PROFILE_FLOPS(16*np*(vmax-a-b-c-d));
					for (int v = 0;v < vmax-a-b-c-d;v++)
					{
						ALIGNED_LOOP(int p = 0;p < np;p++)
						{
							table[p+ainc] = afac[p]*table[p] +
											pfac[p]*table[p+vinc] +
								            d*gfac[p]*table[p-dinc] +
								            a*(s1fac[p]*table[p-ainc] +
											   t1fac[p]*table[p-ainc+vinc]) +
								            b*(s1fac[p]*table[p-binc] +
											   t1fac[p]*table[p-binc+vinc]);
						}
						table += vinc;
					}
					table -= vinc*(vmax-a-b-c-d);
				}
				table += ainc;
			}
			table += binc-ainc*(la+1);
		}
		table += cinc-binc*(lb+1);

        for (c = 1;c <= lc;c++)
        {
            if (c < lc)
            {
                PROFILE_FLOPS(13*np*(vmax-c-d));
                for (int v = 0;v < vmax-c-d;v++)
                {
					ALIGNED_LOOP(int p = 0;p < np;p++)
                	{
						table[p+cinc] = cfac[p]*table[p] +
										qfac[p]*table[p+vinc] +
		            		            c*(s2fac[p]*table[p-cinc] +
		            				       t2fac[p]*table[p-cinc+vinc]) +
		            		            d*(s2fac[p]*table[p-dinc] +
		            				       t2fac[p]*table[p-dinc+vinc]);
                	}
                    table += vinc;
                }
                table -= vinc*(vmax-c-d);
            }

            b = 0;
			if (b < lb)
			{
		        PROFILE_FLOPS(9*np*(vmax-b-c-d));
				for (int v = 0;v < vmax-b-c-d;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+binc] = bfac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            c*gfac[p]*table[p-cinc] +
							            d*gfac[p]*table[p-dinc];
					}
					table += vinc;
				}
				table -= vinc*(vmax-b-c-d);
			}

			a = 0;
			if (a < la)
			{
		        PROFILE_FLOPS(9*np*(vmax-a-b-c-d));
				for (int v = 0;v < vmax-a-b-c-d;v++)
				{
					ALIGNED_LOOP(int p = 0;p < np;p++)
					{
						table[p+ainc] = afac[p]*table[p] +
										pfac[p]*table[p+vinc] +
							            c*gfac[p]*table[p-cinc] +
							            d*gfac[p]*table[p-dinc];
					}
					table += vinc;
				}
				table -= vinc*(vmax-a-b-c-d);
			}
			table += ainc;

			for (a = 1;a <= la;a++)
			{
				if (a < la)
				{
			        PROFILE_FLOPS(14*np*(vmax-a-b-c-d));
					for (int v = 0;v < vmax-a-b-c-d;v++)
					{
						ALIGNED_LOOP(int p = 0;p < np;p++)
						{
							table[p+ainc] = afac[p]*table[p] +
											pfac[p]*table[p+vinc] +
								            c*gfac[p]*table[p-cinc] +
							                d*gfac[p]*table[p-dinc] +
								            a*(s1fac[p]*table[p-ainc] +
											   t1fac[p]*table[p-ainc+vinc]);
						}
						table += vinc;
					}
					table -= vinc*(vmax-a-b-c-d);
				}
				table += ainc;
			}
			table += binc-ainc*(la+1);

            for (b = 1;b <= lb;b++)
            {
                if (b < lb)
                {
                    PROFILE_FLOPS(14*np*(vmax-b-c-d));
                    for (int v = 0;v < vmax-b-c-d;v++)
                    {
						ALIGNED_LOOP(int p = 0;p < np;p++)
                    	{
							table[p+binc] = bfac[p]*table[p] +
											pfac[p]*table[p+vinc] +
			            		            c*gfac[p]*table[p-cinc] +
			            		            d*gfac[p]*table[p-dinc] +
			            		            b*(s1fac[p]*table[p-binc] +
			            				       t1fac[p]*table[p-binc+vinc]);
                    	}
                        table += vinc;
                    }
                    table -= vinc*(vmax-b-c-d);
                }

                a = 0;
				if (a < la)
				{
			        PROFILE_FLOPS(14*np*(vmax-a-b-c-d));
					for (int v = 0;v < vmax-a-b-c-d;v++)
					{
						ALIGNED_LOOP(int p = 0;p < np;p++)
						{
							table[p+ainc] = afac[p]*table[p] +
											pfac[p]*table[p+vinc] +
								            c*gfac[p]*table[p-cinc] +
								            d*gfac[p]*table[p-dinc] +
								            b*(s1fac[p]*table[p-binc] +
											   t1fac[p]*table[p-binc+vinc]);
						}
						table += vinc;
					}
					table -= vinc*(vmax-a-b-c-d);
				}
				table += ainc;

                for (a = 1;a <= la;a++)
                {
                    if (a < la)
                    {
                        PROFILE_FLOPS(19*np*(vmax-a-b-c-d));
                        for (int v = 0;v < vmax-a-b-c-d;v++)
                        {
							ALIGNED_LOOP(int p = 0;p < np;p++)
                        	{
    							table[p+ainc] = afac[p]*table[p] +
    											pfac[p]*table[p+vinc] +
				            		            c*gfac[p]*table[p-cinc] +
				            		            d*gfac[p]*table[p-dinc] +
									            a*(s1fac[p]*table[p-ainc] +
												   t1fac[p]*table[p-ainc+vinc]) +
									            b*(s1fac[p]*table[p-binc] +
												   t1fac[p]*table[p-binc+vinc]);
                        	}
                            table += vinc;
                        }
                        table -= vinc*(vmax-a-b-c-d);
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
