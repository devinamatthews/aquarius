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

#define IDX_EQ(i,r,e,j,s,f) ((i) == (j) && (r) == (s) && (e) == (f))
#define IDX_GE(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >= (f)))))
#define IDX_GT(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >  (f)))))

size_t SLIDE_process_2e_ints(context_t* context, const size_t nprocess, double* integrals, idx4_t* indices, const double cutoff)
{
    shell_t* a = context->a;
    shell_t* b = context->b;
    shell_t* c = context->c;
    shell_t* d = context->d;

    size_t m = 0;
    size_t n = 0;
    for (int l = 0;l < d->nfunc;l++)
    {
        for (int k = 0;k < c->nfunc;k++)
        {
            for (int j = 0;j < b->nfunc;j++)
            {
                for (int i = 0;i < a->nfunc;i++)
                {
                    for (int u = 0;u < d->pos->degeneracy;u++)
                    {
                        for (int t = 0;t < c->pos->degeneracy;t++)
                        {
                            for (int s = 0;s < b->pos->degeneracy;s++)
                            {
                                for (int r = 0;r < a->pos->degeneracy;r++)
                                {
                                    int w = a->irreps[i][r];
                                    int x = b->irreps[j][s];
                                    int y = c->irreps[k][t];
                                    int z = d->irreps[l][u];

                                    if (dirprd[dirprd[w][x]][dirprd[y][z]] != SLIDE_IRREP_TOT_SYM) continue;

                                    for (int h = 0;h < d->ncontr;h++)
                                    {
                                        for (int g = 0;g < c->ncontr;g++)
                                        {
                                            for (int f = 0;f < b->ncontr;f++)
                                            {
                                                for (int e = 0;e < a->ncontr;e++)
                                                {
                                                    if (context->num_processed > m)
                                                    {
                                                        m++;
                                                        continue;
                                                    }

                                                    bool bad = false;

                                                    if (a == b && !IDX_GE(i,r,e,j,s,f)) bad = true;
                                                    if (c == d && !IDX_GE(k,t,g,l,u,h)) bad = true;
                                                    if (a == c && b == d && !(IDX_GT(i,r,e,k,t,g) ||
                                                        (IDX_EQ(i,r,e,k,t,g) && IDX_GE(j,s,f,l,u,h)))) bad = true;

                                                    if (!bad && fabs(context->integrals[m]) > cutoff)
                                                    {
                                                        if (a->spherical)
                                                        {
                                                            switch (ordering)
                                                            {
                                                                case SLIDE_ORDER_ISCF: indices[n].i = SLIDE_INDEX_ISCF(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_ISFC: indices[n].i = SLIDE_INDEX_ISFC(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SICF: indices[n].i = SLIDE_INDEX_SICF(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SIFC: indices[n].i = SLIDE_INDEX_SIFC(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SCIF: indices[n].i = SLIDE_INDEX_SCIF(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SCFI: indices[n].i = SLIDE_INDEX_SCFI(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SFIC: indices[n].i = SLIDE_INDEX_SFIC(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SFCI: indices[n].i = SLIDE_INDEX_SFCI(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                                            }
                                                        }
                                                        else
                                                        {
                                                            switch (ordering)
                                                            {
                                                                case SLIDE_ORDER_ISCF: indices[n].i = SLIDE_INDEX_ISCF(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_ISFC: indices[n].i = SLIDE_INDEX_ISFC(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SICF: indices[n].i = SLIDE_INDEX_SICF(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SIFC: indices[n].i = SLIDE_INDEX_SIFC(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SCIF: indices[n].i = SLIDE_INDEX_SCIF(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SCFI: indices[n].i = SLIDE_INDEX_SCFI(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SFIC: indices[n].i = SLIDE_INDEX_SFIC(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                                                case SLIDE_ORDER_SFCI: indices[n].i = SLIDE_INDEX_SFCI(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                                            }
                                                        }

                                                        if (b->spherical)
                                                        {
                                                            switch (ordering)
                                                            {
                                                                case SLIDE_ORDER_ISCF: indices[n].j = SLIDE_INDEX_ISCF(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_ISFC: indices[n].j = SLIDE_INDEX_ISFC(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SICF: indices[n].j = SLIDE_INDEX_SICF(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SIFC: indices[n].j = SLIDE_INDEX_SIFC(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SCIF: indices[n].j = SLIDE_INDEX_SCIF(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SCFI: indices[n].j = SLIDE_INDEX_SCFI(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SFIC: indices[n].j = SLIDE_INDEX_SFIC(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SFCI: indices[n].j = SLIDE_INDEX_SFCI(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                                            }
                                                        }
                                                        else
                                                        {
                                                            switch (ordering)
                                                            {
                                                                case SLIDE_ORDER_ISCF: indices[n].j = SLIDE_INDEX_ISCF(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_ISFC: indices[n].j = SLIDE_INDEX_ISFC(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SICF: indices[n].j = SLIDE_INDEX_SICF(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SIFC: indices[n].j = SLIDE_INDEX_SIFC(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SCIF: indices[n].j = SLIDE_INDEX_SCIF(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SCFI: indices[n].j = SLIDE_INDEX_SCFI(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SFIC: indices[n].j = SLIDE_INDEX_SFIC(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                                                case SLIDE_ORDER_SFCI: indices[n].j = SLIDE_INDEX_SFCI(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                                            }
                                                        }

                                                        if (c->spherical)
                                                        {
                                                            switch (ordering)
                                                            {
                                                                case SLIDE_ORDER_ISCF: indices[n].k = SLIDE_INDEX_ISCF(c,spherical_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_ISFC: indices[n].k = SLIDE_INDEX_ISFC(c,spherical_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SICF: indices[n].k = SLIDE_INDEX_SICF(c,spherical_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SIFC: indices[n].k = SLIDE_INDEX_SIFC(c,spherical_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SCIF: indices[n].k = SLIDE_INDEX_SCIF(c,spherical_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SCFI: indices[n].k = SLIDE_INDEX_SCFI(c,spherical_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SFIC: indices[n].k = SLIDE_INDEX_SFIC(c,spherical_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SFCI: indices[n].k = SLIDE_INDEX_SFCI(c,spherical_order_inverse[c->L][k],g,t,y); break;
                                                            }
                                                        }
                                                        else
                                                        {
                                                            switch (ordering)
                                                            {
                                                                case SLIDE_ORDER_ISCF: indices[n].k = SLIDE_INDEX_ISCF(c,cartesian_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_ISFC: indices[n].k = SLIDE_INDEX_ISFC(c,cartesian_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SICF: indices[n].k = SLIDE_INDEX_SICF(c,cartesian_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SIFC: indices[n].k = SLIDE_INDEX_SIFC(c,cartesian_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SCIF: indices[n].k = SLIDE_INDEX_SCIF(c,cartesian_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SCFI: indices[n].k = SLIDE_INDEX_SCFI(c,cartesian_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SFIC: indices[n].k = SLIDE_INDEX_SFIC(c,cartesian_order_inverse[c->L][k],g,t,y); break;
                                                                case SLIDE_ORDER_SFCI: indices[n].k = SLIDE_INDEX_SFCI(c,cartesian_order_inverse[c->L][k],g,t,y); break;
                                                            }
                                                        }

                                                        if (d->spherical)
                                                        {
                                                            switch (ordering)
                                                            {
                                                                case SLIDE_ORDER_ISCF: indices[n].l = SLIDE_INDEX_ISCF(d,spherical_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_ISFC: indices[n].l = SLIDE_INDEX_ISFC(d,spherical_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SICF: indices[n].l = SLIDE_INDEX_SICF(d,spherical_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SIFC: indices[n].l = SLIDE_INDEX_SIFC(d,spherical_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SCIF: indices[n].l = SLIDE_INDEX_SCIF(d,spherical_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SCFI: indices[n].l = SLIDE_INDEX_SCFI(d,spherical_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SFIC: indices[n].l = SLIDE_INDEX_SFIC(d,spherical_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SFCI: indices[n].l = SLIDE_INDEX_SFCI(d,spherical_order_inverse[d->L][l],h,u,z); break;
                                                            }
                                                        }
                                                        else
                                                        {
                                                            switch (ordering)
                                                            {
                                                                case SLIDE_ORDER_ISCF: indices[n].l = SLIDE_INDEX_ISCF(d,cartesian_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_ISFC: indices[n].l = SLIDE_INDEX_ISFC(d,cartesian_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SICF: indices[n].l = SLIDE_INDEX_SICF(d,cartesian_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SIFC: indices[n].l = SLIDE_INDEX_SIFC(d,cartesian_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SCIF: indices[n].l = SLIDE_INDEX_SCIF(d,cartesian_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SCFI: indices[n].l = SLIDE_INDEX_SCFI(d,cartesian_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SFIC: indices[n].l = SLIDE_INDEX_SFIC(d,cartesian_order_inverse[d->L][l],h,u,z); break;
                                                                case SLIDE_ORDER_SFCI: indices[n].l = SLIDE_INDEX_SFCI(d,cartesian_order_inverse[d->L][l],h,u,z); break;
                                                            }
                                                        }

                                                        integrals[n++] = context->integrals[m];
                                                    }

                                                    context->num_processed++;
                                                    m++;

                                                    if (n >= nprocess) return n;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return n;
}

size_t SLIDE_process_1e_ints(context_t* context, const size_t nprocess, double* integrals, idx2_t* indices, const double cutoff)
{
    shell_t* a = context->a;
    shell_t* b = context->b;

    size_t m = 0;
    size_t n = 0;
    for (int j = 0;j < b->nfunc;j++)
    {
        for (int i = 0;i < a->nfunc;i++)
        {
            for (int s = 0;s < b->pos->degeneracy;s++)
            {
                for (int r = 0;r < a->pos->degeneracy;r++)
                {
                    int w = a->irreps[i][r];
                    int x = b->irreps[j][s];

                    if (dirprd[w][x] != SLIDE_IRREP_TOT_SYM) continue;

                    for (int f = 0;f < b->ncontr;f++)
                    {
                        for (int e = 0;e < a->ncontr;e++)
                        {
                            if (context->num_processed > m)
                            {
                                m++;
                                continue;
                            }

                            if (fabs(context->integrals[m]) > cutoff && (a != b || IDX_GE(i,r,e,j,s,f)))
                            {
                                if (a->spherical)
                                {
                                    switch (ordering)
                                    {
                                        case SLIDE_ORDER_ISCF: indices[n].i = SLIDE_INDEX_ISCF(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_ISFC: indices[n].i = SLIDE_INDEX_ISFC(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SICF: indices[n].i = SLIDE_INDEX_SICF(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SIFC: indices[n].i = SLIDE_INDEX_SIFC(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SCIF: indices[n].i = SLIDE_INDEX_SCIF(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SCFI: indices[n].i = SLIDE_INDEX_SCFI(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SFIC: indices[n].i = SLIDE_INDEX_SFIC(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SFCI: indices[n].i = SLIDE_INDEX_SFCI(a,spherical_order_inverse[a->L][i],e,r,w); break;
                                    }
                                }
                                else
                                {
                                    switch (ordering)
                                    {
                                        case SLIDE_ORDER_ISCF: indices[n].i = SLIDE_INDEX_ISCF(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_ISFC: indices[n].i = SLIDE_INDEX_ISFC(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SICF: indices[n].i = SLIDE_INDEX_SICF(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SIFC: indices[n].i = SLIDE_INDEX_SIFC(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SCIF: indices[n].i = SLIDE_INDEX_SCIF(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SCFI: indices[n].i = SLIDE_INDEX_SCFI(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SFIC: indices[n].i = SLIDE_INDEX_SFIC(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                        case SLIDE_ORDER_SFCI: indices[n].i = SLIDE_INDEX_SFCI(a,cartesian_order_inverse[a->L][i],e,r,w); break;
                                    }
                                }

                                if (b->spherical)
                                {
                                    switch (ordering)
                                    {
                                        case SLIDE_ORDER_ISCF: indices[n].j = SLIDE_INDEX_ISCF(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_ISFC: indices[n].j = SLIDE_INDEX_ISFC(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SICF: indices[n].j = SLIDE_INDEX_SICF(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SIFC: indices[n].j = SLIDE_INDEX_SIFC(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SCIF: indices[n].j = SLIDE_INDEX_SCIF(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SCFI: indices[n].j = SLIDE_INDEX_SCFI(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SFIC: indices[n].j = SLIDE_INDEX_SFIC(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SFCI: indices[n].j = SLIDE_INDEX_SFCI(b,spherical_order_inverse[b->L][j],f,s,x); break;
                                    }
                                }
                                else
                                {
                                    switch (ordering)
                                    {
                                        case SLIDE_ORDER_ISCF: indices[n].j = SLIDE_INDEX_ISCF(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_ISFC: indices[n].j = SLIDE_INDEX_ISFC(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SICF: indices[n].j = SLIDE_INDEX_SICF(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SIFC: indices[n].j = SLIDE_INDEX_SIFC(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SCIF: indices[n].j = SLIDE_INDEX_SCIF(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SCFI: indices[n].j = SLIDE_INDEX_SCFI(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SFIC: indices[n].j = SLIDE_INDEX_SFIC(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                        case SLIDE_ORDER_SFCI: indices[n].j = SLIDE_INDEX_SFCI(b,cartesian_order_inverse[b->L][j],f,s,x); break;
                                    }
                                }

                                integrals[n++] = context->integrals[m];
                            }

                            context->num_processed++;
                            m++;

                            if (n >= nprocess) return n;
                        }
                    }
                }
            }
        }
    }

    return n;
}

size_t SLIDE_calc_eri(context_t* context, const double alpha, const double beta,
		              const shell_t* a_, const shell_t* b_, const shell_t* c_, const shell_t* d_)
{
	shell_t* a = (shell_t*)a_;
	shell_t* b = (shell_t*)b_;
	shell_t* c = (shell_t*)c_;
	shell_t* d = (shell_t*)d_;

	if (a->L < b->L)
	{
		shell_t* t = a;
		a = b;
		b = t;
	}

	if (c->L < d->L)
	{
		shell_t* t = c;
		c = d;
		d = t;
	}

    int dcrr = dcrindex[a->pos->stabilizer][b->pos->stabilizer];
    int dcrs = dcrindex[c->pos->stabilizer][d->pos->stabilizer];
    int dcrt = dcrindex[pairstabs[a->pos->stabilizer][b->pos->stabilizer]][pairstabs[c->pos->stabilizer][d->pos->stabilizer]];
    int lambdat = dcrdeg[pairstabs[a->pos->stabilizer][b->pos->stabilizer]][pairstabs[c->pos->stabilizer][d->pos->stabilizer]];

    size_t aosize1 = angmom[a->L]->ncart * angmom[b->L]->ncart * angmom[c->L]->ncart * angmom[d->L]->ncart;
    size_t aosize2 = a->nfunc * b->nfunc * c->nfunc * d->nfunc;
    size_t aosize3 = a->nprim * b->nprim * c->nprim * d->nprim;
    size_t aosize4 = a->ncontr * b->ncontr * c->ncontr * d->ncontr;

    size_t sosizemax = aosize2 * aosize4 * a->pos->degeneracy * b->pos->degeneracy * c->pos->degeneracy * d->pos->degeneracy /
                       MIN(MIN(a->pos->degeneracy,b->pos->degeneracy),MIN(c->pos->degeneracy,d->pos->degeneracy));

    if (context->aosize < aosize1*aosize3)
    {
        context->aosize = aosize1*aosize3;
        if (context->aobuf1 != NULL) FREE(context->aobuf1);
        if (context->aobuf1 != NULL) FREE(context->aobuf2);
        context->aobuf1 = MALLOC(double, aosize1*aosize3);
        context->aobuf2 = MALLOC(double, aosize1*aosize3);
        /*
         * It's okay to not free aobuf1 if only aobuf2 cannot be allocated
         */
        if (context->aobuf1 == NULL || context->aobuf2 == NULL) return 0;
    }

    if (context->sosize < sosizemax)
    {
        context->sosize = sosizemax;
        if (context->integrals != NULL) FREE(context->integrals);
        context->integrals = MALLOC(double, sosizemax);
        /*
         * Similarly, aobuf[12] need not be free'd here
         */
        if (context->integrals == NULL) return 0;
    }

    if (beta == 0.0)
    {
        zero(sosizemax, context->integrals, 1);
    }
    else if (beta != 1.0)
    {
        dscal(sosizemax, beta, context->integrals, 1);
    }

    double coef = alpha * (double)order / (double)lambdat;

    for (int i = 0;i < staborder[dcrr];i++)
    {
        int r = stabs[dcrr][i];
        for (int j = 0;j < staborder[dcrs];j++)
        {
            int s = stabs[dcrs][j];
            for (int k = 0;k < staborder[dcrt];k++)
            {
                int t = stabs[dcrt][k];
                int st = dirprd[s][t];

                //int rv = osinv(a->L, b->L, c->L, d->L,
                int rv = oseri(a->L, b->L, c->L, d->L,
                               a->pos->centers[0],
                               b->pos->centers[b->pos->centermap[r]],
                               c->pos->centers[c->pos->centermap[t]],
                               d->pos->centers[d->pos->centermap[st]],
                               a->nprim, b->nprim, c->nprim, d->nprim,
                               a->exponents, b->exponents, c->exponents, d->exponents,
                               context->aobuf2);
                if (rv == -1) return 0;

                prim2contr4r(a, b, c, d, aosize1, context->aobuf2, context->aobuf1);
                cart2spher4r(a, b, c, d, aosize4, context->aobuf1, context->aobuf2);

                transpose(aosize2, aosize4, coef, context->aobuf2, aosize2,
                		                     0.0, context->aobuf1, aosize4);

                ao2so4(a, b, c, d, aosize4, r, t, st, context->aobuf1, context->integrals);
            }
        }
    }

    context->a = (shell_t*)a;
    context->b = (shell_t*)b;
    context->c = (shell_t*)c;
    context->d = (shell_t*)d;
    context->num_integrals = 0;
    context->num_processed = 0;

    for (int l = 0;l < d->nfunc;l++)
    {
        for (int k = 0;k < c->nfunc;k++)
        {
            for (int j = 0;j < b->nfunc;j++)
            {
                for (int i = 0;i < a->nfunc;i++)
                {
                    for (int u = 0;u < d->pos->degeneracy;u++)
                    {
                        int z = d->irreps[l][u];
                        for (int t = 0;t < c->pos->degeneracy;t++)
                        {
                            int y = c->irreps[k][t];
                            int yz = dirprd[y][z];
                            for (int s = 0;s < b->pos->degeneracy;s++)
                            {
                                int x = b->irreps[j][s];
                                int xyz = dirprd[x][yz];
                                for (int r = 0;r < a->pos->degeneracy;r++)
                                {
                                    int w = a->irreps[i][r];
                                    if (dirprd[w][xyz] != SLIDE_IRREP_TOT_SYM) continue;

                                    context->num_integrals += a->ncontr*b->ncontr*c->ncontr*d->ncontr;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return context->num_integrals;
}

size_t SLIDE_calc_ovi(context_t* context, const double alpha, const double beta, const shell_t* a, const shell_t* b)
{
    int i, j, r, s;
    int e, f, m;
    int w, x;
    size_t idx2;
    size_t aosize1, aosize2, aosize3, aosize4, sosizemax;
    int dcrr, lambdar;
    double coef;
    workspace_t work;

    dcrr = dcrindex[a->pos->stabilizer][b->pos->stabilizer];
    lambdar = dcrdeg[a->pos->stabilizer][b->pos->stabilizer];

    aosize1 = angmom[a->L]->ncart * angmom[b->L]->ncart;
    aosize2 = a->nfunc * b->nfunc;
    aosize3 = a->nprim * b->nprim;
    aosize4 = a->ncontr * b->ncontr;

    sosizemax = aosize2 * aosize4 * MAX(a->pos->degeneracy,b->pos->degeneracy);

    if (context->aosize < aosize1*aosize3)
    {
        context->aosize = aosize1*aosize3;
        if (context->aobuf1 != NULL) FREE(context->aobuf1);
        if (context->aobuf1 != NULL) FREE(context->aobuf2);
        context->aobuf1 = MALLOC(double, aosize1*aosize3);
        context->aobuf2 = MALLOC(double, aosize1*aosize3);
        /*
         * It's okay to not free aobuf1 if only aobuf2 cannot be allocated
         */
        if (context->aobuf1 == NULL || context->aobuf2 == NULL) return 0;
    }

    if (context->sosize < sosizemax)
    {
        context->sosize = sosizemax;
        if (context->integrals != NULL) FREE(context->integrals);
        context->integrals = MALLOC(double, sosizemax);
        /*
         * Similarly, aobuf[12] need not be free'd here
         */
        if (context->integrals == NULL) return 0;
    }

    size_t worksize = 3*(a->L+1)*(b->L+1);
    if (reserve_workspaces(&work, worksize) == -1) return 0;

    if (beta == 0.0)
    {
        zero(sosizemax, context->integrals, 1);
    }
    else if (beta != 1.0)
    {
        dscal(sosizemax, beta, context->integrals, 1);
    }

    coef = alpha * (double)order / (double)lambdar;

    for (i = 0;i < staborder[dcrr];i++)
    {
        r = stabs[dcrr][i];

#ifdef ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) default(shared) private(e,f,m,idx2)
#endif //ENABLE_OPENMP
        for (m = 0;m < a->nprim * b->nprim;m++)
        {
            f = m / a->nprim;
            e = m - f * a->nprim;

            idx2 = aosize1 * m;

            oviprim(a->L, b->L, a->pos->centers[0], b->pos->centers[b->pos->centermap[r]],
                    a->exponents[e], b->exponents[f], &context->aobuf2[idx2],
                    active_workspace(work));
        }

        prim2contr2r(a, b, aosize1, context->aobuf2, context->aobuf1);
        cart2spher2r(a, b, aosize4, context->aobuf1, context->aobuf2);

        transpose(aosize2, aosize4, coef, context->aobuf2, aosize2,
        		                     0.0, context->aobuf1, aosize4);

        ao2so2(a, b, aosize4, r, context->aobuf1, context->integrals);
    }

    context->a = (shell_t*)a;
    context->b = (shell_t*)b;
    context->num_integrals = 0;
    context->num_processed = 0;

    for (j = 0;j < b->nfunc;j++)
    {
        for (i = 0;i < a->nfunc;i++)
        {
            for (s = 0;s < b->pos->degeneracy;s++)
            {
                for (r = 0;r < a->pos->degeneracy;r++)
                {
                    w = a->irreps[i][r];
                    x = b->irreps[j][s];

                    if (dirprd[w][x] != SLIDE_IRREP_TOT_SYM) continue;

                    context->num_integrals += a->ncontr*b->ncontr;
                }
            }
        }
    }

    release_workspaces(work);

    return context->num_integrals;
}

size_t SLIDE_calc_kei(context_t* context, const double alpha, const double beta, const shell_t* a, const shell_t* b)
{
    int i, j, r, s;
    int e, f, m;
    int w, x;
    size_t idx2;
    size_t aosize1, aosize2, aosize3, aosize4, sosizemax;
    int dcrr, lambdar;
    double coef;
    workspace_t work1, work2;

    dcrr = dcrindex[a->pos->stabilizer][b->pos->stabilizer];
    lambdar = dcrdeg[a->pos->stabilizer][b->pos->stabilizer];

    aosize1 = angmom[a->L]->ncart * angmom[b->L]->ncart;
    aosize2 = a->nfunc * b->nfunc;
    aosize3 = a->nprim * b->nprim;
    aosize4 = a->ncontr * b->ncontr;

    sosizemax = aosize2 * aosize4 * MAX(a->pos->degeneracy,b->pos->degeneracy);

    if (context->aosize < aosize1*aosize3)
    {
        context->aosize = aosize1*aosize3;
        if (context->aobuf1 != NULL) FREE(context->aobuf1);
        if (context->aobuf1 != NULL) FREE(context->aobuf2);
        context->aobuf1 = MALLOC(double, aosize1*aosize3);
        context->aobuf2 = MALLOC(double, aosize1*aosize3);
        /*
         * It's okay to not free aobuf1 if only aobuf2 cannot be allocated
         */
        if (context->aobuf1 == NULL || context->aobuf2 == NULL) return 0;
    }

    if (context->sosize < sosizemax)
    {
        context->sosize = sosizemax;
        if (context->integrals != NULL) FREE(context->integrals);
        context->integrals = MALLOC(double, sosizemax);
        /*
         * Similarly, aobuf[12] need not be free'd here
         */
        if (context->integrals == NULL) return 0;
    }

    size_t worksize1 = 3*(a->L+2)*(b->L+2);
    if (reserve_workspaces(&work1, worksize1) == -1) return 0;

    size_t worksize2 = 3*(a->L+1)*(b->L+1);
    if (reserve_workspaces(&work2, worksize2) == -1)
	{
        release_workspaces(work1);
    	return 0;
	}

    if (beta == 0.0)
    {
        zero(sosizemax, context->integrals, 1);
    }
    else if (beta != 1.0)
    {
        dscal(sosizemax, beta, context->integrals, 1);
    }

    coef = alpha * (double)order / (double)lambdar;

    for (i = 0;i < staborder[dcrr];i++)
    {
        r = stabs[dcrr][i];

#ifdef ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) default(shared) private(e,f,m,idx2)
#endif //ENABLE_OPENMP
        for (m = 0;m < a->nprim * b->nprim;m++)
        {
            f = m / a->nprim;
            e = m - f * a->nprim;

            idx2 = aosize1 * m;

            keiprim(a->L, b->L, a->pos->centers[0], b->pos->centers[b->pos->centermap[r]],
                    a->exponents[e], b->exponents[f], &context->aobuf2[idx2],
                    active_workspace(work1), active_workspace(work2));
        }

        prim2contr2r(a, b, aosize1, context->aobuf2, context->aobuf1);
        cart2spher2r(a, b, aosize4, context->aobuf1, context->aobuf2);

        transpose(aosize2, aosize4, coef, context->aobuf2, aosize2,
        		                     0.0, context->aobuf1, aosize4);

        ao2so2(a, b, aosize4, r, context->aobuf1, context->integrals);
    }

    context->a = (shell_t*)a;
    context->b = (shell_t*)b;
    context->num_integrals = 0;
    context->num_processed = 0;

    for (j = 0;j < b->nfunc;j++)
    {
        for (i = 0;i < a->nfunc;i++)
        {
            for (s = 0;s < b->pos->degeneracy;s++)
            {
                for (r = 0;r < a->pos->degeneracy;r++)
                {
                    w = a->irreps[i][r];
                    x = b->irreps[j][s];

                    if (dirprd[w][x] != SLIDE_IRREP_TOT_SYM) continue;

                    context->num_integrals += a->ncontr*b->ncontr;
                }
            }
        }
    }

    release_workspaces(work1);
    release_workspaces(work2);

    return context->num_integrals;
}

size_t SLIDE_calc_nai(context_t* context, const double alpha, const double beta, const shell_t* a, const shell_t* b, const center_t** centers, const int ncenters)
{
    int atm, c;
    int i, j, r, s;
    int e, f, m;
    int w, x;
    size_t idx2;
    size_t aosize1, aosize2, aosize3, aosize4, sosizemax;
    int dcrr, lambdar;
    double coef;
    workspace_t work;

    dcrr = dcrindex[a->pos->stabilizer][b->pos->stabilizer];
    lambdar = dcrdeg[a->pos->stabilizer][b->pos->stabilizer];

    aosize1 = angmom[a->L]->ncart * angmom[b->L]->ncart;
    aosize2 = a->nfunc * b->nfunc;
    aosize3 = a->nprim * b->nprim;
    aosize4 = a->ncontr * b->ncontr;

    sosizemax = aosize2 * aosize4 * MAX(a->pos->degeneracy,b->pos->degeneracy);

    if (context->aosize < aosize1*aosize3)
    {
        context->aosize = aosize1*aosize3;
        if (context->aobuf1 != NULL) FREE(context->aobuf1);
        if (context->aobuf1 != NULL) FREE(context->aobuf2);
        context->aobuf1 = MALLOC(double, aosize1*aosize3);
        context->aobuf2 = MALLOC(double, aosize1*aosize3);
        /*
         * It's okay to not free aobuf1 if only aobuf2 cannot be allocated
         */
        if (context->aobuf1 == NULL || context->aobuf2 == NULL) return 0;
    }

    if (context->sosize < sosizemax)
    {
        context->sosize = sosizemax;
        if (context->integrals != NULL) FREE(context->integrals);
        context->integrals = MALLOC(double, sosizemax);
        /*
         * Similarly, aobuf[12] need not be free'd here
         */
        if (context->integrals == NULL) return 0;
    }

    size_t worksize = (a->L+1)*(b->L+1)*(a->L+b->L+1);
    if (reserve_workspaces(&work, worksize) == -1) return 0;

    if (beta == 0.0)
    {
        zero(sosizemax, context->integrals, 1);
    }
    else if (beta != 1.0)
    {
        dscal(sosizemax, beta, context->integrals, 1);
    }

    coef = alpha * (double)order / (double)lambdar;

    for (i = 0;i < staborder[dcrr];i++)
    {
        r = stabs[dcrr][i];

#ifdef ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) default(shared) private(e,f,m,atm,c,idx2)
#endif //ENABLE_OPENMP
        for (m = 0;m < a->nprim * b->nprim;m++)
        {
            f = m / a->nprim;
            e = m - f * a->nprim;

            idx2 = aosize1 * m;
            zero(aosize1, &context->aobuf2[idx2], 1);

            for (atm = 0;atm < ncenters;atm++)
            {
                for (c = 0;c < centers[atm]->degeneracy;c++)
                {
                    naiprim(a->L, b->L, a->pos->centers[0], b->pos->centers[b->pos->centermap[r]], centers[atm]->centers[c],
                            a->exponents[e], b->exponents[f], centers[atm]->element->charge, &context->aobuf2[idx2],
                            active_workspace(work));
                }
            }
        }

        prim2contr2r(a, b, aosize1, context->aobuf2, context->aobuf1);
        cart2spher2r(a, b, aosize4, context->aobuf1, context->aobuf2);

        transpose(aosize2, aosize4, coef, context->aobuf2, aosize2,
        		                     0.0, context->aobuf1, aosize4);

        ao2so2(a, b, aosize4, r, context->aobuf1, context->integrals);
    }

    context->a = (shell_t*)a;
    context->b = (shell_t*)b;
    context->num_integrals = 0;
    context->num_processed = 0;

    for (j = 0;j < b->nfunc;j++)
    {
        for (i = 0;i < a->nfunc;i++)
        {
            for (s = 0;s < b->pos->degeneracy;s++)
            {
                for (r = 0;r < a->pos->degeneracy;r++)
                {
                    w = a->irreps[i][r];
                    x = b->irreps[j][s];

                    if (dirprd[w][x] != SLIDE_IRREP_TOT_SYM) continue;

                    context->num_integrals += a->ncontr*b->ncontr;
                }
            }
        }
    }

    release_workspaces(work);

    return context->num_integrals;
}

size_t SLIDE_calc_moment(context_t* context, const double alpha, const double beta, const shell_t* a, const shell_t* b, const int L, const double *pos)
{
    int i, j, r, s;
    int e, f, m;
    int w, x;
    size_t idx2;
    size_t aosize1, aosize2, aosize3, aosize4, sosizemax;
    int dcrr, lambdar;
    double coef;
    workspace_t work;
    double pos0[3] = {0, 0, 0};

    if (pos == NULL) pos = pos0;

    dcrr = dcrindex[a->pos->stabilizer][b->pos->stabilizer];
    lambdar = dcrdeg[a->pos->stabilizer][b->pos->stabilizer];

    aosize1 = angmom[a->L]->ncart * angmom[b->L]->ncart * angmom[L]->ncart;
    aosize2 = a->nfunc * b->nfunc;
    aosize3 = a->nprim * b->nprim;
    aosize4 = a->ncontr * b->ncontr * angmom[L]->ncart;

    sosizemax = aosize2 * aosize4 * MAX(a->pos->degeneracy,b->pos->degeneracy);

    if (context->aosize < aosize1*aosize3)
    {
        context->aosize = aosize1*aosize3;
        if (context->aobuf1 != NULL) FREE(context->aobuf1);
        if (context->aobuf1 != NULL) FREE(context->aobuf2);
        context->aobuf1 = MALLOC(double, aosize1*aosize3);
        context->aobuf2 = MALLOC(double, aosize1*aosize3);
        /*
         * It's okay to not free aobuf1 if only aobuf2 cannot be allocated
         */
        if (context->aobuf1 == NULL || context->aobuf2 == NULL) return 0;
    }

    if (context->sosize < sosizemax)
    {
        context->sosize = sosizemax;
        if (context->integrals != NULL) FREE(context->integrals);
        context->integrals = MALLOC(double, sosizemax);
        /*
         * Similarly, aobuf[12] need not be free'd here
         */
        if (context->integrals == NULL) return 0;
    }

    size_t worksize = (a->L+1)*(b->L+1)*(L+1);
    if (reserve_workspaces(&work, worksize) == -1) return 0;

    if (beta == 0.0)
    {
        zero(sosizemax, context->integrals, 1);
    }
    else if (beta != 1.0)
    {
        dscal(sosizemax, beta, context->integrals, 1);
    }

    coef = alpha * (double)order / (double)lambdar;

    for (i = 0;i < staborder[dcrr];i++)
    {
        r = stabs[dcrr][i];

#ifdef ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) default(shared) private(e,f,m,atm,c,idx2)
#endif //ENABLE_OPENMP
        for (m = 0;m < a->nprim * b->nprim;m++)
        {
            f = m / a->nprim;
            e = m - f * a->nprim;

            idx2 = aosize1 * m;
            zero(aosize1, &context->aobuf2[idx2], 1);

            momprim(a->L, b->L, a->pos->centers[0], b->pos->centers[b->pos->centermap[r]],
                    a->exponents[e], b->exponents[f], L, pos, &context->aobuf2[idx2],
                    active_workspace(work));
        }

        // cbaqp -> srcba
        prim2contr2r(a, b, aosize1, context->aobuf2, context->aobuf1);
        // srcba -> jisrc
        cart2spher2r(a, b, aosize4, context->aobuf1, context->aobuf2);

        // jisrc -> srcji
        transpose(aosize2, aosize4, coef, context->aobuf2, aosize2, 0.0, context->aobuf1, aosize4);

        ao2so2(a, b, aosize4, r, context->aobuf1, context->integrals);
    }

    context->a = (shell_t*)a;
    context->b = (shell_t*)b;
    context->num_integrals = 0;
    context->num_processed = 0;

    for (j = 0;j < b->nfunc;j++)
    {
        for (i = 0;i < a->nfunc;i++)
        {
            for (s = 0;s < b->pos->degeneracy;s++)
            {
                for (r = 0;r < a->pos->degeneracy;r++)
                {
                    w = a->irreps[i][r];
                    x = b->irreps[j][s];

                    if (dirprd[w][x] != SLIDE_IRREP_TOT_SYM) continue;

                    context->num_integrals += aosize4;
                }
            }
        }
    }

    release_workspaces(work);

    return context->num_integrals;
}
