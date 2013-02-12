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

/*
 * character constants;
 */
const char* angmomlabels[9] = { "s", "p", "d", "f", "g", "h", "i", "k", "l" };
const char* oplabels[8] = { "E  ", "C2z", "C2y", "C2x", "i  ", "sxy", "sxz", "syz" };
const char* grouplabels[8] = { "C1 ", "C2 ", "Cs ", "Ci ", "C2v", "D2 ", "C2h", "D2h" };
char* irreplabels[8];

/*
 * point group-specific data
 * variables are explained in initsymmetry
 */
int group, order;
int ops[8], staborder[8];
int chars[8][8], stabs[8][8], dcrindex[8][8], dcrdeg[8][8];

/*
 * ordering for angular momentum components
 */
int cartesian_order[9][55];
int spherical_order[9][55];
int cartesian_order_inverse[9][55];
int spherical_order_inverse[9][55];

/*
 * parity and spherical transformation data for each angular momentum
 * variables are explained in newangulardata
 */
angular_data_t* angmom[SLIDE_MAX_L+1];

double geomtol;
int ordering;

int binom(const int a, const int b)
{
    int i, j;

    if (b < 0 || b > a) return 0;

    j = 1;
    for (i = 1;i <= MIN(b,a-b);i++)
    {
        j = (j*(a-i+1))/i;
    }

    return j;
}

int fact(const int n)
{
    int i, j;

    j = 1;

    for (i = n;i > 1;i--)
    {
        j = j*i;
    }

    return j;
}

int dfact(const int n)
{
    int i, j;

    j = 1;

    for (i = n;i > 1;i -= 2)
    {
        j = j*i;
    }

    return j;
}

angular_data_t* newangulardata(const int L)
{
    angular_data_t* angdata;
    int nc, ns, x, y, z, i, j, L2, m;

    angdata = MALLOC(angular_data_t, 1);
    if (angdata == NULL)
    {
        return NULL;
    }

    /*
     * number of cartesian of spherical functions
     * if contaminants are used, they are the same
     */
    nc = (L+1)*(L+2)/2;
    ns = 2*L+1;

    angdata->ncart = nc;
    angdata->nspher = ns;

    angdata->cartlabels = MALLOC(char*, nc);
    if (angdata->cartlabels == NULL)
    {
    	FREE(angdata);
        return NULL;
    }
    angdata->spherlabels = MALLOC(char*, nc);
    if (angdata->cartlabels == NULL)
    {
        FREE(angdata->spherlabels);
    	FREE(angdata);
        return NULL;
    }

    angdata->cart2spher = MALLOC(double, nc*nc);
    if (angdata->cart2spher == NULL)
    {
        FREE(angdata->spherlabels);
        FREE(angdata->cartlabels);
    	FREE(angdata);
        return NULL;
    }

    for (i = 0;i < 8;i++)
    {
        angdata->cartparity[i] = MALLOC(int, nc);
        angdata->spherparity[i] = MALLOC(int, nc);

        if (angdata->cartparity[i] == NULL || angdata->spherparity[i] == NULL)
        {
            if (angdata->cartparity[i] != NULL) FREE(angdata->cartparity[i]);
            if (angdata->spherparity[i] != NULL) FREE(angdata->spherparity[i]);

            for (i--;i >= 0;i--)
            {
                FREE(angdata->cartparity[i]);
                FREE(angdata->spherparity[i]);
            }

            FREE(angdata->spherlabels);
            FREE(angdata->cartlabels);
            FREE(angdata->cart2spher);
        	FREE(angdata);

            return NULL;
        }
    }

    for (i = 0;i < nc;i++)
    {
        angdata->cartlabels[i] = MALLOC(char, 5);
        angdata->spherlabels[i] = MALLOC(char, 5);

        if (angdata->cartlabels[i] == NULL || angdata->spherlabels[i] == NULL)
        {

            if (angdata->cartlabels[i] != NULL) FREE(angdata->cartlabels[i]);
            if (angdata->spherlabels[i] != NULL) FREE(angdata->spherlabels[i]);

            for (i--;i >= 0;i--)
            {
            	FREE(angdata->cartlabels[i]);
            	FREE(angdata->spherlabels[i]);
            }

            for (i = 0;i < 8;i++)
            {
            	FREE(angdata->cartparity[i]);
            	FREE(angdata->spherparity[i]);
            }

        	FREE(angdata->spherlabels);
        	FREE(angdata->cartlabels);
        	FREE(angdata->cart2spher);
        	FREE(angdata);

            return NULL;
        }
    }

    /*
     * generate labels like "zz  " or "g211" for lxyz = (0,0,2) and (2,1,1) respectively
     */
    if (L == 0)
    {
        strcpy(angdata->cartlabels[0], "s   ");
    }
    else if (L == 1)
    {
        strcpy(angdata->cartlabels[0], "z   ");
        strcpy(angdata->cartlabels[1], "y   ");
        strcpy(angdata->cartlabels[2], "x   ");
    }
    else if (L == 2)
    {
        strcpy(angdata->cartlabels[0], "zz  ");
        strcpy(angdata->cartlabels[1], "yz  ");
        strcpy(angdata->cartlabels[2], "yy  ");
        strcpy(angdata->cartlabels[3], "xz  ");
        strcpy(angdata->cartlabels[4], "xy  ");
        strcpy(angdata->cartlabels[5], "xx  ");
    }
    else
    {
        i = 0;
        for (x = 0;x <= L;x++)
        {
            for (y = 0;y <= L-x;y++)
            {
                z = L-x-y;
                snprintf(angdata->cartlabels[i], 5, "%s%d%d%d", angmomlabels[L], x, y, z);
                i++;
            }
        }
    }

    /*
     * generate labels like "4f2+" for (n+1)lm = (4,3,+2) etc.;
     * contaminants have l < n, e.g. "4p0 " for (n+1)lm = (4,1,0);
     */
    i = 0;
    for (L2 = L;L2 >= 0;L2 -= 2)
    {
        for (m = L2;m > 0;m--)
        {
            snprintf(angdata->spherlabels[i], 5, "%d%s%d+", L+1, angmomlabels[L2], m);
            i++;
            snprintf(angdata->spherlabels[i], 5, "%d%s%d-", L+1, angmomlabels[L2], m);
            i++;
        }
        snprintf(angdata->spherlabels[i], 5, "%d%s0 ", L+1, angmomlabels[L2]);
        i++;
    }

    /*
     * signs for cartesian functions:;
     *
     * E   - (-1)^(0);
     * C2z - (-1)^(lx+ly);
     * C2y - (-1)^(lx+lz);
     * C2x - (-1)^(ly+lz);
     * i   - (-1)^(lx+ly+lz);
     * sxy - (-1)^(lz);
     * sxz - (-1)^(ly);
     * syz - (-1)^(lx);
     */
    i = 0;
    for (x = 0;x <= L;x++)
    {
        for (y = 0;y <= L-x;y++)
        {
            z = L-x-y;
            angdata->cartparity[0][i] = 1;
            angdata->cartparity[1][i] = 1-2*((x+y)%2);
            angdata->cartparity[2][i] = 1-2*((x+z)%2);
            angdata->cartparity[3][i] = 1-2*((y+z)%2);
            angdata->cartparity[4][i] = 1-2*((x+y+z)%2);
            angdata->cartparity[5][i] = 1-2*(z%2);
            angdata->cartparity[6][i] = 1-2*(y%2);
            angdata->cartparity[7][i] = 1-2*(x%2);
            i++;
        }
    }

    /*
     * signs for spherical functions:;
     *
     * E   - sign(1)*(-1)^(0);
     * C2z - sign(1)*(-1)^(m);
     * C2y - sign(m)*(-1)^(l);
     * C2x - sign(m)*(-1)^(l+m);
     * i   - sign(1)*(-1)^(l);
     * sxy - sign(1)*(-1)^(l+m);
     * sxz - sign(m)*(-1)^(0);
     * syz - sign(m)*(-1)^(m);
     */
    i = 0;
    for (L2 = L;L2 >= 0;L2 -= 2)
    {
        for (m = L2;m > 0;m--)
        {
            angdata->spherparity[0][i] = 1;
            angdata->spherparity[1][i] = 1-2*(m%2);
            angdata->spherparity[2][i] = 1-2*(L2%2);
            angdata->spherparity[3][i] = 1-2*((L2+m)%2);
            angdata->spherparity[4][i] = 1-2*(L2%2);
            angdata->spherparity[5][i] = 1-2*((L2+m)%2);
            angdata->spherparity[6][i] = 1;
            angdata->spherparity[7][i] = 1-2*(m%2);
            i++;
            angdata->spherparity[0][i] = 1;
            angdata->spherparity[1][i] = 1-2*(m%2);
            angdata->spherparity[2][i] = -1+2*(L2%2);
            angdata->spherparity[3][i] = -1+2*((L2+m)%2);
            angdata->spherparity[4][i] = 1-2*(L2%2);
            angdata->spherparity[5][i] = 1-2*((L2+m)%2);
            angdata->spherparity[6][i] = -1;
            angdata->spherparity[7][i] = -1+2*(m%2);
            i++;
        }
        angdata->spherparity[0][i] = 1;
        angdata->spherparity[1][i] = 1;
        angdata->spherparity[2][i] = 1-2*(L2%2);
        angdata->spherparity[3][i] = 1-2*(L2%2);
        angdata->spherparity[4][i] = 1-2*(L2%2);
        angdata->spherparity[5][i] = 1-2*(L2%2);
        angdata->spherparity[6][i] = 1;
        angdata->spherparity[7][i] = 1;
        i++;
    }

    /*
     * generate normalized cartesian-to-spherical transformation matrix;
     */
    j = 0;
    for (x = 0;x <= L;x++)
    {
        for (y = 0;y <= L-x;y++)
        {
            z = L-x-y;

            i = 0;
            for (L2 = L;L2 >= 0;L2 -= 2)
            {
                for (m = L2;m > 0;m--)
                {
                    angdata->cart2spher[i*nc+j] = cartcoef(L2, m, x, y, z);
                    i++;
                    angdata->cart2spher[i*nc+j] = cartcoef(L2, -m, x, y, z);
                    i++;
                }
                angdata->cart2spher[i*nc+j] = cartcoef(L2, 0, x, y, z);
                i++;
            }

            j++;
        }
    }

    return angdata;
}

void freeangulardata(angular_data_t* ang)
{
    int i;

    if (ang != NULL)
    {
        for (i = 0;i < 8;i++)
        {
            if (ang->cartparity[i]!=NULL) FREE(ang->cartparity[i]);
            if (ang->spherparity[i]!=NULL) FREE(ang->spherparity[i]);
        }

        for (i = 0;i < ang->ncart;i++)
        {
            if (ang->cartlabels[i]!=NULL) FREE(ang->cartlabels[i]);
            if (ang->spherlabels[i]!=NULL) FREE(ang->spherlabels[i]);
        }

        if (ang->cartlabels!=NULL) FREE(ang->cartlabels);
        if (ang->spherlabels!=NULL) FREE(ang->spherlabels);
        if (ang->cart2spher!=NULL) FREE(ang->cart2spher);

        FREE(ang);
    }
}

double cartcoef(const int l, const int m, const int lx, const int ly, const int lz)
{
    int j, k, i, am;
    double c, sum, tmp;

    am = abs(m);
    j = lx+ly-am;
    if (MODULO(j, 2) == 1) return 0.0;
    j = j/2;

    c = sqrt((double)(binom(2*lx, lx)*binom(2*ly, ly)*binom(2*lz, lz)*binom(l+am, am))/(double)(binom(2*l, l)*binom(l, am)*binom(l, lx)*binom(l-lx, ly))/(double)(dfact(2*lx-1)*dfact(2*ly-1)*dfact(2*lz-1)))/pow(2, l);
    if (m != 0) c = c*sqrt(2.0);

    if (m >= 0)
    {
        if (MODULO(am - lx, 2) == 1) return 0.0;
        if (MODULO(am - lx, 4) == 2) c = -c;
    }
    else
    {
        if (MODULO(am - lx, 2) == 0) return 0.0;
        if (MODULO(am - lx, 4) == 3) c = -c;
    }

    sum = 0.0;
    for (i = 0;i <= (l-am)/2;i++)
    {
        for (k = 0;k <= j;k++)
        {
            tmp = binom(2*l-2*i, l+am)*binom(l, i)*binom(i, j)*binom(j, k)*binom(am, lx-2*k);
            if (MODULO(i + k, 2)==1) tmp = -tmp;
            sum += tmp;
        }
    }

    return sum*c;
}

void applysymop(const int n, double (*pos)[3], const int op)
{
    int i;

    for (i = 0;i < n;i++)
    {
        switch (op)
        {
            case 1:
                pos[i][0] = -pos[i][0];
                pos[i][1] = -pos[i][1];
                break;
            case 2:
                pos[i][0] = -pos[i][0];
                pos[i][2] = -pos[i][2];
                break;
            case 3:
                pos[i][1] = -pos[i][1];
                pos[i][2] = -pos[i][2];
                break;
            case 4:
                pos[i][0] = -pos[i][0];
                pos[i][1] = -pos[i][1];
                pos[i][2] = -pos[i][2];
                break;
            case 5:
                pos[i][2] = -pos[i][2];
                break;
            case 6:
                pos[i][1] = -pos[i][1];
                break;
            case 7:
                pos[i][0] = -pos[i][0];
                break;
        }
    }
}

bool checksymop(const int n, const double (*pos)[3], const int op)
{
    int i, j;
    double x[3];

    for (i = 0;i < n;i++)
    {
        x[0] = pos[i][0];
        x[1] = pos[i][1];
        x[2] = pos[i][2];
        applysymop(1, &x, op);

        for (j = 0;j < n;j++)
        {
            if (dist(x, pos[j])<geomtol) continue;
        }

        return false;
    }

    return true;
}

int getstab(const double pos[3])
{
    int stab = 0;

    if (fabs(pos[0]) > geomtol) stab += 1;
    if (fabs(pos[1]) > geomtol) stab += 2;
    if (fabs(pos[2]) > geomtol) stab += 4;

    return stab;
}
