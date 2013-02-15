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

int SLIDE_init()
{
    int L, i;
    int order[(SLIDE_MAX_L+1)*(SLIDE_MAX_L+2)/2];

    SLIDE_set_group(SLIDE_GROUP_C1);
    SLIDE_set_ordering(SLIDE_ORDER_SFIC);
    SLIDE_set_geometry_tolerance(1e-6);

    for (i = 0;i < (SLIDE_MAX_L+1)*(SLIDE_MAX_L+2)/2;i++) order[i] = i;

    for (L = 0;L <= SLIDE_MAX_L;L++)
    {
        if (SLIDE_set_cartesian_ordering(L, order) == -1) return -1;
        if (SLIDE_set_spherical_ordering(L, order) == -1) return -1;
        angmom[L] = newangulardata(L);
        if (angmom[L] == NULL) return -1;
    }

    calcfmtable();

    return 0;
}

int SLIDE_finish()
{
    int L;

    for (L = 0;L <= SLIDE_MAX_L;L++)
    {
        freeangulardata(angmom[L]);
    }

    return 0;
}

int SLIDE_set_cartesian_ordering(const int L, const int* ordering_)
{
    int i, ll;
    bool used[55];

    ll = (L+1)*(L+2)/2;

    for (i = 0;i < ll;i++) used[i] = false;

    for (i = 0;i < ll;i++)
    {
        if (used[ordering_[i]] || ordering_[i] < 0 || ordering_[i] > ll) return -1;
        used[ordering_[i]] = true;
        cartesian_order[L][i] = ordering_[i];
        cartesian_order_inverse[L][ordering_[i]] = i;
    }

    return 0;
}

const int* SLIDE_get_cartesian_ordering(const int L)
{
    return cartesian_order[L];
}

int SLIDE_set_spherical_ordering(const int L, const int* ordering_)
{
    int i, ll;
    bool used[19];

    ll = 2*L+1;

    for (i = 0;i < ll;i++) used[i] = false;

    for (i = 0;i < ll;i++)
    {
        if (used[ordering_[i]] || ordering_[i] < 0 || ordering_[i] > ll) return -1;
        used[ordering_[i]] = true;
        spherical_order[L][i] = ordering_[i];
        spherical_order_inverse[L][ordering_[i]] = i;
    }

    return 0;
}

const int* SLIDE_get_spherical_ordering(const int L)
{
    return spherical_order[L];
}

int SLIDE_set_group(int group_)
{
    if (group_ < 0 || group_ > 7) return -1;

    group = group_;

    order = allorders[group];
    memcpy(ops, allops[group], sizeof(int) * 8);
    memcpy(staborder, allstaborders[group], sizeof(int) * 8);
    memcpy(irreplabels, allirreplabels[group], sizeof(char*) * 8);
    memcpy(chars, allchars[group], sizeof(int) * 8 * 8);
    memcpy(stabs, allstabs[group], sizeof(int) * 8 * 8);
    memcpy(dcrindex, alldcrindexs[group], sizeof(int) * 8 * 8);
    memcpy(dcrdeg, alldcrdegs[group], sizeof(int) * 8 * 8);

    return 0;
}

int SLIDE_get_group()
{
    return group;
}

void SLIDE_set_ordering(int ordering_)
{
    ordering = ordering_;
}

int SLIDE_get_ordering()
{
    return ordering;
}

const int* SLIDE_stabilizer(center_t* center)
{
    return stabs[center->stabilizer];
}

int SLIDE_group_order()
{
    return order;
}

const char* SLIDE_group_label(const int group)
{
    return grouplabels[group];
}

const char* SLIDE_L_label(const int L)
{
    return angmomlabels[L];
}

const char* SLIDE_op_label(const int op)
{
    return oplabels[op];
}

const char* SLIDE_irrep_label(const int irrep)
{
    return irreplabels[irrep];
}

const char* SLIDE_func_label(const int L, const bool spherical, const int func)
{
    if (spherical)
    {
        return angmom[L]->spherlabels[spherical_order[L][func]];
    }
    else
    {
        return angmom[L]->cartlabels[cartesian_order[L][func]];
    }
}

const char* SLIDE_shell_func_label(const shell_t* shell, const int func)
{
    if (shell->spherical)
    {
        return angmom[shell->L]->spherlabels[spherical_order[shell->L][func]];
    }
    else
    {
        return angmom[shell->L]->cartlabels[cartesian_order[shell->L][func]];
    }
}

void SLIDE_set_geometry_tolerance(double geomtol_)
{
    geomtol = geomtol_;
}

double SLIDE_get_geometry_tolerance()
{
    return geomtol;
}

double SLIDE_nuclear_repulsion(const center_t** centers, const int ncenters)
{
    int i, j, k, l;
    double nucrep = 0.0;

    for (i = 0;i < ncenters;i++)
    {
        for (j = i;j < ncenters;j++)
        {
            for (k = 0;k < centers[i]->degeneracy;k++)
            {
                if (i == j)
                {
                    for (l = k+1;l < centers[j]->degeneracy;l++)
                    {
                        nucrep += centers[i]->element->charge*centers[j]->element->charge/
                                  dist(centers[i]->centers[k],centers[j]->centers[l]);
                    }
                }
                else
                {
                    for (l = 0;l < centers[j]->degeneracy;l++)
                    {
                        nucrep += centers[i]->element->charge*centers[j]->element->charge/
                                  dist(centers[i]->centers[k],centers[j]->centers[l]);
                    }
                }
            }
        }
    }

    return nucrep;
}
