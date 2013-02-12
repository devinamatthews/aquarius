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
#include <string.h>

void FFUNC(slide_element_symbol,SLIDE_ELEMENT_SYMBOL)(uint64_t* element, char* symbol)
{
    element_t* e = (element_t*)(*element);
    strncpy(symbol, e->symbol, 2);
}

integer FFUNC(slide_element_nucleon,SLIDE_ELEMENT_NUCLEON)(uint64_t* element)
{
    element_t* e = (element_t*)(*element);
    return e->nucleon[e->isotope];
}

integer FFUNC(slide_element_spin,SLIDE_ELEMENT_SPIN)(uint64_t* element)
{
    element_t* e = (element_t*)(*element);
    return e->spin[e->isotope];
}

double FFUNC(slide_element_mass,SLIDE_ELEMENT_MASS)(uint64_t* element)
{
    element_t* e = (element_t*)(*element);
    return e->mass[e->isotope];
}

double FFUNC(slide_element_charge,SLIDE_ELEMENT_CHARGE)(uint64_t* element)
{
    element_t* e = (element_t*)(*element);
    return e->charge;
}

uint64_t FFUNC(slide_copy_element,SLIDE_COPY_ELEMENT)(uint64_t* element, double* charge)
{
    element_t* e = SLIDE_copy_element((element_t*)(*element));

    if (e != NULL)
    {
        e->charge = *charge;
    }

    return (uint64_t)e;
}

void FFUNC(slide_free_element,SLIDE_FREE_ELEMENT)(uint64_t* element)
{
    SLIDE_free_element((element_t*)(*element));
}

uint64_t FFUNC(slide_get_element,SLIDE_GET_ELEMENT)(char* symbol)
{
    return (uint64_t)SLIDE_get_element(symbol);
}

integer FFUNC(slide_center_degeneracy,SLIDE_CENTER_DEGENERACY)(uint64_t* center)
{
    center_t* c = (center_t*)(*center);
    return c->degeneracy;
}

uint64_t FFUNC(slide_center_element,SLIDE_CENTER_ELEMENT)(uint64_t* center)
{
    center_t* c = (center_t*)(*center);
    return (uint64_t)c->element;
}

integer FFUNC(slide_center_stabilizer,SLIDE_CENTER_STABILIZER)(uint64_t* center)
{
    center_t* c = (center_t*)(*center);
    return c->stabilizer;
}

void FFUNC(slide_center,SLIDE_CENTER)(uint64_t* center, integer* degen, double* pos)
{
    center_t* c = (center_t*)(*center);
    memcpy(pos, &c->centers[*degen-1][0], sizeof(double)*3);
}

void FFUNC(slide_center_after_op,SLIDE_CENTER_AFTER_OP)(uint64_t* center, integer* op, double* pos)
{
    center_t* c = (center_t*)(*center);
    memcpy(pos, &c->centers[0][0], sizeof(double)*3);
    applysymop(1, (double(*)[3])pos, *op-1);
}

uint64_t FFUNC(slide_new_center,SLIDE_NEW_CENTER)(double* pos, uint64_t* element)
{
    return (uint64_t)SLIDE_new_center(pos, (element_t*)(*element));
}

void FFUNC(slide_free_center,SLIDE_FREE_CENTER)(uint64_t* center)
{
    SLIDE_free_center((center_t*)(*center));
}

integer FFUNC(slide_shell_l,SLIDE_SHELL_L)(uint64_t* shell)
{
    shell_t* s = (shell_t*)(*shell);
    return s->L;
}

integer FFUNC(slide_shell_nprim,SLIDE_SHELL_NPRIM)(uint64_t* shell)
{
    shell_t* s = (shell_t*)(*shell);
    return s->nprim;
}

integer FFUNC(slide_shell_ncontr,SLIDE_SHELL_NCONTR)(uint64_t* shell)
{
    shell_t* s = (shell_t*)(*shell);
    return s->ncontr;
}

integer FFUNC(slide_shell_nfunc,SLIDE_SHELL_NFUNC)(uint64_t* shell)
{
    shell_t* s = (shell_t*)(*shell);
    return s->nfunc;
}

integer FFUNC(slide_shell_irrep,SLIDE_SHELL_IRREP)(uint64_t* shell, integer* func, integer* degen)
{
    shell_t* s = (shell_t*)(*shell);
    return s->irreps[*func-1][*degen-1]+1;
}

void FFUNC(slide_shell_irreps,SLIDE_SHELL_IRREPS)(uint64_t* shell, integer (*irreps)[8])
{
    int i, r;
    shell_t* s = (shell_t*)(*shell);

    for (i = 0;i < s->nfunc;i++)
    {
        for (r = 0;r < s->pos->degeneracy;r++)
        {
            irreps[i][r] = s->irreps[i][r]+1;
        }
    }
}

void FFUNC(slide_shell_nfunc_irrep,SLIDE_SHELL_NFUNC_IRREP)(uint64_t* shell, integer *nfunc_irrep)
{
    int i;
    shell_t* s = (shell_t*)(*shell);

    for (i = 0;i < order;i++)
    {
        nfunc_irrep[i] = s->nfunc_irrep[i];
    }
}

void FFUNC(slide_shell_idx,SLIDE_SHELL_IDX)(uint64_t* shell, integer *idx)
{
    int i;
    shell_t* s = (shell_t*)(*shell);

    for (i = 0;i < order;i++)
    {
        idx[i] = s->idx[i];
    }
}

void FFUNC(slide_shell_set_idx,SLIDE_SHELL_SET_IDX)(uint64_t* shell, integer *idx)
{
    int i;
    shell_t* s = (shell_t*)(*shell);

    for (i = 0;i < order;i++)
    {
        s->idx[i] = idx[i];
    }
}

logical FFUNC(slide_shell_spherical,SLIDE_SHELL_SPHERICAL)(uint64_t* shell)
{
    shell_t* s = (shell_t*)(*shell);
    return s->spherical;
}

logical FFUNC(slide_shell_contaminants,SLIDE_SHELL_CONTAMINANTS)(uint64_t* shell)
{
    shell_t* s = (shell_t*)(*shell);
    return s->contaminants;
}

uint64_t FFUNC(slide_shell_pos,SLIDE_SHELL_POS)(uint64_t* shell)
{
    shell_t* s = (shell_t*)(*shell);
    return (uint64_t)s->pos;
}

uint64_t FFUNC(slide_new_shell,SLIDE_NEW_SHELL)(uint64_t* pos, integer* L, integer* nprim, integer* ncontr,
                                                   logical* spherical, logical* contaminants, double* exponents,
                                                   double* coefficients, integer* idx)
{
    int i[8], j;
    double *e, *c;

    e = MALLOC(double, *nprim);
    if (e == NULL) return 0;
    c = MALLOC(double, (*nprim)*(*ncontr));
    if (c == NULL)
	{
    	FREE(e);
    	return 0;
	}

    memcpy(e, exponents, sizeof(double)*(*nprim));
    memcpy(c, coefficients, sizeof(double)*(*nprim)*(*ncontr));

    for (j = 0;j < order;j++) i[j] = idx[j];

    return (uint64_t)SLIDE_new_shell((center_t*)(*pos), *L, *nprim, *ncontr, *spherical, *contaminants, e, c, i);
}

uint64_t FFUNC(slide_copy_shell,SLIDE_COPY_SHELL)(uint64_t* shell, logical* spherical, logical* contaminants)
{
    return (uint64_t)SLIDE_copy_shell((shell_t*)(*shell), *spherical, *contaminants);
}

void FFUNC(slide_free_shell,SLIDE_FREE_SHELL)(uint64_t* shell)
{
    SLIDE_free_shell((shell_t*)(*shell));
}

uint64_t FFUNC(slide_context_num_integrals,SLIDE_CONTEXT_NUM_INTEGRALS)(uint64_t* context)
{
    context_t* c = (context_t*)(*context);
    return c->num_integrals;
}

uint64_t FFUNC(slide_context_num_processed,SLIDE_CONTEXT_NUM_PROCESSED)(uint64_t* context)
{
    context_t* c = (context_t*)(*context);
    return c->num_processed;
}

void FFUNC(slide_context_set_num_processed,SLIDE_CONTEXT_SET_NUM_PROCESSED)(uint64_t* context, uint64_t* num_processed)
{
    context_t* c = (context_t*)(*context);
    c->num_processed = *num_processed;
}

uint64_t FFUNC(slide_new_context,SLIDE_NEW_CONTEXT)()
{
    return (uint64_t)SLIDE_new_context();
}

void FFUNC(slide_free_context,SLIDE_FREE_CONTEXT)(uint64_t* context)
{
    SLIDE_free_context((context_t*)(*context));
}

integer FFUNC(slide_set_cartesian_ordering,SLIDE_SET_CARTESIAN_ORDERING)(integer* L, integer* ordering)
{
    int ordering_[(*L+1)*(*L+2)/2];
    for (int i = 0;i < (*L+1)*(*L+2)/2;i++) ordering_[i] = ordering[i]-1;
    return SLIDE_set_cartesian_ordering(*L, ordering_);
}

void FFUNC(slide_get_cartesian_ordering,SLIDE_GET_CARTESIAN_ORDERING)(integer* L, integer* ordering)
{
    const int* ordering_ = SLIDE_get_cartesian_ordering(*L);
    for (int i = 0;i < (*L+1)*(*L+2)/2;i++) ordering[i] = ordering_[i]+1;
}

integer FFUNC(slide_set_spherical_ordering,SLIDE_SET_SPHERICAL_ORDERING)(integer* L, integer* ordering)
{
    int ordering_[2*(*L)+1];
    for (int i = 0;i < 2*(*L)+1;i++) ordering_[i] = ordering[i]-1;
    return SLIDE_set_spherical_ordering(*L, ordering_);
}

void FFUNC(slide_get_spherical_ordering,SLIDE_GET_SPHERICAL_ORDERING)(integer* L, integer* ordering)
{
    const int* ordering_ = SLIDE_get_spherical_ordering(*L);
    for (int i = 0;i < 2*(*L)+1;i++) ordering[i] = ordering_[i]+1;
}

integer FFUNC(slide_set_group,SLIDE_SET_GROUP)(integer* group)
{
    return SLIDE_set_group(*group-1);
}

integer FFUNC(slide_get_group,SLIDE_GET_GROUP)()
{
    return SLIDE_get_group()+1;
}

void FFUNC(slide_set_ordering,SLIDE_SET_ORDERING)(integer* ordering)
{
    return SLIDE_set_ordering(*ordering);
}

integer FFUNC(slide_get_ordering,SLIDE_GET_ORDERING)()
{
    return SLIDE_get_ordering();
}

void FFUNC(slide_set_geometry_tolerance,SLIDE_SET_GEOMETRY_TOLERANCE)(double* geometry_tolerance)
{
    return SLIDE_set_geometry_tolerance(*geometry_tolerance);
}

double FFUNC(slide_get_geometry_tolerance,SLIDE_GET_GEOMETRY_TOLERANCE)()
{
    return SLIDE_get_geometry_tolerance();
}

integer FFUNC(slide_init,SLIDE_INIT)()
{
    return SLIDE_init();
}

integer FFUNC(slide_finish,SLIDE_FINISH)()
{
    return SLIDE_finish();
}

uint64_t FFUNC(slide_process_2e_ints,SLIDE_PROCESS_2E_INTS)(uint64_t* context, integer* nprocess, double* integrals,
                                                               idx4_t* indices, double* cutoff)
{
    return SLIDE_process_2e_ints((context_t*)(*context), *nprocess, integrals, indices, *cutoff);
}

uint64_t FFUNC(slide_process_1e_ints,SLIDE_PROCESS_2E_INTS)(uint64_t* context, integer* nprocess, double* integrals,
                                                               idx2_t* indices, double* cutoff)
{
    return SLIDE_process_1e_ints((context_t*)(*context), *nprocess, integrals, indices, *cutoff);
}

integer FFUNC(slide_calc_eri,SLIDE_CALC_ERI)(uint64_t* context, double* alpha, double* beta, uint64_t* a, uint64_t* b, uint64_t* c, uint64_t* d)
{
    return SLIDE_calc_eri((context_t*)(*context), *alpha, *beta, (shell_t*)(*a), (shell_t*)(*b), (shell_t*)(*c), (shell_t*)(*d));
}

integer FFUNC(slide_calc_ovi,SLIDE_CALC_OVI)(uint64_t* context, double* alpha, double* beta, uint64_t* a, uint64_t* b)
{
    return SLIDE_calc_ovi((context_t*)(*context), *alpha, *beta, (shell_t*)(*a), (shell_t*)(*b));
}

integer FFUNC(slide_calc_kei,SLIDE_CALC_KEI)(uint64_t* context, double* alpha, double* beta, uint64_t* a, uint64_t* b)
{
    return SLIDE_calc_kei((context_t*)(*context), *alpha, *beta, (shell_t*)(*a), (shell_t*)(*b));
}

integer FFUNC(slide_calc_nai,SLIDE_CALC_NAI)(uint64_t* context, double* alpha, double* beta, uint64_t* a, uint64_t* b, uint64_t* centers, integer* ncenters)
{
    return SLIDE_calc_nai((context_t*)(*context), *alpha, *beta, (shell_t*)(*a), (shell_t*)(*b), (const center_t**)centers, *ncenters);
}

integer FFUNC(slide_irrep,SLIDE_IRREP)(uint16_t* idx, integer* ijkl)
{
    return SLIDE_IRREP(*(idx+*ijkl-1))+1;
}

integer FFUNC(slide_orbital,SLIDE_ORBITAL)(uint16_t* idx, integer* ijkl)
{
    return SLIDE_ORBITAL(*(idx+*ijkl-1))+1;
}

integer FFUNC(slide_index_iscf,SLIDE_INDEX_ISCF)(uint64_t* shell, integer* func, integer* contr, integer* degen, integer* irrep)
{
    return SLIDE_INDEX_ISCF((shell_t*)(*shell), *func-1, *contr-1, *degen-1, *irrep-1);
}

integer FFUNC(slide_index_isfc,SLIDE_INDEX_ISFC)(uint64_t* shell, integer* func, integer* contr, integer* degen, integer* irrep)
{
    return SLIDE_INDEX_ISFC((shell_t*)(*shell), *func-1, *contr-1, *degen-1, *irrep-1);
}

integer FFUNC(slide_index_sicf,SLIDE_INDEX_SICF)(uint64_t* shell, integer* func, integer* contr, integer* degen, integer* irrep)
{
    return SLIDE_INDEX_SICF((shell_t*)(*shell), *func-1, *contr-1, *degen-1, *irrep-1);
}

integer FFUNC(slide_index_sifc,SLIDE_INDEX_SIFC)(uint64_t* shell, integer* func, integer* contr, integer* degen, integer* irrep)
{
    return SLIDE_INDEX_SIFC((shell_t*)(*shell), *func-1, *contr-1, *degen-1, *irrep-1);
}

integer FFUNC(slide_index_scif,SLIDE_INDEX_SCIF)(uint64_t* shell, integer* func, integer* contr, integer* degen, integer* irrep)
{
    return SLIDE_INDEX_SCIF((shell_t*)(*shell), *func-1, *contr-1, *degen-1, *irrep-1);
}

integer FFUNC(slide_index_sfic,SLIDE_INDEX_SFIC)(uint64_t* shell, integer* func, integer* contr, integer* degen, integer* irrep)
{
    return SLIDE_INDEX_SFIC((shell_t*)(*shell), *func-1, *contr-1, *degen-1, *irrep-1);
}

integer FFUNC(slide_index_scfi,SLIDE_INDEX_SCFI)(uint64_t* shell, integer* func, integer* contr, integer* degen, integer* irrep)
{
    return SLIDE_INDEX_SCFI((shell_t*)(*shell), *func-1, *contr-1, *degen-1, *irrep-1);
}

integer FFUNC(slide_index_sfci,SLIDE_INDEX_SFCI)(uint64_t* shell, integer* func, integer* contr, integer* degen, integer* irrep)
{
    return SLIDE_INDEX_SFCI((shell_t*)(*shell), *func-1, *contr-1, *degen-1, *irrep-1);
}

integer FFUNC(slide_func_cart,SLIDE_FUNC_CART)(integer* x, integer* y, integer* z)
{
    return SLIDE_FUNC_CART(*x, *y, *z);
}

integer FFUNC(slide_func_spher,SLIDE_FUNC_SPHER)(integer* n, integer* l, integer* m)
{
    return SLIDE_FUNC_SPHER(*n, *l, *m);
}

void FFUNC(slide_group_label,SLIDE_GROUP_LABEL)(integer* group, char* label)
{
    strncpy(label, SLIDE_group_label(*group-1), 3);
}

void FFUNC(slide_l_label,SLIDE_L_LABEL)(integer* L, char* label)
{
    strncpy(label, SLIDE_L_label(*L), 1);
}

void FFUNC(slide_op_label,SLIDE_OP_LABEL)(integer* op, char* label)
{
    strncpy(label, SLIDE_op_label(*op-1), 3);
}

void FFUNC(slide_irrep_label,SLIDE_IRREP_LABEL)(integer* irrep, char* label)
{
    strncpy(label, SLIDE_irrep_label(*irrep-1), 3);
}

void FFUNC(slide_func_label,SLIDE_FUNC_LABEL)(integer* L, logical* spherical, integer* func, char* label)
{
    strncpy(label, SLIDE_func_label(*L, *spherical, *func-1), 4);
}

void FFUNC(slide_shell_func_label,SLIDE_SHELL_FUNC_LABEL)(uint64_t* shell, integer* func, char* label)
{
    strncpy(label, SLIDE_shell_func_label((shell_t*)(*shell), *func-1), 4);
}

void FFUNC(slide_stabilizer,SLIDE_STABILIZER)(uint64_t* center, integer* stab)
{
    int i;
    const int* stab_ = SLIDE_stabilizer((center_t*)(*center));
    for (i = 0;i < order;i++) stab[i] = stab_[i]+1;
}

integer FFUNC(slide_group_order,SLIDE_GROUP_ORDER)()
{
    return SLIDE_group_order();
}

void FFUNC(slide_shell_func_centers,SLIDE_SHELL_FUNC_CENTERS)(uint64_t* shell, integer* func, integer* irrep, integer* proj)
{
    int i;
    int proj_[8];
    SLIDE_shell_func_centers((shell_t*)(*shell), *func-1, *irrep-1, proj_);
    for (i = 0;i < ((shell_t*)(*shell))->pos->degeneracy;i++) proj[i] = proj_[i];
}

double FFUNC(slide_nuclear_repulsion,SLIDE_NUCLEAR_REPULSION)(uint64_t* centers, integer* ncenters)
{
    return SLIDE_nuclear_repulsion((const center_t**)centers, *ncenters);
}

void FFUNC(slide_shell_ao_to_so,SLIDE_SHELL_AO_TO_SO)(uint64_t* shell, int* ao_ordering, double* aoso, int* ld)
{
    SLIDE_shell_ao_to_so((shell_t*)(*shell), *ao_ordering, aoso, *ld);
}
