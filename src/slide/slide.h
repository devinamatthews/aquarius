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

#ifndef _AQUARIUS_SLIDE_H_
#define _AQUARIUS_SLIDE_H_

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#define SLIDE_OP_E          0
#define SLIDE_OP_C2Z        1
#define SLIDE_OP_C2Y        2
#define SLIDE_OP_C2X        3
#define SLIDE_OP_I          4
#define SLIDE_OP_SXY        5
#define SLIDE_OP_SXZ        6
#define SLIDE_OP_SYZ        7

#define SLIDE_GROUP_C1      0
#define SLIDE_GROUP_C2      1
#define SLIDE_GROUP_CS      2
#define SLIDE_GROUP_CI      3
#define SLIDE_GROUP_C2V     4
#define SLIDE_GROUP_D2      5
#define SLIDE_GROUP_C2H     6
#define SLIDE_GROUP_D2H     7

#define SLIDE_IRREP_TOT_SYM 0

#define SLIDE_IRREP_C1_A    0

#define SLIDE_IRREP_C2_A    0
#define SLIDE_IRREP_C2_B    1

#define SLIDE_IRREP_CS_AP   0
#define SLIDE_IRREP_CS_APP  1

#define SLIDE_IRREP_CI_AG   0
#define SLIDE_IRREP_CI_AU   1

#define SLIDE_IRREP_C2V_A1  0
#define SLIDE_IRREP_C2V_A2  1
#define SLIDE_IRREP_C2V_B1  2
#define SLIDE_IRREP_C2V_B2  3

#define SLIDE_IRREP_D2_A    0
#define SLIDE_IRREP_D2_B1   1
#define SLIDE_IRREP_D2_B2   2
#define SLIDE_IRREP_D2_B3   3

#define SLIDE_IRREP_C2H_AG  0
#define SLIDE_IRREP_C2H_AU  1
#define SLIDE_IRREP_C2H_BG  2
#define SLIDE_IRREP_C2H_BU  3

#define SLIDE_IRREP_D2H_AG  0
#define SLIDE_IRREP_D2H_B1G 1
#define SLIDE_IRREP_D2H_B2G 2
#define SLIDE_IRREP_D2H_B3G 3
#define SLIDE_IRREP_D2H_AU  4
#define SLIDE_IRREP_D2H_B1U 5
#define SLIDE_IRREP_D2H_B2U 6
#define SLIDE_IRREP_D2H_B3U 7

/**
 * @defgroup Orderings Orbital ordering schemes
 *
 * @see SLIDE_FUNC_CART
 * @see SLIDE_FUNC_SPHER
 *
 * @{*/
/**
 * Group orbitals fist by irrep and then by shell, with each shell being grouped first by contraction
 */
#define SLIDE_ORDER_ISCF    0
/**
 * Group orbitals fist by irrep and then by shell, with each shell being grouped first by contraction
 */
#define SLIDE_ORDER_ISFC    1
/**
 * Group orbitals fist by shell, with each shell being grouped first by irrep then by contraction
 */
#define SLIDE_ORDER_SICF    2
/**
 * Group orbitals fist by shell, with each shell being grouped first by irrep then by angular momentum
 */
#define SLIDE_ORDER_SIFC    3
/**
 * Group orbitals fist by shell, with each shell being grouped first by contraction then by irrep
 */
#define SLIDE_ORDER_SCIF    4
/**
 * Group orbitals fist by shell, with each shell being grouped first by contraction then by angular momentum
 */
#define SLIDE_ORDER_SCFI    5
/**
 * Group orbitals fist by shell, with each shell being grouped first by angular momentum then by irrep
 */
#define SLIDE_ORDER_SFIC    6
/**
 * Group orbitals fist by shell, with each shell being grouped first by angular momentum then by contraction
 */
#define SLIDE_ORDER_SFCI    7
/**@}*/

#define SLIDE_SHELL_ORDER_FCD   0
#define SLIDE_SHELL_ORDER_FDC   1
#define SLIDE_SHELL_ORDER_CFD   2
#define SLIDE_SHELL_ORDER_CDF   3
#define SLIDE_SHELL_ORDER_DFC   4
#define SLIDE_SHELL_ORDER_DCF   5

/**
 * Maximum total angular momentum supported.
 */
#define SLIDE_MAX_L         8

/**
 * Extract the irrep number of an orbital index (element of an idx4_t or idx2_t).
 */
#define SLIDE_IRREP(idx) (((idx)>>13)&0x7)
/**
 * Extract the orbital number of an orbital index (element of an idx4_t or idx2_t).
 */
#define SLIDE_ORBITAL(idx) ((idx)&0x1FFF)

/**
 * @defgroup IndexMacros Macros for computing function indices.
 *
 * @see SLIDE_FUNC_CART
 * @see SLIDE_FUNC_SPHER
 *
 * @{
 *
 * @param shell shell_t* structure containing the function.
 * @param func Index of the function in angular momentum.
 * @param contr Index of the function in contraction.
 * @param degen Index of the function in symmetry.
 * @param irrep Irrep of this function (same as shell->irreps[func][degen]).
 */
#define SLIDE_INDEX_ISCF(shell,func,contr,degen,irrep) \
        MAKEIDX(irrep, (shell)->idx[irrep] + (contr)*(shell)->nfunc_irrep[irrep] + (shell)->func_irrep[func][degen]);
#define SLIDE_INDEX_SICF(shell,func,contr,degen,irrep) \
        MAKEIDX(irrep, (shell)->idx[irrep] + (contr)*(shell)->nfunc_irrep[irrep] + (shell)->func_irrep[func][degen]);
#define SLIDE_INDEX_ISFC(shell,func,contr,degen,irrep) \
        MAKEIDX(irrep, (shell)->idx[irrep] + (contr) + (shell)->func_irrep[func][degen]*(shell)->ncontr);
#define SLIDE_INDEX_SIFC(shell,func,contr,degen,irrep) \
        MAKEIDX(irrep, (shell)->idx[irrep] + (contr) + (shell)->func_irrep[func][degen]*(shell)->ncontr);
#define SLIDE_INDEX_SCIF(shell,func,contr,degen,irrep) \
        MAKEIDX(irrep, (shell)->idx[irrep] + (contr)*(shell)->nfunc*(shell)->pos->degeneracy + (shell)->func_irrep[func][degen]);
#define SLIDE_INDEX_SCFI(shell,func,contr,degen,irrep) \
        MAKEIDX(irrep, (shell)->idx[irrep] + ((contr)*(shell)->nfunc + (func))*(shell)->pos->degeneracy + (shell)->irrep_pos[func][irrep]);
#define SLIDE_INDEX_SFIC(shell,func,contr,degen,irrep) \
        MAKEIDX(irrep, (shell)->idx[irrep] + (contr) + ((func)*(shell)->pos->degeneracy + (shell)->irrep_pos[func][irrep])*(shell)->ncontr);
#define SLIDE_INDEX_SFCI(shell,func,contr,degen,irrep) \
        MAKEIDX(irrep, (shell)->idx[irrep] + ((contr) + (func)*(shell)->ncontr)*(shell)->pos->degeneracy + (shell)->irrep_pos[func][irrep]);
/**@}*/

/**
 * Compute the index of a function in cartesian angular momentum.
 */
#define SLIDE_FUNC_CART(x,y,z) ((((x)*(3+(x)+2*((y)+(z))))/2) + (y))
/**
 * Compute the index of a function in spherical harmonic angular momentum.
 *
 * Regular spherical harmonics are referenced by n=l, l>=m>=-l. Contaminants may also be referenced by
 * n>l>=0, n-l even.
 */
#define SLIDE_FUNC_SPHER(n,l,m) ((((n)-(l))*((n)+(l)-1))/2 + 2*(n) + ((m) > 0 ? -2*(m) : 2*(m)+1))

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * Describes an element. Elements returned by SLIDE_get_element() are const, but an element which may be modified
 * can be returned by SLIDE_copy_element(). Currently, only the charge may be modified.
 */
typedef struct
{
    double charge;      /**< Charge of this element in au. */
    const int Z;        /**< Atomic number. */
    const int isotope;  /**< Isotope number of this element. This is not e.g. the nucleon number, but an index
                             into the spin, mass, and nucleon arrays. */
    const char* symbol; /**< Symbol of this element, e.g. H or Li. */
    const char* name;   /**< Full name of this element. */
    const int* spin;    /**< Spin of each isotope for this element, given in units of \f$\frac{\hbar}{2}\f$. */
    const double* mass; /**< Mass of each isotope for this element, given in units of amu. */
    const int* nucleon; /**< Nucleon number (P+N) of each isotope for this element. */
} element_t;

/**
 * Describes the position of an atom or orbit of degenerate atoms (see SLIDE_new_center()).
 */
typedef struct
{
    int degeneracy;         /**< The number of atoms in this orbit. */
    int stabilizer;         /**< The set of symmetry operations in the current point group which leave
                                 each atom in this orbit unchanged. */
    element_t* element;     /**< The element of this atom. */
    int centermap[8];       /**< For each symmetry operation in the current point group, the index into
                                 centers[8] to which the operation transforms the representative atom
                                 of this orbit. */
    double centers[8][3];   /**< The positions of each degenerate atom in this orbit. */
} center_t;

/**
 * Contains symmetry, primitive, and general contraction information for a shell of
 * basis functions of a given total angular momentum.
 */
typedef struct
{
    int L;                  /**< Total angular momentum for this shell. */
    int nprim;              /**< Number of primitive Gaussians in this shell. */
    int ncontr;             /**< Number of contracted Gaussians in this shell. */
    int nfunc;              /**< Number of angular momentum components in this shell
                                 (after shpherical harmonic transformation if applied). */
    int nfunc_irrep[8];
    int idx[8];
    bool spherical;         /**< True if this shell should be transformed to the spherical harmonic basis. */
    bool contaminants;      /**< True if contaminants should be retained in the transformation of this
                                 shell to the spherical harmonic basis. */
    int (*func_irrep)[8];   /**< Index of each degenerate component of each function within its irrep */
    int (*irrep_pos)[8];
    int (*irreps)[8];       /**< For each degenerate shell in this orbit (see SLIDE_new_center()),
                                 and array of the irreducible representations in the current point group
                                 for each angular momentum component. */
    double* exponents;      /**< The exponents of the primitive Gaussians. */
    double* coefficients;   /**< The transformation coefficients from primitive to contracted Gaussians,
                                 stored as coefficients[a*nprim + i], where i is a primitive index
                                 and a is a contracted index. */
    center_t* pos;          /**< The (possibly degenerate) center on which this shell is placed. */
} shell_t;

/**
 * The orbital indices of a four-center (two-electron) integral. Each index is composed of an irrep designation in
 * the uppermost three bits, and an orbital number in the remaining bits.
 *
 * @see SLIDE_IRREP()
 * @see SLIDE_INDEX()
 */
typedef struct
{
    uint16_t i;
    uint16_t j;
    uint16_t k;
    uint16_t l;
} idx4_t;

/**
 * The orbital indices of a two-center (one-electron) integral. Each index is composed of an irrep designation in
 * the uppermost three bits, and an orbital number in the remaining bits.
 *
 * @see SLIDE_IRREP()
 * @see SLIDE_INDEX()
 */
typedef struct
{
    uint16_t i;
    uint16_t j;
} idx2_t;

/**
 * A structure containing the scratch arrays needed in computing a block of integrals, and information describing
 * the block of integrals after computation. Arrays need not be allocated before the structure is used in a
 * computation, but the buffer sizes should be 0 in this case (this need not be done manually, though,
 * if SLIDE_new_context() is used).
 */
typedef struct
{
    size_t aosize;          /**< Size of the buffers used in computing the primitive cartesian AO-basis integrals,
                                 and in the transformation from primitive to contracted and cartesian to spherical. */
    size_t sosize;          /**< Size of the buffer used for storing the final transformed SO-basis integrals. */
    double* aobuf1;         /**< Buffer used to hold untransformed integrals */
    double* aobuf2;         /**< Buffer used to hold untransformed integrals */
    double* integrals;      /**< Buffer used to hold the final, transformed integrals */
    shell_t *a, *b, *c, *d; /**< The shells which define this block of integrals (may be only a and b for two-center integrals). */
    size_t num_integrals;   /**< The total number of integrals generated, not excluding zero or near-zero integrals. */
    size_t num_processed;   /**< The number of integrals in this block which have been processed by SLIDE_process_2e_ints() or
                                 SLIDE_process_1e_ints(). May be reset to a lower value to reprocess some integrals. */
} context_t;

/**
 * Allocate and initialize an integral computation context.
 */
context_t* SLIDE_new_context();

/**
 * Free an integral computation context.
 */
void SLIDE_free_context(context_t* context);

int SLIDE_set_cartesian_ordering(const int L, const int* ordering_);

const int* SLIDE_get_cartesian_ordering(const int L);

int SLIDE_set_spherical_ordering(const int L, const int* ordering_);

const int* SLIDE_get_spherical_ordering(const int L);

/**
 * Set the point group for further computations (must be one of SLIDE_GROUP_*). This will invalidate any
 * symmetry designations on shells (see SLIDE_shell_sym()).
 */
int SLIDE_set_group(const int group_);

/**
 * Return the current point group.
 */
int SLIDE_get_group();

/**
 * Set the ordering used to calculate the indices of specific functions (must be one of SLIDE_ORDER_*).
 */
void SLIDE_set_ordering(const int ordering_);

/**
 * Return the current ordering.
 */
int SLIDE_get_ordering();

/**
 * Return the operations comprising the stabilizer of this center
 */
const int* SLIDE_stabilizer(center_t* center);

/**
 * Return the order (number of irreps in) the current group
 */
int SLIDE_group_order();

/**
 * Get a string representing the given group
 */
const char* SLIDE_group_label(const int group);

/**
 * Get a string representing the given angular momentum
 */
const char* SLIDE_L_label(const int L);

/**
 * Get a string representing the given operation
 */
const char* SLIDE_op_label(const int op);

/**
 * Get a string representing the given irrep (for the current group)
 */
const char* SLIDE_irrep_label(const int irrep);

/**
 * Get a string representing the given angular momentum type
 */
const char* SLIDE_func_label(const int L, const bool spherical, const int func);

/**
 * Get a string representing the given angular momentum type
 */
const char* SLIDE_shell_func_label(const shell_t* shell, const int func);

/**
 * Set the tolerance (in au) used to determine the equivalence of two atomic centers. Default is 1e-6 au.
 */
void SLIDE_set_geometry_tolerance(const double geomtol_);

/**
 * Return the tolerance (in au) used to determine the equivalence of two atomic centers.
 */
double SLIDE_get_geometry_tolerance();

/**
 * Generate the orbit of atomic centers from the given representative position. An orbit is a set of positions which
 * is closed under the symmetry elements of the current point group.
 */
center_t* SLIDE_new_center(const double pos[3], const element_t* element);

/**
 * Free a center created by SLIDE_new_center().
 */
void SLIDE_free_center(center_t* center);

/**
 * Calculate the total nuclear repulsion for the given molecule
 */
double SLIDE_nuclear_repulsion(const center_t** centers, const int ncenters);

/**
 * Return the element designated by the given one- or two-letter symbol (e.g. H or Li).
 */
const element_t* SLIDE_get_element(const char* symbol);

/**
 * Return a new copy of the given element structure which may be modified.
 */
element_t* SLIDE_copy_element(const element_t* element);

/**
 * Free an element created by SLIDE_copy_element().
 */
void SLIDE_free_element(element_t* element);

/**
 * Initialize a new shell structure with the given information. The exponents and coefficients arrays will not be copied,
 * and their contents may not be reused (nor should they be freed). All other fields of the structure will be filled in
 * automatically.
 */
shell_t* SLIDE_new_shell(const center_t* pos, const int L, const int nprim, const int ncontr, const bool spherical, const bool contaminants,
                         const double* exponents, const double* coefficients, const int idx[8]);

/**
 * Return a new shell with the same basis functions as the given shell (L, exponents, coefficients, nprim, ncontr, etc.), but with a different
 * spherical harmonic transformation
 */
shell_t* SLIDE_copy_shell(const shell_t*, const bool spherical, const bool contaminants);

/**
 * Free a shell created by SLIDE_new_shell().
 */
void SLIDE_free_shell(shell_t* shell);

/**
 * Determine the symmetry of the final SO functions in the given shell. Usually called from SLIDE_new_shell(); should not
 * be called directly unless the shell has been initialized by hand.
 */
int SLIDE_shell_sym(shell_t* s);

/**
 * Determine the contribution of each of the degenerate centers in the orbit to this function (may be all zeros)
 */
void SLIDE_shell_func_centers(const shell_t* shell, const int func, const int irrep, int* proj);

/**
 * Normalize the transformation coefficients of a shell. Usually called from SLIDE_new_shell(); should not
 * be called directly unless the shell has been initialized by hand.
 *
 * The contraction coefficients are usually already normalized in the given basis set, but not always.
 */
void SLIDE_shell_norm(shell_t* s);

/**
 * Generate the normalized transformation matrix from the AO basis to the SO basis for this shell.
 *
 * The transformation matrix is stored as an ld x nso column-major matrix, with ld >= nao
 */
void SLIDE_shell_ao_to_so(const shell_t* shell, const int primitive_ordering, double* aoso, const int ld);

/**
 * Initialize all structures necessary for computing integrals (besides user-provided structures such as centers and shells).
 * The point group is initialized to C1.
 */
int SLIDE_init();

/**
 * Free all structures associated with integral computation. There is a small amount information in static storage (not
 * per-thread) which cannot be freed.
 */
int SLIDE_finish();

/**
 * Filter a number of integrals from the buffer of the given context into the given buffer, removing integrals whose absolute
 * vaue is <= the given cutoff. Indices for each integral are calculated and stored as well. The context is modified to
 * reflect the integrals processed.
 *
 * @param context The integral computation context to copy integrals from. context->num_processed is updated after processing.
 * @param nprocess The maximum number of integrals to process from the context. Fewer may be processed if the end of the
 *                 buffer is reached.
 * @param integrals The buffer into which to place the integrals, must be at least nprocess elements long.
 * @param indices The buffer into which to place the integral indices, must be at least nprocess elements long.
 * @param cutoff The integral cutoff value. Integrals with absolute value <= cutoff will not be retained. A negative value
 *               will retain all integrals.
 *
 * @return The number of integrals processed. A value less than nprocess means that the end of the integral buffer in the
 *         context has been reached.
 */
size_t SLIDE_process_2e_ints(context_t* context, const size_t nprocess, double* integrals, idx4_t* indices, const double cutoff);

/**
 * Filter a number of integrals from the buffer of the given context into the given buffer, removing integrals whose absolute
 * vaue is <= the given cutoff. Indices for each integral are calculated and stored as well. The context is modified to
 * reflect the integrals processed.
 *
 * @param context The integral computation context to copy integrals from. context->num_processed is updated after processing.
 * @param nprocess The maximum number of integrals to process from the context. Fewer may be processed if the end of the
 *                 buffer is reached.
 * @param integrals The buffer into which to place the integrals, must be at least nprocess elements long.
 * @param indices The buffer into which to place the integral indices, must be at least nprocess elements long.
 * @param cutoff The integral cutoff value. Integrals with absolute value <= cutoff will not be retained. A negative value
 *               will retain all integrals.
 *
 * @return The number of integrals processed. A value less than nprocess means that the end of the integral buffer in the
 *         context has been reached.
 */
size_t SLIDE_process_1e_ints(context_t* context, const size_t nprocess, double* integrals, idx2_t* indices, const double cutoff);

size_t SLIDE_calc_eri(context_t* context, const double alpha, const double beta,
                      const shell_t* a, const shell_t* b, const shell_t* c, const shell_t* d);
size_t SLIDE_calc_ovi(context_t* context, const double alpha, const double beta,
                      const shell_t* a, const shell_t* b);
size_t SLIDE_calc_kei(context_t* context, const double alpha, const double beta,
                      const shell_t* a, const shell_t* b);
size_t SLIDE_calc_nai(context_t* context, const double alpha, const double beta,
                      const shell_t* a, const shell_t* b, const center_t** centers, const int ncenters);

#ifdef __cplusplus
}
#endif

#endif
