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

#ifndef _AQUARIUS_TENSOR_H_
#define _AQUARIUS_TENSOR_H_

#ifdef __cplusplus
#include <vector>
#endif

#include <stdint.h>
#include <stdio.h>

/*
 * Symmetry types
 */
#define TENSOR_NONSYMMETRIC     0
#define TENSOR_SYMMETRIC        1
#define TENSOR_ANTISYMMETRIC    2
#define TENSOR_SYMMETRIC_HOLLOW 3

/*
 * Shorthand notation for symmetry types
 */
#ifndef NS
#define NS  TENSOR_NONSYMMETRIC
#define SY  TENSOR_SYMMETRIC
#define AS  TENSOR_ANTISYMMETRIC
#define SH  TENSOR_SYMMETRIC_HOLLOW
#endif

/*
 * Return codes
 */
#define TENSOR_SUCCESS              0
#define TENSOR_OUT_OF_BOUNDS        -1
#define TENSOR_LENGTH_MISMATCH      -2
#define TENSOR_INDEX_MISMATCH       -3
#define TENSOR_INVALID_NDIM         -4
#define TENSOR_INVALID_LENGTH       -5
#define TENSOR_INVALID_LD           -6
#define TENSOR_LD_TOO_SMALL         -7
#define TENSOR_SYMMETRY_MISMATCH    -8
#define TENSOR_INVALID_SYMMETRY     -9
#define TENSOR_INVALID_START        -10

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * Main interface definitions
 *
 * In all cases, a repeated index label, either within the same tensor or in more than one tensor requires that the edge length be the same
 * in all labeled indices. With the exception of the resym operation and operations supporting it, the symmetry relation must also stay the same.
 * However, even in this case the symmetry relation may not change from AS or SH to SY or vice versa. 0-dimension tensors (scalars) are
 * allowed in all functions. In this case, the arguments len_*, sym_*, and idx_* may be NULL and will not be referenced. The data array for a scalar
 * must NOT be NULL, and should be of size 1. The special case of beta == +/-0.0 will overwrite special floating point values such as NaN and Inf
 * in the output.
 */

/**
 * Binary operations: C = alpha*A*B + beta*C
 */
typedef int (*tensor_func_binary)(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                                                      const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B,
                                  const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* sym_C, const int* idx_C);

/**
 * Multiply two tensors together and sum onto a third
 *
 * This form generalizes contraction and weighting with the unary operations trace, transpose, resym, diagonal, and replicate. Note that
 * the binary contraction operation is similar in form to the unary trace operation, while the binary weighting operation is similar in form to the
 * unary diagonal operation. Any combination of these operations may be performed. Even in the case that only a subset of the elements of C are written
 * to by the multiplication, all elements of C are first scaled by beta. Replication is performed in-place.
 */
int tensor_mult_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                                     const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B,
                 const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* sym_C, const int* idx_C);

/**
 * Contract two tensors into a third
 *
 * The general form for a contraction is ab...ef... * ef...cd... -> ab...cd... where the indices ef... will be summed over.
 * Indices may be transposed in any tensor, but the symmetry relations must not be changed, with the exception that ab... may or may not
 * be (anti)symmetrized with cd.... Any index group may be empty (in the case that ef... is empty, this reduces to an outer product).
 */
int tensor_contract_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                                         const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B,
                     const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* sym_C, const int* idx_C);

/**
 * Weight a tensor by a second and sum onto a third
 *
 * The general form for a weighting is ab...ef... * ef...cd... -> ab...cd...ef... with no indices being summed over.
 * Indices may be transposed in any tensor, but the symmetry relations must not be changed, with the addition that ab... must
 * be (anti)symmetrized with cd.... as ab... and cd... are with ef.... Any index group may be empty
 * (in the case that ef... is empty, this reduces to an outer product).
 */
int tensor_weight_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                                       const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B,
                   const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* sym_C, const int* idx_C);

/**
 * Sum the outer product of two tensors onto a third
 *
 * The general form for an outer product is ab... * cd... -> ab...cd... with no indices being summed over.
 * Indices may be transposed in any tensor, but the symmetry relations must not be changed, with the exception that ab... may or may not
 * be (anti)symmetrized with cd....
 */
int tensor_outer_prod_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                                           const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B,
                       const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* sym_C, const int* idx_C);

typedef int (*tensor_func_binary_dense)(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                                                            const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B,
                                        const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* idx_C);

int tensor_mult_dense_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                                           const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B,
                       const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* idx_C);

int tensor_contract_dense_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                                               const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B,
                           const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* idx_C);

int tensor_weight_dense_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                                             const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B,
                         const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* idx_C);

int tensor_outer_prod_dense_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                                                 const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B,
                             const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* idx_C);

/**
 * Unary operations: B = alpha*A + beta*B
 */
typedef int (*tensor_func_unary)(const double alpha, const double* A, const int ndim_A, const int* lda, const int* len_A, const int* sym_A, const int* idx_A,
                                 const double beta,        double* B, const int ndim_B, const int* ldb, const int* len_B, const int* sym_B, const int* idx_B);

/**
 * sum a tensor (presumably operated on in one or more ways) onto a second
 *
 * This form generalizes all of the unary operations trace, transpose, resym, diagonal, and replicate, which may be performed
 * in any combination. Even in the case that only a subset of the elements of B are written to by the operation, all elements
 * of B are first scaled by beta. Replication is performed in-place.
 */
int tensor_sum_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B);

/**
 * Transpose a tensor and sum onto a second
 *
 * The general form for a transposition operation is ab... -> P(ab...) where P is some permutation. Transposition may change
 * the order in which the elements of the tensor are physically stored, with the exception that (anti)symmetric
 * groups of indices must remain together. The symmetry relations of the indices must not be changed. A transposition among
 * symmetric indices produces no effect, while a transposition among antisymmetric indices may induce a change in sign,
 * but no physical reordering of elements.
 */
int tensor_transpose_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                      const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B);

/**
 * Change the symmetry relations of a tensor and sum onto a second
 *
 * The general form for a transposition operation is ab... -> ab... with a different set of symmetry relations in A and B (this means
 * also that A and B will not be the same size in general). Indices may either be unpacked or (anti)symmetrized by this operation, with
 * both possibilities allowed concurrently on disjoint sets of indices. A set of (anti)symmetric indices may be partially unpacked, leaving
 * two or more indices in packed storage, while one or more nonsymmetric indices may be (anti)symmetrized onto an existing packed set of
 * indices. Even in the case that only a subset of the elements of B are written to by the operation, all elements of B are first scaled
 * by beta. Transposition of the indices is allowed. Each (partial) unpacking includes a permutation factor of n'!/n! where n and n' are
 * number of (anti)symmetric indices before and after the unpacking.
 */
int tensor_resym_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                  const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B);

/**
 * Sum over (semi)diagonal elements of a tensor and sum onto a second
 *
 * The general form for a trace operation is ab...k*l*... -> ab... where k* denotes the index k appearing one or more times, etc. and where
 * the indices kl... will be summed (traced) over. Indices may be transposed (except for (anti)symmetric groups), and multiple appearances
 * of the traced indices kl... need not appear together. Either set of indices may be empty, with the special case that when no indices
 * are traced over, the result is the same as transpose. A trace over two or more antisymmetric indices produces no effect, except that B
 * is still scaled by beta.
 */
int tensor_trace_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                  const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B);

/**
 * Replicate a tensor and sum onto a second
 *
 * The general form for a replication operation is ab... -> ab...c*d*... where c* denotes the index c appearing one or more times.
 * Any indices may be transposed, with the exception that (anti)symmetric groups of indices must stay together. Replication is
 * performed in-place (meaning that no additional scratch space is required). Even in the case that only a subset of the elements of
 * B are written two, all elements of B are first scaled by beta.
 */
int tensor_replicate_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                      const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B);

/**
 * Extract the (semi)diagonal of a tensor and sum onto a second
 *
 * The general form for a diagonal extraction is ab...k*l*... -> ab...k*l*... where k* denotes the index k appearing one or more times
 * (which need not be the same number of times in A and B). Extracting the diagonal of a group of antisymmetric indices produces no effect,
 * except that B is still scaled by beta. Indices may be transposed, except (anti)symmetric groups of indices must stay together.
 */
int tensor_diagonal_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A,
                     const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const int* idx_B);

typedef int (*tensor_func_unary_dense)(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                                       const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B);

int tensor_sum_dense_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                      const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B);

int tensor_transpose_dense_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                            const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B);

int tensor_trace_dense_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                        const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B);

int tensor_replicate_dense_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                            const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B);

int tensor_diagonal_dense_(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A,
                           const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* idx_B);

/*
 * Helper function definitions
 */

int tensor_scale_(const double alpha, double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const int* idx_A);

int tensor_scale_dense_(const double alpha, double* A, const int ndim_A, const int* len_A, const int* lda, const int* idx_A);

int tensor_slice_dense(const double*  A, const int  ndim_A, const int* len_A, const int* lda,
                             double** B,       int* ndim_B,       int* len_B,       int* ldb,
                       const int* start, const int* len);

/**
 * Print a tensor in the form:
 *
 * i_0;k i_1;k ... i_ndim;k A_k
 *
 * where i_j;k are the index string for the kth element of A (A_k) in packed order, for k = 0...size-1.
 */
int tensor_print(FILE* fp, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A);

/**
 * Print a tensor in the form:
 *
 * i_0;k i_1;k ... i_ndim;k A_k
 *
 * where i_j;k are the index string for the kth element of A (A_k) in order, for k = 0...size-1.
 */
int tensor_print_dense(FILE* fp, const double* A, const int ndim_A, const int* len_A, const int* lda);

/**
 * Calculate the number of non-zero (and hence stored) elements in the given tensor.
 */
int64_t tensor_size(const int ndim, const int* len, const int* ld, const int* sym);

/**
 * Calculate the number of elements in the given tensor, assuming a dense (rectangular) layout.
 */
int64_t tensor_size_dense(const int ndim, const int* len, const int* ld);

/**
 * Combine each group of (anti)symmetric indices in the given tensor into a single index of the same total length.
 */
int tensor_densify(int* ndim, int* len, const int* sym);

/**
 * Take a packed (anti)symmetric tensor and expand it to a dense (rectangular) layout, with all redundant and
 * zero elements written explicitly. The result is scaled by a factor prod_k 1/[(l_k)^(n_k)] for k=0...ngroup, where
 * l_k is the length of the indices in the kth (anti)symmetric group of indices and n_k is the number of such indices.
 * This factor produces a normalized result when combined with tensor_symmetrize or tensor_pack.
 */
int tensor_unpack(const double* A, double* B, const int ndim_A, const int* len_A, const int* sym_A);

/**
 * Perform the reverse of tensor_unpack, referencing only non-redundant, non-zero elements of A (although which set of
 * redundant elements will be accessed is unspecified). The result is scaled by a factor prod_k (l_k)^(n_k) for
 *  k=0...ngroup, where l_k is the length of the indices in the kth (anti)symmetric group of indices and n_k is the
 *  number of such indices. The behavior is the same as tensor_symmetrize for the same input, except that fewer elements
 *  must be accessed.
 */
int tensor_pack(const double* A, double* B, const int ndim_A, const int* len_A, const int* sym_A);

/**
 * Perform an explicit (anti)symmetrization of the given dense tensor.
 */
int tensor_symmetrize(const double* A, double* B, const int ndim_A, const int* len_A, const int* sym_A);

#ifdef __cplusplus
}
#endif

/**
 * Macro interface to alloc for char* idx_? arrays
 */
#ifdef __cplusplus

#define CHAR_TO_INT_A \
std::vector<int> idx_A__(ndim_A); \
int* idx_A_ = idx_A__.data(); \
for (int i__ = 0;i__ < ndim_A;i__++) idx_A_[i__] = (int)idx_A[i__];

#define CHAR_TO_INT_AB \
CHAR_TO_INT_A \
std::vector<int> idx_B__(ndim_B); \
int* idx_B_ = idx_B__.data(); \
for (int i__ = 0;i__ < ndim_B;i__++) idx_B_[i__] = (int)idx_B[i__];

#define CHAR_TO_INT_ABC \
CHAR_TO_INT_AB \
std::vector<int> idx_C__(ndim_C); \
int* idx_C_ = idx_C__.data(); \
for (int i__ = 0;i__ < ndim_C;i__++) idx_C_[i__] = (int)idx_C[i__];

#else

#define CHAR_TO_INT_A \
int idx_A_[ndim_A]; \
for (int i__ = 0;i__ < ndim_A;i__++) idx_A_[i__] = (int)idx_A[i__];

#define CHAR_TO_INT_AB \
CHAR_TO_INT_A \
int idx_B_[ndim_B]; \
for (int i__ = 0;i__ < ndim_B;i__++) idx_B_[i__] = (int)idx_B[i__];

#define CHAR_TO_INT_ABC \
CHAR_TO_INT_AB \
int idx_C_[ndim_C]; \
for (int i__ = 0;i__ < ndim_C;i__++) idx_C_[i__] = (int)idx_C[i__];

#endif

static inline int tensor_mult(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const char* idx_A,
                                                  const double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const char* idx_B,
                              const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const int* sym_C, const char* idx_C)
{
    CHAR_TO_INT_ABC
    return tensor_mult_(alpha, A, ndim_A, len_A, lda, sym_A, idx_A_,
                               B, ndim_B, len_B, ldb, sym_B, idx_B_,
                        beta,  C, ndim_C, len_C, ldc, sym_C, idx_C_);
}

static inline int tensor_mult_dense(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const char* idx_A,
                                                        const double* B, const int ndim_B, const int* len_B, const int* ldb, const char* idx_B,
                                    const double beta,        double* C, const int ndim_C, const int* len_C, const int* ldc, const char* idx_C)
{
    CHAR_TO_INT_ABC
    return tensor_mult_dense_(alpha, A, ndim_A, len_A, lda, idx_A_,
                                     B, ndim_B, len_B, ldb, idx_B_,
                              beta,  C, ndim_C, len_C, ldc, idx_C_);
}

static inline int tensor_sum(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const char* idx_A,
                             const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const int* sym_B, const char* idx_B)
{
    CHAR_TO_INT_AB
    return tensor_sum_(alpha, A, ndim_A, len_A, lda, sym_A, idx_A_,
                       beta,  B, ndim_B, len_B, ldb, sym_B, idx_B_);
}

static inline int tensor_sum_dense(const double alpha, const double* A, const int ndim_A, const int* len_A, const int* lda, const char* idx_A,
                                   const double beta,        double* B, const int ndim_B, const int* len_B, const int* ldb, const char* idx_B)
{
    CHAR_TO_INT_AB
    return tensor_sum_dense_(alpha, A, ndim_A, len_A, lda, idx_A_,
                             beta,  B, ndim_B, len_B, ldb, idx_B_);
}

static inline int tensor_scale(const double alpha, double* A, const int ndim_A, const int* len_A, const int* lda, const int* sym_A, const char* idx_A)
{
    CHAR_TO_INT_A
    return tensor_scale_(alpha, A, ndim_A, len_A, lda, sym_A, idx_A_);
}

static inline int tensor_scale_dense(const double alpha, double* A, const int ndim_A, const int* len_A, const int* lda, const char* idx_A)
{
    CHAR_TO_INT_A
    return tensor_scale_dense_(alpha, A, ndim_A, len_A, lda, idx_A_);
}

#endif
