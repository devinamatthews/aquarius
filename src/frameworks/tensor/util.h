#ifndef _AQUARIUS_TENSOR_UTIL_H_
#define _AQUARIUS_TENSOR_UTIL_H_

#define EQ          SH
#define CLUSTER_    0x1
#define UNIT_       0x2

#define LOOKUP3(what,where) ((where) >= ndim_A+ndim_B ? what##_C[(where)-ndim_A-ndim_B]                  \
                                                      : ((where) >= ndim_A ? what##_B[(where)-ndim_A]    \
                                                                           : what##_A[where]))

#define LOOKUP2(what,where) ((where) >= ndim_A ? what##_B[(where)-ndim_A]    \
                                               : what##_A[where])

#define VALIDATE_TENSOR(ndim,len,ld,sym)            \
do                                                  \
{                                                   \
    int ret__ = validate_tensor(ndim,len,ld,sym);   \
    if (ret__ != TENSOR_SUCCESS) return ret__;      \
} while (0)

#ifdef __cplusplus
extern "C"
{
#endif

int tensor_sign(const int ndim, const int* sym, const int* idx1, const int* idx2);

void tensor_info(const int ndim, const int* len, const int* ld, const int* sym,
                 int* group, int* skip, int* stride, size_t* size);

int tensor_iterate(const double alpha, const int ntensor, double * const restrict * restrict data, const size_t* restrict size,
                   const int ndim, const int* restrict len, const int* restrict sym,
                   const int* restrict which, const int* restrict stride, const int* restrict pos,
                   const int* restrict skip);

void dzero(const int n, double* a, const int inca);

void izero(const int n, int* a, const int inca);

int idot(const int n, const int *x, const int incx, const int *y, const int incy);

uint64_t factorial(const int n);

uint64_t binomial(const int a, const int b);

size_t offset_to_next_index(const int current_index, const int group_pos, const int sym, const size_t stride);

size_t offset_between_indices(const int index1, const int index2, const int group_pos, const int sym, const size_t stride);

size_t offset_to_index(const int index, const int group_pos, const int sym, const size_t stride);

int reverse_sequence(const int n, int* s);

void first_packed_indices(const int ndim, const int* len, const int* sym, int* idx);

bool next_dense_indices(const int ndim, const int* len, int* idx);

bool next_packed_indices(const int ndim, const int* len, const int* sym, int* idx);

int relative_sign(const int n, const int* s1, const int* s2);

int next_permutation(const int ndim, const int* sym, int* idx);

void indices_from_labels(const int n, const int* labels, int* used_labels, int* indices);

void reorder_from_indices(const int n, const int* idx1, const int* idx2, int* reorder);

int validate_tensor(const int ndim, const int* len, const int* ld, const int* sym);

#ifndef __cplusplus

void index_connectivity(const int ndim_A, const int* sym_A, const int* idx_A,
                        const int ndim_B, const int* sym_B, const int* idx_B,
                        bool G[ndim_A+ndim_B][ndim_A+ndim_B]);

int connected_components(const int n, const bool G[n][n], int cc[n], int len[n]);

int connected_component(const int n, const bool G[n][n], const int i, bool seen[n], int** cc);

#else
}
#endif

#endif
