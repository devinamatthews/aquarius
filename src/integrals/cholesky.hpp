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

#ifndef _AQUARIUS_INTEGRALS_CHOLESKY_HPP_
#define _AQUARIUS_INTEGRALS_CHOLESKY_HPP_

#include "2eints.hpp"

#include "tensor/symblocked_tensor.hpp"
#include "tensor/dense_tensor.hpp"
#include "input/molecule.hpp"
#include "input/config.hpp"
#include "util/util.h"
#include "util/blas.h"
#include "task/task.hpp"

namespace aquarius
{
namespace integrals
{

template <typename T>
class CholeskyIntegrals : public task::Resource
{
    public:
        const input::Molecule& molecule;

    protected:
        enum Status {TODO, ACTIVE, DONE};
        struct diag_elem_t
        {
            T elem;
            int shelli;
            int shellj;
            int funci;
            int funcj;
            int contri;
            int contrj;
            int idx;
            Status status;

            bool operator<(const diag_elem_t& other) const
            {
                return elem < other.elem;
            }
        };

        const Context& ctx;
        int nvec;
        std::vector<Shell> shells;
        T delta;
        T cond;
        tensor::SymmetryBlockedTensor<T>* L;
        tensor::SymmetryBlockedTensor<T>* D;
        int ndiag;
        int nfunc;
        int nblock;

    public:
        CholeskyIntegrals(const Arena& arena, const Context& ctx, const input::Config& config, const input::Molecule& molecule);

        ~CholeskyIntegrals();

        void test();

        int getRank() const { return nvec; }

        const tensor::SymmetryBlockedTensor<T>& getL() const { return *L; }

        const tensor::SymmetryBlockedTensor<T>& getD() const { return *D; }

    protected:

        void decompose();

        void resortBlock(const int block_size, T* L, diag_elem_t* diag, T* tmp);

        int collectActiveRows(const int block_size, const T* L, diag_elem_t* diag,
                              T* L_active, diag_elem_t* diag_active);

        bool isBlockConverged(const diag_elem_t* diag, int block_size, T& max_elem);

        void sortBlocks(diag_elem_t* diag, int* block_start, int* block_size);

        void getShellOffsets(const Shell& a, const Shell& b, const Shell& c, const Shell& d,
                             size_t& controffa, size_t& funcoffa, size_t& controffb, size_t& funcoffb,
                             size_t& controffc, size_t& funcoffc, size_t& controffd, size_t& funcoffd);

        int getDiagonalBlock(const Shell& a, const Shell& b, diag_elem_t* diag);

        /*
         * perform a modified Cholesky decomposition (LDL^T) on a matrix subblock
         *
         * block_size       - size of the block to update
         * L                - Cholesky std::vectors for the block, dimensions L[block_size][>rank]
         * D                - the diagonal factor of the decomposition, dimensions D[>rank]
         * diag             - diagonal elements of the residual matrix for the block and other accounting information,
         *                    dimensions diag[block_size]
         *
         * return value     - rank after update
         */
        void decomposeBlock(int block_size, T* L_, T* D, diag_elem_t* diag);

        /*
         * update a block of Cholesky vectors
         *
         * rank             - first rank to update
         * rank             - rank after update
         * block_size_i     - size of the block of Cholesky vectors to update
         * L_i              - Cholesky vectors to update, dimensions L_i[block_size_i][rank]
         * diag_i           - diagonal elements of the residual matrix and other accounting information
         *                    for the block to be updated, dimensions diag_i[block_size_i]
         * block_size_j     - size of the block containing the rows leading to this update
         * L_j              - Cholesky vectors for the rows leading to the update, dimensions L_j[rank-rank][rank]
         *                    (only those vectors actually needed, and not the entire block)
         * diag_j           - diagonal elements of the residual matrix for the block leading to the update,
         *                    dimensions diag_j[block_size_j] (the full block)
         * D                - the diagonal factor of the decomposition, dimensions D[rank]
         */
        void updateBlock(int old_rank, int block_size_i, T* L_i_, diag_elem_t* diag_i,
                         int block_size_j, T* L_j_, diag_elem_t* diag_j, T* D);

        T testBlock(const tensor::DenseTensor<T>& block, const Shell& a, const Shell& b,
                                                         const Shell& c, const Shell& d);
};

}
}

#endif
