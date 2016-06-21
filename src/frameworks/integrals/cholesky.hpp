#ifndef _AQUARIUS_FRAMEWORKS_INTEGRALS_CHOLESKY_HPP_
#define _AQUARIUS_FRAMEWORKS_INTEGRALS_CHOLESKY_HPP_

#include "frameworks/util.hpp"
#include "frameworks/molecule.hpp"
#include "frameworks/tensor.hpp"
#include "frameworks/task.hpp"

#include "context.hpp"

namespace aquarius
{
namespace integrals
{

class CholeskyIntegrals
{
    public:
        const molecule::Molecule& molecule;

    protected:
        enum Status {TODO, ACTIVE, DONE};
        struct diag_elem_t
        {
            double elem;
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
        vector<molecule::Shell> shells;
        double delta;
        double cond;
        tensor::Tensor<PGSYMMETRIC> L;
        tensor::Tensor<PGSYMMETRIC> D;
        int ndiag;
        int nfunc;
        int nblock;

    public:
        CholeskyIntegrals(const Context& ctx, task::Config& config, const molecule::Molecule& molecule);

        void test();

        int getRank() const { return nvec; }

        const tensor::Tensor<PGSYMMETRIC>& getL() const { return L; }

        const tensor::Tensor<PGSYMMETRIC>& getD() const { return D; }

    protected:

        void decompose();

        void resortBlock(const int block_size, double* L, diag_elem_t* diag, double* tmp);

        int collectActiveRows(const int block_size, const double* L, diag_elem_t* diag,
                              double* L_active, diag_elem_t* diag_active);

        bool isBlockConverged(const diag_elem_t* diag, int block_size, double& max_elem);

        void sortBlocks(diag_elem_t* diag, int* block_start, int* block_size);

        void getShellOffsets(const molecule::Shell& a, const molecule::Shell& b,
                             const molecule::Shell& c, const molecule::Shell& d,
                             size_t& controffa, size_t& funcoffa, size_t& controffb, size_t& funcoffb,
                             size_t& controffc, size_t& funcoffc, size_t& controffd, size_t& funcoffd);

        int getDiagonalBlock(const molecule::Shell& a, const molecule::Shell& b,
                             diag_elem_t* diag);

        /*
         * perform a modified Cholesky decomposition (LDL^T) on a matrix subblock
         *
         * block_size       - size of the block to update
         * L                - Cholesky vectors for the block, dimensions L[block_size][>rank]
         * D                - the diagonal factor of the decomposition, dimensions D[>rank]
         * diag             - diagonal elements of the residual matrix for the block and other accounting information,
         *                    dimensions diag[block_size]
         *
         * return value     - rank after update
         */
        void decomposeBlock(int block_size, double* L_, double* D, diag_elem_t* diag);

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
        void updateBlock(int old_rank, int block_size_i, double* L_i_, diag_elem_t* diag_i,
                         int block_size_j, double* L_j_, diag_elem_t* diag_j, double* D);

        double testBlock(const tensor::Tensor<BOUNDED>& block,
                         const molecule::Shell& a, const molecule::Shell& b,
                         const molecule::Shell& c, const molecule::Shell& d);
};

}
}

#endif
