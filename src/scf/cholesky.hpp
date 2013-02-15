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

#ifndef _AQUARIUS_SCF_CHOLESKY_HPP_
#define _AQUARIUS_SCF_CHOLESKY_HPP_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <vector>

#include "mpi.h"

#include "tensor.hpp"
#include "dist_tensor.hpp"

#include "slide/slide.hpp"
#include "input/molecule.hpp"
#include "input/config.hpp"
#include "util/util.h"
#include "util/blas.h"

namespace aquarius
{
namespace scf
{

class CholeskyIntegrals
{
    protected:
        enum Status {TODO, ACTIVE, DONE};
        struct diag_elem_t
        {
            double elem;
            int shelli;
            int shellj;
            int func;
            int contr;
            int idx;
            Status status;
        };

        int rank;
        const input::Molecule& molecule;
        std::vector<slide::Shell> shells;
        double delta;
        double cond;
        DistTensor* L;
        DistTensor* D;
        int ndiag;
        int nfunc;
        int nblock;
        slide::Context* context;
        DistWorld* dw;
        MPI::Intracomm comm;

        void decompose();

        void resortBlock(const int block_size, double* L, diag_elem_t* diag, double* tmp);

        int collectActiveRows(const int block_size, const double* L, diag_elem_t* diag,
                              double* L_active, diag_elem_t* diag_active);

        bool isBlockConverged(const diag_elem_t* diag, int block_size, double& max_elem);

        void sortBlocks(diag_elem_t* diag, int* block_start, int* block_size);

        int getDiagonalBlock(const slide::Shell& a, const slide::Shell& b, diag_elem_t* diag);

        void decomposeBlock(int block_size, double* L, double* D, diag_elem_t* diag);

        void updateBlock(int old_rank, int block_size_i, double* L_i, diag_elem_t* diag_i,
                         int block_size_j, double* L_j, diag_elem_t* diag_j, double* D);

        double testBlock(const DenseTensor& block, const slide::Shell& a, const slide::Shell& b,
                         const slide::Shell& c, const slide::Shell& d);

        static bool sort_by_integral(const diag_elem_t& e1, const diag_elem_t& e2) { return e1.elem > e2.elem; }

        //CholeskyIntegrals(const input::Config& config, const input::Molecule& molecule);

    public:
        CholeskyIntegrals(DistWorld* dw, const input::Config& config, const input::Molecule& molecule);

        ~CholeskyIntegrals();

        const input::Molecule& getMolecule() const { return molecule; }

        int getRank() const { return rank; }

        const DistTensor& getL() const { return *L; }

        const DistTensor& getD() const { return *D; }

        void test();
};

}
}

#endif
