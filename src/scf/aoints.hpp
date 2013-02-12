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

#ifndef _AQUARIUS_SCF_AOINTS_HPP_
#define _AQUARIUS_SCF_AOINTS_HPP_

#include <vector>

#include "mpi.h"

#include "slide/slide.hpp"
#include "input/molecule.hpp"

namespace aquarius
{
namespace scf
{

struct integral_t
{
    static MPI::Datatype mpi_type;
    double value;
    idx4_t idx;
    integral_t(double value, idx4_t idx)
    : value(value), idx(idx) {}
};

class AOIntegrals
{
    protected:
        slide::Context context;
        size_t num_ints;
        integral_t *ints;
        MPI::Intracomm comm;
        const input::Molecule& molecule;

        void generateInts();

        void loadBalance();

    public:
        AOIntegrals(const input::Molecule& molecule, const MPI::Comm& comm)
        : molecule(molecule), comm(comm)
        {
            generateInts();
            loadBalance();
        }

        void canonicalize();

        integral_t* getInts() { return ints; }

        const integral_t* getInts() const { return ints; }

        size_t getNumInts() const { return num_ints; }

        const input::Molecule& getMolecule() const { return molecule; }
};

}
}

#endif
