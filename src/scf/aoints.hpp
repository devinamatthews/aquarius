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

#ifndef _AQUARIUS_SCF_AOINTS_HPP_
#define _AQUARIUS_SCF_AOINTS_HPP_

#include <vector>
#include <cassert>

#include "mpi.h"

#include "util/distributed.hpp"
#include "slide/slide.hpp"
#include "input/molecule.hpp"
#include "memory/memory.h"

#define TMP_BUFSIZE 100
#define INTEGRAL_CUTOFF 1e-14

namespace aquarius
{
namespace scf
{

template <typename T>
struct integral_t
{
    static MPI::Datatype mpi_type;
    T value;
    idx4_t idx;
    integral_t(T value, idx4_t idx)
    : value(value), idx(idx)
    {
        if (mpi_type == MPI::DATATYPE_NULL)
        {
            assert(sizeof(short) == sizeof(uint16_t));
            MPI::Datatype _types[] = {MPI_TYPE_<T>::value(), MPI::SHORT};
            int _counts[] = {1, 4};
            MPI::Aint _offsets[] = {0, 8};
            mpi_type = MPI::Datatype::Create_struct(2, _counts, _offsets, _types);
        }
    }
};

template <typename T>
class AOIntegrals : public Distributed<T>
{
    protected:
        slide::Context context;
        size_t num_ints;
        integral_t<T> *ints;
        const input::Molecule& molecule;

    public:
        AOIntegrals(tCTF_World<T>& ctf, const input::Molecule& molecule)
        : Distributed<T>(ctf), molecule(molecule)
        {
            generateInts();
            loadBalance();
        }

    protected:
        void generateInts()
        {
            double tmpval[TMP_BUFSIZE];
            idx4_t tmpidx[TMP_BUFSIZE];

            std::vector< integral_t<T> > tmpints;

            int abcd = 0;
            for (input::Molecule::const_shell_iterator a = molecule.getShellsBegin();a != molecule.getShellsEnd();++a)
            {
                for (input::Molecule::const_shell_iterator b = molecule.getShellsBegin();b != a+1;++b)
                {
                    for (input::Molecule::const_shell_iterator c = molecule.getShellsBegin();c != a+1;++c)
                    {
                        input::Molecule::const_shell_iterator dmax = c+1;
                        if (a == c) dmax = b+1;
                        for (input::Molecule::const_shell_iterator d = molecule.getShellsBegin();d != dmax;++d)
                        {
                            if (abcd%this->nproc == this->rank)
                            {
                                context.calcERI(1.0, 0.0, *a, *b, *c, *d);

                                size_t n;
                                while ((n = context.process2eInts(TMP_BUFSIZE, tmpval, tmpidx, INTEGRAL_CUTOFF)) != 0)
                                {
                                    for (size_t i = 0;i < n;i++)
                                    {
                                        //printf("%d %d %d %d %20.15f\n", tmpidx[i].i, tmpidx[i].j, tmpidx[i].k, tmpidx[i].l, tmpval[i]);
                                        tmpints.push_back(integral_t<T>(tmpval[i], tmpidx[i]));
                                    }
                                }
                            }
                            abcd++;
                        }
                    }
                }
            }

            num_ints = tmpints.size();
            ints = SAFE_MALLOC(integral_t<T>, num_ints);
            std::copy(tmpints.begin(), tmpints.end(), ints);
        }

        void loadBalance()
        {

        }

        void canonicalize()
        {
            for (int i = 0;i < num_ints;++i)
            {
                if (ints[i].idx.i > ints[i].idx.j) std::swap(ints[i].idx.i, ints[i].idx.j);
                if (ints[i].idx.k > ints[i].idx.l) std::swap(ints[i].idx.k, ints[i].idx.l);
                if (ints[i].idx.i > ints[i].idx.k ||
                    (ints[i].idx.i == ints[i].idx.j && ints[i].idx.j > ints[i].idx.l))
                {
                    std::swap(ints[i].idx.i, ints[i].idx.k);
                    std::swap(ints[i].idx.j, ints[i].idx.l);
                }
            }
        }
};

}
}

#endif
