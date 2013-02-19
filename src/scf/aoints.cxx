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

#include <algorithm>

#include "memory/memory.h"

#include "aoints.hpp"

using namespace std;
using namespace MPI;
using namespace aquarius::slide;
using namespace aquarius::input;

namespace aquarius
{
namespace scf
{

#define TMP_BUFSIZE 100
#define INTEGRAL_CUTOFF 1e-14

Datatype integral_t::mpi_type = DATATYPE_NULL;

AOIntegrals::AOIntegrals(tCTF_World<double>& ctf, const Molecule& molecule)
: Distributed<double>(ctf), molecule(molecule)
{
    Datatype _types[] = {DOUBLE, SHORT};
    int _counts[] = {1, 4};
    Aint _offsets[] = {0, 8};
    integral_t::mpi_type = Datatype::Create_struct(2, _counts, _offsets, _types);
    generateInts();
    loadBalance();
}

void AOIntegrals::generateInts()
{
    double tmpval[TMP_BUFSIZE];
    idx4_t tmpidx[TMP_BUFSIZE];

    vector<integral_t> tmpints;

    int abcd = 0;
    for (Molecule::const_shell_iterator a = molecule.getShellsBegin();a != molecule.getShellsEnd();++a)
    {
        for (Molecule::const_shell_iterator b = molecule.getShellsBegin();b != a+1;++b)
        {
            for (Molecule::const_shell_iterator c = molecule.getShellsBegin();c != a+1;++c)
            {
                Molecule::const_shell_iterator dmax = c+1;
                if (a == c) dmax = b+1;
                for (Molecule::const_shell_iterator d = molecule.getShellsBegin();d != dmax;++d)
                {
                    if (abcd%comm.Get_size() == comm.Get_rank())
                    {
                        context.calcERI(1.0, 0.0, *a, *b, *c, *d);

                        size_t n;
                        while ((n = context.process2eInts(TMP_BUFSIZE, tmpval, tmpidx, INTEGRAL_CUTOFF)) != 0)
                        {
                            for (size_t i = 0;i < n;i++)
                            {
                                //printf("%d %d %d %d %20.15f\n", tmpidx[i].i, tmpidx[i].j, tmpidx[i].k, tmpidx[i].l, tmpval[i]);
                                tmpints.push_back(integral_t(tmpval[i], tmpidx[i]));
                            }
                        }
                    }
                    abcd++;
                }
            }
        }
    }

    num_ints = tmpints.size();
    ints = SAFE_MALLOC(integral_t, num_ints);
    copy(tmpints.begin(), tmpints.end(), ints);
}

void AOIntegrals::loadBalance()
{

}

void AOIntegrals::canonicalize()
{
    for (int i = 0;i < num_ints;++i)
    {
        if (ints[i].idx.i > ints[i].idx.j) swap(ints[i].idx.i, ints[i].idx.j);
        if (ints[i].idx.k > ints[i].idx.l) swap(ints[i].idx.k, ints[i].idx.l);
        if (ints[i].idx.i > ints[i].idx.k ||
            (ints[i].idx.i == ints[i].idx.j && ints[i].idx.j > ints[i].idx.l))
        {
            swap(ints[i].idx.i, ints[i].idx.k);
            swap(ints[i].idx.j, ints[i].idx.l);
        }
    }
}

}
}
