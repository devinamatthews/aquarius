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

#include "ovi.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::input;
using namespace aquarius::slide;
using namespace aquarius::op;

template <typename T>
OVI<T>::OVI(Arena<T>& arena, const Molecule& molecule)
: AOOperator<T>(arena, molecule.getNumOrbitals()),
  S(arena, 2, vec(molecule.getNumOrbitals(), molecule.getNumOrbitals()), vec(NS,NS), true)
{
    Context context;

    vector< tkv_pair<T> > pairs;

    int block = 0;
    for (Molecule::const_shell_iterator i = molecule.getShellsBegin();i != molecule.getShellsEnd();++i)
    {
        for (Molecule::const_shell_iterator j = molecule.getShellsBegin();j != molecule.getShellsEnd();++j)
        {
            if (i < j) continue;

            if (block%this->nproc == this->rank)
            {
                context.calcOVI(1.0, 0.0, *i, *j);

                size_t nint = context.getNumIntegrals();
                vector<T> ints(nint);
                vector<idx2_t> idxs(nint);

                size_t nproc = context.process1eInts(nint, ints.data(), idxs.data(), -1.0);
                for (int i = 0;i < nproc;i++)
                {
                    //printf("%d %d %25.15f\n", idxs[i].i+1, idxs[i].j+1, ints[i]);

                    pairs.push_back(tkv_pair<T>(idxs[i].i*this->norb+idxs[i].j, ints[i]));
                    if (idxs[i].i != idxs[i].j)
                    {
                        pairs.push_back(tkv_pair<T>(idxs[i].j*this->norb+idxs[i].i, ints[i]));
                    }
                }
            }

            block++;
        }
    }

    S.writeRemoteData(pairs);
}

INSTANTIATE_SPECIALIZATIONS(OVI);
