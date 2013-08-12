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

#include "kei.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::op;

template <typename U>
KEI<U>::KEI(Arena<U>& arena, const Context& ctx, const Molecule& molecule)
: AOOperator<U>(arena, molecule.getNumOrbitals()),
  T(arena, 2, vec(molecule.getNumOrbitals(), molecule.getNumOrbitals()), vec(NS,NS), true)
{
    vector<vector<int> > idx = Shell::setupIndices(ctx, molecule);
    vector<Shell> shells(molecule.getShellsBegin(), molecule.getShellsEnd());
    vector<tkv_pair<U> > pairs;
    vector<Center> centers;

    for (vector<Atom>::const_iterator i = molecule.getAtomsBegin();i != molecule.getAtomsEnd();++i)
    {
        centers.push_back(i->getCenter());
    }

    int block = 0;
    for (int i = 0;i < shells.size();++i)
    {
        for (int j = 0;j <= i;++j)
        {
            if (block%this->nproc == this->rank)
            {
                OneElectronIntegrals kei(shells[i], shells[j], KEIEvaluator());

                size_t nint = kei.getNumInts();
                vector<U> ints(nint);
                vector<idx2_t> idxs(nint);

                size_t nproc = kei.process(ctx, idx[i], idx[j], nint, ints.data(), idxs.data());
                for (int i = 0;i < nproc;i++)
                {
                    //printf("%d %d %25.15f\n", idxs[i].i+1, idxs[i].j+1, ints[i]);

                    pairs.push_back(tkv_pair<U>(idxs[i].i*molecule.getNumOrbitals()+idxs[i].j, ints[i]));
                    if (idxs[i].i != idxs[i].j)
                    {
                        pairs.push_back(tkv_pair<U>(idxs[i].j*molecule.getNumOrbitals()+idxs[i].i, ints[i]));
                    }
                }
            }

            block++;
        }
    }

    T.writeRemoteData(pairs);
}

INSTANTIATE_SPECIALIZATIONS(KEI);
