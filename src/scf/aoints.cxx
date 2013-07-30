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

#include "aoints.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::scf;
using namespace aquarius::input;
using namespace aquarius::slide;

template <typename T>
AOIntegrals<T>::AOIntegrals(tCTF_World<T>& ctf, const Molecule& molecule)
: Distributed<T>(ctf), molecule(molecule)
{
    generateInts();
    loadBalance();
}

template <typename T>
void AOIntegrals<T>::generateInts()
{
    double tmpval[TMP_BUFSIZE];
    idx4_t tmpidx[TMP_BUFSIZE];

    vector<T> tmpints;
    vector<idx4_t> tmpidxs;

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
                    if (abcd%this->nproc == this->rank)
                    {
                        context.calcERI(1.0, 0.0, *a, *b, *c, *d);

                        size_t n;
                        while ((n = context.process2eInts(TMP_BUFSIZE, tmpval, tmpidx, INTEGRAL_CUTOFF)) != 0)
                        {
                            tmpints.insert(tmpints.end(), tmpval, tmpval+n);
                            tmpidxs.insert(tmpidxs.end(), tmpidx, tmpidx+n);
                        }
                    }
                    abcd++;
                }
            }
        }
    }

    nints = tmpints.size();
    ints = SAFE_MALLOC(T, nints);
    copy(tmpints.begin(), tmpints.end(), ints);
    idxs = SAFE_MALLOC(idx4_t, nints);
    copy(tmpidxs.begin(), tmpidxs.end(), idxs);
}

template <typename T>
void AOIntegrals<T>::loadBalance()
{
    //TODO
}

template <typename T>
void AOIntegrals<T>::canonicalize()
{
    for (int i = 0;i < nints;++i)
    {
        if (idxs[i].i > idxs[i].j) swap(idxs[i].i, idxs[i].j);
        if (idxs[i].k > idxs[i].l) swap(idxs[i].k, idxs[i].l);
        if (idxs[i].i > idxs[i].k ||
            (idxs[i].i == idxs[i].j && idxs[i].j > idxs[i].l))
        {
            swap(idxs[i].i, idxs[i].k);
            swap(idxs[i].j, idxs[i].l);
        }
    }
}

INSTANTIATE_SPECIALIZATIONS(AOIntegrals);
