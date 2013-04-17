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

#include "util/distributed.hpp"
#include "slide/slide.hpp"
#include "input/molecule.hpp"
#include "memory/memory.h"

#include <boost/foreach.hpp>

#define TMP_BUFSIZE 100
#define INTEGRAL_CUTOFF 1e-14

namespace aquarius
{
namespace scf
{

template <typename T>
class AOIntegrals : public Distributed<T>
{
    public:
        const input::Molecule& molecule;

    protected:
        slide::Context context;
        size_t nints;
        T *ints;
        idx4_t *idxs;

    public:
        AOIntegrals(tCTF_World<T>& ctf, const input::Molecule& molecule)
        : Distributed<T>(ctf), molecule(molecule)
        {
            generateInts();
            loadBalance();
        }

        size_t getNumInts() const { return nints; }

        const T* getInts() const { return ints; }

        const idx4_t* getIndices() const { return idxs; }

    protected:
        void generateInts()
        {
            double tmpval[TMP_BUFSIZE];
            idx4_t tmpidx[TMP_BUFSIZE];

            std::vector<T> tmpints;
            std::vector<idx4_t> tmpidxs;

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
            std::copy(tmpints.begin(), tmpints.end(), ints);
            idxs = SAFE_MALLOC(idx4_t, nints);
            std::copy(tmpidxs.begin(), tmpidxs.end(), idxs);
        }

        void loadBalance()
        {
            //TODO
        }

        void canonicalize()
        {
            for (int i = 0;i < nints;++i)
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
