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

#include "multipole.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::tensor;
using namespace aquarius::integrals;
using namespace aquarius::scf;
using namespace aquarius::input;

template <typename T>
Multipole<T>::Multipole(const UHF<T>& uhf, int Lmin_, int Lmax_)
: MOOperator<T>(uhf), CompositeTensor<Multipole<T>,
  OneElectronOperator<T>,T>(Lmax_ == -1 ? (Lmin_+1)*(Lmin_+2)/2 :
          (Lmax_+1)*(Lmax_+2)*(Lmax_+3)/6-Lmin_*(Lmin_+1)*(Lmin_+2)/6),
  Lmin(Lmin_), Lmax(Lmax_ == -1 ? Lmin_ : Lmax_)
{
    Context context;
    const Molecule& m = uhf.getMolecule();
    int N = m.getNumOrbitals();

    DistTensor<T> ao(uhf.arena, 2, vec(N,N), vec(NS,NS), false);

    int xyztot = 0;
    for (int L = Lmin;L <= Lmax;L++)
    {
        vector< vector< tkv_pair<T> > > pairs;

        int ij = 0;
        for (Molecule::const_shell_iterator i = m.getShellsBegin();i != m.getShellsEnd();++i)
        {
            for (Molecule::const_shell_iterator j = i;j != m.getShellsEnd();++j)
            {
                if (i < j) continue;

                if (ij%this->nproc == this->rank)
                {
                    //calc integrals

                    //put into pairs
                }

                ij++;
            }
        }

        for (int xyz = 0;xyz < (L+1)*(L+2)/2;xyz++)
        {
            ao = (T)0;
            ao.writeRemoteData(pairs[xyz]);
            this->tensors[xyztot++].tensor = new OneElectronOperator<T>(uhf, ao);
        }
    }
}

template <typename T>
const OneElectronOperator<T>& Multipole<T>::operator()(int L, int xyz) const
{
    assert(L >= Lmin && L <= Lmax);
    assert(xyz >= 0 && xyz < (L+1)*(L+2)/2);
    return (*this)(L*(L+1)*(L+2)/6-Lmin*(Lmin+1)*(Lmin+2)/6+xyz);
}

template <typename T>
const OneElectronOperator<T>& Multipole<T>::operator()(int x, int y, int z) const
{
    int L = x+y+z;
    assert(L >= Lmin && L <= Lmax);
    int xyz = z+(L+1-x)*(L+2-x)/2;
    return (*this)(L, xyz);
}
