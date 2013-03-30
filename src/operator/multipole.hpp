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

#ifndef _AQUARIUS_OPERATOR_MULTIPOLE_HPP_
#define _AQUARIUS_OPERATOR_MULTIPOLE_HPP_

#include "slide/slide.hpp"
#include "scf/scf.hpp"

#include "1eoperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class Multipole : public MOOperator<T>, public tensor::CompositeTensor<Multipole<T>,OneElectronOperator<T>,T>
{
    protected:
        const int Lmin, Lmax;

    public:
        Multipole(const scf::UHF<T>& uhf, const int Lmin, const int Lmax=Lmin)
        : MOOperator<T>(uhf), tensor::CompositeTensor<Multipole<T>,
              OneElectronOperator<T>,T>((Lmax+1)*(Lmax+2)*(Lmax+3)/6-Lmin*(Lmin+1)*(Lmin+2)/6)
        {
            slide::Context context;
            const input::Molecule& m = uhf.getMolecule();
            int N = m.getNumOrbitals();

            int sizeNN[] = {N, N};
            int shapeNN[] = {NS, NS};

            int xyztot = 0;
            for (int L = Lmin;L <= Lmax;L++)
            {
                std::vector< tensor::DistTensor<T>* > ao;
                std::vector< std::vector< tkv_pair<T> > > pairs;

                for (int xyz = 0;xyz < (L+1)*(L+2)/2;xyz++)
                {
                    ao[xyz] = new tensor::DistTensor<T>(ctf, 2, sizeNN, shapeNN, true);
                }

                int ij = 0;
                for (input::Molecule::const_shell_iterator i = m.getShellsBegin();i != m.getShellsEnd();++i)
                {
                    for (input::Molecule::const_shell_iterator j = i;j != m.getShellsEnd();++j)
                    {
                        //calc integrals

                        //put into pairs
                    }
                }

                for (int xyz = 0;xyz < (L+1)*(L+2)/2;xyz++)
                {
                    ao[xyz]->writeRemoteData(pairs[xyz].size(), pairs[xyz].data());
                    this->tensors_[xyztot++].tensor_ = new OneElectronOperator<T>(uhf, *ao[xyz]);
                }
            }
        }

        const OneElectronOperator<T>& operator()(const int L, const int xyz) const
        {
            return (*this)(L*(L+1)*(L+2)/6+xyz);
        }

        const OneElectronOperator<T>& operator()(const int x, const int y, const int z) const
        {
            const int L = x+y+z;
            const int xyz = z+(L+1-x)*(L+2-x)/2;
            return (*this)(L, xyz);
        }
};

}
}

#endif
