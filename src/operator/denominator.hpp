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

#ifndef _AQUARIUS_OPERATOR_DENOMINATOR_HPP_
#define _AQUARIUS_OPERATOR_DENOMINATOR_HPP_

#include <vector>

#include "1eoperator.hpp"
#include "mooperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class Denominator : public op::MOOperator
{
    protected:
        std::vector<T> dA, da, dI, di;

    public:
        template <typename Derived>
        Denominator(const OneElectronOperatorBase<T,Derived>& F)
        : MOOperator(F),
          dA(F.vrt.nalpha), da(F.vrt.nbeta),
          dI(F.occ.nalpha), di(F.occ.nbeta)
        {
            if (arena.rank == 0)
            {
                int nA = dA.size();
                int na = da.size();
                int nI = dI.size();
                int ni = di.size();

                {
                    std::vector<tkv_pair<T> > pairs(nA);
                    for (int i = 0;i < nA;i++) pairs[i].k = i+i*nA;
                    F.getAB()(std::vec(1,0),std::vec(1,0)).getRemoteData(pairs);
                    for (int i = 0;i < nA;i++) dA[pairs[i].k/nA] = -pairs[i].d;
                }

                {
                    std::vector<tkv_pair<T> > pairs(na);
                    for (int i = 0;i < na;i++) pairs[i].k = i+i*na;
                    F.getAB()(std::vec(0,0),std::vec(0,0)).getRemoteData(pairs);
                    for (int i = 0;i < na;i++) da[pairs[i].k/na] = -pairs[i].d;
                }

                {
                    std::vector<tkv_pair<T> > pairs(nI);
                    for (int i = 0;i < nI;i++) pairs[i].k = i+i*nI;
                    F.getIJ()(std::vec(0,1),std::vec(0,1)).getRemoteData(pairs);
                    for (int i = 0;i < nI;i++) dI[pairs[i].k/nI] = pairs[i].d;
                }

                {
                    std::vector<tkv_pair<T> > pairs(ni);
                    for (int i = 0;i < ni;i++) pairs[i].k = i+i*ni;
                    F.getIJ()(std::vec(0,0),std::vec(0,0)).getRemoteData(pairs);
                    for (int i = 0;i < ni;i++) di[pairs[i].k/ni] = pairs[i].d;
                }
            }
            else
            {
                F.getAB()(std::vec(1,0),std::vec(1,0)).getRemoteData();
                F.getAB()(std::vec(0,0),std::vec(0,0)).getRemoteData();
                F.getIJ()(std::vec(0,1),std::vec(0,1)).getRemoteData();
                F.getIJ()(std::vec(0,0),std::vec(0,0)).getRemoteData();
            }

            arena.Bcast(dA, 0);
            arena.Bcast(da, 0);
            arena.Bcast(dI, 0);
            arena.Bcast(di, 0);
        }

        const std::vector<T>& getDA() const { return dA; }
        const std::vector<T>& getDa() const { return da; }
        const std::vector<T>& getDI() const { return dI; }
        const std::vector<T>& getDi() const { return di; }
};

}
}

#endif
