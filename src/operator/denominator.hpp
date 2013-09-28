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
        std::vector<std::vector<T> > dA, da, dI, di;

    public:
        template <typename Derived>
        Denominator(const OneElectronOperatorBase<T,Derived>& F)
        : MOOperator(F)
        {
            int n = vrt.group.getNumIrreps();

            dA.resize(n);
            da.resize(n);
            dI.resize(n);
            di.resize(n);

            for (int j = 0;j < n;j++)
            {
                dA[j].resize(vrt.nalpha[j]);
                da[j].resize(vrt.nbeta[j]);
                dI[j].resize(occ.nalpha[j]);
                di[j].resize(occ.nbeta[j]);

                std::vector<int> irreps(2,j);

                if (arena.rank == 0)
                {
                    int nA = dA[j].size();
                    int na = da[j].size();
                    int nI = dI[j].size();
                    int ni = di[j].size();

                    {
                        std::vector<tkv_pair<T> > pairs(nA);
                        for (int i = 0;i < nA;i++) pairs[i].k = i+i*nA;
                        F.getAB()(std::vec(1,0),std::vec(1,0))(irreps).getRemoteData(pairs);
                        for (int i = 0;i < nA;i++) dA[j][pairs[i].k/nA] = -pairs[i].d;
                    }

                    {
                        std::vector<tkv_pair<T> > pairs(na);
                        for (int i = 0;i < na;i++) pairs[i].k = i+i*na;
                        F.getAB()(std::vec(0,0),std::vec(0,0))(irreps).getRemoteData(pairs);
                        for (int i = 0;i < na;i++) da[j][pairs[i].k/na] = -pairs[i].d;
                    }

                    {
                        std::vector<tkv_pair<T> > pairs(nI);
                        for (int i = 0;i < nI;i++) pairs[i].k = i+i*nI;
                        F.getIJ()(std::vec(0,1),std::vec(0,1))(irreps).getRemoteData(pairs);
                        for (int i = 0;i < nI;i++) dI[j][pairs[i].k/nI] = pairs[i].d;
                    }

                    {
                        std::vector<tkv_pair<T> > pairs(ni);
                        for (int i = 0;i < ni;i++) pairs[i].k = i+i*ni;
                        F.getIJ()(std::vec(0,0),std::vec(0,0))(irreps).getRemoteData(pairs);
                        for (int i = 0;i < ni;i++) di[j][pairs[i].k/ni] = pairs[i].d;
                    }
                }
                else
                {
                    F.getAB()(std::vec(1,0),std::vec(1,0))(irreps).getRemoteData();
                    F.getAB()(std::vec(0,0),std::vec(0,0))(irreps).getRemoteData();
                    F.getIJ()(std::vec(0,1),std::vec(0,1))(irreps).getRemoteData();
                    F.getIJ()(std::vec(0,0),std::vec(0,0))(irreps).getRemoteData();
                }

                arena.Bcast(dA[j], 0);
                arena.Bcast(da[j], 0);
                arena.Bcast(dI[j], 0);
                arena.Bcast(di[j], 0);
            }
        }

        const std::vector<std::vector<T> >& getDA() const { return dA; }
        const std::vector<std::vector<T> >& getDa() const { return da; }
        const std::vector<std::vector<T> >& getDI() const { return dI; }
        const std::vector<std::vector<T> >& getDi() const { return di; }
};

}
}

#endif
