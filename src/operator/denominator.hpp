#ifndef _AQUARIUS_OPERATOR_DENOMINATOR_HPP_
#define _AQUARIUS_OPERATOR_DENOMINATOR_HPP_

#include "util/global.hpp"

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
        vector<vector<T> > dA, da, dI, di;

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

                vector<int> irreps(2,j);

                if (arena.rank == 0)
                {
                    int nA = dA[j].size();
                    int na = da[j].size();
                    int nI = dI[j].size();
                    int ni = di[j].size();

                    {
                        vector<tkv_pair<T> > pairs(nA);
                        for (int i = 0;i < nA;i++) pairs[i].k = i+i*nA;
                        F.getAB()({1,0},{1,0}).getRemoteData(irreps, pairs);
                        for (int i = 0;i < nA;i++) dA[j][pairs[i].k/nA] = -pairs[i].d;
                    }

                    {
                        vector<tkv_pair<T> > pairs(na);
                        for (int i = 0;i < na;i++) pairs[i].k = i+i*na;
                        F.getAB()({0,0},{0,0}).getRemoteData(irreps, pairs);
                        for (int i = 0;i < na;i++) da[j][pairs[i].k/na] = -pairs[i].d;
                    }

                    {
                        vector<tkv_pair<T> > pairs(nI);
                        for (int i = 0;i < nI;i++) pairs[i].k = i+i*nI;
                        F.getIJ()({0,1},{0,1}).getRemoteData(irreps, pairs);
                        for (int i = 0;i < nI;i++) dI[j][pairs[i].k/nI] = pairs[i].d;
                    }

                    {
                        vector<tkv_pair<T> > pairs(ni);
                        for (int i = 0;i < ni;i++) pairs[i].k = i+i*ni;
                        F.getIJ()({0,0},{0,0}).getRemoteData(irreps, pairs);
                        for (int i = 0;i < ni;i++) di[j][pairs[i].k/ni] = pairs[i].d;
                    }
                }
                else
                {
                    F.getAB()({1,0},{1,0}).getRemoteData(irreps);
                    F.getAB()({0,0},{0,0}).getRemoteData(irreps);
                    F.getIJ()({0,1},{0,1}).getRemoteData(irreps);
                    F.getIJ()({0,0},{0,0}).getRemoteData(irreps);
                }

                arena.Bcast(dA[j], 0);
                arena.Bcast(da[j], 0);
                arena.Bcast(dI[j], 0);
                arena.Bcast(di[j], 0);
            }
        }

        const vector<vector<T> >& getDA() const { return dA; }
        const vector<vector<T> >& getDa() const { return da; }
        const vector<vector<T> >& getDI() const { return dI; }
        const vector<vector<T> >& getDi() const { return di; }
};

}
}

#endif
