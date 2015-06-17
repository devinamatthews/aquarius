#ifndef _AQUARIUS_OPERATOR_COMPLEXDENOMINATOR_HPP_
#define _AQUARIUS_OPERATOR_COMPLEXDENOMINATOR_HPP_

#include "util/global.hpp"

#include "tensor/spinorbital_tensor.hpp"
#include "tensor/symblocked_tensor.hpp"
#include "tensor/ctf_tensor.hpp"

#include "1eoperator.hpp"
#include "mooperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class ComplexDenominator : public op::MOOperator
{
    protected:
        vector<vector<T>> dA, da, dI, di;

        void weight(const vector<int>& alpha_out, const vector<int>& alpha_in,
                    T omega, tensor::SymmetryBlockedTensor<T>& x,
                    T   eta, tensor::SymmetryBlockedTensor<T>& y)
        {
            const symmetry::PointGroup& group = x.getGroup();
            int n = group.getNumIrreps();
            int ndim = x.getDimension();
            vector<int> irreps(ndim);
            vector<symmetry::Representation> reps(ndim+1, x.getRepresentation());

            for (bool done = false;!done && ndim;)
            {
                if (reps[0].isTotallySymmetric())
                {
                    weight(alpha_out, alpha_in, irreps,
                           omega, x(irreps), eta, y(irreps));
                }

                for (int i = 0;i < ndim;i++)
                {
                    irreps[i]++;

                    if (irreps[i] == n)
                    {
                        irreps[i] = 0;
                        if (i == ndim-1) done = true;
                    }
                    else
                    {
                        for (;i >= 0;i--)
                        {
                            reps[i] = reps[i+1]*group.getIrrep(irreps[i]);
                        }
                        break;
                    }
                }
            }
        }

        void weight(const vector<int>& alpha_out, const vector<int>& alpha_in,
                    const vector<int>& irreps,
                    T omega, tensor::CTFTensor<T>& x,
                    T   eta, tensor::CTFTensor<T>& y)
        {

        }

    public:
        template <typename Derived>
        ComplexDenominator(const OneElectronOperatorBase<T,Derived>& F)
        : op::MOOperator(F)
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
                        vector<tkv_pair<T>> pairs(nA);
                        for (int i = 0;i < nA;i++) pairs[i].k = i+i*nA;
                        F.getAB()({1,0},{1,0}).getRemoteData(irreps, pairs);
                        for (int i = 0;i < nA;i++) dA[j][pairs[i].k/nA] = -pairs[i].d;
                    }

                    {
                        vector<tkv_pair<T>> pairs(na);
                        for (int i = 0;i < na;i++) pairs[i].k = i+i*na;
                        F.getAB()({0,0},{0,0}).getRemoteData(irreps, pairs);
                        for (int i = 0;i < na;i++) da[j][pairs[i].k/na] = -pairs[i].d;
                    }

                    {
                        vector<tkv_pair<T>> pairs(nI);
                        for (int i = 0;i < nI;i++) pairs[i].k = i+i*nI;
                        F.getIJ()({0,1},{0,1}).getRemoteData(irreps, pairs);
                        for (int i = 0;i < nI;i++) dI[j][pairs[i].k/nI] = pairs[i].d;
                    }

                    {
                        vector<tkv_pair<T>> pairs(ni);
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

                arena.comm().Bcast(dA[j], 0);
                arena.comm().Bcast(da[j], 0);
                arena.comm().Bcast(dI[j], 0);
                arena.comm().Bcast(di[j], 0);
            }
        }

        void weight(T omega, tensor::SpinorbitalTensor<T>& x,
                    T   eta, tensor::SpinorbitalTensor<T>& y)
        {

        }
};

}
}

#endif
