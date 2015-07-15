#ifndef _AQUARIUS_CC_COMPLEX_DENOMINATOR_HPP_
#define _AQUARIUS_CC_COMPLEX_DENOMINATOR_HPP_

#include "util/global.hpp"

#include "operator/denominator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/deexcitationoperator.hpp"
#include "symmetry/symmetry.hpp"
#include "tensor/ctf_tensor.hpp"
#include "tensor/symblocked_tensor.hpp"
#include "tensor/spinorbital_tensor.hpp"

namespace aquarius
{
namespace cc
{

template <typename T>
class ComplexDenominator : op::Denominator<T>
{
    protected:
        using op::Denominator<T>::dA;
        using op::Denominator<T>::da;
        using op::Denominator<T>::dI;
        using op::Denominator<T>::di;

    public:
        template <typename Derived>
        ComplexDenominator(const op::OneElectronOperatorBase<T,Derived>& F)
        : op::Denominator<T>(F) {}

        template <int np, int nh>
        void weight(op::ExcitationOperator<T,np,nh>& R,
                    op::ExcitationOperator<T,np,nh>& I, const complex<T>& omega) const
        {
            for (int ex = abs(np-nh);ex <= max(np,nh);ex++)
            {
                weight(R(ex), I(ex), omega);
            }
        }

        template <int np, int nh>
        void weight(op::DeexcitationOperator<T,np,nh>& R,
                    op::DeexcitationOperator<T,np,nh>& I, const complex<T>& omega) const
        {
            for (int ex = abs(np-nh);ex <= max(np,nh);ex++)
            {
                weight(R(ex), I(ex), omega);
            }
        }

        void weight(tensor::SpinorbitalTensor<T>& R,
                    tensor::SpinorbitalTensor<T>& I, const complex<T>& omega) const
        {
            using namespace tensor;

            T omega2 = omega.real()*omega.real() +
                       omega.imag()+omega.imag();

            const vector<int>& nout = R.getNumOut();
            const vector<int>& nin = R.getNumIn();
            int spin = R.getSpin();

            const symmetry::PointGroup& group = R.getGroup();
            int n = group.getNumIrreps();

            int nouttot = sum(nout);
            int nintot = sum(nin);
            int ndim = nouttot+nintot;

            for (int alpha_out = 0;alpha_out <= nouttot;alpha_out++)
            {
                int alpha_in = alpha_out + (nintot-nouttot-spin)/2;

                for (int nA = max(0,alpha_out-nout[1]);nA <= min(alpha_out,nout[0]);nA++)
                {
                    int nI = alpha_out-nA;

                    for (int nB = max(0,alpha_in-nin[1]);nB <= min(alpha_in,nin[0]);nB++)
                    {
                        int nJ = alpha_in-nB;

                        ptr_vector<const vector<vector<T>>> dens;
                        for (int i = 0;i <         nA;i++) dens.push_back(&dA);
                        for (int i = 0;i < nout[0]-nA;i++) dens.push_back(&da);
                        for (int i = 0;i <         nI;i++) dens.push_back(&dI);
                        for (int i = 0;i < nout[1]-nI;i++) dens.push_back(&di);
                        for (int i = 0;i <         nB;i++) dens.push_back(&dA);
                        for (int i = 0;i <  nin[0]-nB;i++) dens.push_back(&da);
                        for (int i = 0;i <         nJ;i++) dens.push_back(&dI);
                        for (int i = 0;i <  nin[1]-nJ;i++) dens.push_back(&di);

                        vector<int> sym(ndim);
                        vector<tkv_pair<T>> pairsr;
                        vector<tkv_pair<T>> pairsi;
                        for (bool done = false;!done;)
                        {
                            if (R({nA,nI},{nB,nJ}).exists(sym))
                            {
                                R({nA,nI},{nB,nJ})(sym).getLocalData(pairsr);
                                I({nA,nI},{nB,nJ})(sym).getLocalData(pairsi);
                                assert(pairsr.size() == pairsi.size());

                                sort(pairsr);
                                sort(pairsi);

                                for (int64_t i = 0;i < pairsr.size();i++)
                                {
                                    int64_t k = pairsr[i].k;
                                    assert(k == pairsi[i].k);

                                    T den = T();
                                    for (int i = 0;i < ndim;i++)
                                    {
                                        int len = dens[i][sym[i]].size();
                                        int idx = k%len;
                                        k /= len;
                                        den += dens[i][sym[i]][idx];
                                    }

                                    T num = den+omega.real();
                                    den = num*num+omega.imag()*omega.imag();
                                    T re = pairsr[i].d;
                                    T im = pairsi[i].d;
                                    pairsr[i].d = (num*re+omega.imag()*im)/den;
                                    pairsi[i].d = (num*im-omega.imag()*re)/den;
                                }

                                R({nA,nI},{nB,nJ})(sym).writeRemoteData(pairsr);
                                I({nA,nI},{nB,nJ})(sym).writeRemoteData(pairsi);
                            }

                            for (int i = 0;i < ndim;i++)
                            {
                                sym[i]++;

                                if (sym[i] == n)
                                {
                                    sym[i] = 0;
                                    if (i == ndim-1) done = true;
                                }
                                else break;
                            }
                        }
                    }
                }
            }
        }
};

}
}

#endif
