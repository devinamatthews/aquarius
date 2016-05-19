#ifndef _AQUARIUS_TENSOR_SPINBLOCKED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_SPINBLOCKED_TENSOR_HPP_

#include "../../../frameworks/autocc/autocc.hpp"
#include "../../../frameworks/tensor/tensor.hpp"
#include "../../../frameworks/util/global.hpp"

namespace aquarius
{

template <capability_type C> class SpinBlockedTensor;

namespace detail
{
    template <capability_type C>
    vector<vector<int>> getBlockLengths(const TensorInitializerList<C>& init)
    {
        auto& idxinit = init.template as<INDEXABLE>();
        auto& soinit = init.template as<SPINORBITAL_>();

        vector<vector<int>> len(idxinit.ndim, vector<int>(2));

        int p = 0;
        for (int i = 0;i < soinit.ncreation.size();i++)
        {
            for (int j = 0;j < soinit.ncreation[i];j++, p++)
            {
                len[i][0] = soinit.nalpha[i];
                len[i][1] = soinit.nbeta[i];
            }
        }

        for (int i = 0;i < soinit.nannihilation.size();i++)
        {
            for (int j = 0;j < soinit.nannihilation[i];j++, p++)
            {
                len[i][0] = soinit.nalpha[i];
                len[i][1] = soinit.nbeta[i];
            }
        }
        assert(p == len.size());

        return len;
    }
}

template <capability_type C>
class SpinBlockedTensor : public BlockedTensor<C&~SPINORBITAL_,TensorImplementation<C>>
{
    protected:
        typedef BlockedTensor<C&~SPINORBITAL_,TensorImplementation<C>> Base;
        using typename Base::SubTensor;
        using Base::sliceByBlock;
        using Base::operator();
        using TensorInitializer<INDEXABLE_>::ndim;
        using TensorInitializer<SPINORBITAL_>::nalpha;
        using TensorInitializer<SPINORBITAL_>::nbeta;
        using TensorInitializer<SPINORBITAL_>::nalpha_per_irrep;
        using TensorInitializer<SPINORBITAL_>::nbeta_per_irrep;
        using TensorInitializer<SPINORBITAL_>::ncreation;
        using TensorInitializer<SPINORBITAL_>::nannihilation;
        using TensorInitializer<SPINORBITAL_>::spin;

        template <typename SubFactory, capability_type C_, capability_type C__=C>
        enable_if_t<C__&PGSYMMETRIC_> put(const vector<int>& idx, SubFactory f,
                                          const TensorInitializerList<C_>& subinit)
        {
            vector<int> subsym(ndim, AS);
            int nspaces = nalpha.size();
            int nirrep = this->group.getNumIrreps();
            vector<vector<int>> sublen(ndim, vector<int>(nirrep));

            int dim = 0;
            for (int i = 0;i < nspaces;i++)
            {
                for (int j = 0;j < ncreation[i];j++, dim++)
                {
                    for (int k = 0;k < nirrep;k++)
                    {
                        sublen[dim][k] = (idx[dim] == 0 ?
                                          nalpha_per_irrep[i][k] :
                                          nbeta_per_irrep[i][k]);
                    }
                    if (dim > 0 && idx[i-1] != idx[i]) subsym[dim-1] = NS;
                }
                if (ncreation[i] > 0) subsym[dim-1] = NS;
            }

            for (int i = 0;i < nspaces;i++)
            {
                for (int j = 0;j < nannihilation[i];j++, dim++)
                {
                    for (int k = 0;k < nirrep;k++)
                    {
                        sublen[dim][k] = (idx[dim] == 0 ?
                                          nalpha_per_irrep[i][k] :
                                          nbeta_per_irrep[i][k]);
                    }
                    if (dim > 0 && idx[i-1] != idx[i]) subsym[dim-1] = NS;
                }
                if (nannihilation[i] > 0) subsym[dim-1] = NS;
            }
            assert(dim == ndim);

            Base::put(idx, f(subinit <<
                             TensorInitializer<PGSYMMETRIC|BOUNDED|IPSYMMETRIC>(this->group,
                                                                                sublen,
                                                                                subsym,
                                                                                this->rep)));
        }

        template <typename SubFactory, capability_type C_, capability_type C__=C>
        enable_if_t<!C__&PGSYMMETRIC_> put(const vector<int>& idx, SubFactory f,
                                           const TensorInitializerList<C_>& subinit)
        {
            vector<int> subsym(ndim, AS);
            int nspaces = nalpha.size();
            vector<int> sublen(ndim);

            int dim = 0;
            for (int i = 0;i < nspaces;i++)
            {
                for (int j = 0;j < ncreation[i];j++, dim++)
                {
                    sublen[dim] = (idx[dim] == 0 ? nalpha[i] : nbeta[i]);
                    if (dim > 0 && idx[i-1] != idx[i]) subsym[dim-1] = NS;
                }
                if (ncreation[i] > 0) subsym[dim-1] = NS;
            }

            for (int i = 0;i < nspaces;i++)
            {
                for (int j = 0;j < nannihilation[i];j++, dim++)
                {
                    sublen[dim] = (idx[dim] == 0 ? nalpha[i] : nbeta[i]);
                    if (dim > 0 && idx[i-1] != idx[i]) subsym[dim-1] = NS;
                }
                if (nannihilation[i] > 0) subsym[dim-1] = NS;
            }
            assert(dim == ndim);

            Base::put(idx, f(subinit <<
                             TensorInitializer<BOUNDED|IPSYMMETRIC>(sublen, subsym)));
        }

        template <typename SubFactory>
        void initialize(const INITIALIZER_TYPE(C)& init, SubFactory f)
        {
            INITIALIZER_TYPE(C&~(SPINORBITAL_|PGSYMMETRIC|BOUNDED|IPSYMMETRIC)) subinit(init);

            int nspaces = nalpha.size();
            int ncreationtot = aquarius::sum(ncreation);
            int nannihilationtot = aquarius::sum(nannihilation);
            vector<int> whichout(ncreationtot), whichin(nannihilationtot);
            vector<int> alpha_creation(nspaces), alpha_annihilation(nspaces);

            assert(abs(spin) <= ncreationtot+nannihilationtot);
            assert(abs(spin) >= abs(ncreationtot-nannihilationtot));
            assert(abs(spin)%2 == abs(ncreationtot-nannihilationtot)%2);

            vector<int> idx(ndim);
            for (int alphaout = 0;alphaout <= ncreationtot;alphaout++)
            {
                int alphain = alphaout + (ncreationtot-nannihilationtot-spin)/2;
                if (alphain < 0 || alphain > nannihilationtot) continue;

                fill(whichout.begin(), whichout.end(), 0);

                for (bool doneout = false;!doneout;)
                {
                    fill(alpha_creation.begin(), alpha_creation.end(), 0);

                    for (int i = 0;i < alphaout;i++)
                    {
                        alpha_creation[whichout[i]]++;
                    }

                    fill(whichin.begin(), whichin.end(), 0);

                    for (bool donein = false;!donein;)
                    {
                        fill(alpha_annihilation.begin(), alpha_annihilation.end(), 0);

                        for (int i = 0;i < alphain;i++)
                        {
                            alpha_annihilation[whichin[i]]++;
                        }

                        bool ok = true;
                        for (int i = 0;i < nspaces;i++)
                        {
                            if (alpha_creation[i] > ncreation[i] ||
                                alpha_annihilation[i] > nannihilation[i])
                            {
                                ok = false;
                                break;
                            }
                        }

                        if (ok)
                        {
                            int i = 0;
                            for (int j = 0;j < nspaces;j++)
                            {
                                for (int k = 0;k < alpha_creation[j];k++, i++) idx[i] = 0;
                                for (int k = 0;k <  beta_creation[j];k++, i++) idx[i] = 1;
                            }
                            for (int j = 0;j < nspaces;j++)
                            {
                                for (int k = 0;k < alpha_annihilation[j];k++, i++) idx[i] = 0;
                                for (int k = 0;k <  beta_annihilation[j];k++, i++) idx[i] = 1;
                            }

                            put(idx, f, subinit);
                        }

                        for (int i = 0;i < alphain;i++)
                        {
                            whichin[i]++;
                            if (i < alphain-1)
                            {
                                if (whichin[i] <= whichin[i+1]) break;
                                whichin[i] = 0;
                            }
                            else
                            {
                                if (whichin[i] < nspaces) break;
                                if (i == alphain-1) donein = true;
                            }
                        }

                        if (alphain == 0) donein = true;
                    }

                    for (int i = 0;i < alphaout;i++)
                    {
                        whichout[i]++;
                        if (i < alphaout-1)
                        {
                            if (whichout[i] <= whichout[i+1]) break;
                            whichout[i] = 0;
                        }
                        else
                        {
                            if (whichout[i] < nspaces) break;
                            if (i == alphaout-1) doneout = true;
                        }
                    }

                    if (alphaout == 0) doneout = true;
                }
            }
        }

    public:
        template <typename SubFactory>
        class Factory
        {
            protected:
                SubFactory f;

            public:
                Factory(SubFactory f) : f(f) {}

                Tensor<C> operator()(const INITIALIZER_TYPE(C)& init)
                {
                    return new SpinBlockedTensor(init, f);
                }
        };

        template <typename SubFactory = typename SubTensor::Factory>
        SpinBlockedTensor(const INITIALIZER_TYPE(C)& init, SubFactory f = SubFactory())
        : Base(init, detail::getBlockLengths(init))
        {
            initialize(init, f);
        }

        void sliceBySpin(const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, const TensorImplementation<>& A,
                         const Scalar&  beta,             const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            sliceByBlock(alpha, conja, start_A, A, beta, start_B, length, SubTensor::slice);
        }

        void sliceBySpinAndIrrep(const Scalar& alpha, bool conja, const vector<vector<vector<int>>>& start_A, const TensorImplementation<>& A_,
                                 const Scalar&  beta,             const vector<vector<vector<int>>>& start_B, const vector<vector<vector<int>>>& length)
        {
            sliceByBlock(alpha, conja, start_A, A, beta, start_B, length, SubTensor::sliceByIrrep);
        }

        void getAllKeysBySpin(const vector<int>& spins, KeyVector& keys) const
        {
            (*this)(spins).template impl<BOUNDED_>().getAllKeys(keys);
        }

        void getAllDataBySpin(const vector<int>& spins, KeyValueVector& kv) const
        {
            (*this)(spins).template impl<BOUNDED_>().getAllData(kv);
        }

        void getDataBySpin(const vector<int>& spins, key_type n, key_type* keys, void* values) const
        {
            (*this)(spins).template impl<BOUNDED_>().getData(n, keys, values);
        }

        void setDataBySpin(const vector<int>& spins, key_type n, const key_type* keys, const void* values)
        {
            (*this)(spins).template impl<BOUNDED_>().setData(n, keys, values);
        }

        void addDataBySpin(const vector<int>& spins, key_type n, const Scalar& alpha, const key_type* keys,
                           const void* values, const Scalar& beta)
        {
            (*this)(spins).template impl<BOUNDED_>().addData(n, alpha, keys, values, beta);
        }

        void getLocalKeysBySpin(const vector<int>& spins, KeyVector& keys) const
        {
            (*this)(spins).template impl<DISTRIBUTED_>().getLocalKeys(keys);
        }

        void getLocalDataBySpin(const vector<int>& spins, KeyValueVector& kv) const
        {
            (*this)(spins).template impl<DISTRIBUTED_>().getLocalData(kv);
        }

        void setLocalDataBySpin(const vector<int>& spins, key_type n, const key_type* keys, const void* values)
        {
            (*this)(spins).template impl<DISTRIBUTED_>().setLocalData(n, keys, values);
        }

        void addLocalDataBySpin(const vector<int>& spins, key_type n, const Scalar& alpha, const key_type* keys,
                          const void* values, const Scalar& beta)
        {
            (*this)(spins).template impl<DISTRIBUTED_>().addLocalData(n, alpha, keys, values, beta);
        }

        void getRemoteDataBySpin(const vector<int>& spins, key_type n, key_type* keys, void* values) const
        {
            (*this)(spins).template impl<DISTRIBUTED_>().getRemoteData(n, keys, values);
        }

        void setRemoteDataBySpin(const vector<int>& spins, key_type n, const key_type* keys, const void* values)
        {
            (*this)(spins).template impl<DISTRIBUTED_>().setRemoteData(n, keys, values);
        }

        void addRemoteDataBySpin(const vector<int>& spins, key_type n, const Scalar& alpha, const key_type* keys,
                           const void* values, const Scalar& beta)
        {
            (*this)(spins).template impl<DISTRIBUTED_>().addRemoteData(n, alpha, keys, values, beta);
        }

        void getRemoteDataBySpin(const vector<int>& spins) const
        {
            (*this)(spins).template impl<DISTRIBUTED_>().getRemoteData();
        }

        void setRemoteDataBySpin(const vector<int>& spins)
        {
            (*this)(spins).template impl<DISTRIBUTED_>().setRemoteData();
        }

        void addRemoteDataBySpin(const vector<int>& spins, const Scalar& alpha, const Scalar& beta)
        {
            (*this)(spins).template impl<DISTRIBUTED_>().addRemoteData(alpha, beta);
        }

        void getAllKeysBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyVector& keys) const
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().getAllKeysByIrrep(irreps, keys);
        }

        void getAllDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyValueVector& kv) const
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().getAllDataByIrrep(irreps, kv);
        }

        void getDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, key_type* keys, void* values) const
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().getDataByIrrep(irreps, n, keys, values);
        }

        void setDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const key_type* keys, const void* values)
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().setDataByIrrep(irreps, n, keys, values);
        }

        void addDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                           const void* values, const Scalar& beta)
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().addDataByIrrep(irreps, n, alpha, keys, values, beta);
        }

        void getLocalKeysBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyVector& keys) const
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().getLocalKeysByIrrep(irreps, keys);
        }

        void getLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyValueVector& kv) const
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().getLocalDataByIrrep(irreps, kv);
        }

        void setLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const key_type* keys, const void* values)
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().setLocalDataByIrrep(irreps, n, keys, values);
        }

        void addLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                          const void* values, const Scalar& beta)
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().addLocalDataByIrrep(irreps, n, alpha, keys, values, beta);
        }

        void getRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, key_type* keys, void* values) const
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().getRemoteDataByIrrep(irreps, n, keys, values);
        }

        void setRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const key_type* keys, const void* values)
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().setRemoteDataByIrrep(irreps, n, keys, values);
        }

        void addRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                           const void* values, const Scalar& beta)
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().addRemoteDataByIrrep(irreps, n, alpha, keys, values, beta);
        }

        void getRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps) const
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().getRemoteDataByIrrep(irreps);
        }

        void setRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps)
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().setRemoteDataByIrrep(irreps);
        }

        void addRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const Scalar& alpha, const Scalar& beta)
        {
            (*this)(spins).template impl<PGSYMMETRIC_>().addRemoteDataByIrrep(irreps, alpha, beta);
        }
};

}

#endif
