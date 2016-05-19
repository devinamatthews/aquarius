#ifndef _AQUARIUS_TENSOR_SYMBLOCKED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_SYMBLOCKED_TENSOR_HPP_

#include "../../../frameworks/symmetry/symmetry.hpp"
#include "../../../frameworks/task/task.hpp"
#include "../../../frameworks/tensor/implementations/blocked_tensor.hpp"
#include "../../../frameworks/tensor/tensor.hpp"
#include "../../../frameworks/util/global.hpp"

namespace aquarius
{

/*
 * SpinBlockedTensor can wrap tensor types which are at least BOUNDED and IPSYMMETRIC
 * and optionally DIVISIBLE and/or DISTRIBUTED.
 */
template <capability_type C>
class SymmetryBlockedTensor : public BlockedTensor<C&~PGSYMMETRIC,TensorImplementation<C>>
{
    static_assert(IS_SUPERSET_OF(C,IPSYMMETRIC|BOUNDED|PGSYMMETRIC)&&
                  IS_SUPERSET_OF(PGSYMMETRIC|BOUNDED|IPSYMMETRIC|DIVISIBLE|DISTRIBUTED,C), "");

    protected:
        typedef BlockedTensor<C&~PGSYMMETRIC,TensorImplementation<C>> Base;
        using typename Base::SubTensor;
        using Base::put;
        using Base::sliceByBlock;
        using Base::operator();
        using TensorInitializer<INDEXABLE_>::ndim;
        using TensorInitializer<PGSYMMETRIC_>::group;
        using TensorInitializer<PGSYMMETRIC_>::len_per_irrep;

        template <typename SubFactory, capability_type C_>
        enable_if_t<C_&IPSYMMETRIC_> initialize(const TensorInitializerList<C_>& init, SubFactory f)
        {
            INITIALIZER_TYPE(C&~(PGSYMMETRIC|BOUNDED|IPSYMMETRIC)) subinit(init);

            int n = group.getNumIrreps();
            vector<Representation> irreps;
            for (int i = 0;i < n;i++) irreps.push_back(group.getIrrep(i));

            vector<int> idx(ndim);
            vector<int> sublen(ndim);
            vector<int> subsym(ndim);
            vector<Representation> prod(ndim+1, rep);
            for (bool done = false;!done;t++)
            {
                if (prod[0].isTotallySymmetric())
                {
                    for (int i = 0;i < ndim-1;i++)
                    {
                        sublen[i] = len_per_irrep[i][idx[i]];
                        subsym[i] = (idx[i] == idx[i+1] ? this->sym[i] : NS);
                    }

                    put(idx, f(subinit <<
                               TensorInitializer<BOUNDED|IPSYMMETRIC>(sublen, subsym)));
                }

                for (int i = 0;i < ndim;i++)
                {
                    idx[i]++;

                    if (idx[i] == (this->sym[i] == NS ? n : idx[i+1]+1))
                    {
                        idx[i] = (i == 0 || this->sym[i-1] == NS ? 0 : idx[i-1]);
                        if (i == ndim-1) done = true;
                    }
                    else
                    {
                        for (int j = i;j >= 0;j--)
                        {
                            prod[j] = prod[j+1];
                            prod[j] *= irreps[idx[j]];
                        }
                        break;
                    }
                }

                if (ndim == 0) done = true;
            }
        }

        template <typename SubFactory, capability_type C_>
        enable_if_t<!C_&IPSYMMETRIC_> initialize(const TensorInitializerList<C_>& init, SubFactory f)
        {
            INITIALIZER_TYPE(C&~(PGSYMMETRIC|BOUNDED)) subinit(init);

            int n = group.getNumIrreps();
            vector<Representation> irreps;
            for (int i = 0;i < n;i++) irreps.push_back(group.getIrrep(i));

            vector<int> idx(ndim);
            vector<int> sublen(ndim);
            vector<Representation> prod(ndim+1, rep);
            for (bool done = false;!done;t++)
            {
                if (prod[0].isTotallySymmetric())
                {
                    for (int i = 0;i < ndim-1;i++)
                    {
                        sublen[i] = len_per_irrep[i][idx[i]];
                    }

                    put(idx, f(subinit << TensorInitializer<BOUNDED>(sublen)));
                }

                for (int i = 0;i < ndim;i++)
                {
                    idx[i]++;

                    if (idx[i] == n)
                    {
                        idx[i] = 0;
                        if (i == ndim-1) done = true;
                    }
                    else
                    {
                        for (int j = i;j >= 0;j--)
                        {
                            prod[j] = prod[j+1];
                            prod[j] *= irreps[idx[j]];
                        }
                        break;
                    }
                }

                if (ndim == 0) done = true;
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
                    return Tensor<C>(new SymmetryBlockedTensor(init, f));
                }
        };

        template <typename SubFactory = typename SubTensor::Factory>
        SymmetryBlockedTensor(const INITIALIZER_TYPE(C)& init, SubFactory f = SubFactory())
        : Base(init, init.template as<PGSYMMETRIC_>().len_per_irrep)
        {
            initialize(init, f);
        }

        void sliceByIrrep(const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, const TensorImplementation<>& A,
                          const Scalar&  beta,             const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            sliceByBlock(alpha, conja, start_A, A, beta, start_B, length, SubTensor::slice);
        }

        void getAllKeysByIrrep(const vector<int>& irreps, KeyVector& keys) const
        {
            (*this)(irreps).template impl<BOUNDED_>().getAllKeys(keys);
        }

        void getAllDataByIrrep(const vector<int>& irreps, KeyValueVector& kv) const
        {
            (*this)(irreps).template impl<BOUNDED_>().getAllData(kv);
        }

        void getDataByIrrep(const vector<int>& irreps, key_type n, key_type* keys, void* values) const
        {
            (*this)(irreps).template impl<BOUNDED_>().getData(n, keys, values);
        }

        void setDataByIrrep(const vector<int>& irreps, key_type n, const key_type* keys, const void* values)
        {
            (*this)(irreps).template impl<BOUNDED_>().setData(n, keys, values);
        }

        void addDataByIrrep(const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                           const void* values, const Scalar& beta)
        {
            (*this)(irreps).template impl<BOUNDED_>().addData(n, alpha, keys, values, beta);
        }

        void getLocalKeysByIrrep(const vector<int>& irreps, KeyVector& keys) const
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().getLocalKeys(keys);
        }

        void getLocalDataByIrrep(const vector<int>& irreps, KeyValueVector& kv) const
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().getLocalData(kv);
        }

        void setLocalDataByIrrep(const vector<int>& irreps, key_type n, const key_type* keys, const void* values)
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().setLocalData(n, keys, values);
        }

        void addLocalDataByIrrep(const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                          const void* values, const Scalar& beta)
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().addLocalData(n, alpha, keys, values, beta);
        }

        void getRemoteDataByIrrep(const vector<int>& irreps, key_type n, key_type* keys, void* values) const
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().getRemoteData(n, keys, values);
        }

        void setRemoteDataByIrrep(const vector<int>& irreps, key_type n, const key_type* keys, const void* values)
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().setRemoteData(n, keys, values);
        }

        void addRemoteDataByIrrep(const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                           const void* values, const Scalar& beta)
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().addRemoteData(n, alpha, keys, values, beta);
        }

        void getRemoteDataByIrrep(const vector<int>& irreps) const
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().getRemoteData();
        }

        void setRemoteDataByIrrep(const vector<int>& irreps)
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().setRemoteData();
        }

        void addRemoteDataByIrrep(const vector<int>& irreps, const Scalar& alpha, const Scalar& beta)
        {
            (*this)(irreps).template impl<DISTRIBUTED_>().addRemoteData(alpha, beta);
        }
};

}

#endif
