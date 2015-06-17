#ifndef _AQUARIUS_TENSOR_BLOCKED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_BLOCKED_TENSOR_HPP_

#include "util/global.hpp"

#include "tensor/tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <capability_type C, class Base>
class BlockedTensor : public Base
{
    static_assert(IS_SUPERSET_OF(C,BOUNDED), "");

    protected:
        vector<vector<int>> blocklen;
        unique_vector<Tensor<C>> tensors;

        void sortIntoBlocks(key_type n, const key_type* keys, const void* values,
                            vector<pair<vector<key_type>,vector<char>>>& kvs) const
        {
            int sz = this->R.size;

            int nb = 1;
            for (int d = 0;d < this->ndim;d++)
            {
                nb *= this->blocklen[d].size();
            }

            kvs.clear();
            kvs.resize(nb);

            if (nb == 1)
            {
                kvs.clear();
                kvs.resize(1);
                kvs[0].first.assign(keys, keys+n);
                kvs[0].second.assign((const char*)values, (const char*)values+n*sz);
                return;
            }

            for (key_type ki = 0;ki < n;ki++)
            {
                key_type k = keys[ki];
                key_type bk = 0;
                key_type bks = 1;
                int b = 0;
                int bs = 1;
                for (int d = 0;d < this->ndim;d++)
                {
                    key_type i = k%this->len[d];
                    k          = k/this->len[d];

                    int bi = 0;
                    for (;bi < this->blocklen[d].size() && i >= this->blocklen[d][bi];
                          bi++, i -= this->blocklen[d][bi]);
                    assert(bi < this->blocklen[d].size());

                    b += bs*bi;
                    bs *= this->blocklen[d].size();
                    bk += bks*i;
                    bks *= this->blocklen[d][bi];
                }
                assert(b < kvs.size());
                kvs[b].first.push_back(bk);
                kvs[b].second.insert(kvs[b].second.end(),
                                     (const char*)values+ki*sz,
                                     (const char*)values+(ki+1)*sz);
            }
        }

        const Tensor<C>* get(const vector<int>& idx) const
        {
            assert(idx.size() == blocklen.size());

            int off = 0;
            int stride = 1;
            for (int i = 0;i < idx.size();i++)
            {
                assert(idx[i] >= 0 && idx[i] < blocklen[i].size());
                off += stride*idx[i];
                stride *= blocklen[i].size();
            }

            return tensors.ptr(off).get();
        }

        void put(const vector<int>& idx, Tensor<C>* t)
        {
            assert(idx.size() == blocklen.size());

            int off = 0;
            int stride = 1;
            for (int i = 0;i < idx.size();i++)
            {
                assert(idx[i] >= 0 && idx[i] < blocklen[i].size());
                off += stride*idx[i];
                stride *= blocklen[i].size();
            }

            tensors.ptr(off).reset(t);
        }

    public:
        template <capability_type C_>
        BlockedTensor(const TensorInitializerList<C_>& ilist,
                      const vector<vector<int>>& blocklen)
        : Base(ilist), blocklen(blocklen)
        {
            size_t n = 1;
            for (auto& v : blocklen) n *= v.size();
            tensors.resize(n);
        }

        const vector<vector<int>>& getBlockLengths() const
        {
            return this->blocklen;
        }

        bool exists(const vector<int>& idx) const
        {
            return get(idx);
        }

        Tensor<C>& operator()(const vector<int>& idx)
        {
            return *const_cast<Tensor<C>*>(const_cast<const BlockedTensor&>(*this).get(idx));
        }

        const Tensor<C>& operator()(const vector<int>& idx) const
        {
            return *get(idx);
        }

        void getAllKeys(KeyVector& keys) const
        {
            vector<int> b(this->ndim);
            KeyVector subkeys;
            for (bool done = false;!done && this->ndim;)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().getAllKeys(subkeys);
                    keys += subkeys;
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getAllData(KeyValueVector& kv) const
        {
            vector<int> b(this->ndim);
            KeyValueVector subkv(kv.field());
            for (bool done = false;!done && this->ndim;)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().getAllData(subkv);
                    kv += subkv;
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getData(key_type n, key_type* keys, void* values) const
        {
            int sz = this->R.size;

            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(this->ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && this->ndim;bi++)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().getData(kvs[bi].first.size(),
                                       kvs[bi].first.data(),
                                       kvs[bi].second.data());
                }

                copy(kvs[bi].first.begin(), kvs[bi].first.end(), keys+i);
                copy(kvs[bi].second.begin(), kvs[bi].second.end(), (char*)values+i*sz);
                i += kvs[bi].first.size();

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void setData(key_type n, const key_type* keys, const void* values)
        {
            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(this->ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && this->ndim;bi++)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().setData(kvs[bi].first.size(),
                                       kvs[bi].first.data(),
                                       kvs[bi].second.data());
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void addData(key_type n, const Scalar& alpha, const key_type* keys,
                     const void* values, const Scalar& beta)
        {
            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(this->ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && this->ndim;bi++)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().addData(kvs[bi].first.size(), alpha,
                                       kvs[bi].first.data(),
                                       kvs[bi].second.data(), beta);
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getLocalKeys(KeyVector& keys) const
        {
            int sz = this->R.size;

            vector<int> b(this->ndim);
            KeyVector subkeys;
            for (bool done = false;!done && this->ndim;)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().getLocalKeys(subkeys);
                    keys += subkeys;
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getLocalData(KeyValueVector& kv) const
        {
            int sz = this->R.size;

            vector<int> b(this->ndim);
            KeyValueVector subkv(kv.field());
            for (bool done = false;!done && this->ndim;)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().getLocalData(subkv);
                    kv += subkv;
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void setLocalData(key_type n, const key_type* keys, const void* values)
        {
            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(this->ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && this->ndim;bi++)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().setLocalData(kvs[bi].first.size(),
                                            kvs[bi].first.data(),
                                            kvs[bi].second.data());
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void addLocalData(key_type n, const Scalar& alpha, const key_type* keys,
                          const void* values, const Scalar& beta)
        {
            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(this->ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && this->ndim;bi++)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().addLocalData(kvs[bi].first.size(), alpha,
                                            kvs[bi].first.data(),
                                            kvs[bi].second.data(), beta);
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getRemoteData(key_type n, key_type* keys, void* values) const
        {
            int sz = this->R.size;

            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(this->ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && this->ndim;bi++)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().getRemoteData(kvs[bi].first.size(),
                                             kvs[bi].first.data(),
                                             kvs[bi].second.data());
                }

                copy(kvs[bi].first.begin(), kvs[bi].first.end(), keys+i);
                copy(kvs[bi].second.begin(), kvs[bi].second.end(), (char*)values+i*sz);
                i += kvs[bi].first.size();

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void setRemoteData(key_type n, const key_type* keys, const void* values)
        {
            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(this->ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && this->ndim;bi++)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().setRemoteData(kvs[bi].first.size(),
                                             kvs[bi].first.data(),
                                             kvs[bi].second.data());
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void addRemoteData(key_type n, const Scalar& alpha, const key_type* keys,
                           const void* values, const Scalar& beta)
        {
            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(this->ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && this->ndim;bi++)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().addRemoteData(kvs[bi].first.size(), alpha,
                                             kvs[bi].first.data(),
                                             kvs[bi].second.data(), beta);
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getRemoteData() const
        {
            vector<int> b(this->ndim);
            for (bool done = false;!done && this->ndim;)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().getRemoteData();
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void setRemoteData()
        {
            vector<int> b(this->ndim);
            for (bool done = false;!done && this->ndim;)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().setRemoteData();
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void addRemoteData(const Scalar& alpha, const Scalar& beta)
        {
            vector<int> b(this->ndim);
            for (bool done = false;!done && this->ndim;)
            {
                if (this->exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().addRemoteData(alpha, beta);
                }

                for (int d = 0;d < this->ndim;d++)
                {
                    b[d]++;
                    if (b[d] == this->blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == this->ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void slice(const Scalar& alpha, bool conja, const vector<int>& start_A, const TensorImplementation<>& A,
                   const Scalar&  beta,             const vector<int>& start_B, const vector<int>& length)
        {
            assert(0);
        }
};

}
}

#endif
