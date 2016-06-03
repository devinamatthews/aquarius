#ifndef _AQUARIUS_TENSOR_IMPLEMENTATIONS_BLOCKED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_IMPLEMENTATIONS_BLOCKED_TENSOR_HPP_

#include "frameworks/util.hpp"
#include "frameworks/tensor.hpp"

namespace aquarius
{
namespace tensor
{

namespace detail
{
    inline bool checkAlignment(int idx, const string& idx_A, const vector<int>& sym_A,
                                        const string& idx_B, const vector<int>& sym_B)
    {
        if (idx == idx_A.size()-1 || sym_A[idx] == NS) return false;

        for (auto i = idx_B.begin();;)
        {
            i = find(i, idx_B.end(), idx_A[idx]);
            if (i == idx_B.end()) return true;
            if (i == idx_B.end()--) return false;

            int j = i-idx_B.begin();

            if (sym_A[idx] != sym_B[j] ||
                idx_A[idx+1] != idx_B[j+1]) return false;
        }

        return true;
    }

    inline vector<int> getStrides(const string& indices, const vector<int>& len, const string& idx_A)
    {
        vector<int> strides(indices.size());
        vector<int> stride_A(idx_A.size());

        if (!stride_A.empty()) stride_A[0] = 1;
        for (int i = 1;i < idx_A.size();i++)
        {
            stride_A[i] = stride_A[i-1]*len[i-1];
        }

        for (int i = 0;i < indices.size();i++)
        {
            for (int j = 0;j < idx_A.size();j++)
            {
                if (indices[i] == idx_A[j])
                {
                    strides[i] += stride_A[j];
                }
            }
        }

        return strides;
    }

    template<typename RAIterator>
    int relative_sign(RAIterator s1b, RAIterator s1e, RAIterator s2b, RAIterator s2e)
    {
        int sz = s1e-s1b;
        assert(sz == (int)(s2e-s2b));
        int i, k;
        int sign = 1;
        std::vector<bool> seen(sz);

        for (i = 0;i < sz;i++) seen[i] = false;

        for (i = 0;i < sz;i++)
        {
            if (seen[i]) continue;
            int j = i;
            while (true)
            {
                for (k = 0;k < sz && (!(s1b[k] == s2b[j]) || seen[k]);k++);
                assert(k < sz);
                j = k;
                seen[j] = true;
                if (j == i) break;
                sign = -sign;
            }
        }

        return sign;
    }

    template <typename Container>
    int relative_sign(Container&& a, Container&& b)
    {
        return relative_sign(a.begin(), a.end(), b.begin(), b.end());
    }

    struct index_locator_
    {
        int sort;
        int idx;
        int pos_A;
        int pos_B;
        int pos_C;

        index_locator_(int sort, int idx, int pos_A, int pos_B, int pos_C)
        : sort(sort), idx(idx), pos_A(pos_A), pos_B(pos_B), pos_C(pos_C) {}

        static bool sortA(const index_locator_& a, const index_locator_& b)
        {
            return a.pos_A < b.pos_A;
        }

        static bool sortB(const index_locator_& a, const index_locator_& b)
        {
            return a.pos_B < b.pos_B;
        }

        static bool sortC(const index_locator_& a, const index_locator_& b)
        {
            return a.pos_C < b.pos_C;
        }

        bool operator==(int idx)
        {
            return this->idx == idx;
        }
    };

    template <typename T>
    int align_symmetric_indices(int order_A, T& idx_A, const int* sym_A,
                                int order_B, T& idx_B, const int* sym_B)
    {
        int fact = 1;

        std::vector<index_locator_> indices;

        for (int i = 0;i < order_A;i++)
        {
            int i_in_B; for (i_in_B = 0;i_in_B < order_B && idx_A[i] != idx_B[i_in_B];i_in_B++);
            if (i_in_B == order_B) continue;

            indices.push_back(index_locator_(0, idx_A[i], i, i_in_B, 0));
        }

        while (!indices.empty())
        {
            std::vector<index_locator_> group;
            group.push_back(indices[0]);
            group.back().sort = 0;
            indices.erase(indices.begin());

            int s = 1;
            for (std::vector<index_locator_>::iterator it = indices.begin();;)
            {
                if (it == indices.end()) break;

                if ((group[0].pos_A == -1 && it->pos_A != -1) ||
                    (group[0].pos_A != -1 && it->pos_A == -1) ||
                    (group[0].pos_B == -1 && it->pos_B != -1) ||
                    (group[0].pos_B != -1 && it->pos_B == -1))
                {
                    ++it;
                    continue;
                }

                bool sym_in_A = false;
                for (int k = group[0].pos_A-1;k >= 0 && sym_A[k] != NS;k--)
                {
                    if (idx_A[k] == it->idx)
                    {
                        sym_in_A = true;
                        break;
                    }
                }
                for (int k = group[0].pos_A+1;k < order_A && sym_A[k-1] != NS;k++)
                {
                    if (idx_A[k] == it->idx)
                    {
                        sym_in_A = true;
                        break;
                    }
                }
                if (!sym_in_A)
                {
                    ++it;
                    continue;
                }

                bool sym_in_B = false;
                for (int k = group[0].pos_B-1;k >= 0 && sym_B[k] != NS;k--)
                {
                    if (idx_B[k] == it->idx)
                    {
                        sym_in_B = true;
                        break;
                    }
                }
                for (int k = group[0].pos_B+1;k < order_B && sym_B[k-1] != NS;k++)
                {
                    if (idx_B[k] == it->idx)
                    {
                        sym_in_B = true;
                        break;
                    }
                }
                if (!sym_in_B)
                {
                    ++it;
                    continue;
                }

                group.push_back(*it);
                group.back().sort = s++;
                it = indices.erase(it);
            }

            if (group.size() <= 1) continue;

            std::vector<int> order_A, order_B;

            for (int i = 0;i < (int)group.size();i++)
                order_A.push_back(group[i].sort);

            std::sort(group.begin(), group.end(), index_locator_::sortB);
            for (int i = 0;i < (int)group.size();i++)
            {
                order_B.push_back(group[i].sort);
                idx_B[group[group[i].sort].pos_B] = group[i].idx;
            }
            if (sym_B[group[0].pos_B] == AS)
                fact *= relative_sign(order_A, order_B);
        }

        return fact;
    }

    template <typename T>
    int align_symmetric_indices(int order_A, T& idx_A, const int* sym_A,
                                int order_B, T& idx_B, const int* sym_B,
                                int order_C, T& idx_C, const int* sym_C)
    {
        int fact = 1;

        std::vector<index_locator_> indices;

        for (int i = 0;i < order_A;i++)
        {
            int i_in_B; for (i_in_B = 0;i_in_B < order_B && idx_A[i] != idx_B[i_in_B];i_in_B++);
            if (i_in_B == order_B) i_in_B = -1;

            int i_in_C; for (i_in_C = 0;i_in_C < order_C && idx_A[i] != idx_C[i_in_C];i_in_C++);
            if (i_in_C == order_C) i_in_C = -1;

            if (i_in_B == -1 && i_in_C == -1) continue;

            indices.push_back(index_locator_(0, idx_A[i], i, i_in_B, i_in_C));
        }

        for (int i = 0;i < order_B;i++)
        {
            int i_in_A; for (i_in_A = 0;i_in_A < order_A && idx_B[i] != idx_A[i_in_A];i_in_A++);
            if (i_in_A == order_A) i_in_A = -1;

            int i_in_C; for (i_in_C = 0;i_in_C < order_C && idx_B[i] != idx_C[i_in_C];i_in_C++);
            if (i_in_C == order_C) i_in_C = -1;

            if (i_in_A != -1 || i_in_C == -1) continue;

            indices.push_back(index_locator_(0, idx_B[i], i_in_A, i, i_in_C));
        }

        while (!indices.empty())
        {
            std::vector<index_locator_> group;
            group.push_back(indices[0]);
            group.back().sort = 0;
            indices.erase(indices.begin());

            int s = 1;
            for (std::vector<index_locator_>::iterator it = indices.begin();;)
            {
                if (it == indices.end()) break;

                if ((group[0].pos_A == -1 && it->pos_A != -1) ||
                    (group[0].pos_A != -1 && it->pos_A == -1) ||
                    (group[0].pos_B == -1 && it->pos_B != -1) ||
                    (group[0].pos_B != -1 && it->pos_B == -1) ||
                    (group[0].pos_C == -1 && it->pos_C != -1) ||
                    (group[0].pos_C != -1 && it->pos_C == -1))
                {
                    ++it;
                    continue;
                }

                if (group[0].pos_A != -1)
                {
                    bool sym_in_A = false;
                    for (int k = group[0].pos_A-1;k >= 0 && sym_A[k] != NS;k--)
                    {
                        if (idx_A[k] == it->idx)
                        {
                            sym_in_A = true;
                            break;
                        }
                    }
                    for (int k = group[0].pos_A+1;k < order_A && sym_A[k-1] != NS;k++)
                    {
                        if (idx_A[k] == it->idx)
                        {
                            sym_in_A = true;
                            break;
                        }
                    }
                    if (!sym_in_A)
                    {
                        ++it;
                        continue;
                    }
                }

                if (group[0].pos_B != -1)
                {
                    bool sym_in_B = false;
                    for (int k = group[0].pos_B-1;k >= 0 && sym_B[k] != NS;k--)
                    {
                        if (idx_B[k] == it->idx)
                        {
                            sym_in_B = true;
                            break;
                        }
                    }
                    for (int k = group[0].pos_B+1;k < order_B && sym_B[k-1] != NS;k++)
                    {
                        if (idx_B[k] == it->idx)
                        {
                            sym_in_B = true;
                            break;
                        }
                    }
                    if (!sym_in_B)
                    {
                        ++it;
                        continue;
                    }
                }

                if (group[0].pos_C != -1)
                {
                    bool sym_in_C = false;
                    for (int k = group[0].pos_C-1;k >= 0 && sym_C[k] != NS;k--)
                    {
                        if (idx_C[k] == it->idx)
                        {
                            sym_in_C = true;
                            break;
                        }
                    }
                    for (int k = group[0].pos_C+1;k < order_C && sym_C[k-1] != NS;k++)
                    {
                        if (idx_C[k] == it->idx)
                        {
                            sym_in_C = true;
                            break;
                        }
                    }
                    if (!sym_in_C)
                    {
                        ++it;
                        continue;
                    }
                }

                group.push_back(*it);
                group.back().sort = s++;
                it = indices.erase(it);
            }

            if (group.size() <= 1) continue;

            std::vector<int> order_A, order_B, order_C;

            if (group[0].pos_A != -1)
            {
                for (int i = 0;i < (int)group.size();i++)
                    order_A.push_back(group[i].sort);

                if (group[0].pos_B != -1)
                {
                    std::sort(group.begin(), group.end(), index_locator_::sortB);
                    for (int i = 0;i < (int)group.size();i++)
                    {
                        order_B.push_back(group[i].sort);
                        idx_B[group[group[i].sort].pos_B] = group[i].idx;
                    }
                    if (sym_B[group[0].pos_B] == AS)
                        fact *= relative_sign(order_A, order_B);
                }

                if (group[0].pos_C != -1)
                {
                    std::sort(group.begin(), group.end(), index_locator_::sortC);
                    for (int i = 0;i < (int)group.size();i++)
                    {
                        order_C.push_back(group[i].sort);
                        idx_C[group[group[i].sort].pos_C] = group[i].idx;
                    }
                    if (sym_C[group[0].pos_C] == AS)
                        fact *= relative_sign(order_A, order_C);
                }
            }
            else
            {
                for (int i = 0;i < (int)group.size();i++)
                    order_B.push_back(group[i].sort);

                std::sort(group.begin(), group.end(), index_locator_::sortC);
                for (int i = 0;i < (int)group.size();i++)
                {
                    order_C.push_back(group[i].sort);
                    idx_C[group[group[i].sort].pos_C] = group[i].idx;
                }
                if (sym_C[group[0].pos_C] == AS)
                    fact *= relative_sign(order_B, order_C);
            }
        }

        return fact;
    }

    template <typename T>
    int overcounting_factor(int order_A, const T& idx_A, const int* sym_A,
                            int order_B, const T& idx_B, const int* sym_B,
                            int order_C, const T& idx_C, const int* sym_C)
    {
        int fact = 1;

        for (int i = 0;i < order_A;i++)
        {
            int j;
            for (j = 0;j < order_B && idx_A[i] != idx_B[j];j++);
            if (j == order_B) continue;

            int k;
            for (k = 0;k < order_C && idx_A[i] != idx_C[k];k++);
            if (k != order_C) continue;

            int ninarow = 1;
            while (i < order_A &&
                   j < order_B &&
                   sym_A[i] != NS &&
                   sym_B[j] != NS &&
                   idx_A[i] == idx_B[j])
            {
                ninarow++;
                i++;
                j++;
            }
            if (i < order_A &&
                j < order_B &&
                idx_A[i] != idx_B[j]) ninarow--;

            if (ninarow >= 2){
              if (sym_A[i-ninarow+1]!=SY)
                for (;ninarow > 1;ninarow--) fact *= ninarow;
            }
        }

        return fact;
    }

    template <typename T>
    int overcounting_factor(int order_A, const T& idx_A, const int* sym_A,
                            int order_B, const T& idx_B, const int* sym_B)
    {
        int fact;
        int ninarow;
        fact = 1.0;

        for (int i = 0;i < order_A;i++)
        {
            int j;
            ninarow = 0;
            for (j = 0;j < order_B && idx_A[i] != idx_B[j];j++);
            if (j>=order_B){
                ninarow = 1;
                while (sym_A[i] != NS)
                {
                    i++;
                    for (j = 0;j < order_B && idx_A[i] != idx_B[j];j++);
                    if (j>=order_B) ninarow++;
                }
            }
            if (ninarow >= 2){
                if (sym_A[i-ninarow+1]==AS) return 0.0;
                if (sym_A[i-ninarow+1]==SY) {
                    /*printf("CTF error: sum over SY index pair currently not functional, ABORTING\n");
                    assert(0);*/
                }
                if (sym_A[i-ninarow+1]!=SY)
                  for (;ninarow > 1;ninarow--) fact *= ninarow;
            }
        }
        return fact;
    }
}

template <capability_type C, class Base>
class BlockedTensorBase : public Base
{
    static_assert(IS_SUPERSET_OF(C,BOUNDED), "");

    public:
        typedef Tensor<C> SubTensor;

    protected:
        vector<int> alloced;
        vector<int> stride;
        vector<vector<int>> blocklen;
        unique_vector<SubTensor> tensors;
        using TensorInitializer_::R;
        using TensorInitializer<INDEXABLE_>::ndim;
        using TensorInitializer<BOUNDED_>::len;

        void sortIntoBlocks(key_type n, const key_type* keys, const void* values,
                            vector<pair<vector<key_type>,vector<char>>>& kvs) const
        {
            int sz = R.size;

            int nb = 1;
            for (int d = 0;d < ndim;d++)
            {
                nb *= blocklen[d].size();
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
                for (int d = 0;d < ndim;d++)
                {
                    key_type i = k%len[d];
                    k          = k/len[d];

                    int bi = 0;
                    for (;bi < blocklen[d].size() && i >= blocklen[d][bi];
                          bi++, i -= blocklen[d][bi]);
                    assert(bi < blocklen[d].size());

                    b += bs*bi;
                    bs *= blocklen[d].size();
                    bk += bks*i;
                    bks *= blocklen[d][bi];
                }
                assert(b < kvs.size());
                kvs[b].first.push_back(bk);
                kvs[b].second.insert(kvs[b].second.end(),
                                     (const char*)values+ki*sz,
                                     (const char*)values+(ki+1)*sz);
            }
        }

        const SubTensor* get(const vector<int>& idx) const
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

        template <typename T, typename U>
        void sliceByBlock(const Scalar& alpha, bool conja, const vector<vector<T>>& start_A, const TensorImplementation<>& A_,
                          const Scalar&  beta,             const vector<vector<T>>& start_B, const vector<vector<T>>& length,
                          U slice_func)
        {
            auto& A = static_cast<const BlockedTensorBase&>(A_);

            int off_A = 0;
            vector<int> idx(ndim, 0);
            vector<T> start_A_sub(ndim);
            vector<T> start_B_sub(ndim);
            vector<T> len_sub(ndim);
            for (bool done = false;!done;)
            {
                if (  tensors.ptr(off_A) &&   alloced[off_A] &&
                    A.tensors.ptr(off_A) && A.alloced[off_A])
                {
                    for (int i = 0;i < ndim;i++)
                    {
                        start_A_sub[i] = start_A[i][idx[i]];
                        start_B_sub[i] = start_B[i][idx[i]];
                        len_sub[i] = blocklen[i][idx[i]];
                    }

                    tensors[off_A].*slice_func(alpha, conja, start_A_sub, A.tensors[off_A],
                                                beta,        start_B_sub, len_sub);
                }

                for (int i = 0;i < ndim;i++)
                {
                    idx[i]++;
                    off_A += stride[i];

                    int n = blocklen[i].size();
                    if (idx[i] == n)
                    {
                        idx[i] = 0;
                        off_A -= stride[i]*n;
                        if (i == ndim-1) done = true;
                    }
                    else break;
                }

                if (ndim == 0) done = true;
            }
        }

    public:
        template <capability_type C_>
        BlockedTensorBase(const TensorInitializerList<C_>& ilist,
                          const vector<vector<int>>& blocklen)
        : Base(ilist), stride(blocklen.size(), 1), blocklen(blocklen)
        {
            int ndim = blocklen.size();
            for (int i = 1;i < ndim;i++)
                stride[i] = stride[i-1]*blocklen[i-1].size();
            tensors.resize((ndim == 0 ? 1 : stride[ndim-1]*blocklen[ndim-1].size()));
        }

        const vector<vector<int>>& getBlockLengths() const
        {
            return blocklen;
        }

        bool exists(const vector<int>& idx) const
        {
            return get(idx);
        }

        SubTensor& operator()(const vector<int>& idx)
        {
            return *const_cast<SubTensor*>(const_cast<const BlockedTensorBase&>(*this).get(idx));
        }

        const SubTensor& operator()(const vector<int>& idx) const
        {
            return *get(idx);
        }

        void getAllKeys(KeyVector& keys) const
        {
            vector<int> b(ndim);
            KeyVector subkeys;
            for (bool done = false;!done && ndim;)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().getAllKeys(subkeys);
                    keys += subkeys;
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getAllData(KeyValueVector& kv) const
        {
            vector<int> b(ndim);
            KeyValueVector subkv(kv.field());
            for (bool done = false;!done && ndim;)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().getAllData(subkv);
                    kv += subkv;
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getData(key_type n, key_type* keys, void* values) const
        {
            int sz = R.size;

            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && ndim;bi++)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().getData(kvs[bi].first.size(),
                                       kvs[bi].first.data(),
                                       kvs[bi].second.data());
                }

                copy(kvs[bi].first.begin(), kvs[bi].first.end(), keys+i);
                copy(kvs[bi].second.begin(), kvs[bi].second.end(), (char*)values+i*sz);
                i += kvs[bi].first.size();

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void setData(key_type n, const key_type* keys, const void* values)
        {
            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && ndim;bi++)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().setData(kvs[bi].first.size(),
                                       kvs[bi].first.data(),
                                       kvs[bi].second.data());
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
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

            vector<int> b(ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && ndim;bi++)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<BOUNDED_>().addData(kvs[bi].first.size(), alpha,
                                       kvs[bi].first.data(),
                                       kvs[bi].second.data(), beta);
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getLocalKeys(KeyVector& keys) const
        {
            int sz = R.size;

            vector<int> b(ndim);
            KeyVector subkeys;
            for (bool done = false;!done && ndim;)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().getLocalKeys(subkeys);
                    keys += subkeys;
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getLocalData(KeyValueVector& kv) const
        {
            int sz = R.size;

            vector<int> b(ndim);
            KeyValueVector subkv(kv.field());
            for (bool done = false;!done && ndim;)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().getLocalData(subkv);
                    kv += subkv;
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void setLocalData(key_type n, const key_type* keys, const void* values)
        {
            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && ndim;bi++)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().setLocalData(kvs[bi].first.size(),
                                            kvs[bi].first.data(),
                                            kvs[bi].second.data());
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
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

            vector<int> b(ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && ndim;bi++)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().addLocalData(kvs[bi].first.size(), alpha,
                                            kvs[bi].first.data(),
                                            kvs[bi].second.data(), beta);
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getRemoteData(key_type n, key_type* keys, void* values) const
        {
            int sz = R.size;

            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && ndim;bi++)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().getRemoteData(kvs[bi].first.size(),
                                             kvs[bi].first.data(),
                                             kvs[bi].second.data());
                }

                copy(kvs[bi].first.begin(), kvs[bi].first.end(), keys+i);
                copy(kvs[bi].second.begin(), kvs[bi].second.end(), (char*)values+i*sz);
                i += kvs[bi].first.size();

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void setRemoteData(key_type n, const key_type* keys, const void* values)
        {
            vector<pair<vector<key_type>,vector<char>>> kvs;
            sortIntoBlocks(n, keys, values, kvs);

            vector<int> b(ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && ndim;bi++)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().setRemoteData(kvs[bi].first.size(),
                                             kvs[bi].first.data(),
                                             kvs[bi].second.data());
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
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

            vector<int> b(ndim);
            int bi = 0;
            key_type i = 0;
            for (bool done = false;!done && ndim;bi++)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().addRemoteData(kvs[bi].first.size(), alpha,
                                             kvs[bi].first.data(),
                                             kvs[bi].second.data(), beta);
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void getRemoteData() const
        {
            vector<int> b(ndim);
            for (bool done = false;!done && ndim;)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().getRemoteData();
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void setRemoteData()
        {
            vector<int> b(ndim);
            for (bool done = false;!done && ndim;)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().setRemoteData();
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void addRemoteData(const Scalar& alpha, const Scalar& beta)
        {
            vector<int> b(ndim);
            for (bool done = false;!done && ndim;)
            {
                if (exists(b))
                {
                    (*this)(b).template impl<DISTRIBUTED_>().addRemoteData(alpha, beta);
                }

                for (int d = 0;d < ndim;d++)
                {
                    b[d]++;
                    if (b[d] == blocklen[d].size())
                    {
                        b[d] = 0;
                        if (d == ndim) done = true;
                    }
                    else break;
                }
            }
        }

        void scale(const Scalar& alpha, const string& idx_A)
        {
            string inds = uniqued(idx_A);
            vector<int> stride_A = detail::getStrides(inds, stride, idx_A);

            int m = inds.size();
            int off_A = 0;
            vector<int> idx(m, 0);
            for (bool done = false;!done;)
            {
                if (alloced[off_A])
                {
                    tensors[off_A][idx_A] *= alpha;
                }

                for (int i = 0;i < m;i++)
                {
                    idx[i]++;
                    off_A += stride_A[i];

                    int n = blocklen[i].size();
                    if (idx[i] == n)
                    {
                        idx[i] = 0;
                        off_A -= stride_A[i]*n;
                        if (i == m-1) done = true;
                    }
                    else break;
                }

                if (m == 0) done = true;
            }
        }

        void div(const Scalar& alpha, bool conja, const TensorImplementation<>& A_,
                                      bool conjb, const TensorImplementation<>& B_,
                 const Scalar& beta)
        {
            auto& A = static_cast<const BlockedTensorBase&>(A_);
            auto& B = static_cast<const BlockedTensorBase&>(B_);

            int off_A = 0;
            vector<int> idx(ndim, 0);
            for (bool done = false;!done;)
            {
                if (tensors.ptr(off_A) && alloced[off_A])
                {
                    beta*tensors[off_A] += alpha*A.tensors[off_A]/B.tensors[off_A];
                }

                for (int i = 0;i < ndim;i++)
                {
                    idx[i]++;
                    off_A += stride[i];

                    int n = blocklen[i].size();
                    if (idx[i] == n)
                    {
                        idx[i] = 0;
                        off_A -= stride[i]*n;
                        if (i == ndim-1) done = true;
                    }
                    else break;
                }

                if (ndim == 0) done = true;
            }
        }

        void invert(const Scalar& alpha, bool conja, const TensorImplementation<>& A_,
                    const Scalar& beta)
        {
            auto& A = static_cast<const BlockedTensorBase&>(A_);

            int off_A = 0;
            vector<int> idx(ndim, 0);
            for (bool done = false;!done;)
            {
                if (tensors.ptr(off_A) && alloced[off_A])
                {
                    beta*tensors[off_A] += alpha/A.tensors[off_A];
                }

                for (int i = 0;i < ndim;i++)
                {
                    idx[i]++;
                    off_A += stride[i];

                    int n = blocklen[i].size();
                    if (idx[i] == n)
                    {
                        idx[i] = 0;
                        off_A -= stride[i]*n;
                        if (i == ndim-1) done = true;
                    }
                    else break;
                }

                if (ndim == 0) done = true;
            }
        }
};

template <capability_type C, class Base, typename=void>
class BlockedTensor;

template <capability_type C, class Base>
class BlockedTensor<C, Base, enable_if_t<C&IPSYMMETRIC_>> : public BlockedTensorBase<C, Base>
{
    protected:
        vector<double> factor;
        vector<vector<int>> reorder;
        using BlockedTensorBase<C, Base>::alloced;
        using BlockedTensorBase<C, Base>::stride;
        using BlockedTensorBase<C, Base>::blocklen;
        using BlockedTensorBase<C, Base>::tensors;
        using TensorInitializer<INDEXABLE_>::ndim;
        using TensorInitializer<IPSYMMETRIC_>::sym;
        using typename BlockedTensorBase<C, Base>::SubTensor;

        void put(const vector<int>& idx, Tensor<C>&& t)
        {
            assert(idx.size() == blocklen.size());

            vector<int> r = range(ndim);
            vector<int> r0 = r;
            const vector<int>& subsym = t.getSymmetry();

            for (bool done = false, first = true;!done;first = false)
            {
                int off = 0;
                for (int i = 0;i < idx.size();i++)
                {
                    assert(idx[i] >= 0 && idx[i] < blocklen[i].size());
                    off += stride[i]*idx[i];
                }

                //TODO: wrong!!!

                if (first) alloced[off] = true;
                tensors.ptr(off).reset(new SubTensor(t));
                reorder[off] = r;
                factor[off] = detail::relative_sign(r, r0);

                for (int i = 0;;)
                {
                    int j = i+1; for (;j < ndim && subsym[j] != NS;j++);

                    if (j == i+1 || !next_permutation(r.begin()+i, r.begin()+j))
                    {
                        if (j == ndim) done = true;
                    }
                    else break;

                    i = j;
                }
            }
        }

    public:
        template <capability_type C_>
        BlockedTensor(const TensorInitializerList<C_>& ilist,
                      const vector<vector<int>>& blocklen)
        : BlockedTensorBase<C, Base>(ilist, blocklen)
        {
            factor.resize(tensors.size());
            reorder.resize(tensors.size());
        }

        void mult(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idx_A,
                                       bool conjb, const TensorImplementation<>& B_, const string& idx_B,
                  const Scalar&  beta,                                               const string& idx_C)
        {
            auto& A = static_cast<const BlockedTensor&>(A_);
            auto& B = static_cast<const BlockedTensor&>(B_);

            string idx_A_(idx_A);
            string idx_B_(idx_B);
            string idx_C_(idx_C);
            string idx_A__(idx_A);
            string idx_B__(idx_B);
            string idx_C__(idx_C);

            double f1 = detail::align_symmetric_indices(A.ndim, idx_A_, A.sym.data(),
                                                        B.ndim, idx_B_, B.sym.data(),
                                                          ndim, idx_C_,   sym.data());
            f1 *= detail::overcounting_factor(A.ndim, idx_A_, A.sym.data(),
                                              B.ndim, idx_B_, B.sym.data(),
                                                ndim, idx_C_,   sym.data());

            vector<char> inds; inds.reserve(A.ndim+B.ndim+ndim);
            vector<int> syms; syms.reserve(A.ndim+B.ndim+ndim);
            vector<int> nb; nb.reserve(A.ndim+B.ndim+ndim);

            for (int i = 0;i < A.ndim;i++)
            {
                if (!contains(inds, idx_A_[i]))
                {
                    nb.push_back(A.blocklen[i].size());
                    inds.push_back(idx_A_[i]);
                    syms.push_back(detail::checkAlignment(i, idx_A_, A.sym, idx_B_, B.sym) &&
                                   detail::checkAlignment(i, idx_A_, A.sym, idx_C_,   sym) ? AS : NS);
                }
            }
            for (int i = 0;i < B.ndim;i++)
            {
                if (!contains(inds, idx_B_[i]))
                {
                    nb.push_back(B.blocklen[i].size());
                    inds.push_back(idx_B_[i]);
                    syms.push_back(detail::checkAlignment(i, idx_B_, B.sym, idx_C_, sym) ? AS : NS);
                }
            }
            for (int i = 0;i < ndim;i++)
            {
                if (!contains(inds, idx_C_[i]))
                {
                    nb.push_back(blocklen[i].size());
                    inds.push_back(idx_C_[i]);
                    syms.push_back(sym[i] != NS ? AS: NS);
                }
            }

            vector<int> stride_A = detail::getStrides(inds, A.stride, idx_A_);
            vector<int> stride_B = detail::getStrides(inds, B.stride, idx_B_);
            vector<int> stride_C = detail::getStrides(inds,   stride, idx_C_);

            vector<Scalar> beta_(tensors.size(), beta);

            int m = inds.size();
            int off_A = 0;
            int off_B = 0;
            int off_C = 0;
            vector<int> idx(m, 0);
            for (bool done = false;!done;)
            {
                if ((A.tensors.ptr(off_A) && B.tensors.ptr(off_B) && tensors.ptr(off_C)) &&
                    (A.alloced[off_A] || B.alloced[off_B] || alloced[off_C]))
                {
                    double f2 = detail::overcounting_factor(A.ndim, idx_A_, A.tensors[off_A].getSymmetry().data(),
                                                            B.ndim, idx_B_, B.tensors[off_B].getSymmetry().data(),
                                                              ndim, idx_C_,   tensors[off_C].getSymmetry().data());
                    double f3 = A.factor[off_A]*B.factor[off_B]*factor[off_C];

                    for (int i = 0;i < A.ndim;i++) idx_A__[i] = idx_A_[A.reorder[off_A][i]];
                    for (int i = 0;i < B.ndim;i++) idx_B__[i] = idx_B_[B.reorder[off_B][i]];
                    for (int i = 0;i <   ndim;i++) idx_C__[i] = idx_C_[  reorder[off_C][i]];

                    assert(off_A >= 0 && off_A < A.tensors.size());
                    assert(off_B >= 0 && off_B < B.tensors.size());
                    assert(off_C >= 0 && off_C <   tensors.size());
                    beta_[off_C]*tensors[off_C][idx_C__] +=
                        (alpha*f1*f3/f2)*aquarius::conj(A.tensors[off_A][idx_A__], conja)*
                                         aquarius::conj(B.tensors[off_B][idx_B__], conjb);

                    beta_[off_C] = 1.0;
                }

                for (int i = 0;i < m;i++)
                {
                    idx[i]++;
                    off_A += stride_A[i];
                    off_B += stride_B[i];
                    off_C += stride_C[i];

                    if (idx[i] == (syms[i] == NS ? nb[i] : idx[i+1]+1))
                    {
                        int lower = (i == 0 || syms[i-1] == NS ? 0 : idx[i-1]);
                        off_A -= stride_A[i]*(idx[i]-lower);
                        off_B -= stride_B[i]*(idx[i]-lower);
                        off_C -= stride_C[i]*(idx[i]-lower);
                        idx[i] = lower;
                        if (i == m-1) done = true;
                    }
                    else break;
                }

                if (m == 0) done = true;
            }
        }

        void sum(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idx_A,
                 const Scalar&  beta,                                               const string& idx_B)
        {
            auto& A = static_cast<const BlockedTensor&>(A_);

            string idx_A_(idx_A);
            string idx_B_(idx_B);
            string idx_A__(idx_A);
            string idx_B__(idx_B);

            double f = detail::align_symmetric_indices(A.ndim, idx_A_, A.sym.data(),
                                                         ndim, idx_B_,   sym.data());

            vector<char> inds; inds.reserve(A.ndim+ndim);
            vector< int> syms; syms.reserve(A.ndim+ndim);
            vector< int> nb; nb.reserve(A.ndim+ndim);

            for (int i = 0;i < A.ndim;i++)
            {
                if (!contains(inds, idx_A_[i]))
                {
                    nb.push_back(A.blocklen[i].size());
                    inds.push_back(idx_A_[i]);
                    syms.push_back(detail::checkAlignment(i, idx_A_, A.sym, idx_B_, sym) ? AS : NS);
                }
            }
            for (int i = 0;i < ndim;i++)
            {
                if (!contains(inds, idx_B_[i]))
                {
                    nb.push_back(blocklen[i].size());
                    inds.push_back(idx_B_[i]);
                    syms.push_back(sym[i] != NS ? AS : NS);
                }
            }

            vector<int> stride_A = detail::getStrides(inds, A.stride, idx_A_);
            vector<int> stride_B = detail::getStrides(inds,   stride, idx_B_);

            vector<Scalar> beta_(tensors.size(), beta);

            int m = inds.size();
            int off_A = 0;
            int off_B = 0;
            vector<int> idx(m, 0);
            for (bool done = false;!done;)
            {
                if ((A.tensors.ptr(off_A) && tensors.ptr(off_B)) &&
                    (A.alloced[off_A] || alloced[off_B]))
                {
                    double f3 = A.factor[off_A]*factor[off_B];

                    for (int i = 0;i < A.ndim;i++) idx_A__[i] = idx_A_[A.reorder[off_A][i]];
                    for (int i = 0;i <   ndim;i++) idx_B__[i] = idx_B_[  reorder[off_B][i]];

                    assert(off_A >= 0 && off_A < A.tensors.size());
                    assert(off_B >= 0 && off_B <   tensors.size());
                    beta_[off_B]*tensors[off_B][idx_B__] +=
                        (alpha*f*f3)*aquarius::conj(A.tensors[off_A][idx_A__], conja);

                    beta_[off_B] = 1.0;
                }

                for (int i = 0;i < m;i++)
                {
                    idx[i]++;
                    off_A += stride_A[i];
                    off_B += stride_B[i];

                    if (idx[i] == (syms[i] == NS ? nb[i] : idx[i+1]+1))
                    {
                        int lower = (i == 0 || syms[i-1] == NS ? 0 : idx[i-1]);
                        off_A -= stride_A[i]*(idx[i]-lower);
                        off_B -= stride_B[i]*(idx[i]-lower);
                        idx[i] = lower;
                        if (i == m-1) done = true;
                    }
                    else break;
                }

                if (m == 0) done = true;
            }
        }

        Scalar dot(bool conja, const TensorImplementation<>& A_, const string& idx_A,
                   bool conjb,                                   const string& idx_B) const
        {
            auto& A = static_cast<const BlockedTensor&>(A_);

            string idx_A_(idx_A);
            string idx_B_(idx_B);
            string idx_A__(idx_A);
            string idx_B__(idx_B);

            double f1 = detail::align_symmetric_indices(A.ndim, idx_A_, A.sym.data(),
                                                          ndim, idx_B_,   sym.data());
            f1 *= detail::overcounting_factor(A.ndim,   idx_A_, A.sym.data(),
                                                ndim,   idx_B_,   sym.data(),
                                                   0, string(),         NULL);

            vector<char> inds; inds.reserve(A.ndim+ndim);
            vector< int> syms; syms.reserve(A.ndim+ndim);
            vector< int> nb; nb.reserve(A.ndim+ndim);

            for (int i = 0;i < A.ndim;i++)
            {
                if (!contains(inds, idx_A_[i]))
                {
                    nb.push_back(A.blocklen[i].size());
                    inds.push_back(idx_A_[i]);
                    syms.push_back(detail::checkAlignment(i, idx_A_, A.sym, idx_B_, sym) ? AS : NS);
                }
            }
            for (int i = 0;i < ndim;i++)
            {
                if (!contains(inds, idx_B_[i]))
                {
                    nb.push_back(blocklen[i].size());
                    inds.push_back(idx_B_[i]);
                    syms.push_back(sym[i] != NS ? AS : NS);
                }
            }

            vector<int> stride_A = detail::getStrides(inds, A.stride, idx_A_);
            vector<int> stride_B = detail::getStrides(inds,   stride, idx_B_);

            Scalar total(this->F);

            int m = inds.size();
            int off_A = 0;
            int off_B = 0;
            vector<int> idx(m, 0);
            for (bool done = false;!done;)
            {
                if ((A.tensors.ptr(off_A) && tensors.ptr(off_B)) &&
                    (A.alloced[off_A] || alloced[off_B]))
                {
                    double f2 = detail::overcounting_factor(A.ndim,   idx_A_, tensors[off_A].getSymmetry().data(),
                                                              ndim,   idx_B_, tensors[off_B].getSymmetry().data(),
                                                                 0, string(),                                NULL);
                    double f3 = A.factor[off_A]*factor[off_B];

                    for (int i = 0;i < A.ndim;i++) idx_A__[i] = idx_A_[A.reorder[off_A][i]];
                    for (int i = 0;i <   ndim;i++) idx_B__[i] = idx_B_[  reorder[off_B][i]];

                    assert(off_A >= 0 && off_A < A.tensors.size());
                    assert(off_B >= 0 && off_B <   tensors.size());
                    total += (f1*f3/f2)*scalar(aquarius::conj(A.tensors[off_A][idx_A__], conja)*
                                               aquarius::conj(  tensors[off_B][idx_B__], conjb));
                }

                for (int i = 0;i < m;i++)
                {
                    idx[i]++;
                    off_A += stride_A[i];
                    off_B += stride_B[i];

                    if (idx[i] == (syms[i] == NS ? nb[i] : idx[i+1]+1))
                    {
                        int lower = (i == 0 || syms[i-1] == NS ? 0 : idx[i-1]);
                        off_A -= stride_A[i]*(idx[i]-lower);
                        off_B -= stride_B[i]*(idx[i]-lower);
                        idx[i] = lower;
                        if (i == m-1) done = true;
                    }
                    else break;
                }

                if (m == 0) done = true;
            }

            return total;
        }

        Scalar norm(int p) const
        {
            Scalar nrm(this->F);

            int off_A = 0;
            vector<int> idx(ndim, 0);
            for (bool done = false;!done;)
            {
                if (tensors.ptr(off_A) && alloced[off_A])
                {
                    double factor = 1;
                    const vector<int>& subsym = tensors[off_A].getSymmetry();
                    for (int i = 0;i < ndim;)
                    {
                        int j; for (j = i;sym[j] != NS;j++) continue; j++;

                        int m = j-i;
                        for (int k = i;k < j;)
                        {
                            int l; for (l = k;subsym[l] != NS;l++) continue; l++;
                            int o = l-k;
                            factor *= binom(m,o);
                            m -= o;
                            k = l;
                        }

                        i = j;
                    }

                    Scalar subnrm = tensors[off_A].norm(p);

                    if (p == 2)
                    {
                        nrm += factor*subnrm*subnrm;
                    }
                    else if (p == 0)
                    {
                        nrm = max(nrm,subnrm);
                    }
                    else if (p == 1)
                    {
                        nrm += factor*subnrm;
                    }
                }

                for (int i = 0;i < ndim;i++)
                {
                    idx[i]++;
                    off_A += stride[i];

                    int n = blocklen[i].size();
                    if (idx[i] == n)
                    {
                        idx[i] = 0;
                        off_A -= stride[i]*n;
                        if (i == ndim-1) done = true;
                    }
                    else break;
                }

                if (ndim == 0) done = true;
            }

            if (p == 2) nrm = sqrt(nrm);

            return nrm;
        }
};

template <capability_type C, class Base>
class BlockedTensor<C, Base, enable_if_t<!C&IPSYMMETRIC_>> : public BlockedTensorBase<C, Base>
{
    protected:
        using BlockedTensorBase<C, Base>::alloced;
        using BlockedTensorBase<C, Base>::stride;
        using BlockedTensorBase<C, Base>::blocklen;
        using BlockedTensorBase<C, Base>::tensors;
        using TensorInitializer<INDEXABLE_>::ndim;
        using typename BlockedTensorBase<C, Base>::SubTensor;

        void put(const vector<int>& idx, Tensor<C>&& t)
        {
            assert(idx.size() == blocklen.size());

            int off = 0;
            for (int i = 0;i < idx.size();i++)
            {
                assert(idx[i] >= 0 && idx[i] < blocklen[i].size());
                off += stride[i]*idx[i];
            }

            alloced[off] = true;
            tensors.ptr(off).reset(new SubTensor(t));
        }

    public:
        template <capability_type C_>
        BlockedTensor(const TensorInitializerList<C_>& ilist,
                      const vector<vector<int>>& blocklen)
        : BlockedTensorBase<C, Base>(ilist, blocklen) {}

        void mult(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idx_A,
                                       bool conjb, const TensorImplementation<>& B_, const string& idx_B,
                  const Scalar&  beta,                                               const string& idx_C)
        {
            auto& A = static_cast<const BlockedTensor&>(A_);
            auto& B = static_cast<const BlockedTensor&>(B_);

            vector<char> inds; inds.reserve(A.ndim+B.ndim+ndim);
            vector<int> nb; nb.reserve(A.ndim+B.ndim+ndim);

            for (int i = 0;i < A.ndim;i++)
            {
                if (!contains(inds, idx_A[i]))
                {
                    nb.push_back(A.blocklen[i].size());
                    inds.push_back(idx_A[i]);
                }
            }
            for (int i = 0;i < B.ndim;i++)
            {
                if (!contains(inds, idx_B[i]))
                {
                    nb.push_back(B.blocklen[i].size());
                    inds.push_back(idx_B[i]);
                }
            }
            for (int i = 0;i < ndim;i++)
            {
                if (!contains(inds, idx_C[i]))
                {
                    nb.push_back(blocklen[i].size());
                    inds.push_back(idx_C[i]);
                }
            }

            vector<int> stride_A = detail::getStrides(inds, A.stride, idx_A);
            vector<int> stride_B = detail::getStrides(inds, B.stride, idx_B);
            vector<int> stride_C = detail::getStrides(inds,   stride, idx_C);

            vector<Scalar> beta_(tensors.size(), beta);

            int m = inds.size();
            int off_A = 0;
            int off_B = 0;
            int off_C = 0;
            vector<int> idx(m, 0);
            for (bool done = false;!done;)
            {
                if ((A.tensors.ptr(off_A) && B.tensors.ptr(off_B) && tensors.ptr(off_C)) &&
                    (A.alloced[off_A] || B.alloced[off_B] || alloced[off_C]))
                {
                    assert(off_A >= 0 && off_A < A.tensors.size());
                    assert(off_B >= 0 && off_B < B.tensors.size());
                    assert(off_C >= 0 && off_C <   tensors.size());
                    beta_[off_C]*tensors[off_C][idx_C] +=
                        alpha*aquarius::conj(A.tensors[off_A][idx_A], conja)*
                              aquarius::conj(B.tensors[off_B][idx_B], conjb);

                    beta_[off_C] = 1.0;
                }

                for (int i = 0;i < m;i++)
                {
                    idx[i]++;
                    off_A += stride_A[i];
                    off_B += stride_B[i];
                    off_C += stride_C[i];

                    if (idx[i] == nb[i])
                    {
                        off_A -= stride_A[i]*idx[i];
                        off_B -= stride_B[i]*idx[i];
                        off_C -= stride_C[i]*idx[i];
                        idx[i] = 0;
                        if (i == m-1) done = true;
                    }
                    else break;
                }

                if (m == 0) done = true;
            }
        }

        void sum(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idx_A,
                 const Scalar&  beta,                                               const string& idx_B)
        {
            auto& A = static_cast<const BlockedTensor&>(A_);

            vector<char> inds; inds.reserve(A.ndim+ndim);
            vector< int> nb; nb.reserve(A.ndim+ndim);

            for (int i = 0;i < A.ndim;i++)
            {
                if (!contains(inds, idx_A[i]))
                {
                    nb.push_back(A.blocklen[i].size());
                    inds.push_back(idx_A[i]);
                }
            }
            for (int i = 0;i < ndim;i++)
            {
                if (!contains(inds, idx_B[i]))
                {
                    nb.push_back(blocklen[i].size());
                    inds.push_back(idx_B[i]);
                }
            }

            vector<int> stride_A = detail::getStrides(inds, A.stride, idx_A);
            vector<int> stride_B = detail::getStrides(inds,   stride, idx_B);

            vector<Scalar> beta_(tensors.size(), beta);

            int m = inds.size();
            int off_A = 0;
            int off_B = 0;
            vector<int> idx(m, 0);
            for (bool done = false;!done;)
            {
                if ((A.tensors.ptr(off_A) && tensors.ptr(off_B)) &&
                    (A.alloced[off_A] || alloced[off_B]))
                {
                    assert(off_A >= 0 && off_A < A.tensors.size());
                    assert(off_B >= 0 && off_B <   tensors.size());
                    beta_[off_B]*tensors[off_B][idx_B] +=
                        alpha*aquarius::conj(A.tensors[off_A][idx_A], conja);

                    beta_[off_B] = 1.0;
                }

                for (int i = 0;i < m;i++)
                {
                    idx[i]++;
                    off_A += stride_A[i];
                    off_B += stride_B[i];

                    if (idx[i] == nb[i])
                    {
                        off_A -= stride_A[i]*idx[i];
                        off_B -= stride_B[i]*idx[i];
                        idx[i] = 0;
                        if (i == m-1) done = true;
                    }
                    else break;
                }

                if (m == 0) done = true;
            }
        }

        Scalar dot(bool conja, const TensorImplementation<>& A_, const string& idx_A,
                   bool conjb,                                   const string& idx_B) const
        {
            auto& A = static_cast<const BlockedTensor&>(A_);

            vector<char> inds; inds.reserve(A.ndim+ndim);
            vector< int> nb; nb.reserve(A.ndim+ndim);

            for (int i = 0;i < A.ndim;i++)
            {
                if (!contains(inds, idx_A[i]))
                {
                    nb.push_back(A.blocklen[i].size());
                    inds.push_back(idx_A[i]);
                }
            }
            for (int i = 0;i < ndim;i++)
            {
                if (!contains(inds, idx_B[i]))
                {
                    nb.push_back(blocklen[i].size());
                    inds.push_back(idx_B[i]);
                }
            }

            vector<int> stride_A = detail::getStrides(inds, A.stride, idx_A);
            vector<int> stride_B = detail::getStrides(inds,   stride, idx_B);

            Scalar total(this->F);

            int m = inds.size();
            int off_A = 0;
            int off_B = 0;
            vector<int> idx(m, 0);
            for (bool done = false;!done;)
            {
                if ((A.tensors.ptr(off_A) && tensors.ptr(off_B)) &&
                    (A.alloced[off_A] || alloced[off_B]))
                {
                    assert(off_A >= 0 && off_A < A.tensors.size());
                    assert(off_B >= 0 && off_B <   tensors.size());
                    total += scalar(aquarius::conj(A.tensors[off_A][idx_A], conja)*
                                    aquarius::conj(  tensors[off_B][idx_B], conjb));
                }

                for (int i = 0;i < m;i++)
                {
                    idx[i]++;
                    off_A += stride_A[i];
                    off_B += stride_B[i];

                    if (idx[i] == nb[i])
                    {
                        off_A -= stride_A[i]*idx[i];
                        off_B -= stride_B[i]*idx[i];
                        idx[i] = 0;
                        if (i == m-1) done = true;
                    }
                    else break;
                }

                if (m == 0) done = true;
            }

            return total;
        }

        Scalar norm(int p) const
        {
            Scalar nrm(this->F);

            int off_A = 0;
            vector<int> idx(ndim, 0);
            for (bool done = false;!done;)
            {
                if (tensors.ptr(off_A) && alloced[off_A])
                {
                    Scalar subnrm = tensors[off_A].norm(p);

                    if (p == 2)
                    {
                        nrm += subnrm*subnrm;
                    }
                    else if (p == 0)
                    {
                        nrm = max(nrm,subnrm);
                    }
                    else if (p == 1)
                    {
                        nrm += subnrm;
                    }
                }

                for (int i = 0;i < ndim;i++)
                {
                    idx[i]++;
                    off_A += stride[i];

                    int n = blocklen[i].size();
                    if (idx[i] == n)
                    {
                        idx[i] = 0;
                        off_A -= stride[i]*n;
                        if (i == ndim-1) done = true;
                    }
                    else break;
                }

                if (ndim == 0) done = true;
            }

            if (p == 2) nrm = sqrt(nrm);

            return nrm;
        }
};

}
}

#endif
