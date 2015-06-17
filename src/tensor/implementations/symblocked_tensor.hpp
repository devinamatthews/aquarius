#ifndef _AQUARIUS_TENSOR_SYMBLOCKED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_SYMBLOCKED_TENSOR_HPP_

#include "util/global.hpp"

#include "symmetry/symmetry.hpp"
#include "task/task.hpp"
#include "tensor/tensor.hpp"

#include "blocked_tensor.hpp"

namespace aquarius
{
namespace tensor
{

namespace detail
{
    inline vector<int> getStrides(const string& indices, int ndim, int len, const string& idx_A)
    {
        vector<int> strides(indices.size(), 0);
        vector<int> stride_A(ndim);

        if (ndim > 0) stride_A[0] = 1;
        for (int i = 1;i < ndim;i++)
        {
            stride_A[i] = stride_A[i-1]*len;
        }

        for (int i = 0;i < indices.size();i++)
        {
            for (int j = 0;j < ndim;j++)
            {
                if (indices[i] == idx_A[j])
                {
                    strides[i] += stride_A[j];
                }
            }
        }

        return strides;
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
        vector<int> alloced;
        vector<double> factor;
        vector<vector<int>> reorder;
        using BlockedTensor<C&~PGSYMMETRIC,TensorImplementation<C>>::tensors;

        template <typename SubFactory>
        void initialize(const INITIALIZER_TYPE(C)& init, SubFactory f)
        {
            INITIALIZER_TYPE(C&~(PGSYMMETRIC|BOUNDED|IPSYMMETRIC)) subinit(init);

            int n = this->group.getNumIrreps();
            vector<symmetry::Representation> irreps;
            for (int i = 0;i < n;i++) irreps.push_back(this->group.getIrrep(i));

            int ntensors = 1;
            vector<int> sublen(this->ndim);

            for (int i = 0;i < this->ndim;i++)
            {
                ntensors *= n;
                sublen[i] = this->len_per_irrep[i][0];
            }

            alloced.resize(ntensors);
            factor.resize(ntensors, 1.0);
            reorder.resize(ntensors);

            int t = 0;
            vector<int> idx(this->ndim);
            vector<symmetry::Representation> prod(this->ndim+1, this->rep);
            for (bool done = false;!done;t++)
            {
                if (prod[0].isTotallySymmetric())
                {
                    reorder[t].resize(this->ndim);
                    vector<int> subsym(this->sym);

                    bool ok = true;
                    for (int i = 0;i < this->ndim-1;i++)
                    {
                        reorder[t][i] = i;
                        if (this->sym[i] != NS)
                        {
                            if (idx[i] < idx[i+1]) subsym[i] = NS;
                            else if (idx[i] > idx[i+1]) ok = false;
                        }
                    }
                    if (this->ndim > 0) reorder[t][this->ndim-1] = this->ndim-1;

                    assert(t < ntensors);
                    if (ok)
                    {
                        tensors.ptr(t).reset(new Tensor<C&~PGSYMMETRIC>(f(subinit <<
                                               TensorInitializer<BOUNDED|IPSYMMETRIC>(sublen, subsym))));
                        alloced[t] = true;
                    }
                }

                for (int i = 0;i < this->ndim;i++)
                {
                    idx[i] = (idx[i] == n-1 ? 0 : idx[i]+1);
                    sublen[i] = this->len_per_irrep[i][idx[i]];

                    if (idx[i] != 0)
                    {
                        for (int j = i;j >= 0;j--)
                        {
                            prod[j] = prod[j+1];
                            prod[j] *= irreps[idx[j]];
                        }
                        break;
                    }
                    else if (i == this->ndim-1)
                    {
                        done = true;
                    }
                }

                if (this->ndim == 0) done = true;
            }

            vector<int> norder = range(this->ndim);

            t = 0;
            idx.assign(this->ndim, 0);
            prod.assign(this->ndim+1, this->rep);
            for (bool done = false;!done;t++)
            {
                if (prod[0].isTotallySymmetric() && !tensors.ptr(t))
                {
                    vector<int> idxreal(idx);
                    for (int i = 0;i < this->ndim;)
                    {
                        int j; for (j = i;j < this->ndim && this->sym[j] != NS;j++); j++;
                        cosort(   idxreal.begin()+i,    idxreal.begin()+j,
                               reorder[t].begin()+i, reorder[t].begin()+j);
                        if (this->sym[i] == AS)
                        {
                            factor[t] *= relative_sign(    norder.begin()+i,     norder.begin()+j,
                                                       reorder[t].begin()+i, reorder[t].begin()+j);
                        }
                        i = j;
                    }

                    vector<int> invorder(this->ndim);
                    for (int i = 0;i < this->ndim;i++)
                    {
                        int j; for (j = 0;j < this->ndim && reorder[t][j] != i;j++);
                        invorder[i] = j;
                    }
                    //swap(invorder, reorder[t]);

                    int treal = 0;
                    int stride = 1;
                    for (int i = 0;i < this->ndim;i++)
                    {
                        treal += stride*idxreal[i];
                        stride *= n;
                    }

                    assert(t < ntensors && treal < ntensors);
                    assert(tensors.ptr(treal));
                    tensors.ptr(t).reset(new Tensor<C&~PGSYMMETRIC>(tensors[treal]));
                }

                for (int i = 0;i < this->ndim;i++)
                {
                    idx[i] = (idx[i] == n-1 ? 0 : idx[i]+1);
                    sublen[i] = this->len_per_irrep[i][idx[i]];

                    if (idx[i] != 0)
                    {
                        for (int j = i;j >= 0;j--)
                        {
                            prod[j] = prod[j+1];
                            prod[j] *= irreps[idx[j]];
                        }
                        break;
                    }
                    else if (i == this->ndim-1)
                    {
                        done = true;
                    }
                }

                if (this->ndim == 0) done = true;
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

        SymmetryBlockedTensor(const INITIALIZER_TYPE(C)& init)
        : BlockedTensor<C&~PGSYMMETRIC,TensorImplementation<C>>(init, init.template as<PGSYMMETRIC>().len_per_irrep)
        {
            initialize(init, typename Tensor<C&~PGSYMMETRIC>::Factory());
        }

        template <typename SubFactory>
        SymmetryBlockedTensor(const INITIALIZER_TYPE(C)& init, SubFactory f)
        : BlockedTensor<C&~PGSYMMETRIC,TensorImplementation<C>>(init, init.template as<PGSYMMETRIC>().len_per_irrep)
        {
            initialize(init, f);
        }

        void mult(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idx_A,
                                       bool conjb, const TensorImplementation<>& B_, const string& idx_B,
                  const Scalar&  beta,                                               const string& idx_C)
        {
            auto& A = static_cast<const SymmetryBlockedTensor&>(A_);
            auto& B = static_cast<const SymmetryBlockedTensor&>(B_);

            assert(this->group == A.group);
            assert(this->group == B.group);

            int n = this->group.getNumIrreps();

            string idx_A_(idx_A);
            string idx_B_(idx_B);
            string idx_C_(idx_C);
            string idx_A__(idx_A);
            string idx_B__(idx_B);
            string idx_C__(idx_C);

            double f1 = detail::align_symmetric_indices(A.ndim, idx_A_, A.sym.data(),
                                                        B.ndim, idx_B_, B.sym.data(),
                                                          this->ndim, idx_C_,   this->sym.data());
            f1 *= detail::overcounting_factor(A.ndim, idx_A_, A.sym.data(),
                                              B.ndim, idx_B_, B.sym.data(),
                                                this->ndim, idx_C_,   this->sym.data());

            vector<int> stride_idx_A = detail::getStrides(idx_A_, A.ndim, n, idx_A_);
            vector<int> stride_idx_B = detail::getStrides(idx_B_, B.ndim, n, idx_B_);
            vector<int> stride_idx_C = detail::getStrides(idx_C_,   this->ndim, n, idx_C_);

            vector<char> inds(A.ndim+B.ndim+this->ndim);
            vector<int> syms(A.ndim+B.ndim+this->ndim, NS);
            vector<int> stride_A(A.ndim+B.ndim+this->ndim, 0);
            vector<int> stride_B(A.ndim+B.ndim+this->ndim, 0);
            vector<int> stride_C(A.ndim+B.ndim+this->ndim, 0);

            int m = 0;
            for (int i = 0;i < A.ndim;i++)
            {
                bool found = false;
                for (int j = 0;j < m;j++)
                {
                    if (inds[j] == idx_A_[i])
                    {
                        found = true;
                        stride_A[j] += stride_idx_A[i];
                    }
                }
                if (!found)
                {
                    stride_A[m] += stride_idx_A[i];
                    inds[m++] = idx_A_[i];
                    if (i == A.ndim-1) continue;

                    bool in_B = false, aligned_in_B = true;
                    for (int i_in_B = 0;;i_in_B++)
                    {
                        for (;i_in_B < B.ndim && idx_A_[i] != idx_B_[i_in_B];i_in_B++);
                        if (i_in_B == B.ndim) break;

                        in_B = true;

                        if (A.sym[i] == NS || A.sym[i] != B.sym[i_in_B] ||
                            idx_A_[i+1] != idx_B_[i_in_B+1])
                        {
                            aligned_in_B = false;
                            break;
                        }
                    }

                    bool in_C = false, aligned_in_C = true;
                    for (int i_in_C = 0;;i_in_C++)
                    {
                        for (;i_in_C < this->ndim && idx_A_[i] != idx_C_[i_in_C];i_in_C++);
                        if (i_in_C == this->ndim) break;

                        in_C = true;

                        if (A.sym[i] == NS || A.sym[i] != this->sym[i_in_C] ||
                            idx_A_[i+1] != idx_C_[i_in_C+1])
                        {
                            aligned_in_C = false;
                            break;
                        }
                    }

                    if (!(in_B && !aligned_in_B) &&
                        !(in_C && !aligned_in_C) &&
                        A.sym[i] != NS) syms[m-1] = AS;
                }
            }
            for (int i = 0;i < B.ndim;i++)
            {
                bool found = false;
                for (int j = 0;j < m;j++)
                {
                    if (inds[j] == idx_B_[i])
                    {
                        found = true;
                        stride_B[j] += stride_idx_B[i];
                    }
                }
                if (!found)
                {
                    stride_B[m] += stride_idx_B[i];
                    inds[m++] = idx_B_[i];
                    if (i == B.ndim-1) continue;

                    bool in_C = false, aligned_in_C = true;
                    for (int i_in_C = 0;;i_in_C++)
                    {
                        for (;i_in_C < this->ndim && idx_B_[i] != idx_C_[i_in_C];i_in_C++);
                        if (i_in_C == this->ndim) break;

                        in_C = true;

                        if (B.sym[i] == NS || B.sym[i] != this->sym[i_in_C] ||
                            idx_B_[i+1] != idx_C_[i_in_C+1])
                        {
                            aligned_in_C = false;
                            break;
                        }
                    }

                    if (!(in_C && !aligned_in_C) &&
                        B.sym[i] != NS) syms[m-1] = AS;
                }
            }
            for (int i = 0;i < this->ndim;i++)
            {
                bool found = false;
                for (int j = 0;j < m;j++)
                {
                    if (inds[j] == idx_C_[i])
                    {
                        found = true;
                        stride_C[j] += stride_idx_C[i];
                    }
                }
                if (!found)
                {
                    stride_C[m] += stride_idx_C[i];
                    inds[m++] = idx_C_[i];
                    if (i == this->ndim-1) continue;
                    if (this->sym[i] != NS) syms[m-1] = AS;
                }
            }
            inds.resize(m);
            syms.resize(m);
            stride_A.resize(m);
            stride_B.resize(m);
            stride_C.resize(m);

            vector<Scalar> beta_(tensors.size(), beta);

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
                                                              this->ndim, idx_C_,   tensors[off_C].getSymmetry().data());
                    double f3 = A.factor[off_A]*B.factor[off_B]*factor[off_C];

                    for (int i = 0;i < A.ndim;i++) idx_A__[i] = idx_A_[A.reorder[off_A][i]];
                    for (int i = 0;i < B.ndim;i++) idx_B__[i] = idx_B_[B.reorder[off_B][i]];
                    for (int i = 0;i <   this->ndim;i++) idx_C__[i] = idx_C_[  reorder[off_C][i]];

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

                    if (idx[i] >= (syms[i] == NS ? n : idx[i+1]+1))
                    {
                        int lower = (i == 0 || syms[i-1] == NS ? 0 : idx[i-1]);
                        off_A -= stride_A[i]*(idx[i]-lower);
                        off_B -= stride_B[i]*(idx[i]-lower);
                        off_C -= stride_C[i]*(idx[i]-lower);
                        idx[i] = lower;
                        if (i == m-1) done = true;
                    }
                    else
                    {
                        break;
                    }
                }

                if (m == 0) done = true;
            }
        }

        void sum(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idx_A,
                 const Scalar&  beta,                                               const string& idx_B)
        {
            auto& A = static_cast<const SymmetryBlockedTensor&>(A_);

            assert(this->group == A.group);

            int n = this->group.getNumIrreps();

            string idx_A_(idx_A);
            string idx_B_(idx_B);
            string idx_A__(idx_A);
            string idx_B__(idx_B);

            double f = detail::align_symmetric_indices(A.ndim, idx_A_, A.sym.data(),
                                                 this->ndim, idx_B_,   this->sym.data());

            string inds_A = idx_A_;
            string inds_B = idx_B_;

            vector<int> stride_idx_A = detail::getStrides(idx_A_, A.ndim, n, idx_A_);
            vector<int> stride_idx_B = detail::getStrides(idx_B_,   this->ndim, n, idx_B_);

            vector<char> inds(A.ndim+this->ndim);
            vector<int> syms(A.ndim+this->ndim, NS);
            vector<int> stride_A(A.ndim+this->ndim, 0);
            vector<int> stride_B(A.ndim+this->ndim, 0);

            int m = 0;
            for (int i = 0;i < A.ndim;i++)
            {
                bool found = false;
                for (int j = 0;j < m;j++)
                {
                    if (inds[j] == idx_A_[i])
                    {
                        found = true;
                        stride_A[j] += stride_idx_A[i];
                    }
                }
                if (!found)
                {
                    stride_A[m] += stride_idx_A[i];
                    inds[m++] = idx_A_[i];
                    if (i == A.ndim-1) continue;

                    bool in_B = false, aligned_in_B = true;
                    for (int i_in_B = 0;;i_in_B++)
                    {
                        for (;i_in_B < this->ndim && idx_A_[i] != idx_B_[i_in_B];i_in_B++);
                        if (i_in_B == this->ndim) break;

                        in_B = true;

                        if (A.sym[i] == NS || A.sym[i] != this->sym[i_in_B] ||
                            idx_A_[i+1] != idx_B_[i_in_B+1])
                        {
                            aligned_in_B = false;
                            break;
                        }
                    }

                    if (!(in_B && !aligned_in_B) &&
                        A.sym[i] != NS) syms[m-1] = AS;
                }
            }
            for (int i = 0;i < this->ndim;i++)
            {
                bool found = false;
                for (int j = 0;j < m;j++)
                {
                    if (inds[j] == idx_B_[i])
                    {
                        found = true;
                        stride_B[j] += stride_idx_B[i];
                    }
                }
                if (!found)
                {
                    stride_B[m] += stride_idx_B[i];
                    inds[m++] = idx_B_[i];
                    if (i == this->ndim-1) continue;
                    if (this->sym[i] != NS) syms[m-1] = AS;
                }
            }
            inds.resize(m);
            stride_A.resize(m);
            stride_B.resize(m);

            vector<Scalar> beta_(tensors.size(), beta);

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
                    for (int i = 0;i <   this->ndim;i++) idx_B__[i] = idx_B_[  reorder[off_B][i]];

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

                    if (idx[i] == (syms[i] == NS ? n : idx[i+1]+1))
                    {
                        int lower = (i == 0 || syms[i-1] == NS ? 0 : idx[i-1]);
                        off_A -= stride_A[i]*(idx[i]-lower);
                        off_B -= stride_B[i]*(idx[i]-lower);
                        idx[i] = lower;
                        if (i == m-1) done = true;
                    }
                    else
                    {
                        break;
                    }
                }

                if (m == 0) done = true;
            }
        }

        void scale(const Scalar& alpha, const string& idx_A)
        {
            int n = this->group.getNumIrreps();

            string inds_A(idx_A);

            uniq(inds_A);

            int nA = inds_A.size();

            vector<int> stride_A = detail::getStrides(inds_A, this->ndim, n, idx_A);

            int off_A = 0;
            vector<int> iA(nA, 0);
            for (bool doneA = false;!doneA;)
            {
                if (alloced[off_A])
                {
                    tensors[off_A][idx_A] *= alpha;
                }

                for (int i = 0;i < nA;i++)
                {
                    iA[i]++;
                    off_A += stride_A[i];

                    if (iA[i] == n)
                    {
                        iA[i] = 0;
                        off_A -= stride_A[i]*n;
                        if (i == nA-1) doneA = true;
                    }
                    else
                    {
                        break;
                    }
                }

                if (nA == 0) doneA = true;
            }
        }

        Scalar dot(bool conja, const TensorImplementation<>& A_, const string& idx_A,
                   bool conjb,                                   const string& idx_B) const
        {
            auto& A = static_cast<const SymmetryBlockedTensor&>(A_);

            assert(this->group == A.group);

            int n = this->group.getNumIrreps();

            string idx_A_(idx_A);
            string idx_B_(idx_B);
            string idx_A__(idx_A);
            string idx_B__(idx_B);

            double f1 = detail::align_symmetric_indices(A.ndim, idx_A_, A.sym.data(),
                                                        this->ndim, idx_B_, this->sym.data());
            f1 *= detail::overcounting_factor(A.ndim, idx_A_, A.sym.data(),
                                              this->ndim, idx_B_, this->sym.data(),
                                                   0,     string(),         NULL);

            vector<int> stride_idx_A = detail::getStrides(idx_A_, A.ndim, n, idx_A_);
            vector<int> stride_idx_B = detail::getStrides(idx_B_, this->ndim, n, idx_B_);

            vector<char> inds(A.ndim+this->ndim);
            vector<int> syms(A.ndim+this->ndim, NS);
            vector<int> stride_A(A.ndim+this->ndim, 0);
            vector<int> stride_B(A.ndim+this->ndim, 0);

            int m = 0;
            for (int i = 0;i < A.ndim;i++)
            {
                bool found = false;
                for (int j = 0;j < m;j++)
                {
                    if (inds[j] == idx_A_[i])
                    {
                        found = true;
                        stride_A[j] += stride_idx_A[i];
                    }
                }
                if (!found)
                {
                    stride_A[m] += stride_idx_A[i];
                    inds[m++] = idx_A_[i];
                    if (i == A.ndim-1) continue;

                    bool in_B = false, aligned_in_B = true;
                    for (int i_in_B = 0;;i_in_B++)
                    {
                        for (;i_in_B < this->ndim && idx_A_[i] != idx_B_[i_in_B];i_in_B++);
                        if (i_in_B == this->ndim) break;

                        in_B = true;

                        if (A.sym[i] == NS || A.sym[i] != this->sym[i_in_B] ||
                            idx_A_[i+1] != idx_B_[i_in_B+1])
                        {
                            aligned_in_B = false;
                            break;
                        }
                    }

                    if (!(in_B && !aligned_in_B) &&
                        A.sym[i] != NS) syms[m-1] = AS;
                }
            }
            for (int i = 0;i < this->ndim;i++)
            {
                bool found = false;
                for (int j = 0;j < m;j++)
                {
                    if (inds[j] == idx_B_[i])
                    {
                        found = true;
                        stride_B[j] += stride_idx_B[i];
                    }
                }
                if (!found)
                {
                    stride_B[m] += stride_idx_B[i];
                    inds[m++] = idx_B_[i];
                    if (i == this->ndim-1) continue;
                    if (this->sym[i] != NS) syms[m-1] = AS;
                }
            }
            inds.resize(m);
            syms.resize(m);
            stride_A.resize(m);
            stride_B.resize(m);

            Scalar total;

            int off_A = 0;
            int off_B = 0;
            vector<int> idx(m, 0);
            for (bool done = false;!done;)
            {
                if ((A.tensors.ptr(off_A) && tensors.ptr(off_B)) &&
                    (A.alloced[off_A] || alloced[off_B]))
                {
                    double f2 = detail::overcounting_factor(A.ndim, idx_A_, tensors[off_A].getSymmetry().data(),
                                                            this->ndim, idx_B_, tensors[off_B].getSymmetry().data(),
                                                                 0,     string(),                                  NULL);
                    double f3 = A.factor[off_A]*factor[off_B];

                    for (int i = 0;i < A.ndim;i++) idx_A__[i] = idx_A_[A.reorder[off_A][i]];
                    for (int i = 0;i < this->ndim;i++) idx_B__[i] = idx_B_[reorder[off_B][i]];

                    assert(off_A >= 0 && off_A < A.tensors.size());
                    assert(off_B >= 0 && off_B < tensors.size());
                    total += (f1*f3/f2)*scalar(aquarius::conj(A.tensors[off_A][idx_A__], conja)*
                                               aquarius::conj(tensors[off_B][idx_B__], conjb));
                }

                for (int i = 0;i < m;i++)
                {
                    idx[i]++;
                    off_A += stride_A[i];
                    off_B += stride_B[i];

                    if (idx[i] >= (syms[i] == NS ? n : idx[i+1]+1))
                    {
                        int lower = (i == 0 || syms[i-1] == NS ? 0 : idx[i-1]);
                        off_A -= stride_A[i]*(idx[i]-lower);
                        off_B -= stride_B[i]*(idx[i]-lower);
                        idx[i] = lower;
                        if (i == m-1) done = true;
                    }
                    else
                    {
                        break;
                    }
                }

                if (m == 0) done = true;
            }

            return total;
        }

        Scalar norm(int p) const
        {
            Scalar nrm;

            int n = this->group.getNumIrreps();

            vector<int> stride(this->ndim,1);
            for (int i = 1;i < this->ndim;i++) stride[i] = stride[i-1]*n;

            int off_A = 0;
            vector<int> iA(this->ndim, 0);
            for (bool doneA = false;!doneA;)
            {
                if (tensors.ptr(off_A) && alloced[off_A])
                {
                    double factor = 1;
                    const vector<int>& subsym = tensors[off_A].getSymmetry();
                    for (int i = 0;i < this->ndim;)
                    {
                        int j; for (j = i;this->sym[j] != NS;j++); j++;

                        int m = j-i;
                        for (int k = i;k < j;)
                        {
                            int l; for (l = k;subsym[l] != NS;l++); l++;
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

                for (int i = 0;i < this->ndim;i++)
                {
                    iA[i]++;
                    off_A += stride[i];

                    if (iA[i] == n)
                    {
                        iA[i] = 0;
                        off_A -= stride[i]*n;
                        if (i == this->ndim-1) doneA = true;
                    }
                    else
                    {
                        break;
                    }
                }

                if (this->ndim == 0) doneA = true;
            }

            if (p == 2) nrm = sqrt(nrm);

            return nrm;
        }

        void slice(const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, const TensorImplementation<>& A_,
                   const Scalar&  beta,             const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            auto& A = static_cast<const SymmetryBlockedTensor&>(A_);

            assert(this->ndim == A.ndim);

            int n = this->group.getNumIrreps();

            vector<vector<int>> end_A(this->ndim);
            vector<vector<int>> end_B(this->ndim);

            for (int i = 0;i < this->ndim;i++)
            {
                end_A[i].resize(n);
                end_B[i].resize(n);

                assert(this->sym[i] == A.sym[i] && this->sym[i] == NS);
                for (int j = 0;j < n;j++)
                {
                    end_A[i][j] = start_A[i][j]+this->len_per_irrep[i][j];
                    end_B[i][j] = start_B[i][j]+this->len_per_irrep[i][j];
                    assert(start_A[i][j] >= 0);
                    assert(start_B[i][j] >= 0);
                    assert(end_A[i][j] <= A.len_per_irrep[i][j]);
                    assert(end_B[i][j] <= this->len_per_irrep[i][j]);
                }
            }

            vector<int> stride(this->ndim,1);
            for (int i = 1;i < this->ndim;i++) stride[i] = stride[i-1]*n;

            int off_A = 0;
            vector<int> iA(this->ndim, 0);
            vector<int> start_A_sub(this->ndim);
            vector<int> start_B_sub(this->ndim);
            vector<int> len_sub(this->ndim);
            for (bool doneA = false;!doneA;)
            {
                if (  tensors.ptr(off_A) &&   alloced[off_A] &&
                    A.tensors.ptr(off_A) && A.alloced[off_A])
                {
                    for (int i = 0;i < this->ndim;i++)
                    {
                        start_A_sub[i] = start_A[i][iA[i]];
                        start_B_sub[i] = start_B[i][iA[i]];
                        len_sub[i] = this->len_per_irrep[i][iA[i]];
                    }

                    tensors[off_A].slice(alpha, conja, start_A_sub, A.tensors[off_A],
                                          beta,        start_B_sub, len_sub);
                }

                for (int i = 0;i < this->ndim;i++)
                {
                    iA[i]++;
                    off_A += stride[i];

                    if (iA[i] == n)
                    {
                        iA[i] = 0;
                        off_A -= stride[i]*n;
                        if (i == this->ndim-1) doneA = true;
                    }
                    else
                    {
                        break;
                    }
                }

                if (this->ndim == 0) doneA = true;
            }
        }

        void div(const Scalar& alpha, bool conja, const TensorImplementation<>& A_,
                                      bool conjb, const TensorImplementation<>& B_,
                 const Scalar& beta)
        {
            auto& A = static_cast<const SymmetryBlockedTensor&>(A_);
            auto& B = static_cast<const SymmetryBlockedTensor&>(B_);

            assert(this->group == A.group);
            assert(this->group == B.group);
            assert(this->sym == A.sym);
            assert(this->sym == B.sym);
            assert(this->len_per_irrep == A.len_per_irrep);
            assert(this->len_per_irrep == B.len_per_irrep);

            int n = this->group.getNumIrreps();

            vector<int> stride(this->ndim,1);
            for (int i = 1;i < this->ndim;i++) stride[i] = stride[i-1]*n;

            int off_A = 0;
            vector<int> iA(this->ndim, 0);
            for (bool doneA = false;!doneA;)
            {
                if (tensors.ptr(off_A) && alloced[off_A])
                {
                    beta*tensors[off_A] += alpha*A.tensors[off_A]/B.tensors[off_A];
                }

                for (int i = 0;i < this->ndim;i++)
                {
                    iA[i]++;
                    off_A += stride[i];

                    if (iA[i] == n)
                    {
                        iA[i] = 0;
                        off_A -= stride[i]*n;
                        if (i == this->ndim-1) doneA = true;
                    }
                    else
                    {
                        break;
                    }
                }

                if (this->ndim == 0) doneA = true;
            }
        }

        void invert(const Scalar& alpha, bool conja, const TensorImplementation<>& A_,
                    const Scalar& beta)
        {
            auto& A = static_cast<const SymmetryBlockedTensor&>(A_);

            assert(this->group == A.group);
            assert(this->sym == A.sym);
            assert(this->len_per_irrep == A.len_per_irrep);

            int n = this->group.getNumIrreps();

            vector<int> stride(this->ndim,1);
            for (int i = 1;i < this->ndim;i++) stride[i] = stride[i-1]*n;

            int off_A = 0;
            vector<int> iA(this->ndim, 0);
            for (bool doneA = false;!doneA;)
            {
                if (tensors.ptr(off_A) && alloced[off_A])
                {
                    beta*tensors[off_A] += alpha/A.tensors[off_A];
                }

                for (int i = 0;i < this->ndim;i++)
                {
                    iA[i]++;
                    off_A += stride[i];

                    if (iA[i] == n)
                    {
                        iA[i] = 0;
                        off_A -= stride[i]*n;
                        if (i == this->ndim-1) doneA = true;
                    }
                    else
                    {
                        break;
                    }
                }

                if (this->ndim == 0) doneA = true;
            }
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
}

#endif
