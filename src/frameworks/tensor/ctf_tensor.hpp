#ifndef _AQUARIUS_TENSOR_CTF_TENSOR_HPP_
#define _AQUARIUS_TENSOR_CTF_TENSOR_HPP_

#include "util/global.hpp"

#include "task/task.hpp"

#include "indexable_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <typename T>
class CTFTensor : public IndexableTensor< CTFTensor<T>,T >, public Distributed
{
    INHERIT_FROM_INDEXABLE_TENSOR(CTFTensor<T>,T)

    protected:
        tCTF_Tensor<T>* dt;
        vector<int> len;
        vector<int> sym;
        static map<const tCTF_World<T>*,pair<int,CTFTensor<T>*>> scalars;

        void allocate();

        void free();

        void register_scalar();

        void unregister_scalar();

        CTFTensor<T>& scalar() const;

        static void first_packed_indices(int ndim, const int* len, const int* sym, int* idx)
        {
            int i;

            if (ndim > 0) idx[0] = 0;
            for (i = 0;i < ndim - 1;i++)
            {
                switch (sym[i])
                {
                    case AS:
                    case SH:
                        idx[i+1] = idx[i] + 1;
                        break;
                    case SY:
                        idx[i+1] = idx[i];
                        break;
                    case NS:
                        idx[i+1] = 0;
                        break;
                }
            }
        }

        static bool next_packed_indices(int ndim, const int* len, const int* sym, int* idx)
        {
            int i;

            for (i = 0;i < ndim;i++)
            {
                if (i == ndim - 1)
                {
                    if (idx[i] >= len[i]-1)
                    {
                        return false;
                    }
                    else
                    {
                        idx[i]++;
                        break;
                    }
                }
                else
                {
                    if ((sym[i] == SY && idx[i] >= idx[i+1]) ||
                        ((sym[i] == AS || sym[i] == SH) && idx[i] >= idx[i+1]-1) ||
                        (idx[i] >= len[i]-1))
                    {
                        if (i == 0)
                        {
                            idx[i] = 0;
                        }
                        else if (sym[i-1] == NS)
                        {
                            idx[i] = 0;
                        }
                        else if (sym[i-1] == SY)
                        {
                            idx[i] = idx[i-1];
                        }
                        else // AS and SH
                        {
                            idx[i] = idx[i-1] + 1;
                        }
                    }
                    else
                    {
                        idx[i]++;
                        break;
                    }
                }
            }

            return (ndim > 0 ? true : false);
        }

    public:
        CTFTensor(const string& name, const Arena& arena, T scalar = (T)0);

        CTFTensor(const string& name, const CTFTensor<T>& A, T scalar);

        CTFTensor(const CTFTensor<T>& A, bool copy=true, bool zero=false);

        CTFTensor(const string& name, const CTFTensor<T>& A, bool copy=true, bool zero=false);

        CTFTensor(const string& name, CTFTensor<T>* A);

        CTFTensor(const string& name, const CTFTensor<T>& A, const vector<int>& start_A, const vector<int>& len_A);

        CTFTensor(const string& name, const Arena& arena, int ndim, const vector<int>& len, const vector<int>& sym,
                   bool zero=true);

        ~CTFTensor();

        void resize(int ndim, const vector<int>& len, const vector<int>& sym, bool zero);

        const vector<int>& getLengths() const { return len; }

        const vector<int>& getSymmetry() const { return sym; }

        T* getRawData(int64_t& size);

        const T* getRawData(int64_t& size) const;

        template <typename Container>
        void getLocalData(Container& pairs) const
        {
            int64_t npair;
            tkv_pair<T> *data;
            dt->read_local(&npair, &data);
            pairs.assign(data, data+npair);
            if (npair > 0) ::free(data);
        }

        template <typename Container>
        void getRemoteData(Container& pairs) const
        {
            dt->read(pairs.size(), pairs.data());
        }

        void getRemoteData() const
        {
            dt->read(0, NULL);
        }

        template <typename Container>
        void writeRemoteData(const Container& pairs)
        {
            dt->write(pairs.size(), pairs.data());
        }

        void writeRemoteData()
        {
            dt->write(0, NULL);
        }

        template <typename Container>
        void writeRemoteData(double alpha, double beta, const Container& pairs)
        {
            dt->write(pairs.size(), alpha, beta, pairs.data());
        }

        void writeRemoteData(double alpha, double beta)
        {
            dt->write(0, alpha, beta, NULL);
        }

        template <typename Container>
        void getAllData(Container& vals) const
        {
            int64_t npair;
            if (this->arena.rank == 0)
            {
                getAllData(vals, 0);
                npair = vals.size();
            }
            else
            {
                getAllData(0);
            }
            this->arena.comm().Bcast(&npair, 1, 0);
            if (this->arena.rank != 0) vals.resize(npair);
            this->arena.comm().Bcast(vals, 0);
        }

        template <typename Container>
        void getAllData(Container& vals, int rank) const
        {
            assert(this->arena.rank == rank);

            for (int i = 0;i < ndim;i++)
            {
                if (len[i] == 0)
                {
                    vals.clear();
                    return;
                }
            }

            vector<tkv_pair<T>> pairs;
            vector<int> idx(ndim, 0);

            first_packed_indices(ndim, len.data(), sym.data(), idx.data());

            do
            {
                int64_t key = 0, stride = 1;
                for (int i = 0;i < ndim;i++)
                {
                    key += idx[i]*stride;
                    stride *= len[i];
                }
                pairs.push_back(tkv_pair<T>(key, (T)0));
            }
            while (next_packed_indices(ndim, len.data(), sym.data(), idx.data()));

            dt->read(pairs.size(), pairs.data());

            sort(pairs.begin(), pairs.end());
            size_t npair = pairs.size();
            vals.resize(npair);

            for (size_t i = 0;i < npair;i++)
            {
                vals[i] = pairs[i].d;
            }
        }

        void getAllData(int rank) const
        {
            assert(this->arena.rank != rank);
            dt->read(0, NULL);
        }

        void slice(T alpha, bool conja, const CTFTensor<T>& A,
                   const vector<int>& start_A, T beta);

        void slice(T alpha, bool conja, const CTFTensor<T>& A,
                   T beta, const vector<int>& start_B);

        void slice(T alpha, bool conja, const CTFTensor<T>& A, const vector<int>& start_A,
                   T  beta,                                     const vector<int>& start_B,
                                                                const vector<int>& len);

        void div(T alpha, bool conja, const CTFTensor<T>& A,
                          bool conjb, const CTFTensor<T>& B, T beta);

        void invert(T alpha, bool conja, const CTFTensor<T>& A, T beta);

        void weight(const vector<const vector<T>*>& d, double shift = 0);

        void print(FILE* fp, double cutoff = -1.0) const;

        void compare(FILE* fp, const CTFTensor<T>& other, double cutoff = 0.0) const;

        real_type_t<T> norm(int p) const;

        void mult(T alpha, bool conja, const CTFTensor<T>& A, const string& idx_A,
                           bool conjb, const CTFTensor<T>& B, const string& idx_B,
                  T  beta,                                     const string& idx_C);

        void sum(T alpha, T beta);

        void sum(T alpha, bool conja, const CTFTensor<T>& A, const string& idx_A,
                 T  beta,                                     const string& idx_B);

        void scale(T alpha, const string& idx_A);

        T dot(bool conja, const CTFTensor<T>& A, const string& idx_A,
              bool conjb,                         const string& idx_B) const;
};

}
}

#endif
