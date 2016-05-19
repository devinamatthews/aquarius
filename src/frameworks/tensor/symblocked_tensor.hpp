#ifndef _AQUARIUS_TENSOR_SYMBLOCKED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_SYMBLOCKED_TENSOR_HPP_

#include "util/global.hpp"

#include "symmetry/symmetry.hpp"
#include "task/task.hpp"

#include "composite_tensor.hpp"
#include "ctf_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <class T>
class SymmetryBlockedTensor : public IndexableCompositeTensor<SymmetryBlockedTensor<T>,CTFTensor<T>,T>,
                              public Distributed
{
    INHERIT_FROM_INDEXABLE_COMPOSITE_TENSOR(SymmetryBlockedTensor<T>,CTFTensor<T>,T)

    protected:
        const symmetry::PointGroup& group;
        symmetry::Representation rep;
        vector<vector<int>> len;
        vector<int> sym;
        vector<double> factor;
        vector<vector<int>> reorder;
        static map<const tCTF_World<T>*,map<const symmetry::PointGroup*,pair<int,SymmetryBlockedTensor<T>*>>> scalars;

        static vector<int> getStrides(const string& indices, int ndim,
                                      int len, const string& idx_A);

        void allocate(bool zero);

        void register_scalar();

        void unregister_scalar();

        SymmetryBlockedTensor<T>& scalar() const;

    public:
        CTFTensor<T>& operator()(const vector<int>& irreps);

        const CTFTensor<T>& operator()(const vector<int>& irreps) const;

        SymmetryBlockedTensor(const SymmetryBlockedTensor<T>& other);

        SymmetryBlockedTensor(SymmetryBlockedTensor<T>&& other);

        SymmetryBlockedTensor(const string& name, const SymmetryBlockedTensor<T>& other);

        SymmetryBlockedTensor(const string& name, SymmetryBlockedTensor<T>&& other);

        SymmetryBlockedTensor(const string& name, const SymmetryBlockedTensor<T>& other, T scalar);

        SymmetryBlockedTensor(const string& name, const SymmetryBlockedTensor<T>& A,
                              const vector<vector<int>>& start_A,
                              const vector<vector<int>>& len_A);

        SymmetryBlockedTensor(const string& name, const Arena& arena, const symmetry::PointGroup& group,
                              int ndim, const vector<vector<int>>& len,
                              const vector<int>& sym, bool zero=true);

        SymmetryBlockedTensor(const string& name, const Arena& arena, const symmetry::PointGroup& group,
                              const symmetry::Representation& rep, int ndim, const vector<vector<int>>& len,
                              const vector<int>& sym, bool zero=true);

        ~SymmetryBlockedTensor();

        const symmetry::PointGroup& getGroup() const { return group; }

        const vector<vector<int>>& getLengths() const { return len; }

        const vector<int>& getSymmetry() const { return sym; }

        bool exists(const vector<int>& irreps) const;

        T* getRawData(const vector<int>& irreps, int64_t& size)
        {
            return (*this)(irreps).getRawData(size);
        }

        const T* getRawData(const vector<int>& irreps, int64_t& size) const
        {
            return (*this)(irreps).getRawData(size);
        }

        template <typename Container>
        void getLocalData(const vector<int>& irreps, Container& pairs) const
        {
            (*this)(irreps).getLocalData(pairs);
        }

        template <typename Container>
        void getRemoteData(const vector<int>& irreps, Container& pairs) const
        {
            (*this)(irreps).getRemoteData(pairs);
        }

        void getRemoteData(const vector<int>& irreps) const
        {
            (*this)(irreps).getRemoteData();
        }

        template <typename Container>
        void writeRemoteData(const vector<int>& irreps, const Container& pairs)
        {
            (*this)(irreps).writeRemoteData(pairs);
        }

        void writeRemoteData(const vector<int>& irreps)
        {
            (*this)(irreps).writeRemoteData();
        }

        template <typename Container>
        void writeRemoteData(const vector<int>& irreps, double alpha, double beta, const Container& pairs)
        {
            (*this)(irreps).writeRemoteData(alpha, beta, pairs);
        }

        void writeRemoteData(const vector<int>& irreps, double alpha, double beta)
        {
            (*this)(irreps).writeRemoteData(alpha, beta);
        }

        template <typename Container>
        void getAllData(const vector<int>& irreps, Container& vals) const
        {
            (*this)(irreps).getAllData(vals);
        }

        template <typename Container>
        void getAllData(const vector<int>& irreps, Container& vals, int rank) const
        {
            (*this)(irreps).getAllData(vals, rank);
        }

        void getAllData(const vector<int>& irreps, int rank) const
        {
            (*this)(irreps).getAllData(rank);
        }

        void slice(T alpha, bool conja, const SymmetryBlockedTensor<T>& A,
                   const vector<vector<int>>& start_A, T beta);

        void slice(T alpha, bool conja, const SymmetryBlockedTensor<T>& A,
                   T beta, const vector<vector<int>>& start_B);

        void slice(T alpha, bool conja, const SymmetryBlockedTensor<T>& A,
                                        const vector<vector<int>>& start_A,
                   T  beta,             const vector<vector<int>>& start_B,
                                        const vector<vector<int>>& len);

        virtual void mult(T alpha, bool conja, const SymmetryBlockedTensor<T>& A, const string& idx_A,
                                   bool conjb, const SymmetryBlockedTensor<T>& B, const string& idx_B,
                          T beta,                                                 const string& idx_C);

        virtual void sum(T alpha, bool conja, const SymmetryBlockedTensor<T>& A, const string& idx_A,
                         T beta,                                                 const string& idx_B);

        virtual void scale(T alpha, const string& idx_A);

        virtual T dot(bool conja, const SymmetryBlockedTensor<T>& A, const string& idx_A,
                      bool conjb,                                    const string& idx_B) const;

        void weight(const vector<const vector<vector<T>>*>& d,
                    double shift = 0);

        real_type_t<T> norm(int p) const;
};

}
}

#endif
