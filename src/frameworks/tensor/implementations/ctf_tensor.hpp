#ifndef _AQUARIUS_TENSOR_IMPLEMENTATIONS_CTF_TENSOR_HPP_
#define _AQUARIUS_TENSOR_IMPLEMENTATIONS_CTF_TENSOR_HPP_

#include "ctf.hpp"

#include "frameworks/util.hpp"
#include "frameworks/tensor.hpp"

namespace aquarius
{
namespace tensor
{

class CTFTensor : public TensorImplementation<BOUNDED|IPSYMMETRIC|INDEXABLE|DISTRIBUTED|DIVISIBLE>
{
    public:
        CTFTensor(const INITIALIZER_TYPE(BOUNDED|IPSYMMETRIC|INDEXABLE|DISTRIBUTED|DIVISIBLE)& init);

        ~CTFTensor();

        void mult(const Scalar& alpha, bool conja, const TensorImplementation<>& A, const string& idxA,
                                       bool conjb, const TensorImplementation<>& B, const string& idxB,
                  const Scalar& beta,                                               const string& idxC);

        void sum(const Scalar& alpha, bool conja, const TensorImplementation<>& A, const string& idxA,
                 const Scalar& beta,                                               const string& idxB);

        void scale(const Scalar& alpha, const string& idxA);

        Scalar dot(bool conja, const TensorImplementation<>& A, const string& idxA,
                   bool conjb,                                  const string& idxB) const;

        void div(const Scalar& alpha, bool conja, const TensorImplementation<>& A,
                                      bool conjb, const TensorImplementation<>& B,
                 const Scalar& beta);

        void invert(const Scalar& alpha, bool conja, const TensorImplementation<>& A,
                    const Scalar& beta);

        const vector<key_type>& getKeyStrides() const;

        void getAllKeys(KeyVector& keys) const;

        void getAllData(KeyValueVector& kv) const;

        void getLocalKeys(KeyVector& keys) const;

        void getLocalData(KeyValueVector& kv) const;

        void setLocalData(key_type n, const key_type* keys, const void* values)
        {
            setRemoteData(n, keys, values);
        }

        void addLocalData(key_type n, const Scalar& alpha, const key_type* keys,
                          const void* values, const Scalar& beta)
        {
            addRemoteData(n, alpha, keys, values, beta);
        }

        void getRemoteData(key_type n, key_type* keys, void* values) const;

        void getRemoteData() const;

        void setRemoteData(key_type n, const key_type* keys, const void* values);

        void setRemoteData();

        void addRemoteData(key_type n, const Scalar& alpha, const key_type* keys,
                           const void* values, const Scalar& beta);

        void addRemoteData(const Scalar& alpha, const Scalar& beta);

        Scalar norm(int p) const;

        void slice(const Scalar& alpha, bool conja, const vector<int>& start_A, const TensorImplementation<>& A,
                   const Scalar&  beta,             const vector<int>& start_B, const vector<int>& len);

    protected:
        unique_ptr<CTF_int::tensor> ctf;
        vector<key_type> strides;
        //static map<pair<int,const CTF::World*>, pair<int,void*> > scalars;

        //template <typename T> CTF::Tensor<T>& scalar() const
        //{
        //    assert(F == Field(T()));
        //    auto it = scalars.find(make_pair(F.type(), &arena.ctf()));
        //    assert(it != scalars.end());
        //    return *static_cast<CTF::Tensor<T>*>(it->second.second);
        //}
};

}
}

#endif
