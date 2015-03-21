#ifndef _AQUARIUS_COMPOSITE_TENSOR_HPP_
#define _AQUARIUS_COMPOSITE_TENSOR_HPP_

#include "util/global.hpp"

#include "indexable_tensor.hpp"

namespace aquarius
{
namespace tensor
{

#define INHERIT_FROM_COMPOSITE_TENSOR(Derived,Base,T) \
    protected: \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::tensors; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::addTensor; \
    public: \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::mult; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::div; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::sum; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::invert; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::dot; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::exists; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::operator(); \
        using aquarius::tensor::Tensor< Derived,T >::getDerived; \
        using aquarius::tensor::Tensor< Derived,T >::operator=; \
        using aquarius::tensor::Tensor< Derived,T >::operator+=; \
        using aquarius::tensor::Tensor< Derived,T >::operator-=; \
        using aquarius::tensor::Tensor< Derived,T >::operator*=; \
        using aquarius::tensor::Tensor< Derived,T >::operator/=; \
        using aquarius::tensor::Tensor< Derived,T >::operator*; \
        using aquarius::tensor::Tensor< Derived,T >::operator/; \
        Derived & operator=(const Derived & other) \
        { \
            sum((T)1, false, other, (T)0); \
            return *this; \
        } \
    private:

#define INHERIT_FROM_INDEXABLE_COMPOSITE_TENSOR(Derived,Base,T) \
    protected: \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::tensors; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::addTensor; \
        using aquarius::tensor::IndexableTensorBase< Derived, T >::ndim; \
    public: \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::mult; \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::sum; \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::div; \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::invert; \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::scale; \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::dot; \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::operator=; \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::operator+=; \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::operator-=; \
        using aquarius::tensor::IndexableTensorBase< Derived, T >::operator[]; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::operator(); \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::exists; \
        using aquarius::tensor::Tensor< Derived,T >::getDerived; \
        using aquarius::tensor::Tensor< Derived,T >::operator*=; \
        using aquarius::tensor::Tensor< Derived,T >::operator/=; \
        using aquarius::tensor::Tensor< Derived,T >::operator*; \
        using aquarius::tensor::Tensor< Derived,T >::operator/; \
        Derived & operator=(const Derived & other) \
        { \
            sum((T)1, false, other, (T)0); \
            return *this; \
        } \
    private:

template <class Derived, class Base, class T>
class CompositeTensor : public Tensor<Derived,T>
{
    protected:
        struct TensorRef
        {
            Base* tensor;
            bool isAlloced;
            int ref;
            TensorRef(Base* tensor_=NULL, bool isAlloced=false, int ref=-1)
            : tensor(tensor_), isAlloced(isAlloced), ref(ref) {}
            bool operator==(const Base* other) const { return tensor == other; }
            bool operator!=(const Base* other) const { return tensor != other; }
        };

        vector<TensorRef> tensors;

        Base& addTensor(Base* new_tensor, bool isAlloced=true)
        {
            tensors.push_back(TensorRef(new_tensor, isAlloced));
            return *new_tensor;
        }

        Base& addTensor(Base& new_tensor, bool isAlloced=false)
        {
            tensors.push_back(TensorRef(&new_tensor, isAlloced));
            return new_tensor;
        }

        Base& addTensor(int ref)
        {
            assert(ref >= -1 && ref < tensors.size());
            tensors.push_back(TensorRef(tensors[ref].tensor, false, ref));
            return *tensors[ref].tensor;
        }

    public:
        CompositeTensor(const CompositeTensor<Derived,Base,T>& other)
        : Tensor<Derived,T>(other.name), tensors(other.tensors)
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].ref == -1)
                {
                    tensors[i].tensor = new Base(*tensors[i].tensor);
                }
            }
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i].ref != -1)
                {
                    tensors[i].tensor = tensors[tensors[i].ref].tensor;
                }
            }
        }

        CompositeTensor(const string& name, const CompositeTensor<Derived,Base,T>& other)
        : Tensor<Derived,T>(name), tensors(other.tensors)
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].ref == -1)
                {
                    tensors[i].tensor = new Base(name, *tensors[i].tensor);
                }
            }
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i].ref != -1)
                {
                    tensors[i].tensor = tensors[tensors[i].ref].tensor;
                }
            }
        }

        CompositeTensor(const string& name, int ntensors = 0)
        : Tensor<Derived,T>(name), tensors(ntensors) {}

        virtual ~CompositeTensor()
        {
            for (int i = tensors.size()-1;i >= 0;i--)
            {
                if (tensors[i].isAlloced)
                {
                    delete tensors[i].tensor;
                }
            }
        }

        int getNumTensors() const { return tensors.size(); }

        bool exists(int idx) const
        {
            return tensors[idx] != NULL;
        }

        /**********************************************************************
         *
         * Subtensor indexing
         *
         *********************************************************************/
        Base& operator()(int idx)
        {
            if (tensors[idx] == NULL)
                throw logic_error("tensor component does not exist");
            return *tensors[idx].tensor;
        }

        const Base& operator()(int idx) const
        {
            if (tensors[idx] == NULL)
                throw logic_error("tensor component does not exist");
            return *tensors[idx].tensor;
        }

        /**********************************************************************
         *
         * Implementation of Tensor stubs
         *
         *********************************************************************/
        void mult(const T alpha)
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].ref == -1)
                {
                    *tensors[i].tensor *= alpha;
                }
            }
        }

        void mult(const T alpha, bool conja, const Derived& A,
                                 bool conjb, const Derived& B, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors.size() ||
                tensors.size() != B.tensors.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].ref == -1 && A.exists(i) && B.exists(i))
                {
                    beta*(*tensors[i].tensor) += alpha*A(i)*B(i);
                }
            }
        }

        void div(const T alpha, bool conja, const Derived& A,
                                bool conjb, const Derived& B, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors.size() ||
                tensors.size() != B.tensors.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].ref == -1 && A.exists(i) && B.exists(i))
                {
                    beta*(*tensors[i].tensor) += alpha*A(i)/B(i);
                }
            }
        }

        void sum(const T alpha, const T beta)
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].ref == -1)
                {
                    beta*(*tensors[i].tensor) += alpha;
                }
            }
        }

        void sum(const T alpha, bool conja, const Derived& A, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].ref == -1 && A.exists(i))
                {
                    beta*(*tensors[i].tensor) += alpha*A(i);
                }
            }
        }

        void invert(const T alpha, bool conja, const Derived& A, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].ref == -1 && A.exists(i))
                {
                    beta*(*tensors[i].tensor) += alpha/A(i);
                }
            }
        }

        T dot(bool conja, const Derived& A, bool conjb) const
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            T s = (T)0;

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].ref == -1 && A.exists(i))
                {
                    s += tensors[i].tensor->dot(conja, A(i), conjb);
                }
            }

            return s;
        }
};

template <class Derived, class Base, class T>
class IndexableCompositeTensor : public IndexableTensorBase<Derived,T>, public CompositeTensor<Derived,Base,T>
{
    INHERIT_FROM_TENSOR(Derived,T)

    protected:
        using IndexableTensorBase<Derived,T>::ndim;

    public:
        using IndexableTensorBase< Derived, T >::operator=;
        using IndexableTensorBase< Derived, T >::operator+=;
        using IndexableTensorBase< Derived, T >::operator-=;
        //using CompositeTensor<Derived,Base,T>::div;
        //using CompositeTensor<Derived,Base,T>::invert;
        using IndexableTensorBase<Derived,T>::scale;
        using IndexableTensorBase<Derived,T>::dot;
        using IndexableTensorBase<Derived,T>::mult;
        using IndexableTensorBase<Derived,T>::sum;
        using IndexableTensorBase<Derived,T>::implicit;

    public:
        IndexableCompositeTensor(const Derived& other)
        : IndexableTensorBase<Derived,T>(other), CompositeTensor<Derived,Base,T>(other) {}

        IndexableCompositeTensor(const string& name, const Derived& other)
        : IndexableTensorBase<Derived,T>(other), CompositeTensor<Derived,Base,T>(name, other) {}

        IndexableCompositeTensor(const string& name, int ndim=0, int ntensors=0)

        : IndexableTensorBase<Derived,T>(ndim), CompositeTensor<Derived,Base,T>(name, ntensors) {}

        virtual ~IndexableCompositeTensor() {}

        void mult(const T alpha)
        {
            scale(alpha);
        }

        void mult(const T alpha, bool conja, const Derived& A,
                                 bool conjb, const Derived& B,
                  const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim != A.ndim || ndim != B_.ndim) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            mult(alpha, conja, A, A.implicit(),
                        conjb, B, B.implicit(),
                  beta,             implicit());
        }

        virtual Derived& scalar() const = 0;

        void sum(const T alpha, const T beta)
        {
            CompositeTensor<Derived,Base,T>::sum(alpha, beta);
        }

        void sum(const T alpha, bool conja, const Derived& A, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim != A.ndim) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            sum(alpha, conja, A, A.implicit(),
                 beta,             implicit());
        }

        void div(const T alpha, bool conja, const Derived& A,
                                bool conjb, const Derived& B, const T beta)
        {
            CompositeTensor<Derived,Base,T>::div(alpha, conja, A, conjb, B, beta);
        }

        void invert(const T alpha, bool conja, const Derived& A, const T beta)
        {
            CompositeTensor<Derived,Base,T>::invert(alpha, conja, A, beta);
        }

        void scale(const T alpha)
        {
            scale(alpha, implicit());
        }

        T dot(bool conja, const Derived& A, bool conjb) const
        {
            #ifdef VALIDATE_INPUTS
            if (ndim != A.ndim) throw InvalidNdimError();
            #endif //VALIDATE_INPUTS

            return dot(conja, A, A.implicit(),
                       conjb,      implicit());
        }
};

}
}

#endif
