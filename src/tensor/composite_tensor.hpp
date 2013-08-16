/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_COMPOSITE_TENSOR_HPP_
#define _AQUARIUS_COMPOSITE_TENSOR_HPP_

#include <vector>
#include <string>
#include <algorithm>

#include "stl_ext/stl_ext.hpp"

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
            TensorRef() : tensor(NULL), isAlloced(true) {}
            TensorRef(Base* tensor_, bool isAlloced=true)
            : tensor(tensor_), isAlloced(isAlloced) {}
            bool operator==(const Base* other) const { return tensor == other; }
            bool operator!=(const Base* other) const { return tensor != other; }
        };

        std::vector<TensorRef> tensors;

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

    public:
        void set_name(const char * name_){
            int i;
            for (i=0; i<tensors_.size(); i++){
                tensors_[i].tensor_->set_name(name_);
            }
        }


        CompositeTensor(const CompositeTensor<Derived,Base,T>& other)
        : tensors(other.tensors)
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL)
                {
                    tensors[i] = TensorRef(new Base(*tensors[i].tensor));
                }
            }
        }

        CompositeTensor(const int ntensors = 0)
        : tensors(ntensors) {}

        virtual ~CompositeTensor()
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL && tensors[i].isAlloced)
                {
                    delete tensors[i].tensor;
                }
            }
        }

        int getNumTensors() const { return tensors.size(); }

        bool componentExists(const int idx) const
        {
            return tensors[idx] != NULL;
        }

        /**********************************************************************
         *
         * Subtensor indexing
         *
         *********************************************************************/
        Base& operator()(const int idx)
        {
            if (tensors[idx] == NULL)
                throw std::logic_error("tensor component does not exist");
            return *tensors[idx].tensor;
        }

        const Base& operator()(const int idx) const
        {
            if (tensors[idx] == NULL)
                throw std::logic_error("tensor component does not exist");
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
                if (tensors[i] != NULL)
                {
                    *tensors[i].tensor *= alpha;
                }
            }
        }

        void mult(const T alpha, bool conja, const Derived& A,
                                 bool conjb, const Derived& B, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors_.size() ||
                tensors.size() != B.tensors_.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL &&
                    A.componentExists(i) &&
                    B.componentExists(i))
                {
                    beta*(*tensors[i].tensor) += alpha*(A(i))*
                                                         (B(i));
                }
            }
        }

        void div(const T alpha, bool conja, const Derived& A,
                                bool conjb, const Derived& B, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors_.size() ||
                tensors.size() != B.tensors_.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL &&
                    A.componentExists(i) &&
                    B.componentExists(i))
                {
                    beta*(*tensors[i].tensor) += alpha*(A(i))/
                                                         (B(i));
                }
            }
        }

        void sum(const T alpha, const T beta)
        {
            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL)
                {
                    beta*(*tensors[i].tensor) += alpha;
                }
            }
        }

        void sum(const T alpha, bool conja, const Derived& A, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors_.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL &&
                    A.componentExists(i))
                {
                    beta*(*tensors[i].tensor) += alpha*(A(i));
                }
            }
        }

        void invert(const T alpha, bool conja, const Derived& A, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors_.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL &&
                    A.componentExists(i))
                {
                    beta*(*tensors[i].tensor) += alpha/(A(i));
                }
            }
        }

        T dot(bool conja, const Derived& A, bool conjb) const
        {
            #ifdef VALIDATE_INPUTS
            if (tensors.size() != A.tensors_.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            T s = (T)0;

            for (int i = 0;i < tensors.size();i++)
            {
                if (tensors[i] != NULL &&
                    A.componentExists(i))
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
        IndexableCompositeTensor(const int ndim=0, const int ntensors=0)
        : IndexableTensorBase<Derived,T>(ndim), CompositeTensor<Derived,Base,T>(ntensors) {}

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

        void sum(const T alpha, const T beta)
        {
            Derived tensor(static_cast<const Derived&>(*this), alpha);
            beta*(*this)[this->implicit()] = tensor[""];
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
