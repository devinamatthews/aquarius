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

#include "indexabletensor.hpp"

namespace aquarius
{
namespace tensor
{

#define INHERIT_FROM_COMPOSITE_TENSOR(Derived,Base,T) \
    protected: \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::tensors_; \
    public: \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::mult; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::div; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::sum; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::invert; \
    INHERIT_FROM_TENSOR(CONCAT(Derived),T)

#define INHERIT_FROM_INDEXABLE_COMPOSITE_TENSOR(Derived,Base,T) \
    protected: \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::tensors_; \
        using aquarius::tensor::IndexableTensor< Derived, T >::ndim_; \
    public: \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::mult; \
        using aquarius::tensor::IndexableCompositeTensor< Derived, Base, T >::sum; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::div; \
        using aquarius::tensor::CompositeTensor< Derived, Base, T >::invert; \
        using aquarius::tensor::IndexableTensor< Derived, T >::mult; \
        using aquarius::tensor::IndexableTensor< Derived, T >::sum; \
        using aquarius::tensor::IndexableTensor< Derived, T >::scale; \
    INHERIT_FROM_TENSOR(CONCAT(Derived),T)

template <class Derived, class Base, class T>
class CompositeTensor : virtual public Tensor<Derived,T>
{
    protected:
        struct TensorRef
        {
            Base* tensor_;
            bool isAlloced;
            TensorRef() : tensor_(NULL), isAlloced(true) {}
            TensorRef(Base* tensor_, bool isAlloced=true)
            : tensor_(tensor_), isAlloced(isAlloced) {}
            bool operator==(const Base* other) const { return tensor_ == other; }
            bool operator!=(const Base* other) const { return tensor_ != other; }
        };

        std::vector<TensorRef> tensors_;

    public:
        CompositeTensor(const CompositeTensor<Derived,Base,T>& other)
        : Tensor<Derived,T>(*this), tensors_(other.tensors_)
        {
            for (int i = 0;i < tensors_.size();i++)
            {
                if (tensors_[i] != NULL)
                {
                    tensors_[i] = TensorRef(new Base(*tensors_[i].tensor_));
                }
            }
        }

        CompositeTensor(const int ntensors)
        : Tensor<Derived,T>(*this), tensors_(ntensors) {}

        virtual ~CompositeTensor()
        {
            for (int i = 0;i < tensors_.size();i++)
            {
                if (tensors_[i] != NULL && tensors_[i].isAlloced)
                {
                    delete tensors_[i].tensor_;
                }
            }
        }

        int getNumTensors() const { return tensors_.size(); }

        bool componentExists(const int idx) const
        {
            return tensors_[idx] != NULL;
        }

        /**********************************************************************
         *
         * Subtensor indexing
         *
         *********************************************************************/
        Base& operator()(const int idx)
        {
            if (tensors_[idx] == NULL)
                throw std::logic_error("tensor component does not exist");
            return *tensors_[idx].tensor_;
        }

        const Base& operator()(const int idx) const
        {
            if (tensors_[idx] == NULL)
                throw std::logic_error("tensor component does not exist");
            return *tensors_[idx].tensor_;
        }

        /**********************************************************************
         *
         * Implementation of Tensor stubs
         *
         *********************************************************************/
        void mult(const T alpha)
        {
            for (int i = 0;i < tensors_.size();i++)
            {
                if (tensors_[i] != NULL)
                {
                    *tensors_[i].tensor_ *= alpha;
                }
            }
        }

        void mult(const T alpha, bool conja, const Tensor<Derived,T>& A,
                                 bool conjb, const Tensor<Derived,T>& B, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors_.size() != A.getDerived().tensors_.size() ||
                tensors_.size() != B.getDerived().tensors_.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors_.size();i++)
            {
                if (tensors_[i] != NULL &&
                    A.getDerived().tensors_[i] != NULL &&
                    B.getDerived().tensors_[i] != NULL)
                {
                    beta*(*tensors_[i].tensor_) += alpha*(*A.getDerived().tensors_[i].tensor_)*
                                                         (*B.getDerived().tensors_[i].tensor_);
                }
            }
        }

        void div(const T alpha, bool conja, const Tensor<Derived,T>& A,
                                bool conjb, const Tensor<Derived,T>& B, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors_.size() != A.getDerived().tensors_.size() ||
                tensors_.size() != B.getDerived().tensors_.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors_.size();i++)
            {
                if (tensors_[i] != NULL &&
                    A.getDerived().tensors_[i] != NULL &&
                    B.getDerived().tensors_[i] != NULL)
                {
                    beta*(*tensors_[i].tensor_) += alpha*(*A.getDerived().tensors_[i].tensor_)/
                                                         (*B.getDerived().tensors_[i].tensor_);
                }
            }
        }

        void sum(const T alpha, const T beta)
        {
            for (int i = 0;i < tensors_.size();i++)
            {
                if (tensors_[i] != NULL)
                {
                    beta*(*tensors_[i].tensor_) += alpha;
                }
            }
        }

        void sum(const T alpha, bool conja, const Tensor<Derived,T>& A, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors_.size() != A.getDerived().tensors_.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors_.size();i++)
            {
                if (tensors_[i] != NULL &&
                    A.getDerived().tensors_[i] != NULL)
                {
                    beta*(*tensors_[i].tensor_) += alpha*(*A.getDerived().tensors_[i].tensor_);
                }
            }
        }

        void invert(const T alpha, bool conja, const Tensor<Derived,T>& A, const T beta)
        {
            #ifdef VALIDATE_INPUTS
            if (tensors_.size() != A.getDerived().tensors_.size()) throw LengthMismatchError();
            #endif //VALIDATE_INPUTS

            for (int i = 0;i < tensors_.size();i++)
            {
                if (tensors_[i] != NULL &&
                    A.getDerived().tensors_[i] != NULL)
                {
                    beta*(*tensors_[i].tensor_) += alpha/(*A.getDerived().tensors_[i].tensor_);
                }
            }
        }
};

template <class Derived, class Base, class T>
class IndexableCompositeTensor : public IndexableTensor<Derived,T>, public CompositeTensor<Derived,Base,T>
{
    public:
        IndexableCompositeTensor(const IndexableCompositeTensor<Derived,Base,T>& other)
        : Tensor<Derived,T>(*this), IndexableTensor<Derived,T>(other), CompositeTensor<Derived,Base,T>(other) {}

        IndexableCompositeTensor(const int ndim=0, const int ntensors=0)
        : Tensor<Derived,T>(*this), IndexableTensor<Derived,T>(ndim), CompositeTensor<Derived,Base,T>(ntensors) {}

        virtual ~IndexableCompositeTensor() {}

        void mult(const T alpha)
        {
            IndexableTensor<Derived,T>::mult(alpha);
        }

        void mult(const T alpha, bool conja, const Tensor<Derived,T>& A,
                                 bool conjb, const Tensor<Derived,T>& B, const T beta)
        {
            IndexableTensor<Derived,T>::mult(alpha, conja, A, conjb, B, beta);
        }

        void sum(const T alpha, const T beta)
        {
            IndexableTensor<Derived,T>::sum(alpha, beta);
        }

        void sum(const T alpha, bool conja, const Tensor<Derived,T>& A, const T beta)
        {
            IndexableTensor<Derived,T>::sum(alpha, conja, A, beta);
        }
};

template <class Derived, class T>
struct Scalar<TensorMult<Derived,T>, typename std::enable_if<(sizeof(&Derived::getNumTensors)>0) && !std::is_base_of<IndexableTensor<Derived,T>,Derived>::value>::type>
{
    static T value(const TensorMult<Derived,T>& other)
    {
        T s = (T)0;

        for (int i = 0;i < other.A_.tensor_.getNumTensors();i++)
        {
            if (other.A_.tensor_.componentExists(i) &&
                other.B_.tensor_.componentExists(i))
            {
                s += scalar(other.A_.tensor_(i)*
                            other.B_.tensor_(i));
            }
        }

        return s;
    }
};

}
}

#endif
