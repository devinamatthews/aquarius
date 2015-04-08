#ifndef _AQUARIUS_TENSOR_HPP_
#define _AQUARIUS_TENSOR_HPP_

#include "util/global.hpp"

namespace aquarius
{
namespace tensor
{

template <class Derived, class T> class Tensor;
template <class Derived, class T> class ScaledTensor;
template <class Derived, class T> class InvertedTensor;
template <class Derived, class T> class TensorMult;
template <class Derived, class T> class TensorDiv;

class TensorError;
class OutOfBoundsError;
class LengthMismatchError;
class IndexMismatchError;
class InvalidNdimError;
class InvalidLengthError;
class InvalidLdError;
class LdTooSmallError;
class SymmetryMismatchError;
class InvalidSymmetryError;
class InvalidStartError;

#define INHERIT_FROM_TENSOR(Derived,T) \
    public: \
        using aquarius::tensor::Tensor< Derived,T >::getDerived; \
        using aquarius::tensor::Tensor< Derived,T >::operator=; \
        using aquarius::tensor::Tensor< Derived,T >::operator+=; \
        using aquarius::tensor::Tensor< Derived,T >::operator-=; \
        using aquarius::tensor::Tensor< Derived,T >::operator*=; \
        using aquarius::tensor::Tensor< Derived,T >::operator/=; \
        using aquarius::tensor::Tensor< Derived,T >::operator*; \
        using aquarius::tensor::Tensor< Derived,T >::operator/; \
        using aquarius::tensor::Tensor< Derived,T >::name; \
        Derived & operator=(const Derived & other) \
        { \
            sum((T)1, false, other, (T)0); \
            return *this; \
        } \
    private:

template <class Derived, typename T>
class Tensor
{
    public:
        typedef T dtype;
        const string name;
        
        Tensor(const string& name) : name(name) {}

        virtual ~Tensor() {}

        Derived& getDerived() { return static_cast<Derived&>(*this); }

        const Derived& getDerived() const { return static_cast<const Derived&>(*this); }

        /**********************************************************************
         *
         * Operators with scalars
         *
         *********************************************************************/
        Derived& operator=(const T val)
        {
            sum(val, (T)0);
            return getDerived();
        }

        Derived& operator+=(const T val)
        {
            sum(val, (T)1);
            return getDerived();
        }

        Derived& operator-=(const T val)
        {
            sum(-val, (T)1);
            return getDerived();
        }

        Derived& operator*=(const T val)
        {
            mult(val);
            return getDerived();
        }

        Derived& operator/=(const T val)
        {
            mult(1.0/val);
            return getDerived();
        }

        /**********************************************************************
         *
         * Binary operations (multiplication and division)
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator=(const TensorMult<cvDerived,T>& other)
        {
            mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)0);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator+=(const TensorMult<cvDerived,T>& other)
        {
            mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)1);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator-=(const TensorMult<cvDerived,T>& other)
        {
            mult(-other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)1);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator=(const TensorDiv<cvDerived,T>& other)
        {
            div(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)0);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator+=(const TensorDiv<cvDerived,T>& other)
        {
            div(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)1);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator-=(const TensorDiv<cvDerived,T>& other)
        {
            div(-other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)1);
            return getDerived();
        }

        /**********************************************************************
         *
         * Unary operations (assignment, summation, multiplication, and division)
         *
         *********************************************************************/

        /*
         * Automatically generated assignment operator exists in spite of
         * (and precedes) template assignment operator, so a non-template
         * one must be given
         */
        Derived& operator=(const Derived& other)
        {
            sum((T)1, false, other, (T)0);
            return getDerived();
        }

        template <typename cvDerived> if_exists_t<typename cvDerived::dtype, Derived&>
        //ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator=(cvDerived& other)
        {
            sum((T)1, false, other, (T)0);
            return getDerived();
        }

        template <typename cvDerived> if_exists_t<typename cvDerived::dtype, Derived&>
        //ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator+=(cvDerived& other)
        {
            sum((T)1, false, other, (T)1);
            return getDerived();
        }

        template <typename cvDerived> if_exists_t<typename cvDerived::dtype, Derived&>
        //ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator-=(cvDerived& other)
        {
            sum((T)(-1), false, other, (T)1);
            return getDerived();
        }

        template <typename cvDerived> if_exists_t<typename cvDerived::dtype, Derived&>
        //ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator*=(cvDerived& other)
        {
            mult((T)1, false, getDerived(), false, other, (T)0);
            return getDerived();
        }

        template <typename cvDerived> if_exists_t<typename cvDerived::dtype, Derived&>
        //ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator/=(cvDerived& other)
        {
            div((T)1, false, getDerived(), false, other, (T)0);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator=(const ScaledTensor<cvDerived,T>& other)
        {
            sum(other.factor_, other.conj_, other.tensor_, (T)0);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator+=(const ScaledTensor<cvDerived,T>& other)
        {
            sum(other.factor_, other.conj_, other.tensor_, (T)1);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator-=(const ScaledTensor<cvDerived,T>& other)
        {
            sum(-other.factor_, other.conj_, other.tensor_, (T)1);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator*=(const ScaledTensor<cvDerived,T>& other)
        {
            mult(other.factor_, false, getDerived(), other.conj_, other.tensor_, (T)0);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator/=(const ScaledTensor<cvDerived,T>& other)
        {
            div((T)1/other.factor_, false, getDerived(), other.conj_, other.tensor_, (T)0);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator=(const InvertedTensor<cvDerived,T>& other)
        {
            invert(other.factor_, other.conj_, other.tensor_, (T)0);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator+=(const InvertedTensor<cvDerived,T>& other)
        {
            invert(other.factor_, other.conj_, other.tensor_, (T)1);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator-=(const InvertedTensor<cvDerived,T>& other)
        {
            invert(-other.factor_, other.conj_, other.tensor_, (T)0);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator*=(const InvertedTensor<cvDerived,T>& other)
        {
            div(other.factor_, false, getDerived(), other.conj_, other.tensor_, (T)0);
            return getDerived();
        }

        ENABLE_IF_SAME(Derived,cvDerived,Derived&)
        operator/=(const InvertedTensor<cvDerived,T>& other)
        {
            mult((T)1/other.factor_, false, getDerived(), other.conj_, other.tensor_, (T)0);
            return getDerived();
        }

        /**********************************************************************
         *
         * Intermediate operations
         *
         *********************************************************************/
        friend ScaledTensor<Derived,T> operator*(const T factor, Derived& other)
        {
            return ScaledTensor<Derived,T>(other.getDerived(), factor);
        }

        friend ScaledTensor<const Derived,T> operator*(const T factor, const Derived& other)
        {
            return ScaledTensor<const Derived,T>(other.getDerived(), factor);
        }

        ScaledTensor<Derived,T> operator*(const T factor)
        {
            return ScaledTensor<Derived,T>(getDerived(), factor);
        }

        ScaledTensor<const Derived,T> operator*(const T factor) const
        {
            return ScaledTensor<const Derived,T>(getDerived(), factor);
        }

        friend InvertedTensor<const Derived,T> operator/(const T factor, const Derived& other)
        {
            return InvertedTensor<const Derived,T>(other.getDerived(), factor);
        }

        ScaledTensor<Derived,T> operator/(const T factor)
        {
            return ScaledTensor<Derived,T>(getDerived(), (T)1/factor);
        }

        ScaledTensor<const Derived,T> operator/(const T factor) const
        {
            return ScaledTensor<const Derived,T>(getDerived(), (T)1/factor);
        }

        ScaledTensor<Derived,T> operator-()
        {
            return ScaledTensor<Derived,T>(getDerived(), (T)(-1));
        }

        ScaledTensor<const Derived,T> operator-() const
        {
            return ScaledTensor<const Derived,T>(getDerived(), (T)(-1));
        }

        friend ScaledTensor<const Derived,T> conj(const Derived& t)
        {
            return ScaledTensor<const Derived,T>(t.getDerived(), (T)1, true);
        }

        template <typename cvDerived> if_exists_t<typename cvDerived::dtype, TensorMult<Derived,T>>
        //ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorMult<Derived,T>))
        operator*(const cvDerived& other) const
        {
            return TensorMult<Derived,T>(ScaledTensor<const Derived,T>(getDerived(), (T)1),
                                         ScaledTensor<const Derived,T>(other.getDerived(), (T)1));
        }

        template <typename cvDerived> if_exists_t<typename cvDerived::dtype, TensorDiv<Derived,T>>
        //ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorDiv<Derived,T>))
        operator/(const cvDerived& other) const
        {
            return TensorDiv<Derived,T>(ScaledTensor<const Derived,T>(getDerived(), (T)1),
                                        ScaledTensor<const Derived,T>(other.getDerived(), (T)1));
        }

        /**********************************************************************
         *
         * Stubs
         *
         *********************************************************************/

        /*
         * this = alpha*this + beta*A*B
         */
        virtual void mult(const T alpha, bool conja, const Derived& A,
                                         bool conjb, const Derived& B, const T beta) = 0;

        /*
         * this = alpha*this
         */
        virtual void mult(const T alpha) = 0;

        /*
         * this = alpha*this + beta*A/B
         */
        virtual void div(const T alpha, bool conja, const Derived& A,
                                        bool conjb, const Derived& B, const T beta) = 0;

        /*
         * this = alpha*this + beta*A
         */
        virtual void sum(const T alpha, bool conja, const Derived& A, const T beta) = 0;

        /*
         * this = alpha*this + beta
         */
        virtual void sum(const T alpha, const T beta) = 0;

        /*
         * this = alpha*this + beta/A
         */
        virtual void invert(const T alpha, bool conja, const Derived& A, const T beta) = 0;

        /*
         * scalar = A*this
         */
        virtual T dot(bool conja, const Derived& A, bool conjb) const = 0;
};

template <class Derived, typename T>
class ScaledTensor
{
    public:
        Derived& tensor_;
        T factor_;
        bool conj_;

        template <typename cvDerived>
        ScaledTensor(const ScaledTensor<cvDerived,T>& other)
        : tensor_(other.tensor_), factor_(other.factor_), conj_(other.conj_) {}

        ScaledTensor(Derived& tensor, const T factor, const bool conj_=false)
        : tensor_(tensor), factor_(factor), conj_(conj_) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        ScaledTensor<Derived,T> operator-() const
        {
            ScaledTensor<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return ret;
        }

        friend ScaledTensor<const Derived,T> conj(const ScaledTensor<Derived,T>& st)
        {
            ScaledTensor<Derived,T> ret(st);
            ret.conj_ = !ret.conj_;
            return ret;
        }

        /**********************************************************************
         *
         * Unary tensor operations
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const cvDerived& other)
        {
            tensor_.sum((T)1, false, other, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const cvDerived& other)
        {
            tensor_.sum((T)1, false, other, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const cvDerived& other)
        {
            tensor_.sum((T)(-1), false, other, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator*=(const cvDerived& other)
        {
            tensor_.mult(factor_, false, tensor_, false, other, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator/=(const cvDerived& other)
        {
            tensor_.div(factor_, false, tensor_, false, other, (T)0);
            return *this;
        }

        ScaledTensor<Derived,T>& operator=(const ScaledTensor<Derived,T>& other)
        {
            tensor_.sum(other.factor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.sum(other.factor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.sum(other.factor_, other.conj_, other.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.sum(-other.factor_, other.conj_, other.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator*=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.mult(factor_*other.factor_, false, tensor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator/=(const ScaledTensor<cvDerived,T>& other)
        {
            tensor_.div(factor_/other.factor_, false, tensor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.invert(other.factor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.invert(other.factor_, other.conj_, other.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.invert(-other.factor_, other.conj_, other.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator*=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.div(factor_*other.factor_, false, tensor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator/=(const InvertedTensor<cvDerived,T>& other)
        {
            tensor_.mult(factor_/other.factor_, false, tensor_, other.conj_, other.tensor_, (T)0);
            return *this;
        }

        /**********************************************************************
         *
         * Binary tensor operations
         *
         *********************************************************************/
        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const TensorMult<cvDerived,T>& other)
        {
            tensor_.mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const TensorMult<cvDerived,T>& other)
        {
            tensor_.mult(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const TensorMult<cvDerived,T>& other)
        {
            tensor_.mult(-other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator=(const TensorDiv<cvDerived,T>& other)
        {
            tensor_.div(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, (T)0);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator+=(const TensorDiv<cvDerived,T>& other)
        {
            tensor_.div(other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(ScaledTensor<Derived,T>&))
        operator-=(const TensorDiv<cvDerived,T>& other)
        {
            tensor_.div(-other.factor_, other.A_.conj_, other.A_.tensor_, other.B_.conj_, other.B_.tensor_, factor_);
            return *this;
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorMult<Derived,T>))
        operator*(const ScaledTensor<cvDerived,T>& other) const
        {
            return TensorMult<Derived,T>(*this, other);
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorMult<Derived,T>))
        operator*(const cvDerived& other) const
        {
            return TensorMult<Derived,T>(*this, ScaledTensor<const Derived,T>(other.getDerived(), (T)1));
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorDiv<Derived,T>))
        operator/(const ScaledTensor<cvDerived,T>& other) const
        {
            return TensorDiv<Derived,T>(*this, other);
        }

        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorDiv<Derived,T>))
        operator/(const cvDerived& other) const
        {
            return TensorDiv<Derived,T>(*this, ScaledTensor<const Derived,T>(other.getDerived(), (T)1));
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        ScaledTensor<Derived,T> operator*(const T factor) const
        {
            ScaledTensor<Derived,T> it(*this);
            it.factor_ *= factor;
            return it;
        }

        friend ScaledTensor<Derived,T> operator*(const T factor, const ScaledTensor<Derived,T>& other)
        {
            return other*factor;
        }

        ScaledTensor<Derived,T> operator/(const T factor) const
        {
            ScaledTensor<Derived,T> it(*this);
            it.factor_ /= factor;
            return it;
        }

        friend InvertedTensor<Derived,T> operator/(const T factor, const ScaledTensor<Derived,T>& other)
        {
            return InvertedTensor<Derived,T>(other.tensor_, factor/other.factor_);
        }

        ScaledTensor<Derived,T>& operator=(const T val)
        {
            tensor_.sum(val, (T)0);
            return *this;
        }

        ScaledTensor<Derived,T>& operator+=(const T val)
        {
            tensor_.sum(val, factor_);
            return *this;
        }

        ScaledTensor<Derived,T>& operator-=(const T val)
        {
            tensor_.sum(-val, factor_);
            return *this;
        }

        ScaledTensor<Derived,T>& operator*=(const T val)
        {
            tensor_.mult(val);
            return *this;
        }

        ScaledTensor<Derived,T>& operator/=(const T val)
        {
            tensor_.mult((T)1/val);
            return *this;
        }
};

template <class Derived1, class Derived2, class T>
//enable_if_t<is_same<const Derived1, const Derived2>::value,TensorMult<Derived1,T>>
TensorMult<Derived1,T>
operator*(const Derived1& t1, const ScaledTensor<Derived2,T>& t2)
{
    return TensorMult<Derived1,T>(ScaledTensor<const Derived1,T>(t1.getDerived(), (T)1), t2);
}

template <class Derived1, class Derived2, class T>
//enable_if_t<is_same<const Derived1, const Derived2>::value,TensorDiv<Derived1,T>>
TensorDiv<Derived1,T>
operator/(const Derived1& t1, const ScaledTensor<Derived2,T>& t2)
{
    return TensorDiv<Derived1,T>(ScaledTensor<const Derived1,T>(t1.getDerived(), (T)1), t2);
}

template <class Derived, typename T>
class InvertedTensor
{
    private:
        const InvertedTensor& operator=(const InvertedTensor<Derived,T>& other);

    public:
        Derived& tensor_;
        T factor_;
        bool conj_;

        InvertedTensor(Derived& tensor, const T factor, const bool conj_=false)
        : tensor_(tensor), factor_(factor), conj_(conj_) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        InvertedTensor<Derived,T> operator-() const
        {
            InvertedTensor<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return *this;
        }

        friend InvertedTensor<Derived,T> conj(const InvertedTensor<Derived,T>& tm)
        {
            InvertedTensor<Derived,T> ret(tm);
            ret.conj_ = !ret.conj_;
            return ret;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        InvertedTensor<Derived,T> operator*(const T factor) const
        {
            InvertedTensor<Derived,T> ret(*this);
            ret.factor_ *= factor;
            return ret;
        }

        InvertedTensor<Derived,T> operator/(const T factor) const
        {
            InvertedTensor<Derived,T> ret(*this);
            ret.factor_ /= factor;
            return ret;
        }

        friend InvertedTensor<Derived,T> operator*(const T factor, const InvertedTensor<Derived,T>& other)
        {
            return other*factor;
        }
};

template <class Derived, typename T>
class TensorMult
{
    private:
        const TensorMult& operator=(const TensorMult<Derived,T>& other);

    public:
        ScaledTensor<const Derived,T> A_;
        ScaledTensor<const Derived,T> B_;
        T factor_;

        template <class Derived1, class Derived2>
        TensorMult(const ScaledTensor<Derived1,T>& A, const ScaledTensor<Derived2,T>& B)
        : A_(A), B_(B), factor_(A.factor_*B.factor_) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        TensorMult<Derived,T> operator-() const
        {
            TensorMult<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return ret;
        }

        friend TensorMult<Derived,T> conj(const TensorMult<Derived,T>& tm)
        {
            TensorMult<Derived,T> ret(tm);
            ret.A_.conj_ = !ret.A_.conj_;
            ret.B_.conj_ = !ret.B_.conj_;
            return ret;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        TensorMult<Derived,T> operator*(const T factor) const
        {
            TensorMult<Derived,T> ret(*this);
            ret.factor_ *= factor;
            return ret;
        }

        TensorMult<Derived,T> operator/(const T factor) const
        {
            TensorMult<Derived,T> ret(*this);
            ret.factor_ /= factor;
            return ret;
        }

        friend TensorMult<Derived,T> operator*(const T factor, const TensorMult<Derived,T>& other)
        {
            return other*factor;
        }
};

template <class Derived, typename T>
class TensorDiv
{
    private:
        const TensorDiv& operator=(const TensorDiv<Derived,T>& other);

    public:
        ScaledTensor<const Derived,T> A_;
        ScaledTensor<const Derived,T> B_;
        T factor_;

        template <class Derived1, class Derived2>
        TensorDiv(const ScaledTensor<Derived1,T>& A, const ScaledTensor<Derived2,T>& B)
        : A_(A), B_(B), factor_(A.factor_/B.factor_) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/
        TensorDiv<Derived,T> operator-() const
        {
            TensorDiv<Derived,T> ret(*this);
            ret.factor_ = -ret.factor_;
            return ret;
        }

        friend TensorDiv<Derived,T> conj(const TensorDiv<Derived,T>& tm)
        {
            TensorDiv<Derived,T> ret(tm);
            ret.A_.conj_ = !ret.A_.conj_;
            ret.B_.conj_ = !ret.B_.conj_;
            return ret;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/
        TensorDiv<Derived,T> operator*(const T factor) const
        {
            TensorDiv<Derived,T> ret(*this);
            ret.factor_ *= factor;
            return ret;
        }

        TensorDiv<Derived,T> operator/(const T factor) const
        {
            TensorDiv<Derived,T> ret(*this);
            ret.factor_ /= factor;
            return ret;
        }

        friend TensorDiv<Derived,T> operator*(const T factor, const TensorDiv<Derived,T>& other)
        {
            return other*factor;
        }
};

class TensorError : public exception
{
    public:
        virtual const char* what() const throw() = 0;
};

class OutOfBoundsError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "out-of-bounds read or write"; }
};

class LengthMismatchError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "length mismatch error"; }
};

class IndexMismatchError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "index mismatch error"; }
};

class InvalidNdimError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid number of dimensions"; }
};

class InvalidLengthError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid length"; }
};

class InvalidLdError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid leading dimension"; }
};

class LdTooSmallError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "leading dimension is too small"; }
};

class SymmetryMismatchError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "symmetry mismatch error"; }
};

class InvalidSymmetryError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid symmetry value"; }
};

class InvalidStartError : public TensorError
{
    public:
        virtual const char* what() const throw() { return "invalid start value"; }
};

}

template <class Derived, typename T>
T scalar(const tensor::TensorMult<Derived,T>& tm)
{
    return tm.factor_*tm.B_.tensor_.dot(tm.A_.conj_, tm.A_.tensor_, tm.B_.conj_);
}

}

#endif
