#ifndef _AQUARIUS_FRAMEWORKS_TENSOR_HPP_
#define _AQUARIUS_FRAMEWORKS_TENSOR_HPP_

#include "frameworks/util.hpp"

#include "tensor/ring.hpp"

namespace aquarius
{

typedef uint64_t capability_type;

enum Capability : capability_type
{
    CONST_       = 0x8000000000000000u,
    DIVISIBLE_   = 0x4000000000000000u,
    INDEXABLE_   = 0x2000000000000000u,
    BOUNDED_     = 0x1000000000000000u,
    DISTRIBUTED_ = 0x0800000000000000u,
    IPSYMMETRIC_ = 0x0400000000000000u,
    PGSYMMETRIC_ = 0x0200000000000000u,
    SPINORBITAL_ = 0x0100000000000000u,

    DIVISIBLE    = DIVISIBLE_,
    INDEXABLE    = INDEXABLE_,
    BOUNDED      = (BOUNDED_|INDEXABLE),
    DISTRIBUTED  = DISTRIBUTED_,
    IPSYMMETRIC  = (IPSYMMETRIC_|INDEXABLE),
    PGSYMMETRIC  = PGSYMMETRIC_,
    SPINORBITAL  = (SPINORBITAL_|IPSYMMETRIC|BOUNDED),
};

/*******************************************************************************
 *
 * Macro Helpers
 *
 ******************************************************************************/

#define IS_POWER_OF_TWO(x) ((x) > 0 && ((x)&(x-1)) == 0)
#define LOW_BIT(x) ((x)^((x)&((x)-1)))

#define ARE_DISTINCT(a,b) (((a)&(b))==0)
#define IS_SUPERSET_OF(a,b) (((a)&(b))==(b))
#define IS_STRICT_SUPERSET_OF(a,b) (IS_SUPERSET_OF(a,b)&&((a)!=(b)))

#define TREE_BASE(TreeBase,Base,Root,cargs,pargs,bargs,...) \
template <capability_type C, capability_type C_=C, capability_type C__=0> \
struct TreeBase : public TreeBase<C, (C_-C__), LOW_BIT(C_-C__)>, Base<C__> \
{ TreeBase cargs : TreeBase<C, (C_-C__), LOW_BIT(C_-C__)> pargs , Base<C__> bargs __VA_ARGS__ }; \
\
template <capability_type C> \
struct TreeBase<C, 0, 0> : public Root<C> { TreeBase cargs : Root<C> pargs {} }; \
\
template <capability_type C, capability_type C_> \
struct TreeBase<C, C_, 0> : public TreeBase<C, C_, LOW_BIT(C_)> \
{ TreeBase cargs : TreeBase<C, C_, LOW_BIT(C_)> pargs {} };

#define INLINE_BASE(InlineBase,Base,Root,cargs,pargs,...) \
template <capability_type C, capability_type C_=C, capability_type C__=0> \
struct InlineBase : public Base<C, C__, InlineBase<C, (C_-C__), LOW_BIT(C_-C__)>> \
{ InlineBase cargs : Base<C, C__, InlineBase<C, (C_-C__), LOW_BIT(C_-C__)>> pargs __VA_ARGS__ }; \
\
template <capability_type C> \
struct InlineBase<C, 0, 0> : public Root<C> { InlineBase cargs : Root<C> pargs {} }; \
\
template <capability_type C, capability_type C_> \
struct InlineBase<C, C_, 0> : public InlineBase<C, C_, LOW_BIT(C_)> \
{ InlineBase cargs : InlineBase<C, C_, LOW_BIT(C_)> pargs {} };

#define DONT_NEED_INITIALIZATION (CONST_|DIVISIBLE)
#define INITIALIZER_TYPE(C) TensorInitializerList<((C)&(~DONT_NEED_INITIALIZATION))>

#define TENSOR_INTERFACE(C_) \
template <> class TensorInterface<C_> : private VirtualBase

#define TENSOR_DEFINITION(C_) \
template <capability_type C, class Base> \
class TensorDefinition<C, C_, Base> : public Base

#define TENSOR_WRAPPER(C_) \
template <capability_type C, class Base> \
class TensorWrapper<C, C_, Base> : public Base

/*******************************************************************************
 *
 * Forward Declarations
 *
 ******************************************************************************/

namespace tensor
{

class ConstScaledTensor;
class ScaledTensor;
class TensorMult;

template <capability_type C=0> class Tensor;
template <capability_type C=0> class ConstTensor;

template <capability_type C=0, typename=void> class TensorImplementation;

template <capability_type C=0> class TensorInitializer {};
template <capability_type C> class TensorInitializerList;

template <capability_type C> class TensorInterface {};

template <capability_type C, capability_type C_, class Base> class TensorDefinition : public Base
{
    public:
        TensorDefinition(const INITIALIZER_TYPE(C)& ilist) : Base(ilist) {}
};

template <capability_type C, capability_type C_, class Base> class TensorWrapper : public Base {};

struct VirtualBase
{
    virtual ~VirtualBase() {}
};

/**************************************************************************************************
 *
 * Tensor Initialization
 *
 *************************************************************************************************/

/*
 * Class to hold initialization data for the base Tensor class.
 */
class TensorInitializer_
{
    public:
        const string name;
        const Field F;
        const Ring R;

        TensorInitializer_(const string& name, Field F, Ring R)
        : name(name), F(F), R(R) {}
};

/*
 * Generic initializer type which stores any number of specialized initializers.
 * This type is constructed from TensorInitializer by concatenation with the << operator.
 */
template <capability_type C=0> class TensorInitializerList
: public map<capability_type, any>
{
    template <capability_type C_> friend class TensorInitializerList;

    protected:
        template <capability_type C_>
        enable_if_t<IS_POWER_OF_TWO(C_)>
        addInitializer(const TensorInitializer<C_>& init)
        {
            (*this)[C_] = init;
        }

        template <capability_type C_,
                  capability_type C__=C_,
                  capability_type C___=0>
        struct AddInitializer : AddInitializer<C_, (C__-C___), LOW_BIT(C__-C___)>
        {
            AddInitializer(TensorInitializerList& parent, const TensorInitializerList<C_>& ilist)
            : AddInitializer<C_, (C__-C___), LOW_BIT(C__-C___)>(parent, ilist)
            {
                parent.addInitializer(ilist.at(C___).template get<TensorInitializer<C___>>());
            }

            AddInitializer(TensorInitializerList& parent, const TensorImplementation<C_>& impl)
            : AddInitializer<C_, (C__-C___), LOW_BIT(C__-C___)>(parent, impl)
            {
                parent.addInitializer(static_cast<const TensorInitializer<C___>&>(impl));
            }
        };

        template <capability_type C_>
        struct AddInitializer<C_, 0, 0>
        {
            AddInitializer(TensorInitializerList& parent, const TensorInitializerList<C_>& ilist)
            {
                auto it1 = parent.find(0);
                auto it2 = ilist.find(0);

                if (it2 != ilist.end())
                {
                    assert(it1 == parent.end());
                    parent[0] = ilist.at(0);
                }
            }

            AddInitializer(TensorInitializerList& parent, const TensorImplementation<C_>& impl)
            {
                parent[0] = static_cast<const TensorInitializer_&>(impl);
            }
        };

        template <capability_type C_, capability_type C__>
        struct AddInitializer<C_, C__, 0> : AddInitializer<C_, C__, LOW_BIT(C__)>
        {
            AddInitializer(TensorInitializerList& parent, const TensorInitializerList<C_>& ilist)
            : AddInitializer<C_, C__, LOW_BIT(C__)>(parent, ilist) {}

            AddInitializer(TensorInitializerList& parent, const TensorImplementation<C_>& impl)
            : AddInitializer<C_, C__, LOW_BIT(C__)>(parent, impl) {}
        };

        template <capability_type C_>
        void addInitializer(const TensorInitializerList<C_>& ilist)
        {
            AddInitializer<C_>(*this, ilist);
        }

    public:
        TensorInitializerList() {}

        /*
         * Allow initialization of a new initializer list from another which has a
         * subset of capabilities.
         */
        template <capability_type C_>
        TensorInitializerList(const TensorInitializerList<C_>& other, enable_if_t<IS_STRICT_SUPERSET_OF(C,C_)>* = nullptr)
        : map<capability_type, any>(other) {}

        /*
         * Allow initialization of a new initializer list from another which has
         * the same set or a superset of capabilities.
         */
        template <capability_type C_>
        TensorInitializerList(const TensorInitializerList<C_>& other, enable_if_t<IS_SUPERSET_OF(C_,C)>* = nullptr)
        {
            for (auto& i : other)
            {
                if (!i.first || !ARE_DISTINCT(i.first,C)) (*this)[i.first] = i.second;
            }
        }

        /*
         * Allow initialization of a new initializer list from an existing
         * implementation which has the same or a superset of capabilities.
         */
        template <capability_type C_, typename=enable_if_t<IS_SUPERSET_OF(C_,C)>>
        TensorInitializerList(const TensorImplementation<C_>& other)
        {
            AddInitializer<C_>(*this, other);
            for (capability_type C__ = CONST_;C__ > 0;C__ >>= 1)
            {
                if (ARE_DISTINCT(C__,C)) erase(C__);
            }
        }

        /*
         * Additional initializers may be added to the list by the << operator,
         * although the new initializer may not already be present.
         *
         * Checking for IS_POWER_OF_TWO(C_) ensures that the new initializar represents
         * only one capability, and that this specialization is not selected for something
         * like class TensorInitializer<A|B> : TensorInitializerList<A|B>.
         */
        template <capability_type C_>
        enable_if_t<IS_POWER_OF_TWO(C_)&&ARE_DISTINCT(C,C_),TensorInitializerList<C|C_>>
        operator<<(const TensorInitializer<C_>& init) const
        {
            /*
             * Can't use ret.swap(*this) since *this is const
             * (a non-member function taking rvalue-refs would be ideal).
             */
            TensorInitializerList<C|C_> ret(*this);
            ret.addInitializer(init);
            return ret;
        }

        /*
         * Initializer lists may also be concatenated providing they
         * provide distinct capabilities.
         */
        template <capability_type C_>
        enable_if_t<ARE_DISTINCT(C,C_),TensorInitializerList<C|C_>>
        operator<<(const TensorInitializerList<C_>& ilist) const
        {
            TensorInitializerList<C|C_> ret(*this);
            typename TensorInitializerList<C|C_>::template AddInitializer<C_>(ret, ilist);
            return ret;
        }

        template <capability_type C_=0>
        TensorInitializer<C_>& as()
        {
            return this->at(C_).template get<TensorInitializer<C_>>();
        }

        template <capability_type C_=0>
        const TensorInitializer<C_>& as() const
        {
            return const_cast<TensorInitializerList<C>&>(*this).as<C_>();
        }
};

/*
 * Convenience function since << should be commutative for this purpose.
 */
template <capability_type C1, capability_type C2>
enable_if_t<IS_POWER_OF_TWO(C1),TensorInitializerList<C1|C2>>
operator<<(const TensorInitializer<C1>& init, const TensorInitializerList<C2>& ilist)
{
    return ilist << init;
}

/*
 * Convenience function.
 */
template <capability_type C1, capability_type C2>
enable_if_t<IS_POWER_OF_TWO(C1)&&IS_POWER_OF_TWO(C2),TensorInitializerList<C1|C2>>
operator<<(const TensorInitializer<C1>& init1, const TensorInitializer<C2>& init2)
{
    return TensorInitializerList<>() << init1 << init2;
}

/*
 * Thin wrapper which just creates a TensorInitializer_. Since this derives from
 * TensorInitializerList<> with a public constructor, this is the only way to
 * generate an initializer list, ensuring this information is present (alas, it
 * does not guarantee that it is not provided multiple times, in which case the behavior
 * is undefined).
 */
template <> class TensorInitializer<> : public TensorInitializerList<>
{
    public:
        TensorInitializer() : TensorInitializerList<>()
        {
            Field F(Field::DOUBLE);
            (*this)[0] = make_any<TensorInitializer_>("?", F, F);
        }

        TensorInitializer(Field F) : TensorInitializerList<>()
        {
            (*this)[0] = make_any<TensorInitializer_>("?", F, F);
        }

        TensorInitializer(Field F, Ring R) : TensorInitializerList<>()
        {
            (*this)[0] = make_any<TensorInitializer_>("?", F, R);
        }

        TensorInitializer(const string& name) : TensorInitializerList<>()
        {
            Field F(Field::DOUBLE);
            (*this)[0] = make_any<TensorInitializer_>(name, F, F);
        }

        TensorInitializer(const string& name, Field F) : TensorInitializerList<>()
        {
            (*this)[0] = make_any<TensorInitializer_>(name, F, F);
        }

        TensorInitializer(const string& name, Field F, Ring R) : TensorInitializerList<>()
        {
            (*this)[0] = make_any<TensorInitializer_>(name, F, R);
        }
};

/**************************************************************************************************
 *
 * Tensor Implementation Base Specialization
 *
 *************************************************************************************************/

/*
 * No-capability specialization; this provides the operations universal to all
 * tensor types
 *
 * This also serves as a type-erased base class for return values and
 * capability-agnostic storage
 */
template <>
class TensorImplementation<> : public TensorInitializer_
{
    template <capability_type C_> friend class ConstTensor;

    protected:
        map<capability_type,void*> ptr;

        virtual Tensor<> construct(const map<capability_type,any>& ilist) = 0;

    public:
        const capability_type C;

        template <capability_type C__>
        TensorImplementation(capability_type C_, const TensorInitializerList<C__>& ilist)
        : TensorInitializer_(ilist.at(0).template get<TensorInitializer_>()),
          C(C_)
        {
            assert(R == F); //for now
        }

        virtual ~TensorImplementation() {}

        template <capability_type C>
        enable_if_t<IS_POWER_OF_TWO(C),TensorInterface<C>&> as()
        {
            return const_cast<TensorInterface<C>&>(const_cast<const TensorImplementation<>&>(*this).as<C>());
        }

        template <capability_type C>
        enable_if_t<IS_POWER_OF_TWO(C),const TensorInterface<C>&> as() const
        {
            auto i = ptr.find(C);
            assert(i != ptr.end());
            return *static_cast<const TensorInterface<C>*>(i->second);
        }

        /*
         * this = a
         */
        TensorImplementation<>& operator=(const Scalar& a);

        /*
         * this = beta*this + alpha*A*B
         */
        virtual void mult(const Scalar& alpha, bool conja, const TensorImplementation<>& A,
                                               bool conjb, const TensorImplementation<>& B, const Scalar& beta) = 0;

        /*
         * this = beta*this + alpha*A
         */
        virtual void sum(const Scalar& alpha, bool conja, const TensorImplementation<>& A, const Scalar& beta) = 0;

        /*
         * this = alpha*this
         */
        virtual void scale(const Scalar& alpha) = 0;

        /*
         * scalar = A*this
         */
        virtual Scalar dot(bool conja, const TensorImplementation<>& A, bool conjb) const = 0;
};

/**************************************************************************************************
 *
 * Intermediate Wrappers
 *
 *************************************************************************************************/

class TensorMult
{
    private:
        const TensorMult& operator=(const TensorMult& other);

    public:
        const TensorImplementation<>& A;
        const TensorImplementation<>& B;
        bool conja;
        bool conjb;
        Scalar factor;

        TensorMult(const TensorImplementation<>& A, const TensorImplementation<>& B,
                   bool conja, bool conjb, const Scalar& factor)
        : A(A), B(B), conja(conja), conjb(conjb), factor(factor) {}

        template <capability_type C>
        TensorMult(const ConstTensor<C>& A, const TensorImplementation<>& B,
                   bool conja, bool conjb, const Scalar& factor)
        : A(A.impl()), B(B), conja(conja), conjb(conjb), factor(factor) {}

        template <capability_type C>
        TensorMult(const TensorImplementation<>& A, const ConstTensor<C>& B,
                   bool conja, bool conjb, const Scalar& factor)
        : A(A), B(B.impl()), conja(conja), conjb(conjb), factor(factor) {}

        template <capability_type C1, capability_type C2>
        TensorMult(const ConstTensor<C1>& A, const ConstTensor<C2>& B,
                   bool conja, bool conjb, const Scalar& factor)
        : A(A.impl()), B(B.impl()), conja(conja), conjb(conjb), factor(factor) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/

        TensorMult operator-() const
        {
            TensorMult ret(*this);
            ret.factor = -ret.factor;
            return ret;
        }

        TensorMult conj() const
        {
            TensorMult ret(*this);
            ret.conja = !ret.conja;
            ret.conjb = !ret.conjb;
            return ret;
        }

        friend TensorMult conj(const TensorMult& tm)
        {
            return tm.conj();
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/

        template <typename T, typename=enable_if_arithmetic_t<T>>
        TensorMult operator*(T factor) const
        {
            TensorMult ret(*this);
            ret.factor *= factor;
            return ret;
        }

        TensorMult operator*(const Scalar& factor) const
        {
            TensorMult ret(*this);
            ret.factor *= factor;
            return ret;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        TensorMult operator/(T factor) const
        {
            TensorMult ret(*this);
            ret.factor /= factor;
            return ret;
        }

        TensorMult operator/(const Scalar& factor) const
        {
            TensorMult ret(*this);
            ret.factor /= factor;
            return ret;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        friend TensorMult operator*(T factor, const TensorMult& other)
        {
            return other*factor;
        }

        friend TensorMult operator*(const Scalar& factor, const TensorMult& other)
        {
            return other*factor;
        }
};

class ConstScaledTensor
{
    public:
        TensorImplementation<>& tensor;
        Scalar factor;
        bool conj_;

        ConstScaledTensor(const TensorImplementation<>& tensor, const Scalar& factor, bool conj=false)
        : tensor(const_cast<TensorImplementation<>&>(tensor)), factor(factor), conj_(conj) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/

        ConstScaledTensor operator-() const
        {
            ConstScaledTensor ret(*this);
            ret.factor = -ret.factor;
            return ret;
        }

        ConstScaledTensor conj() const
        {
            ConstScaledTensor ret(*this);
            ret.conj_ = !ret.conj_;
            return ret;
        }

        friend ConstScaledTensor conj(const ConstScaledTensor& st)
        {
            return st.conj();
        }

        /**********************************************************************
         *
         * Binary tensor operations
         *
         *********************************************************************/

        TensorMult operator*(const ConstScaledTensor& other) const
        {
            return TensorMult(tensor, other.tensor, conj_, other.conj_, factor*other.factor);
        }

        template <capability_type C>
        TensorMult operator*(const ConstTensor<C>& other) const
        {
            return TensorMult(tensor, other, conj_, false, factor);
        }

        template <capability_type C>
        friend TensorMult operator*(const ConstTensor<C>& t, const ConstScaledTensor& other)
        {
            return TensorMult(t, other.tensor, false, other.conj_, other.factor);
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ConstScaledTensor operator*(T factor)
        {
            ConstScaledTensor it(*this);
            it.factor *= factor;
            return it;
        }

        ConstScaledTensor operator*(const Scalar& factor) const
        {
            ConstScaledTensor it(*this);
            it.factor *= factor;
            return it;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        friend ConstScaledTensor operator*(T factor, const ConstScaledTensor& other)
        {
            return other*factor;
        }

        friend ConstScaledTensor operator*(const Scalar& factor, const ConstScaledTensor& other)
        {
            return other*factor;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ConstScaledTensor operator/(T factor) const
        {
            ConstScaledTensor it(*this);
            it.factor /= factor;
            return it;
        }

        ConstScaledTensor operator/(const Scalar& factor) const
        {
            ConstScaledTensor it(*this);
            it.factor /= factor;
            return it;
        }
};

class ScaledTensor : public ConstScaledTensor
{
    public:
        ScaledTensor(const ConstScaledTensor& other)
        : ConstScaledTensor(other) {}

        ScaledTensor(const ScaledTensor& other)
        : ConstScaledTensor(other) {}

        ScaledTensor(const TensorImplementation<>& tensor, const Scalar& factor, bool conj=false)
        : ConstScaledTensor(tensor, factor, conj) {}

        /**********************************************************************
         *
         * Unary negation, conjugation
         *
         *********************************************************************/

        using ConstScaledTensor::operator-;

        ScaledTensor operator-()
        {
            ScaledTensor ret(*this);
            ret.factor = -ret.factor;
            return ret;
        }

        /**********************************************************************
         *
         * Unary tensor operations
         *
         *********************************************************************/

        template <capability_type C>
        ScaledTensor& operator+=(const ConstTensor<C>& other)
        {
            tensor.sum(Scalar(1), conj_, other.impl(), factor);
            return *this;
        }

        template <capability_type C>
        ScaledTensor& operator-=(const ConstTensor<C>& other)
        {
            tensor.sum(Scalar(-1), conj_, other.impl(), factor);
            return *this;
        }

        template <capability_type C>
        ScaledTensor& operator*=(const ConstTensor<C>& other)
        {
            tensor.mult(factor, conj_, tensor, conj_, other.impl(), Scalar());
            return *this;
        }

        ScaledTensor& operator=(const ScaledTensor& other)
        {
            tensor.sum(other.factor, conj_ != other.conj_, other.tensor, Scalar());
            return *this;
        }

        ScaledTensor& operator=(const ConstScaledTensor& other)
        {
            tensor.sum(other.factor, conj_ != other.conj_, other.tensor, Scalar());
            return *this;
        }

        ScaledTensor& operator+=(const ConstScaledTensor& other)
        {
            tensor.sum(other.factor, conj_ != other.conj_, other.tensor, factor);
            return *this;
        }

        ScaledTensor& operator-=(const ConstScaledTensor& other)
        {
            tensor.sum(-other.factor, conj_ != other.conj_, other.tensor, factor);
            return *this;
        }

        ScaledTensor& operator*=(const ConstScaledTensor& other)
        {
            tensor.mult(factor*other.factor, conj_, tensor, conj_ != other.conj_, other.tensor, Scalar());
            return *this;
        }

        /**********************************************************************
         *
         * Binary tensor operations
         *
         *********************************************************************/

        ScaledTensor& operator=(const TensorMult& other)
        {
            tensor.mult(other.factor, conj_ != other.conja, other.A, conj_ != other.conjb, other.B, Scalar());
            return *this;
        }

        ScaledTensor& operator+=(const TensorMult& other)
        {
            tensor.mult(other.factor, conj_ != other.conja, other.A, conj_ != other.conjb, other.B, factor);
            return *this;
        }

        ScaledTensor& operator-=(const TensorMult& other)
        {
            tensor.mult(-other.factor, conj_ != other.conja, other.A, conj_ != other.conjb, other.B, factor);
            return *this;
        }

        /**********************************************************************
         *
         * Operations with scalars
         *
         *********************************************************************/

        using ConstScaledTensor::operator*;

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ScaledTensor operator*(T factor)
        {
            ScaledTensor it(*this);
            it.factor *= factor;
            return it;
        }

        ScaledTensor operator*(const Scalar& factor)
        {
            ScaledTensor it(*this);
            it.factor *= factor;
            return it;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        friend ScaledTensor operator*(T factor, ScaledTensor& other)
        {
            return other*factor;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        friend ScaledTensor operator*(T factor, ScaledTensor&& other)
        {
            return other*factor;
        }

        friend ScaledTensor operator*(const Scalar& factor, ScaledTensor& other)
        {
            return other*factor;
        }

        friend ScaledTensor operator*(const Scalar& factor, ScaledTensor&& other)
        {
            return other*factor;
        }

        using ConstScaledTensor::operator/;

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ScaledTensor operator/(T factor)
        {
            ScaledTensor it(*this);
            it.factor /= factor;
            return it;
        }

        ScaledTensor operator/(const Scalar& factor)
        {
            ScaledTensor it(*this);
            it.factor /= factor;
            return it;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ScaledTensor& operator*=(T val)
        {
            tensor.scale(val*factor);
            return *this;
        }

        ScaledTensor& operator*=(const Scalar& val)
        {
            tensor.scale(val*factor);
            return *this;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ScaledTensor& operator/=(T val)
        {
            tensor.scale(factor/val);
            return *this;
        }

        ScaledTensor& operator/=(const Scalar& val)
        {
            tensor.scale(factor/val);
            return *this;
        }
};

/*******************************************************************************
 *
 * Automatic Inheritance Classes
 *
 ******************************************************************************/

template <capability_type C>
struct TensorImplRoot : TensorImplementation<C&(~C)>
{
    TensorImplRoot(const INITIALIZER_TYPE(C)& ilist)
    : TensorImplementation<C&(~C)>(C, ilist) {}
};

template <capability_type C, capability_type C_>
enable_if_t<!!(C_&DONT_NEED_INITIALIZATION),TensorInitializer<C_>>
getInitializer(const INITIALIZER_TYPE(C)& ilist)
{
    return TensorInitializer<C_>();
}

template <capability_type C, capability_type C_>
enable_if_t<!(C_&DONT_NEED_INITIALIZATION),const TensorInitializer<C_>&>
getInitializer(const INITIALIZER_TYPE(C)& ilist)
{
    return ilist.template as<C_>();
}

/*
 * Tree inheritance chain to include all of the initialization data (which
 * is used later on by the implementation).
 */
TREE_BASE
(
    TensorInitBase,
    TensorInitializer,
    TensorImplRoot,
    (const INITIALIZER_TYPE(C)& ilist),
    (ilist),
    (getInitializer<C,C__>(ilist)),
    {}
)

template <capability_type C>
struct TensorConstructor : TensorInitBase<C>
{
    TensorConstructor(const INITIALIZER_TYPE(C)& ilist)
    : TensorInitBase<C>(ilist) {}

    Tensor<(C&~C)> construct(const map<capability_type,any>& m)
    {
        INITIALIZER_TYPE(C) ilist(static_cast<const TensorImplementation<C>&>(*this));

        for (auto& p : m)
        {
            auto i = ilist.find(p.first);
            assert(i != ilist.end());
            i->second = p.second;
        }

        return Tensor<C>::construct(ilist);
    }
};

/*
 * Tree inheritance chain to include all of the implemented interfaces
 * and to collect appropriate pointers for later casting.
 */
TREE_BASE
(
    TensorIfaceBase,
    TensorInterface,
    TensorConstructor,
    (const INITIALIZER_TYPE(C)& ilist),
    (ilist),
    (),
    {
        this->ptr[C__] = static_cast<void*>(static_cast<TensorInterface<C__>*>(this));
    }
)

/*
 * In-line inheritance chain to include all TensorImpls needed for C
 */
INLINE_BASE
(
    TensorImplBase,
    TensorDefinition,
    TensorIfaceBase,
    (const INITIALIZER_TYPE(C)& ilist),
    (ilist),
    {}
)

template <capability_type C_>
class TensorRoot
{
    friend class ConstScaledTensor;
    friend class ScaledTensor;
    friend class TensorMult;
    template <capability_type C, capability_type C__, class Base> friend class TensorDefinition;

    protected:
        /*
         * Use shared_ptr to ensure clean-up
         */
        shared_ptr<TensorImplementation<>> ptr;

    public:
        const capability_type C;

        TensorImplementation<>& impl()
        {
            return *ptr;
        }

        const TensorImplementation<>& impl() const
        {
            return *ptr;
        }

        template <capability_type C__>
        TensorInterface<C__>& impl()
        {
            return impl().template as<C__>();
        }

        template <capability_type C__>
        const TensorInterface<C__>& impl() const
        {
            return impl().template as<C__>();
        }

        TensorRoot() : C(C_) {}

        operator const ConstTensor<C_^CONST_>&() const
        {
            return static_cast<ConstTensor<C_^CONST_>&>(*this);
        }

        template <capability_type C__=C_, typename=enable_if_t<!(C__&CONST_)>>
        operator Tensor<C_>&()
        {
            return static_cast<Tensor<C_>&>(*this);
        }
};

INLINE_BASE
(
    TensorWrapperBase,
    TensorWrapper,
    TensorRoot,
    (),
    (),
    {}
)

/*******************************************************************************
 *
 * Generic Tensor Implementation
 *
 ******************************************************************************/

/*
 * Base class for all final implementations; provides automatic derivation from each of
 * the needed TensorImpl classes, safe destruction for type-erasure within Tensor,
 * and automatic initialization of TensorImpl base classes from an initialization list.
 */
template <capability_type C>
class TensorImplementation<C, enable_if_t<(C>0)>> : public TensorImplBase<C>
{
    public:
        /*
         * Constructor need only initiate automatic initialization chain
         */
        TensorImplementation(const INITIALIZER_TYPE(C)& ilist)
        : TensorImplBase<C>(ilist)
        {
        	assert(ilist.find(0) != ilist.end());
        }
};

/*******************************************************************************
 *
 * Tensor Wrapper Base Specialization
 *
 ******************************************************************************/

namespace detail
{
    template <capability_type C, typename=void> struct Construct;
}

template <capability_type C>
class ConstTensor : public TensorWrapperBase<C^CONST_>
{
    template <capability_type C_> friend class ConstTensor;

    public:
        /*
         * Initialize from an actual implementation. It will be destroyed when
         * the last wrapper referencing it dies.
         */
        ConstTensor(TensorImplementation<>* t)
        {
            assert(IS_SUPERSET_OF(t->C, C&~CONST_));
            this->ptr.reset(t);
            const_cast<capability_type&>(this->C) = t->C;
        }

        /*
         * Re-wrap another Tensor wrapper. Capabilities are checked at
         * run-time to enable "remembering" lost capabilties.
         */
        template <capability_type C_>
        ConstTensor(const ConstTensor<C_>& t)
        {
            assert(IS_SUPERSET_OF(t.C, C&~CONST_));
            this->ptr = t.ptr;
            const_cast<capability_type&>(this->C) = t.C;
        }

        /*
         * Copy tensor wrapper.
         */
        ConstTensor(const ConstTensor& t)
        {
            this->ptr = t.ptr;
            const_cast<capability_type&>(this->C) = t.C;
        }

        /*
         * Assume the identity of the other tensor.
         */
        ConstTensor(ConstTensor&& t)
        {
            this->ptr = move(t.ptr);
            const_cast<capability_type&>(this->C) = t.C;
        }

        /*
         * Assume the identity of the other tensor.
         */
        ConstTensor& operator=(ConstTensor&& t)
        {
            this->ptr = move(t.ptr);
            const_cast<capability_type&>(this->C) = t.C;
            return *this;
        }

        const string& getName()
        {
            return this->impl().name;
        }

        const Field& getField()
        {
            return this->impl().F;
        }

        const Ring& getRing()
        {
            return this->impl().R;
        }

        bool is(capability_type C_) const
        {
            return IS_SUPERSET_OF(this->C, C_);
        }

        template <capability_type C_>
        ConstTensor<C_> as() const
        {
            return ConstTensor<C_>(*this);
        }

        static Tensor<C&(~CONST_)> construct(const INITIALIZER_TYPE(C)& init)
        {
            return detail::Construct<C&(~CONST_)>::construct(init);
        }

        static Tensor<C&(~CONST_)> construct(const string& name, const INITIALIZER_TYPE(C)& init)
        {
            return detail::Construct<C&(~CONST_)>::construct(TensorInitializer<>(name) << init);
        }

        struct Factory
        {
            Tensor<C&(~CONST_)> operator()(const INITIALIZER_TYPE(C)& init)
            {
                return construct(init);
            }
        };

        Tensor<C&(~CONST_)> construct() const
        {
            return this->ptr->construct(TensorInitializer<>(this->ptr->name, this->ptr->F, this->ptr->R));
        }

        template <capability_type C_>
        enable_if_t<IS_POWER_OF_TWO(C_),Tensor<C&(~CONST_)>>
        construct(const TensorInitializer<C_>& init) const
        {
            return this->ptr->construct(TensorInitializer<>(this->ptr->name, this->ptr->F, this->ptr->R) << init);
        }

        template <capability_type C_>
        enable_if_t<!IS_POWER_OF_TWO(C_),Tensor<C&(~CONST_)>>
        construct(const TensorInitializerList<C_>& ilist) const
        {
            return this->ptr->construct(TensorInitializer<>(this->ptr->name, this->ptr->F, this->ptr->R) << ilist);
        }

        Tensor<C&(~CONST_)> construct(const string& name) const
        {
            return this->ptr->construct(TensorInitializer<>(name, this->ptr->F, this->ptr->R));
        }

        template <capability_type C_>
        enable_if_t<IS_POWER_OF_TWO(C_),Tensor<C&(~CONST_)>>
        construct(const string& name, const TensorInitializer<C_>& init) const
        {
            return this->ptr->construct(TensorInitializer<>(name, this->ptr->F, this->ptr->R) << init);
        }

        template <capability_type C_>
        enable_if_t<!IS_POWER_OF_TWO(C_),Tensor<C&(~CONST_)>>
        construct(const string& name, const TensorInitializerList<C_>& ilist) const
        {
            return this->ptr->construct(TensorInitializer<>(name, this->ptr->F, this->ptr->R) << ilist);
        }

        /**********************************************************************
         *
         * Intermediate operations
         *
         *********************************************************************/

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ConstScaledTensor operator*(T factor) const
        {
            return ConstScaledTensor(this->impl(), factor);
        }

        ConstScaledTensor operator*(const Scalar& factor) const
        {
            return ConstScaledTensor(this->impl(), factor);
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        friend ConstScaledTensor operator*(T factor, const ConstTensor& other)
        {
            return other*factor;
        }

        friend ConstScaledTensor operator*(const Scalar& factor, const ConstTensor& other)
        {
            return other*factor;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ConstScaledTensor operator/(T factor) const
        {
            return (*this)*(1/factor);
        }

        ConstScaledTensor operator/(const Scalar& factor) const
        {
            return (*this)*(1/factor);
        }

        ConstScaledTensor operator-() const
        {
            return (*this)*(-1);
        }

        ConstScaledTensor conj() const
        {
            return ConstScaledTensor(this->impl(), 1, true);
        }

        friend ConstScaledTensor conj(const ConstTensor& other)
        {
            return other.conj();
        }

        template <capability_type C_>
        TensorMult operator*(const ConstTensor<C_>& other) const
        {
            return TensorMult(this->impl(), other.impl(), false, false, 1);
        }
};

template <capability_type C>
class Tensor : public ConstTensor<C|CONST_>
{
    public:
        /*
         * Initialize from an actual implementation. It will be destroyed when
         * the last wrapper referencing it dies.
         */
        Tensor(TensorImplementation<>* t) : ConstTensor<C|CONST_>(t) {}

        /*
         * Re-wrap another Tensor wrapper. Capabilities are checked at
         * run-time to enable "remembering" lost capabilties.
         */
        template <capability_type C_>
        Tensor(Tensor<C_>& t) : ConstTensor<C|CONST_>(t) {}

        template <capability_type C_>
        Tensor(Tensor<C_>&& t) : ConstTensor<C|CONST_>(t) {}

        /*
         * Copy tensor wrapper.
         */
        Tensor(Tensor& t) : ConstTensor<C|CONST_>(t) {}

        /*
         * Assume the identity of the other tensor.
         */
        Tensor(Tensor&& t) : ConstTensor<C|CONST_>(move(t)) {}

        /*
         * Assume the identity of the other tensor.
         */
        Tensor& operator=(Tensor&& t)
        {
            ConstTensor<C|CONST_>::operator=(move(t));
            return *this;
        }

        template <capability_type C_>
        Tensor<C_> as()
        {
            return Tensor<C_>(*this);
        }

        using typename ConstTensor<C|CONST_>::Factory;

        using ConstTensor<C|CONST_>::construct;

        /**********************************************************************
         *
         * Intermediate operations
         *
         *********************************************************************/

        using ConstTensor<C|CONST_>::operator*;

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ScaledTensor operator*(T factor)
        {
            return ScaledTensor(this->impl(), Scalar(factor));
        }

        ScaledTensor operator*(const Scalar& factor)
        {
            return ScaledTensor(this->impl(), factor);
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        friend ScaledTensor operator*(T factor, Tensor& other)
        {
            return other*factor;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        friend ScaledTensor operator*(T factor, Tensor&& other)
        {
            return other*factor;
        }

        friend ScaledTensor operator*(const Scalar& factor, Tensor& other)
        {
            return other*factor;
        }

        friend ScaledTensor operator*(const Scalar& factor, Tensor&& other)
        {
            return other*factor;
        }

        using ConstTensor<C|CONST_>::operator/;

        template <typename T, typename=enable_if_arithmetic_t<T>>
        ScaledTensor operator/(T factor)
        {
            return (*this)*(1/Scalar(factor));
        }

        ScaledTensor operator/(const Scalar& factor)
        {
            return (*this)*(1/factor);
        }

        ScaledTensor operator-()
        {
            return (*this)*(-1);
        }

        ScaledTensor conj() const
        {
            return ScaledTensor(this->impl(), 1, true);
        }

        friend ScaledTensor conj(const Tensor& other)
        {
            return other.conj();
        }

        /**********************************************************************
         *
         * Operators with scalars
         *
         *********************************************************************/

        template <typename T, typename=enable_if_arithmetic_t<T>>
        Tensor& operator*=(T val)
        {
            this->impl().scale(Scalar(val));
            return *this;
        }

        Tensor& operator*=(const Scalar& val)
        {
            this->impl().scale(val);
            return *this;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        Tensor& operator/=(T val)
        {
            return (*this) *= (1/Scalar(val));
        }

        Tensor& operator/=(const Scalar& val)
        {
            return (*this) *= (1/val);
        }

        /**********************************************************************
         *
         * Binary operations (multiplication and division)
         *
         *********************************************************************/

        Tensor& operator=(const TensorMult& other)
        {
            this->impl().mult(other.factor, other.conja, other.A, other.conjb, other.B, Scalar());
            return *this;
        }

        Tensor& operator+=(const TensorMult& other)
        {
            this->impl().mult(other.factor, other.conja, other.A, other.conjb, other.B, Scalar(1));
            return *this;
        }

        Tensor& operator-=(const TensorMult& other)
        {
            this->impl().mult(-other.factor, other.conja, other.A, other.conjb, other.B, Scalar(1));
            return *this;
        }

        /**********************************************************************
         *
         * Unary operations (assignment, summation, scaling, and inversion)
         *
         *********************************************************************/

        Tensor& operator=(const Tensor& other)
        {
            this->impl().sum(Scalar(1), false, other.impl(), Scalar());
            return *this;
        }

        template <capability_type C_>
        Tensor& operator=(const ConstTensor<C_>& other)
        {
            this->impl().sum(Scalar(1), false, other.impl(), Scalar());
            return *this;
        }

        template <capability_type C_>
        Tensor& operator+=(const ConstTensor<C_>& other)
        {
            this->impl().sum(Scalar(1), false, other.impl(), Scalar(1));
            return *this;
        }

        template <capability_type C_>
        Tensor& operator-=(const ConstTensor<C_>& other)
        {
            this->impl().sum(Scalar(-1), false, other.impl(), Scalar(1));
            return *this;
        }

        template <capability_type C_>
        Tensor& operator*=(const ConstTensor<C_>& other)
        {
            this->impl().mult(Scalar(1), false, this->impl(), false, other.impl(), Scalar());
            return *this;
        }

        Tensor& operator=(const ConstScaledTensor& other)
        {
            this->impl().sum(other.factor, other.conj_, other.tensor, Scalar());
            return *this;
        }

        Tensor& operator+=(const ConstScaledTensor& other)
        {
            this->impl().sum(other.factor, other.conj_, other.tensor, Scalar(1));
            return *this;
        }

        Tensor& operator-=(const ConstScaledTensor& other)
        {
            this->impl().sum(-other.factor, other.conj_, other.tensor, Scalar(1));
            return *this;
        }

        Tensor& operator*=(const ConstScaledTensor& other)
        {
            this->impl().mult(other.factor, false, this->impl(), other.conj_, other.tensor, Scalar());
            return *this;
        }
};

inline Scalar scalar(const TensorMult& tm)
{
    return tm.factor*tm.B.dot(tm.conja, tm.A, tm.conjb);
}

}
}

#include "tensor/divisible.hpp"
#include "tensor/indexable.hpp"
#include "tensor/ipsymmetric.hpp"
#include "tensor/bounded.hpp"
#include "tensor/distributed.hpp"
#include "tensor/pgsymmetric.hpp"
#include "tensor/spinorbital.hpp"

#include "tensor/implementations/symblocked_tensor.hpp"
#include "tensor/implementations/spinblocked_tensor.hpp"
#include "tensor/implementations/ctf_tensor.hpp"

namespace aquarius
{
namespace tensor
{
namespace detail
{

/*******************************************************************************
 *
 * Generic Tensor Construction
 *
 ******************************************************************************/

template <capability_type C>
struct Construct<C, enable_if_t<IS_SUPERSET_OF(C,SPINORBITAL)>>
{
    static Tensor<C> construct(const INITIALIZER_TYPE(C)& ilist)
    {
        return new SpinBlockedTensor<C>(ilist);
    }
};

template <capability_type C>
struct Construct<C, enable_if_t<!IS_SUPERSET_OF(C,SPINORBITAL)&&
                                 IS_SUPERSET_OF(C,PGSYMMETRIC)&&
                                 IS_SUPERSET_OF(C,IPSYMMETRIC)>>
{
    static Tensor<C> construct(const INITIALIZER_TYPE(C)& ilist)
    {
        return new SymmetryBlockedTensor<C>(ilist);
    }
};

template <capability_type C>
struct Construct<C, enable_if_t<!IS_SUPERSET_OF(C,SPINORBITAL)&&
                                 IS_SUPERSET_OF(C,PGSYMMETRIC)&&
                                !IS_SUPERSET_OF(C,IPSYMMETRIC)>>
{
    static Tensor<C> construct(const INITIALIZER_TYPE(C)& ilist)
    {
        int ndim = ilist.template as<INDEXABLE>().ndim;

        return new SymmetryBlockedTensor<C|IPSYMMETRIC>(ilist <<
            TensorInitializer<IPSYMMETRIC_>(vector<int>(ndim, NS)));
    }
};

template <capability_type C>
struct Construct<C, enable_if_t<!IS_SUPERSET_OF(C,SPINORBITAL)&&
                                !IS_SUPERSET_OF(C,PGSYMMETRIC)&&
                                 IS_SUPERSET_OF(C,IPSYMMETRIC)&&
                                 IS_SUPERSET_OF(C,DISTRIBUTED)>>
{
    static Tensor<C> construct(const INITIALIZER_TYPE(C)& ilist)
    {
        return new CTFTensor(ilist);
    }
};

template <capability_type C>
struct Construct<C, enable_if_t<!IS_SUPERSET_OF(C,SPINORBITAL)&&
                                !IS_SUPERSET_OF(C,PGSYMMETRIC)&&
                                 IS_SUPERSET_OF(C,IPSYMMETRIC)&&
                                !IS_SUPERSET_OF(C,DISTRIBUTED)>>
{
    static Tensor<C> construct(const INITIALIZER_TYPE(C)& ilist)
    {
        Arena arena;

        return new CTFTensor(ilist << TensorInitializer<DISTRIBUTED>(arena));
    }
};

template <capability_type C>
struct Construct<C, enable_if_t<!IS_SUPERSET_OF(C,SPINORBITAL)&&
                                !IS_SUPERSET_OF(C,PGSYMMETRIC)&&
                                !IS_SUPERSET_OF(C,IPSYMMETRIC)&&
                                 IS_SUPERSET_OF(C,BOUNDED)&&
                                 IS_SUPERSET_OF(C,DISTRIBUTED)>>
{
    static Tensor<C> construct(const INITIALIZER_TYPE(C)& ilist)
    {
        int ndim = ilist.template as<INDEXABLE>().ndim;

        return new CTFTensor(ilist <<
                             TensorInitializer<IPSYMMETRIC_>(vector<int>(ndim, NS)));
    }
};

template <capability_type C>
struct Construct<C, enable_if_t<!IS_SUPERSET_OF(C,SPINORBITAL)&&
                                !IS_SUPERSET_OF(C,PGSYMMETRIC)&&
                                !IS_SUPERSET_OF(C,IPSYMMETRIC)&&
                                 IS_SUPERSET_OF(C,BOUNDED)&&
                                !IS_SUPERSET_OF(C,DISTRIBUTED)>>
{
    static Tensor<C> construct(const INITIALIZER_TYPE(C)& ilist)
    {
        Arena arena;
        int ndim = ilist.template as<INDEXABLE>().ndim;

        return new CTFTensor(ilist <<
                             TensorInitializer<DISTRIBUTED>(arena) <<
                             TensorInitializer<IPSYMMETRIC_>(vector<int>(ndim, NS)));
    }
};

}
}
}

#endif
