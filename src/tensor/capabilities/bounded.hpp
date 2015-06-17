#ifndef _AQUARIUS_TENSOR_BOUNDED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_BOUNDED_TENSOR_HPP_

#include "util/global.hpp"

#include "tensor/tensor.hpp"

namespace aquarius
{
namespace tensor
{

typedef int64_t key_type; //must be the same as CTF::Pair<T>::k

class KeyValueVector
{
    public:
        typedef vector<key_type>::size_type size_type;
        typedef vector<key_type>::difference_type difference_type;

        KeyValueVector(const KeyValueVector& other)
        : F(other.F), keys_(other.keys_), values_(NULL)
        {
            switch (F.type)
            {
                case Field::SINGLE:
                    values_ = static_cast<void*>(new vector<      float>(
                        *static_cast<const vector<      float>*>(other.values_)));
                    break;
                case Field::DOUBLE:
                    values_ = static_cast<void*>(new vector<     double>(
                        *static_cast<const vector<     double>*>(other.values_)));
                    break;
                case Field::LDOUBLE:
                    values_ = static_cast<void*>(new vector<long double>(
                        *static_cast<const vector<long double>*>(other.values_)));
                    break;
                case Field::SCOMPLEX:
                    values_ = static_cast<void*>(new vector<   scomplex>(
                        *static_cast<const vector<   scomplex>*>(other.values_)));
                    break;
                case Field::DCOMPLEX:
                    values_ = static_cast<void*>(new vector<   dcomplex>(
                        *static_cast<const vector<   dcomplex>*>(other.values_)));
                    break;
                case Field::LDCOMPLEX:
                    values_ = static_cast<void*>(new vector<  ldcomplex>(
                        *static_cast<const vector<  ldcomplex>*>(other.values_)));
                    break;
            }
        }

        explicit KeyValueVector(Field F, size_type size = 0)
        : F(F), keys_(size), values_(NULL)
        {
            switch (F.type)
            {
                case Field::SINGLE:
                    values_ = static_cast<void*>(new vector<      float>(size));
                    break;
                case Field::DOUBLE:
                    values_ = static_cast<void*>(new vector<     double>(size));
                    break;
                case Field::LDOUBLE:
                    values_ = static_cast<void*>(new vector<long double>(size));
                    break;
                case Field::SCOMPLEX:
                    values_ = static_cast<void*>(new vector<   scomplex>(size));
                    break;
                case Field::DCOMPLEX:
                    values_ = static_cast<void*>(new vector<   dcomplex>(size));
                    break;
                case Field::LDCOMPLEX:
                    values_ = static_cast<void*>(new vector<  ldcomplex>(size));
                    break;
            }
        }

        ~KeyValueVector()
        {
            switch (F.type)
            {
                case Field::SINGLE:
                    delete static_cast<vector<      float>*>(values_);
                    break;
                case Field::DOUBLE:
                    delete static_cast<vector<     double>*>(values_);
                    break;
                case Field::LDOUBLE:
                    delete static_cast<vector<long double>*>(values_);
                    break;
                case Field::SCOMPLEX:
                    delete static_cast<vector<   scomplex>*>(values_);
                    break;
                case Field::DCOMPLEX:
                    delete static_cast<vector<   dcomplex>*>(values_);
                    break;
                case Field::LDCOMPLEX:
                    delete static_cast<vector<  ldcomplex>*>(values_);
                    break;
            }
        }

        KeyValueVector& operator=(KeyValueVector other)
        {
            swap(other);
            return *this;
        }

        template <typename T> typename enable_if<is_field<T>::value,T&>::type
        value(size_type i)
        {
            return const_cast<T&>(const_cast<const KeyValueVector&>(*this).value<T>(i));
        }

        template <typename T> typename enable_if<is_field<T>::value,const T&>::type
        value(size_type i) const
        {
            assert(F == Field(T()));
            return static_cast<const vector<T>*>(values_)->operator[](i);
        }

        Scalar value(size_type i) const
        {
            switch (F.type)
            {
                case Field::SINGLE:
                    return Scalar(static_cast<vector<      float>*>(values_)->operator[](i));
                    break;
                case Field::DOUBLE:
                    return Scalar(static_cast<vector<     double>*>(values_)->operator[](i));
                    break;
                case Field::LDOUBLE:
                    return Scalar(static_cast<vector<long double>*>(values_)->operator[](i));
                    break;
                case Field::SCOMPLEX:
                    return Scalar(static_cast<vector<   scomplex>*>(values_)->operator[](i));
                    break;
                case Field::DCOMPLEX:
                    return Scalar(static_cast<vector<   dcomplex>*>(values_)->operator[](i));
                    break;
                case Field::LDCOMPLEX:
                    return Scalar(static_cast<vector<  ldcomplex>*>(values_)->operator[](i));
                    break;
            }

            return Scalar(0.0);
        }

        void value(size_type i, const Scalar& s)
        {
            assert(s.field() == F);

            switch (F.type)
            {
                case Field::SINGLE:
                    static_cast<vector<      float>*>(values_)->operator[](i) = s.to<      float>();
                    break;
                case Field::DOUBLE:
                    static_cast<vector<     double>*>(values_)->operator[](i) = s.to<     double>();
                    break;
                case Field::LDOUBLE:
                    static_cast<vector<long double>*>(values_)->operator[](i) = s.to<long double>();
                    break;
                case Field::SCOMPLEX:
                    static_cast<vector<   scomplex>*>(values_)->operator[](i) = s.to<   scomplex>();
                    break;
                case Field::DCOMPLEX:
                    static_cast<vector<   dcomplex>*>(values_)->operator[](i) = s.to<   dcomplex>();
                    break;
                case Field::LDCOMPLEX:
                    static_cast<vector<  ldcomplex>*>(values_)->operator[](i) = s.to<  ldcomplex>();
                    break;
            }
        }

        const Field& field() const { return F; }

        key_type& key(size_type i) { return keys_[i]; }

        const key_type& key(size_type i) const { return keys_[i]; }

        size_type size() const { return keys_.size(); }

        bool empty() const { return keys_.empty(); }

        void resize(size_type i)
        {
            keys_.resize(i);

            switch (F.type)
            {
                case Field::SINGLE:
                    static_cast<vector<      float>*>(values_)->resize(i);
                    break;
                case Field::DOUBLE:
                    static_cast<vector<     double>*>(values_)->resize(i);
                    break;
                case Field::LDOUBLE:
                    static_cast<vector<long double>*>(values_)->resize(i);
                    break;
                case Field::SCOMPLEX:
                    static_cast<vector<   scomplex>*>(values_)->resize(i);
                    break;
                case Field::DCOMPLEX:
                    static_cast<vector<   dcomplex>*>(values_)->resize(i);
                    break;
                case Field::LDCOMPLEX:
                    static_cast<vector<  ldcomplex>*>(values_)->resize(i);
                    break;
            }
        }

        void reserve(size_type i)
        {
            keys_.reserve(i);

            switch (F.type)
            {
                case Field::SINGLE:
                    static_cast<vector<      float>*>(values_)->reserve(i);
                    break;
                case Field::DOUBLE:
                    static_cast<vector<     double>*>(values_)->reserve(i);
                    break;
                case Field::LDOUBLE:
                    static_cast<vector<long double>*>(values_)->reserve(i);
                    break;
                case Field::SCOMPLEX:
                    static_cast<vector<   scomplex>*>(values_)->reserve(i);
                    break;
                case Field::DCOMPLEX:
                    static_cast<vector<   dcomplex>*>(values_)->reserve(i);
                    break;
                case Field::LDCOMPLEX:
                    static_cast<vector<  ldcomplex>*>(values_)->reserve(i);
                    break;
            }
        }

        vector<key_type>& keys() { return keys_; }

        const vector<key_type>& keys() const { return keys_; }

        template <typename T> typename enable_if<is_void<T>::value,const void*>::type
        data() const
        {
            switch (F.type)
            {
                case Field::SINGLE:
                    return static_cast<const void*>(static_cast<const vector<      float>*>(values_)->data());
                    break;
                case Field::DOUBLE:
                    return static_cast<const void*>(static_cast<const vector<     double>*>(values_)->data());
                    break;
                case Field::LDOUBLE:
                    return static_cast<const void*>(static_cast<const vector<long double>*>(values_)->data());
                    break;
                case Field::SCOMPLEX:
                    return static_cast<const void*>(static_cast<const vector<   scomplex>*>(values_)->data());
                    break;
                case Field::DCOMPLEX:
                    return static_cast<const void*>(static_cast<const vector<   dcomplex>*>(values_)->data());
                    break;
                case Field::LDCOMPLEX:
                    return static_cast<const void*>(static_cast<const vector<  ldcomplex>*>(values_)->data());
                    break;
            }
            return NULL;
        }

        template <typename T> typename enable_if<is_void<T>::value,void*>::type
        data()
        {
            return const_cast<void*>(const_cast<const KeyValueVector&>(*this).data<void>());
        }

        template <typename T> typename enable_if<is_field<T>::value,const T*>::type
        data() const
        {
            assert(F == Field(T()));
            return static_cast<const vector<T>*>(values_)->data();
        }

        template <typename T> typename enable_if<is_field<T>::value,T*>::type
        data()
        {
            return const_cast<T*>(const_cast<const KeyValueVector&>(*this).data<T>());
        }

        template <typename T> typename enable_if<is_field<T>::value,const vector<T>&>::type
        values() const
        {
            assert(F == Field(T()));
            return *static_cast<const vector<T>*>(values_);
        }

        template <typename T> typename enable_if<is_field<T>::value,vector<T>&>::type
        values()
        {
            return const_cast<vector<T>&>(const_cast<const KeyValueVector&>(*this).values<T>());
        }

        void push_back(key_type k)
        {
            keys_.push_back(k);

            switch (F.type)
            {
                case Field::SINGLE:
                    static_cast<vector<      float>*>(values_)->push_back(0.0f);
                    break;
                case Field::DOUBLE:
                    static_cast<vector<     double>*>(values_)->push_back(0.0);
                    break;
                case Field::LDOUBLE:
                    static_cast<vector<long double>*>(values_)->push_back(0.0l);
                    break;
                case Field::SCOMPLEX:
                    static_cast<vector<   scomplex>*>(values_)->push_back(0.0f);
                    break;
                case Field::DCOMPLEX:
                    static_cast<vector<   dcomplex>*>(values_)->push_back(0.0);
                    break;
                case Field::LDCOMPLEX:
                    static_cast<vector<  ldcomplex>*>(values_)->push_back(0.0l);
                    break;
            }
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        push_back(key_type k, T v)
        {
            assert(F == Field(v));
            keys_.push_back(k);
            static_cast<vector<T>*>(values_)->push_back(v);
        }

        void push_back(key_type k, const Scalar& v)
        {
            assert(F == v.field());
            keys_.push_back(k);

            switch (F.type)
            {
                case Field::SINGLE:
                    static_cast<vector<      float>*>(values_)->push_back(v.to<      float>());
                    break;
                case Field::DOUBLE:
                    static_cast<vector<     double>*>(values_)->push_back(v.to<     double>());
                    break;
                case Field::LDOUBLE:
                    static_cast<vector<long double>*>(values_)->push_back(v.to<long double>());
                    break;
                case Field::SCOMPLEX:
                    static_cast<vector<   scomplex>*>(values_)->push_back(v.to<   scomplex>());
                    break;
                case Field::DCOMPLEX:
                    static_cast<vector<   dcomplex>*>(values_)->push_back(v.to<   dcomplex>());
                    break;
                case Field::LDCOMPLEX:
                    static_cast<vector<  ldcomplex>*>(values_)->push_back(v.to<  ldcomplex>());
                    break;
            }
        }

        void pop_back()
        {
            keys_.pop_back();

            switch (F.type)
            {
                case Field::SINGLE:
                    static_cast<vector<      float>*>(values_)->pop_back();
                    break;
                case Field::DOUBLE:
                    static_cast<vector<     double>*>(values_)->pop_back();
                    break;
                case Field::LDOUBLE:
                    static_cast<vector<long double>*>(values_)->pop_back();
                    break;
                case Field::SCOMPLEX:
                    static_cast<vector<   scomplex>*>(values_)->pop_back();
                    break;
                case Field::DCOMPLEX:
                    static_cast<vector<   dcomplex>*>(values_)->pop_back();
                    break;
                case Field::LDCOMPLEX:
                    static_cast<vector<  ldcomplex>*>(values_)->pop_back();
                    break;
            }
        }

        KeyValueVector& operator+=(const KeyValueVector& kv)
        {
            assert(F == kv.field());
            keys_.insert(keys_.end(), kv.keys_.begin(), kv.keys_.end());

            switch (F.type)
            {
                case Field::SINGLE:
                {
                    auto& v1 = *static_cast<      vector<float>*>(   values_);
                    auto& v2 = *static_cast<const vector<float>*>(kv.values_);
                    v1.insert(v1.end(), v2.begin(), v2.end());
                }
                break;
                case Field::DOUBLE:
                {
                    auto& v1 = *static_cast<      vector<double>*>(   values_);
                    auto& v2 = *static_cast<const vector<double>*>(kv.values_);
                    v1.insert(v1.end(), v2.begin(), v2.end());
                }
                break;
                case Field::LDOUBLE:
                {
                    auto& v1 = *static_cast<      vector<long double>*>(   values_);
                    auto& v2 = *static_cast<const vector<long double>*>(kv.values_);
                    v1.insert(v1.end(), v2.begin(), v2.end());
                }
                break;
                case Field::SCOMPLEX:
                {
                    auto& v1 = *static_cast<      vector<scomplex>*>(   values_);
                    auto& v2 = *static_cast<const vector<scomplex>*>(kv.values_);
                    v1.insert(v1.end(), v2.begin(), v2.end());
                }
                break;
                case Field::DCOMPLEX:
                {
                    auto& v1 = *static_cast<      vector<dcomplex>*>(   values_);
                    auto& v2 = *static_cast<const vector<dcomplex>*>(kv.values_);
                    v1.insert(v1.end(), v2.begin(), v2.end());
                }
                break;
                case Field::LDCOMPLEX:
                {
                    auto& v1 = *static_cast<      vector<ldcomplex>*>(   values_);
                    auto& v2 = *static_cast<const vector<ldcomplex>*>(kv.values_);
                    v1.insert(v1.end(), v2.begin(), v2.end());
                }
                break;
            }

            return *this;
        }

        friend KeyValueVector operator+(KeyValueVector kv1, const KeyValueVector& kv2)
        {
            kv1 += kv2;
            return kv1;
        }

        void swap(KeyValueVector& other)
        {
            using aquarius::swap;
            assert(F == other.F);
            swap(keys_, other.keys_);
            swap(values_, other.values_);
        }

        friend void swap(KeyValueVector& kv1, KeyValueVector& kv2)
        {
            kv1.swap(kv2);
        }

        void assign(const KeyValueVector& other)
        {
            *this = other;
        }

        void clear()
        {
            keys_.clear();

            switch (F.type)
            {
                case Field::SINGLE:
                    static_cast<vector<      float>*>(values_)->clear();
                    break;
                case Field::DOUBLE:
                    static_cast<vector<     double>*>(values_)->clear();
                    break;
                case Field::LDOUBLE:
                    static_cast<vector<long double>*>(values_)->clear();
                    break;
                case Field::SCOMPLEX:
                    static_cast<vector<   scomplex>*>(values_)->clear();
                    break;
                case Field::DCOMPLEX:
                    static_cast<vector<   dcomplex>*>(values_)->clear();
                    break;
                case Field::LDCOMPLEX:
                    static_cast<vector<  ldcomplex>*>(values_)->clear();
                    break;
            }
        }

    protected:
        Field F;
        vector<key_type> keys_;
        void *values_;
};

typedef vector<key_type> KeyVector;

template <> class TensorInitializer<BOUNDED_> : public Destructible
{
    public:
        vector<int> len;

        TensorInitializer(const vector<int>& len)
        : len(len)
        {
            for (int i = 0;i < len.size();i++) assert(len[i] >= 0);
        }
};

template <> class TensorInitializer<BOUNDED> : public TensorInitializerList<BOUNDED>
{
    public:
        TensorInitializer(const vector<int>& len)
        {
            addInitializer(TensorInitializer<INDEXABLE>(len.size()));
            addInitializer(TensorInitializer<BOUNDED_>(len));
        }
};

template <> class TensorInitializer<BOUNDED|IPSYMMETRIC> : public TensorInitializerList<BOUNDED|IPSYMMETRIC>
{
    public:
        TensorInitializer(const vector<int>& len, const vector<int>& sym)
        {
            assert(len.size() == sym.size());
            addInitializer(TensorInitializer<INDEXABLE>(len.size()));
            addInitializer(TensorInitializer<BOUNDED_>(len));
            addInitializer(TensorInitializer<IPSYMMETRIC_>(sym));
        }
};

TENSOR_DEFINITION(BOUNDED_)
{
    public:
        TensorDefinition(const INITIALIZER_TYPE(C)& ilist) : Base(ilist) {}

        const vector<int>& getLengths() const
        {
            return this->len;
        }
};

TENSOR_INTERFACE(BOUNDED_)
{
    public:
        virtual const vector<int>& getLengths() const = 0;

        virtual void getAllKeys(KeyVector& keys) const = 0;

        virtual void getAllData(KeyValueVector& kv) const = 0;

        virtual void getData(key_type n, key_type* keys, void* values) const = 0;

        virtual void setData(key_type n, const key_type* keys, const void* values) = 0;

        virtual void addData(key_type n, const Scalar& alpha, const key_type* keys,
                             const void* values, const Scalar& beta) = 0;

        virtual Scalar norm(int p) const = 0;

        virtual void slice(const Scalar& alpha, bool conja, const vector<int>& start_A, const TensorImplementation<>& A,
                           const Scalar&  beta,             const vector<int>& start_B, const vector<int>& length) = 0;
};

TENSOR_WRAPPER(BOUNDED_)
{
    public:
        const vector<int>& getLengths() const
        {
            return this->template impl<BOUNDED_>().getLengths();
        }

        key_type getKey(const vector<int>& idx) const
        {
            const vector<int>& len = getLengths();

            int n = len.size();
            key_type key = 0;

            if (areKeyStridesPowersOfTwo())
            {
                const vector<int>& stride_bits = getKeyStrideBits();

                for (int i = 0;i < n;i++)
                {
                    assert(idx[i] >= 0 && idx[i] < len[i]);
                    key |= key_type(len[i])<<stride_bits[i];
                }
            }
            else
            {
                const vector<key_type>& strides = getKeyStrides();

                for (int i = 0;i < n;i++)
                {
                    assert(idx[i] >= 0 && idx[i] < len[i]);
                    key += key_type(len[i])*strides[i];
                }
            }

            return key;
        }

        bool areKeyStridesPowersOfTwo() const
        {
            return this->template impl<BOUNDED_>().areKeyStridesPowersOfTwo();
        }

        const vector<int>& getKeyStrideBits() const
        {
            return this->template impl<BOUNDED_>().getKeyStrideBits();
        }

        const vector<key_type>& getKeyStrideMasks() const
        {
            return this->template impl<BOUNDED_>().getKeyStrideMasks();
        }

        vector<int> getIndices(key_type key) const
        {
            const vector<int>& len = getLengths();

            int n = len.size();
            vector<int> idx(n);

            if (areKeyStridesPowersOfTwo())
            {
                const vector<int>& stride_bits = getKeyStrideBits();
                const vector<key_type>& stride_mask = getKeyStrideMasks();

                for (int i = 0;i < n;i++)
                {
                    idx[i] = (key>>stride_bits[i])&stride_mask[i];
                }
            }
            else
            {
                const vector<int>& strides = getKeyStrides();

                for (int i = 0;i < n;i++)
                {
                    idx[i] = (key/strides[i])%len[i];
                }
            }

            return idx;
        }

        const vector<key_type>& getKeyStrides() const
        {
            return this->template impl<BOUNDED_>().getKeyStrides();
        }

        KeyVector getAllKeys() const
        {
            KeyVector keys;
            this->template impl<BOUNDED_>().getAllKeys(keys);
            return keys;
        }

        void getAllKeys(KeyVector& keys) const
        {
            this->template impl<BOUNDED_>().getAllKeys(keys);
        }

        KeyValueVector getAllData() const
        {
            KeyValueVector kv(this->impl().F);
            this->template impl<BOUNDED_>().getAllData(kv);
            return kv;
        }

        void getAllData(KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<BOUNDED_>().getAllData(kv);
        }

        void getData(KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<BOUNDED_>().getData(kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        getData(KeyVector& keys, vector<T>& values) const
        {
            assert(this->impl().F == Field(T()));
            values.resize(keys.size());
            this->template impl<BOUNDED_>().getData(keys.size(), keys.data(), static_cast<void*>(values.data()));
        }

        Scalar norm(int p) const
        {
            return this->template impl<BOUNDED_>().norm(p);
        }

        void setData(const KeyValueVector& kv)
        {
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->impl().F == kv.field());
            this->template impl<BOUNDED_>().setData(kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        setData(const KeyVector& keys, const vector<T>& values)
        {
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<BOUNDED_>().setData(keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        void addData(const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->impl().F == kv.field());
            this->template impl<BOUNDED_>().addData(kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        addData(const Scalar& alpha, const KeyVector& keys,
                const vector<T>& values, const Scalar& beta)
        {
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<BOUNDED_>().addData(keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_>
        void slice(const Scalar& alpha, bool conja, const vector<int>& start_A, const ConstTensor<C_>& A,
                   const Scalar&  beta,             const vector<int>& start_B, const vector<int>& length)
        {
            static_assert(!(C&CONST_), "The operand must not be const.");
            static_assert(C_&BOUNDED_, "The operands must be BOUNDED.");
            this->template impl<BOUNDED_>().slice(alpha, conja, start_A, A.impl(),
                                                   beta,        start_B, length);
        }

        template <capability_type C_>
        void slice(const Scalar& alpha, const vector<int>& start_A, const ConstTensor<C_>& A,
                   const Scalar&  beta, const vector<int>& start_B, const vector<int>& length)
        {
            slice(alpha, false, A, start_A, beta, start_B, length);
        }

        template <capability_type C_>
        void slice(bool conja, const vector<int>& start_A, const ConstTensor<C_>& A,
                               const vector<int>& start_B, const vector<int>& length)
        {
            slice(1, conja, A, start_A, 0, start_B, length);
        }

        template <capability_type C_>
        void slice(const vector<int>& start_A, const ConstTensor<C_>& A,
                   const vector<int>& start_B, const vector<int>& length)
        {
            slice(1, false, A, start_A, 0, start_B, length);
        }

        template <capability_type C_>
        void slice(const Scalar& alpha, bool conja, const ConstTensor<C_>& A,
                   const Scalar&  beta, const vector<int>& start_B)
        {
            slice(alpha, conja, A, vector<int>(A.getDimension()), beta, start_B, A.getLengths());
        }

        template <capability_type C_>
        void slice(const Scalar& alpha, const ConstTensor<C_>& A,
                   const Scalar&  beta, const vector<int>& start_B)
        {
            slice(alpha, false, A, vector<int>(A.getDimension()), beta, start_B, A.getLengths());
        }

        template <capability_type C_>
        void slice(bool conja, const ConstTensor<C_>& A, const vector<int>& start_B)
        {
            slice(1, conja, A, vector<int>(A.getDimension()), 0, start_B, A.getLengths());
        }

        template <capability_type C_>
        void slice(const ConstTensor<C_>& A, const vector<int>& start_B)
        {
            slice(1, false, A, vector<int>(A.getDimension()), 0, start_B, A.getLengths());
        }

        template <capability_type C_>
        void slice(const Scalar& alpha, bool conja, const vector<int>& start_A, const ConstTensor<C_>& A,
                   const Scalar&  beta)
        {
            slice(alpha, conja, A, start_A, beta, vector<int>(this->getDimension()), this->getLengths());
        }

        template <capability_type C_>
        void slice(const Scalar& alpha, const vector<int>& start_A, const ConstTensor<C_>& A,
                   const Scalar&  beta)
        {
            slice(alpha, false, A, start_A, beta, vector<int>(this->getDimension()), this->getLengths());
        }

        template <capability_type C_>
        void slice(bool conja, const vector<int>& start_A, const ConstTensor<C_>& A)
        {
            slice(1, conja, A, start_A, 0, vector<int>(this->getDimension()), this->getLengths());
        }

        template <capability_type C_>
        void slice(const vector<int>& start_A, const ConstTensor<C_>& A)
        {
            slice(1, false, A, start_A, 0, vector<int>(this->getDimension()), this->getLengths());
        }
};

}
}

#endif
