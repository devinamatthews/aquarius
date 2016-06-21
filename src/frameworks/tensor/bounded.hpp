#ifndef _AQUARIUS_FRAMEWORKS_TENSOR_BOUNDED_HPP_
#define _AQUARIUS_FRAMEWORKS_TENSOR_BOUNDED_HPP_

#include "frameworks/util.hpp"
#include "frameworks/tensor.hpp"

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
        : F(other.F), keys_(other.keys_)
        {
            switch (F.type())
            {
                case Field::SINGLE:   new (& fvalues_) vector<   float>( other.fvalues_); break;
                case Field::DOUBLE:   new (& dvalues_) vector<  double>( other.dvalues_); break;
                case Field::SCOMPLEX: new (&fcvalues_) vector<scomplex>(other.fcvalues_); break;
                case Field::DCOMPLEX: new (&dcvalues_) vector<dcomplex>(other.dcvalues_); break;
            }
        }

        explicit KeyValueVector(size_type size = 0)
        : F(Field::DOUBLE), keys_(size)
        {
            new (&dvalues_) vector<double>(size);
        }

        explicit KeyValueVector(Field F, size_type size = 0)
        : F(F), keys_(size)
        {
            switch (F.type())
            {
                case Field::SINGLE:   new (& fvalues_) vector<   float>(size); break;
                case Field::DOUBLE:   new (& dvalues_) vector<  double>(size); break;
                case Field::SCOMPLEX: new (&fcvalues_) vector<scomplex>(size); break;
                case Field::DCOMPLEX: new (&dcvalues_) vector<dcomplex>(size); break;
            }
        }

        ~KeyValueVector()
        {
            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.~vector(); break;
                case Field::DOUBLE:    dvalues_.~vector(); break;
                case Field::SCOMPLEX: fcvalues_.~vector(); break;
                case Field::DCOMPLEX: dcvalues_.~vector(); break;
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

        template <typename T> enable_if_field_t<T,const T&>
        value(size_type i) const
        {
            assert(F == field_type<T>::value);
            return static_cast<const vector<T>*>(&fvalues_)->operator[](i);
        }

        Scalar value(size_type i) const
        {
            switch (F.type())
            {
                case Field::SINGLE:   return Scalar( fvalues_[i]); break;
                case Field::DOUBLE:   return Scalar( dvalues_[i]); break;
                case Field::SCOMPLEX: return Scalar(fcvalues_[i]); break;
                case Field::DCOMPLEX: return Scalar(dcvalues_[i]); break;
            }

            return Scalar();
        }

        template <typename T> enable_if_field_t<T>
        value(size_type i, T s)
        {
            switch (F.type())
            {
                case Field::SINGLE:    fvalues_[i] = (   float)s; break;
                case Field::DOUBLE:    dvalues_[i] = (  double)s; break;
                case Field::SCOMPLEX: fcvalues_[i] = (scomplex)s; break;
                case Field::DCOMPLEX: dcvalues_[i] = (dcomplex)s; break;
            }
        }

        void value(size_type i, const Scalar& s)
        {
            switch (F.type())
            {
                case Field::SINGLE:    fvalues_[i] = (   float)s; break;
                case Field::DOUBLE:    dvalues_[i] = (  double)s; break;
                case Field::SCOMPLEX: fcvalues_[i] = (scomplex)s; break;
                case Field::DCOMPLEX: dcvalues_[i] = (dcomplex)s; break;
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

            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.resize(i); break;
                case Field::DOUBLE:    dvalues_.resize(i); break;
                case Field::SCOMPLEX: fcvalues_.resize(i); break;
                case Field::DCOMPLEX: dcvalues_.resize(i); break;
            }
        }

        void reserve(size_type i)
        {
            keys_.reserve(i);

            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.reserve(i); break;
                case Field::DOUBLE:    dvalues_.reserve(i); break;
                case Field::SCOMPLEX: fcvalues_.reserve(i); break;
                case Field::DCOMPLEX: dcvalues_.reserve(i); break;
            }
        }

        vector<key_type>& keys() { return keys_; }

        const vector<key_type>& keys() const { return keys_; }

        template <typename T> enable_if_same_t<T,void,const void*>
        data() const
        {
            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.data(); break;
                case Field::DOUBLE:    dvalues_.data(); break;
                case Field::SCOMPLEX: fcvalues_.data(); break;
                case Field::DCOMPLEX: dcvalues_.data(); break;
            }

            return NULL;
        }

        template <typename T> enable_if_same_t<T,void,void*>
        data()
        {
            return const_cast<void*>(const_cast<const KeyValueVector&>(*this).data<void>());
        }

        template <typename T> enable_if_field_t<T,const T*>
        data() const
        {
            return values<T>().data();
        }

        template <typename T> enable_if_field_t<T,T*>
        data()
        {
            return const_cast<T*>(const_cast<const KeyValueVector&>(*this).data<T>());
        }

        template <typename T> enable_if_field_t<T,const vector<T>&>
        values() const
        {
            assert(F == Field(T()));
            return *reinterpret_cast<const vector<T>*>(&fvalues_);
        }

        template <typename T> enable_if_field_t<T,vector<T>&>
        values()
        {
            return const_cast<vector<T>&>(const_cast<const KeyValueVector&>(*this).values<T>());
        }

        void push_back(key_type k)
        {
            keys_.push_back(k);

            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.push_back(0.0f); break;
                case Field::DOUBLE:    dvalues_.push_back(0.0 ); break;
                case Field::SCOMPLEX: fcvalues_.push_back(0.0f); break;
                case Field::DCOMPLEX: dcvalues_.push_back(0.0 ); break;
            }
        }

        template <typename T> enable_if_field_t<T>
        push_back(key_type k, T v)
        {
            keys_.push_back(k);

            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.push_back(v); break;
                case Field::DOUBLE:    dvalues_.push_back(v); break;
                case Field::SCOMPLEX: fcvalues_.push_back(v); break;
                case Field::DCOMPLEX: dcvalues_.push_back(v); break;
            }
        }

        void push_back(key_type k, const Scalar& v)
        {
            keys_.push_back(k);

            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.push_back((   float)v); break;
                case Field::DOUBLE:    dvalues_.push_back((  double)v); break;
                case Field::SCOMPLEX: fcvalues_.push_back((scomplex)v); break;
                case Field::DCOMPLEX: dcvalues_.push_back((dcomplex)v); break;
            }
        }

        void pop_back()
        {
            keys_.pop_back();

            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.pop_back(); break;
                case Field::DOUBLE:    dvalues_.pop_back(); break;
                case Field::SCOMPLEX: fcvalues_.pop_back(); break;
                case Field::DCOMPLEX: dcvalues_.pop_back(); break;
            }
        }

        KeyValueVector& operator+=(const KeyValueVector& kv)
        {
            assert(F == kv.field());
            keys_.insert(keys_.end(), kv.keys_.begin(), kv.keys_.end());

            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.insert( fvalues_.end(), kv. fvalues_.begin(), kv. fvalues_.end()); break;
                case Field::DOUBLE:    dvalues_.insert( dvalues_.end(), kv. dvalues_.begin(), kv. dvalues_.end()); break;
                case Field::SCOMPLEX: fcvalues_.insert(fcvalues_.end(), kv.fcvalues_.begin(), kv.fcvalues_.end()); break;
                case Field::DCOMPLEX: dcvalues_.insert(dcvalues_.end(), kv.dcvalues_.begin(), kv.dcvalues_.end()); break;
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
            assert(F == other.F);
            keys_.swap(other.keys_);
            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.swap( other.fvalues_); break;
                case Field::DOUBLE:    dvalues_.swap( other.dvalues_); break;
                case Field::SCOMPLEX: fcvalues_.swap(other.fcvalues_); break;
                case Field::DCOMPLEX: dcvalues_.swap(other.dcvalues_); break;
            }
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

            switch (F.type())
            {
                case Field::SINGLE:    fvalues_.clear(); break;
                case Field::DOUBLE:    dvalues_.clear(); break;
                case Field::SCOMPLEX: fcvalues_.clear(); break;
                case Field::DCOMPLEX: dcvalues_.clear(); break;
            }
        }

    protected:
        Field F;
        vector<key_type> keys_;
        union
        {
            vector<   float>  fvalues_;
            vector<  double>  dvalues_;
            vector<scomplex> fcvalues_;
            vector<dcomplex> dcvalues_;
        };
};

typedef vector<key_type> KeyVector;

template <> class TensorInitializer<BOUNDED_>
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

        template <typename T> enable_if_field_t<T>
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

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void setData(const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<BOUNDED_>().setData(kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        enable_if_field_t<T>
        setData(const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<BOUNDED_>().setData(keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void addData(const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<BOUNDED_>().addData(kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        enable_if_field_t<T>
        addData(const Scalar& alpha, const KeyVector& keys,
                const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<BOUNDED_>().addData(keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(const Scalar& alpha, bool conja, const vector<int>& start_A, ConstTensor<BOUNDED> A,
                   const Scalar&  beta,             const vector<int>& start_B, const vector<int>& length)
        {
            this->template impl<BOUNDED_>().slice(alpha, conja, start_A, A.impl(),
                                                   beta,        start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(const Scalar& alpha, const vector<int>& start_A, ConstTensor<BOUNDED> A,
                   const Scalar&  beta, const vector<int>& start_B, const vector<int>& length)
        {
            slice(alpha, false, A, start_A, beta, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(bool conja, const vector<int>& start_A, ConstTensor<BOUNDED> A,
                               const vector<int>& start_B, const vector<int>& length)
        {
            slice(Scalar(1), conja, A, start_A, Scalar(0), start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(const vector<int>& start_A, ConstTensor<BOUNDED> A,
                   const vector<int>& start_B, const vector<int>& length)
        {
            slice(Scalar(1), false, A, start_A, Scalar(0), start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(const Scalar& alpha, bool conja, ConstTensor<BOUNDED> A,
                   const Scalar&  beta, const vector<int>& start_B)
        {
            slice(alpha, conja, A, vector<int>(A.getDimension()), beta, start_B, A.getLengths());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(const Scalar& alpha, ConstTensor<BOUNDED> A,
                   const Scalar&  beta, const vector<int>& start_B)
        {
            slice(alpha, false, A, vector<int>(A.getDimension()), beta, start_B, A.getLengths());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(bool conja, ConstTensor<BOUNDED> A, const vector<int>& start_B)
        {
            slice(Scalar(1), conja, A, vector<int>(A.getDimension()), Scalar(0), start_B, A.getLengths());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(ConstTensor<BOUNDED> A, const vector<int>& start_B)
        {
            slice(Scalar(1), false, A, vector<int>(A.getDimension()), Scalar(0), start_B, A.getLengths());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(const Scalar& alpha, bool conja, const vector<int>& start_A, ConstTensor<BOUNDED> A,
                   const Scalar&  beta)
        {
            slice(alpha, conja, A, start_A, beta, vector<int>(this->getDimension()), this->getLengths());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(const Scalar& alpha, const vector<int>& start_A, ConstTensor<BOUNDED> A,
                   const Scalar&  beta)
        {
            slice(alpha, false, A, start_A, beta, vector<int>(this->getDimension()), this->getLengths());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(bool conja, const vector<int>& start_A, ConstTensor<BOUNDED> A)
        {
            slice(Scalar(1), conja, A, start_A, Scalar(0), vector<int>(this->getDimension()), this->getLengths());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void slice(const vector<int>& start_A, ConstTensor<BOUNDED> A)
        {
            slice(Scalar(1), false, A, start_A, Scalar(0), vector<int>(this->getDimension()), this->getLengths());
        }
};

}
}

#endif
