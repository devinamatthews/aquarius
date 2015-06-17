#ifndef _AQUARIUS_TENSOR_DISTRIBUTED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_DISTRIBUTED_TENSOR_HPP_

#include "util/global.hpp"

#include "tensor/tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <> class TensorInitializer<DISTRIBUTED> : public Destructible
{
    public:
        mutable Arena arena;

        TensorInitializer(const Arena& arena)
        : arena(arena) {}
};

TENSOR_DEFINITION(DISTRIBUTED)
{
    public:
        TensorDefinition(const INITIALIZER_TYPE(C)& ilist) : Base(ilist) {}

        const Arena& getArena() const { return this->arena; }

        void getData(key_type n, key_type* keys, void* values) const
        {
            if (this->arena.rank == 0)
            {
                this->template as<DISTRIBUTED>().getRemoteData(n, keys, values);
            }
            else
            {
                this->template as<DISTRIBUTED>().getRemoteData();
            }
            this->arena.comm().Bcast(values, n, 0, Field::MPI_TYPE(this->F));
        }

        void setData(key_type n, const key_type* keys, const void* values)
        {
            if (this->arena.rank == 0)
            {
                this->template as<DISTRIBUTED>().setRemoteData(n, keys, values);
            }
            else
            {
                this->template as<DISTRIBUTED>().setRemoteData();
            }
        }

        void addData(key_type n, const Scalar& alpha, const key_type* keys,
                     const void* values, const Scalar& beta)
        {
            if (this->arena.rank == 0)
            {
                this->template as<DISTRIBUTED>().addRemoteData(n, alpha, keys, values, beta);
            }
            else
            {
                this->template as<DISTRIBUTED>().addRemoteData(alpha, beta);
            }
        }
};

TENSOR_INTERFACE(DISTRIBUTED)
{
    public:
        virtual const Arena& getArena() const = 0;

        virtual void getLocalKeys(KeyVector& keys) const { assert(0); }

        virtual void getLocalData(KeyValueVector& kv) const { assert(0); }

        virtual void setLocalData(key_type n, const key_type* keys, const void* values) { assert(0); }

        virtual void addLocalData(key_type n, const Scalar& alpha, const key_type* keys,
                                  const void* values, const Scalar& beta) { assert(0); }

        virtual void getRemoteData(key_type n, key_type* keys, void* values) const { assert(0); }

        virtual void setRemoteData(key_type n, const key_type* keys, const void* values) { assert(0); }

        virtual void addRemoteData(key_type n, const Scalar& alpha, const key_type* keys,
                                   const void* values, const Scalar& beta) { assert(0); }

        virtual void getRemoteData() const { assert(0); }

        virtual void setRemoteData() { assert(0); }

        virtual void addRemoteData(const Scalar& alpha, const Scalar& beta) { assert(0); }
};

TENSOR_WRAPPER(DISTRIBUTED)
{
    public:
        const Arena& getArena() const
        {
            return this->template impl<DISTRIBUTED>().getArena();
        }

        KeyVector getLocalKeys() const
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            KeyVector keys;
            getLocalKeys(keys);
            return keys;
        }

        void getLocalKeys(KeyVector& keys) const
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            this->template impl<DISTRIBUTED>().getLocalKeys(keys);
        }

        KeyValueVector getLocalData() const
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            KeyValueVector kv(this->impl().F);
            getLocalData(kv);
            return kv;
        }

        void getLocalData(KeyValueVector& kv) const
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            assert(this->impl().F == kv.field());
            this->template impl<DISTRIBUTED>().getLocalData(kv);
        }

        void setLocalData(const KeyValueVector& kv)
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->impl().F == kv.field());
            this->template impl<DISTRIBUTED>().setLocalData(kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        setLocalData(const KeyVector& keys, const vector<T>& values)
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->base().impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<DISTRIBUTED>().setLocalData(keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        void addLocalData(const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->base().impl().F == kv.field());
            this->template impl<DISTRIBUTED>().addLocalData(kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        addLocalData(const Scalar& alpha, const KeyVector& keys,
                     const vector<T>& values, const Scalar& beta)
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->base().impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<DISTRIBUTED>().addLocalData(keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        void getRemoteData() const
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            this->template impl<DISTRIBUTED>().getRemoteData();
        }

        void setRemoteData()
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            this->template impl<DISTRIBUTED>().setRemoteData();
        }

        void addRemoteData(const Scalar& alpha, const Scalar& beta)
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            this->template impl<DISTRIBUTED>().addRemoteData(alpha, beta);
        }

        void getRemoteData(KeyValueVector& kv) const
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            assert(this->base().impl().F == kv.field());
            this->template impl<DISTRIBUTED>().getRemoteData(kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        getRemoteData(KeyVector& keys, vector<T>& values) const
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            assert(this->base().impl().F == Field(T()));
            values.resize(keys.size());
            this->template impl<DISTRIBUTED>().getRemoteData(keys.size(), keys.data(), static_cast<void*>(values.data()));
        }

        void setRemoteData(const KeyValueVector& kv)
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->base().impl().F == kv.field());
            this->template impl<DISTRIBUTED>().setRemoteData(kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        setRemoteData(const KeyVector& keys, const vector<T>& values)
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->base().impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<DISTRIBUTED>().setRemoteData(keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        void addRemoteData(const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->base().impl().F == kv.field());
            this->template impl<DISTRIBUTED>().addRemoteData(kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        addRemoteData(const Scalar& alpha, const KeyVector& keys,
                      const vector<T>& values, const Scalar& beta)
        {
            static_assert(C&BOUNDED_, "The operand must be BOUNDED.");
            static_assert(!(C&CONST_), "The operand must not be const.");
            assert(this->base().impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<DISTRIBUTED>().addRemoteData(keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }
};

}
}

#endif
