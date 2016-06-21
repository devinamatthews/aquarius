#ifndef _AQUARIUS_FRAMEWORKS_TENSOR_DISTRIBUTED_HPP_
#define _AQUARIUS_FRAMEWORKS_TENSOR_DISTRIBUTED_HPP_

#include "frameworks/util.hpp"
#include "frameworks/tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <> class TensorInitializer<DISTRIBUTED>
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

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        KeyVector getLocalKeys() const
        {
            KeyVector keys;
            getLocalKeys(keys);
            return keys;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        void getLocalKeys(KeyVector& keys) const
        {
            this->template impl<DISTRIBUTED>().getLocalKeys(keys);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        KeyValueVector getLocalData() const
        {
            KeyValueVector kv(this->impl().F);
            getLocalData(kv);
            return kv;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        void getLocalData(KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<DISTRIBUTED>().getLocalData(kv);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void setLocalData(const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<DISTRIBUTED>().setLocalData(kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        enable_if_field_t<T>
        setLocalData(const KeyVector& keys, const vector<T>& values)
        {
            assert(this->base().impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<DISTRIBUTED>().setLocalData(keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void addLocalData(const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->base().impl().F == kv.field());
            this->template impl<DISTRIBUTED>().addLocalData(kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        enable_if_field_t<T>
        addLocalData(const Scalar& alpha, const KeyVector& keys,
                     const vector<T>& values, const Scalar& beta)
        {
            assert(this->base().impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<DISTRIBUTED>().addLocalData(keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        void getRemoteData() const
        {
            this->template impl<DISTRIBUTED>().getRemoteData();
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void setRemoteData()
        {
            this->template impl<DISTRIBUTED>().setRemoteData();
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void addRemoteData(const Scalar& alpha, const Scalar& beta)
        {
            this->template impl<DISTRIBUTED>().addRemoteData(alpha, beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        void getRemoteData(KeyValueVector& kv) const
        {
            assert(this->base().impl().F == kv.field());
            this->template impl<DISTRIBUTED>().getRemoteData(kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        enable_if_field_t<T>
        getRemoteData(KeyVector& keys, vector<T>& values) const
        {
            assert(this->base().impl().F == Field(T()));
            values.resize(keys.size());
            this->template impl<DISTRIBUTED>().getRemoteData(keys.size(), keys.data(), static_cast<void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void setRemoteData(const KeyValueVector& kv)
        {
            assert(this->base().impl().F == kv.field());
            this->template impl<DISTRIBUTED>().setRemoteData(kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        enable_if_field_t<T>
        setRemoteData(const KeyVector& keys, const vector<T>& values)
        {
            assert(this->base().impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<DISTRIBUTED>().setRemoteData(keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void addRemoteData(const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->base().impl().F == kv.field());
            this->template impl<DISTRIBUTED>().addRemoteData(kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }
        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        enable_if_field_t<T>
        addRemoteData(const Scalar& alpha, const KeyVector& keys,
                      const vector<T>& values, const Scalar& beta)
        {
            assert(this->base().impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<DISTRIBUTED>().addRemoteData(keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }
};

}
}

#endif
