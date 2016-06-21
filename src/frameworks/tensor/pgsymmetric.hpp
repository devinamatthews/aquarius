#ifndef _AQUARIUS_FRAMEWORKS_TENSOR_PGSYMMETRIC_HPP_
#define _AQUARIUS_FRAMEWORKS_TENSOR_PGSYMMETRIC_HPP_

#include "frameworks/util.hpp"
#include "frameworks/symmetry.hpp"
#include "frameworks/tensor.hpp"

namespace aquarius
{
namespace tensor
{

namespace detail
{
    inline vector<int> getNumTotal(const vector<vector<int>>& n_per_irrep)
    {
        int m = n_per_irrep.size();
        vector<int> n_total(m);

        for (int i = 0;i < m;i++)
        {
            int n = n_per_irrep[i].size();
            for (int j = 0;j < n;j++)
            {
                n_total[i] += n_per_irrep[i][j];
            }
        }

        return n_total;
    }
}

template <> class TensorInitializer<PGSYMMETRIC_>
{
    public:
        const symmetry::PointGroup& group;
        symmetry::Representation rep;
        vector<vector<int>> len_per_irrep;

        TensorInitializer(const symmetry::PointGroup& group)
        : group(group), rep(group.totallySymmetricIrrep()) {}

        TensorInitializer(const symmetry::PointGroup& group,
                          const symmetry::Representation& rep)
        : group(group), rep(rep) {}

        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& len_per_irrep)
        : group(group), rep(group.totallySymmetricIrrep()), len_per_irrep(len_per_irrep) {}

        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& len_per_irrep,
                          const symmetry::Representation& rep)
        : group(group), rep(rep), len_per_irrep(len_per_irrep) {}
};

template <> class TensorInitializer<PGSYMMETRIC|BOUNDED>
: public TensorInitializerList<PGSYMMETRIC|BOUNDED>
{
    public:
        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& len_per_irrep)
        {
            assert(len_per_irrep.size() == 0 ||
                   len_per_irrep[0].size() == group.getNumIrreps());
            addInitializer(TensorInitializer<PGSYMMETRIC>(group, len_per_irrep));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getNumTotal(len_per_irrep)));
            addInitializer(TensorInitializer<INDEXABLE>(len_per_irrep.size()));
        }

        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& len_per_irrep,
                          const symmetry::Representation& rep)
        {
            assert(len_per_irrep.size() == 0 ||
                   len_per_irrep[0].size() == group.getNumIrreps());
            addInitializer(TensorInitializer<PGSYMMETRIC>(group, len_per_irrep, rep));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getNumTotal(len_per_irrep)));
            addInitializer(TensorInitializer<INDEXABLE>(len_per_irrep.size()));
        }
};

template <> class TensorInitializer<PGSYMMETRIC|BOUNDED|IPSYMMETRIC>
: public TensorInitializerList<PGSYMMETRIC|BOUNDED|IPSYMMETRIC>
{
    public:
        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& len_per_irrep,
                          const vector<int>& sym)
        {
            assert(len_per_irrep.size() == sym.size());
            assert(len_per_irrep.size() == 0 ||
                   len_per_irrep[0].size() == group.getNumIrreps());
            addInitializer(TensorInitializer<PGSYMMETRIC>(group, len_per_irrep));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getNumTotal(len_per_irrep)));
            addInitializer(TensorInitializer<INDEXABLE>(len_per_irrep.size()));
            addInitializer(TensorInitializer<IPSYMMETRIC_>(sym));
        }

        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& len_per_irrep,
                          const vector<int>& sym,
                          const symmetry::Representation& rep)
        {
            assert(len_per_irrep.size() == sym.size());
            assert(len_per_irrep.size() == 0 ||
                   len_per_irrep[0].size() == group.getNumIrreps());
            addInitializer(TensorInitializer<PGSYMMETRIC>(group, len_per_irrep, rep));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getNumTotal(len_per_irrep)));
            addInitializer(TensorInitializer<INDEXABLE>(len_per_irrep.size()));
            addInitializer(TensorInitializer<IPSYMMETRIC_>(sym));
        }
};

TENSOR_INTERFACE(PGSYMMETRIC_)
{
    public:
        virtual const symmetry::PointGroup& getPointGroup() const = 0;

        virtual const symmetry::Representation& getRepresentation() const = 0;

        virtual const vector<vector<int>>& getLengthsPerIrrep() const = 0;

        virtual void getAllKeysByIrrep(const vector<int>& irreps, KeyVector& keys) const { assert(0); };

        virtual void getAllDataByIrrep(const vector<int>& irreps, KeyValueVector& kv) const { assert(0); };

        virtual void getDataByIrrep(const vector<int>& irreps, key_type n, key_type* keys, void* values) const { assert(0); };

        virtual void setDataByIrrep(const vector<int>& irreps, key_type n, const key_type* keys, const void* values) { assert(0); };

        virtual void addDataByIrrep(const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                                    const void* values, const Scalar& beta) { assert(0); };

        virtual void getLocalKeysByIrrep(const vector<int>& irreps, KeyVector& keys) const { assert(0); }

        virtual void getLocalDataByIrrep(const vector<int>& irreps, KeyValueVector& kv) const { assert(0); }

        virtual void setLocalDataByIrrep(const vector<int>& irreps, key_type n, const key_type* keys, const void* values) { assert(0); }

        virtual void addLocalDataByIrrep(const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                                  const void* values, const Scalar& beta) { assert(0); }

        virtual void getRemoteDataByIrrep(const vector<int>& irreps, key_type n, key_type* keys, void* values) const { assert(0); }

        virtual void setRemoteDataByIrrep(const vector<int>& irreps, key_type n, const key_type* keys, const void* values) { assert(0); }

        virtual void addRemoteDataByIrrep(const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                                   const void* values, const Scalar& beta) { assert(0); }

        virtual void getRemoteDataByIrrep(const vector<int>& irreps) const { assert(0); }

        virtual void setRemoteDataByIrrep(const vector<int>& irreps) { assert(0); }

        virtual void addRemoteDataByIrrep(const vector<int>& irreps, const Scalar& alpha, const Scalar& beta) { assert(0); }

        virtual void sliceByIrrep(const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, const TensorImplementation<>& A,
                                  const Scalar&  beta,             const vector<vector<int>>& start_B, const vector<vector<int>>& length) = 0;
};

TENSOR_WRAPPER(PGSYMMETRIC_)
{
    public:
        const symmetry::PointGroup& getPointGroup() const
        {
            return this->template impl<PGSYMMETRIC>().getPointGroup();
        }

        const symmetry::Representation& getRepresentation() const
        {
            return this->template impl<PGSYMMETRIC>().getRepresentation();
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        const vector<vector<int>>& getLengthsPerIrrep() const
        {
            return this->template impl<PGSYMMETRIC>().getLengthsPerIrrep();
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        KeyVector getAllKeysByIrrep(const vector<int>& irreps) const
        {
            KeyVector keys;
            this->template impl<PGSYMMETRIC>().getAllKeysByIrrep(irreps, keys);
            return keys;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        void getAllKeysByIrrep(const vector<int>& irreps, KeyVector& keys) const
        {
            this->template impl<PGSYMMETRIC>().getAllKeysByIrrep(irreps, keys);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        KeyValueVector getAllDataByIrrep(const vector<int>& irreps) const
        {
            KeyValueVector kv(this->impl().F);
            this->template impl<PGSYMMETRIC>().getAllDataByIrrep(irreps, kv);
            return kv;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        void getAllDataByIrrep(const vector<int>& irreps, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().getAllDataByIrrep(irreps, kv);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        void getDataByIrrep(const vector<int>& irreps, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().getDataByIrrep(irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED)>>
        enable_if_field_t<T>
        getDataByIrrep(const vector<int>& irreps, KeyVector& keys, vector<T>& values) const
        {
            assert(this->impl().F == Field(T()));
            values.resize(keys.size());
            this->template impl<PGSYMMETRIC>().getDataByIrrep(irreps, keys.size(), keys.data(), static_cast<void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void setDataByIrrep(const vector<int>& irreps, const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().setDataByIrrep(irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        enable_if_field_t<T>
        setDataByIrrep(const vector<int>& irreps, const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<PGSYMMETRIC>().setDataByIrrep(irreps, keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void addDataByIrrep(const vector<int>& irreps, const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().addDataByIrrep(irreps, kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        enable_if_field_t<T>
        addDataByIrrep(const vector<int>& irreps, const Scalar& alpha, const KeyVector& keys,
                const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<PGSYMMETRIC>().addDataByIrrep(irreps, keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        KeyVector getLocalKeysByIrrep(const vector<int>& irreps) const
        {
            KeyVector keys;
            getLocalKeysByIrrep(irreps, keys);
            return keys;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void getLocalKeysByIrrep(const vector<int>& irreps, KeyVector& keys) const
        {
            this->template impl<PGSYMMETRIC>().getLocalKeysByIrrep(irreps, keys);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        KeyValueVector getLocalDataByIrrep(const vector<int>& irreps) const
        {
            KeyValueVector kv(this->impl().F);
            getLocalDataByIrrep(irreps, kv);
            return kv;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void getLocalDataByIrrep(const vector<int>& irreps, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().getLocalDataByIrrep(irreps, kv);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void setLocalDataByIrrep(const vector<int>& irreps, const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().setLocalDataByIrrep(irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        enable_if_field_t<T>
        setLocalDataByIrrep(const vector<int>& irreps, const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<PGSYMMETRIC>().setLocalDataByIrrep(irreps, keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void addLocalDataByIrrep(const vector<int>& irreps, const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().addLocalDataByIrrep(irreps, kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        enable_if_field_t<T>
        addLocalDataByIrrep(const vector<int>& irreps, const Scalar& alpha, const KeyVector& keys,
                     const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<PGSYMMETRIC>().addLocalDataByIrrep(irreps, keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void getRemoteDataByIrrep(const vector<int>& irreps) const
        {
            this->template impl<DISTRIBUTED>().getRemoteDataByIrrep(irreps);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void setRemoteDataByIrrep(const vector<int>& irreps)
        {
            this->template impl<DISTRIBUTED>().setRemoteDataByIrrep(irreps);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void addRemoteDataByIrrep(const vector<int>& irreps, const Scalar& alpha, const Scalar& beta)
        {
            this->template impl<PGSYMMETRIC>().addRemoteDataByIrrep(irreps, alpha, beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void getRemoteDataByIrrep(const vector<int>& irreps, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().getRemoteDataByIrrep(irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        enable_if_field_t<T>
        getRemoteDataByIrrep(const vector<int>& irreps, KeyVector& keys, vector<T>& values) const
        {
            assert(this->impl().F == Field(T()));
            values.resize(keys.size());
            this->template impl<PGSYMMETRIC>().getRemoteDataByIrrep(irreps, keys.size(), keys.data(), static_cast<void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void setRemoteDataByIrrep(const vector<int>& irreps, const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().setRemoteDataByIrrep(irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        enable_if_field_t<T>
        setRemoteDataByIrrep(const vector<int>& irreps, const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<PGSYMMETRIC>().setRemoteDataByIrrep(irreps, keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        void addRemoteDataByIrrep(const vector<int>& irreps, const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<PGSYMMETRIC>().addRemoteDataByIrrep(irreps, kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED|DISTRIBUTED)>>
        enable_if_field_t<T>
        addRemoteDataByIrrep(const vector<int>& irreps, const Scalar& alpha, const KeyVector& keys,
                      const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<PGSYMMETRIC>().addRemoteDataByIrrep(irreps, keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep (const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, ConstTensor<BOUNDED|PGSYMMETRIC> A,
                           const Scalar&  beta,             const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            this->template impl<PGSYMMETRIC>().sliceByIrrep(alpha, conja, start_A, A.impl(),
                                                             beta,        start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(const Scalar& alpha, const vector<vector<int>>& start_A, ConstTensor<BOUNDED|PGSYMMETRIC> A,
                          const Scalar&  beta, const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            sliceByIrrep(alpha, false, start_A, A, beta, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(bool conja, const vector<vector<int>>& start_A, ConstTensor<BOUNDED|PGSYMMETRIC> A,
                                      const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            sliceByIrrep(1, conja, start_A, A, 0, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(const vector<vector<int>>& start_A, ConstTensor<BOUNDED|PGSYMMETRIC> A,
                          const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            sliceByIrrep(1, false, start_A, A, 0, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(const Scalar& alpha, bool conja, ConstTensor<BOUNDED|PGSYMMETRIC> A,
                          const Scalar&  beta,             const vector<vector<int>>& start_B)
        {
            sliceByIrrep(alpha, conja, vector<vector<int>>(this->getDimension(), vector<int>(this->getPointGroup().getNumIrreps())), A, beta, start_B, A.getLengthsPerIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(const Scalar& alpha, ConstTensor<BOUNDED|PGSYMMETRIC> A,
                          const Scalar&  beta, const vector<vector<int>>& start_B)
        {
            sliceByIrrep(alpha, false, vector<vector<int>>(this->getDimension(), vector<int>(this->getPointGroup().getNumIrreps())), A, beta, start_B, A.getLengthsPerIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(bool conja, ConstTensor<BOUNDED|PGSYMMETRIC> A, const vector<vector<int>>& start_B)
        {
            sliceByIrrep(1, conja, vector<vector<int>>(this->getDimension(), vector<int>(this->getPointGroup().getNumIrreps())), A, 0, start_B, A.getLengthsPerIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(ConstTensor<BOUNDED|PGSYMMETRIC> A, const vector<vector<int>>& start_B)
        {
            sliceByIrrep(1, false, vector<vector<int>>(this->getDimension(), vector<int>(this->getPointGroup().getNumIrreps())), A, 0, start_B, A.getLengthsPerIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, ConstTensor<BOUNDED|PGSYMMETRIC> A, const Scalar&  beta)
        {
            sliceByIrrep(alpha, conja, start_A, A, beta, vector<vector<int>>(this->getDimension(), vector<int>(this->getPointGroup().getNumIrreps())), this->getLengthsPerIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(const Scalar& alpha, const vector<vector<int>>& start_A, ConstTensor<BOUNDED|PGSYMMETRIC> A, const Scalar&  beta)
        {
            sliceByIrrep(alpha, false, start_A, A, beta, vector<vector<int>>(this->getDimension(), vector<int>(this->getPointGroup().getNumIrreps())), this->getLengthsPerIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(bool conja, const vector<vector<int>>& start_A, ConstTensor<BOUNDED|PGSYMMETRIC> A)
        {
            sliceByIrrep(1, conja, start_A, A, 0, vector<vector<int>>(this->getDimension(), vector<int>(this->getPointGroup().getNumIrreps())), this->getLengthsPerIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,BOUNDED)>>
        void sliceByIrrep(const vector<vector<int>>& start_A, ConstTensor<BOUNDED|PGSYMMETRIC> A)
        {
            sliceByIrrep(1, false, start_A, A, 0, vector<vector<int>>(this->getDimension(), vector<int>(this->getPointGroup().getNumIrreps())), this->getLengthsPerIrrep());
        }
};

TENSOR_DEFINITION(PGSYMMETRIC_)
{
    public:
        TensorDefinition(const INITIALIZER_TYPE(C)& ilist) : Base(ilist) {}

        const symmetry::PointGroup& getPointGroup() const
        {
            return this->group;
        }

        const symmetry::Representation& getRepresentation() const
        {
            return this->rep;
        }

        const vector<vector<int>>& getLengthsPerIrrep() const
        {
            return this->len_per_irrep;
        }

        void slice(const Scalar& alpha, bool conja, const vector<int>& start_A, const TensorImplementation<>& A,
                   const Scalar&  beta,             const vector<int>& start_B, const vector<int>& length)
        {
            int n = this->group.getNumIrreps();
            int na = A.as<INDEXABLE_>().getDimension();
            int nb = this->ndim;
            const vector<vector<int>>& la = A.as<PGSYMMETRIC_>().getLengthsPerIrrep();
            const vector<vector<int>>& lb = this->len_per_irrep;

            vector<vector<int>> start_A_irr(na, vector<int>(n));
            vector<vector<int>> start_B_irr(nb, vector<int>(n));
            vector<vector<int>> length_irr(nb, vector<int>(n));

            for (int i = 0;i < na;i++)
            {
                int s = start_A[i];
                for (int j = 0;j < n;j++)
                {
                    if (s < la[i][j] && s >= 0) start_A_irr[i][j] = s;
                    s -= la[i][j];
                }
            }

            for (int i = 0;i < nb;i++)
            {
                int s = start_B[i];
                int num = 0;
                for (int j = 0;j < n;j++)
                {
                    if (s < lb[i][j] && s >= 0)
                    {
                        start_B_irr[i][j] = s;
                        num = length[i];
                    }
                    s -= lb[i][j];
                    length_irr[i][j] = min(num,lb[i][j]);
                    num -= length_irr[i][j];
                }
            }

            this->sliceByIrrep(alpha, conja, start_A_irr, A, beta, start_B_irr, length_irr);
        }

        void getDataByIrrep(const vector<int>& irreps, key_type n, key_type* keys, void* values) const
        {
            const Arena& arena = this->template as<DISTRIBUTED_>().getArena();

            if (arena.rank == 0)
            {
                this->getRemoteDataByIrrep(irreps, n, keys, values);
            }
            else
            {
                this->getRemoteDataByIrrep(irreps);
            }
            arena.comm().Bcast(values, n, 0, Field::MPI_TYPE(this->F));
        }

        void setDataByIrrep(const vector<int>& irreps, key_type n, const key_type* keys, const void* values)
        {
            auto& dist = this->template as<DISTRIBUTED_>();

            if (dist.getArena().rank == 0)
            {
                this->setRemoteDataByIrrep(irreps, n, keys, values);
            }
            else
            {
                this->setRemoteDataByIrrep(irreps);
            }
        }

        void addDataByIrrep(const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                     const void* values, const Scalar& beta)
        {
            auto& dist = this->template as<DISTRIBUTED_>();

            if (dist.getArena().rank == 0)
            {
                this->addRemoteDataByIrrep(irreps, n, alpha, keys, values, beta);
            }
            else
            {
                this->addRemoteDataByIrrep(irreps, alpha, beta);
            }
        }
};

}
}

#endif
