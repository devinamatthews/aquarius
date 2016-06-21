#ifndef _AQUARIUS_FRAMEWORKS_TENSOR_SPINORBITAL_HPP_
#define _AQUARIUS_FRAMEWORKS_TENSOR_SPINORBITAL_HPP_

#include "frameworks/util.hpp"
#include "frameworks/tensor.hpp"

namespace aquarius
{
namespace tensor
{

namespace detail
{
    inline vector<int> getSymmetry(const vector<int>& ncreation,
                                   const vector<int>& nannihilation)
    {
        vector<int> sym;

        int m = ncreation.size();
        assert(m == nannihilation.size());

        int n = 0;
        for (int i = 0;i < m;i++) n += ncreation[i];
        for (int i = 0;i < m;i++) n += nannihilation[i];

        sym.resize(n, AS);

        int j = -1;
        for (int i = 0;i < m;i++)
        {
            assert(ncreation[i] >= 0);
            j += ncreation[i];
            if (j >= 0) sym[j] = NS;
        }
        for (int i = 0;i < m;i++)
        {
            assert(nannihilation[i] >= 0);
            j += nannihilation[i];
            if (j >= 0) sym[j] = NS;
        }

        return sym;
    }

    inline vector<int> getLengths(const vector<int>& nalpha,
                                  const vector<int>& nbeta,
                                  const vector<int>& ncreation,
                                  const vector<int>& nannihilation)
    {
        vector<int> len;

        int m = nalpha.size();
        assert(m == nbeta.size());
        assert(m == ncreation.size());
        assert(m == nannihilation.size());

        int n = 0;
        for (int i = 0;i < m;i++) n += ncreation[i];
        for (int i = 0;i < m;i++) n += nannihilation[i];

        len.resize(n, -1);

        int j = 0;
        for (int i = 0;i < m;i++)
        {
            assert(ncreation[i] >= 0);
            n = nalpha[i]+nbeta[i];
            for (int k = 0;k < ncreation[i];k++, j++) len[j] = n;
        }
        for (int i = 0;i < m;i++)
        {
            assert(nannihilation[i] >= 0);
            n = nalpha[i]+nbeta[i];
            for (int k = 0;k < nannihilation[i];k++, j++) len[j] = n;
        }

        return len;
    }
}

template <> class TensorInitializer<SPINORBITAL_>
{
    public:
        vector<int> nalpha;
        vector<int> nbeta;
        vector<vector<int>> nalpha_per_irrep;
        vector<vector<int>> nbeta_per_irrep;
		vector<int> ncreation;
	    vector<int> nannihilation;
	    int spin;

        TensorInitializer(const vector<int>& nalpha,
                          const vector<int>& nbeta,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          int spin)
        : nalpha(nalpha), nbeta(nbeta), ncreation(ncreation),
          nannihilation(nannihilation), spin(spin) {}

        TensorInitializer(const vector<vector<int>>& nalpha_per_irrep,
                          const vector<vector<int>>& nbeta_per_irrep,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          int spin)
        : nalpha(detail::getNumTotal(nalpha_per_irrep)), nbeta(detail::getNumTotal(nbeta_per_irrep)),
          nalpha_per_irrep(nalpha_per_irrep), nbeta_per_irrep(nbeta_per_irrep),
          ncreation(ncreation), nannihilation(nannihilation), spin(spin) {}
};

template <> class TensorInitializer<SPINORBITAL>
: public TensorInitializerList<SPINORBITAL>
{
    public:
        TensorInitializer(const vector<int>& nalpha,
                          const vector<int>& nbeta,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          int spin = 0)
        {
            addInitializer(TensorInitializer<IPSYMMETRIC_>(detail::getSymmetry(ncreation,
                                                                               nannihilation)));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getLengths(nalpha, nbeta,
                                                                          ncreation,
                                                                          nannihilation)));
            addInitializer(TensorInitializer<SPINORBITAL_>(nalpha, nbeta, ncreation,
                                                           nannihilation, spin));
        }

        TensorInitializer(const vector<int>& nalpha,
                          const vector<int>& nbeta)
        {
            addInitializer(TensorInitializer<IPSYMMETRIC_>({}));
            addInitializer(TensorInitializer<BOUNDED_>({}));
            addInitializer(TensorInitializer<SPINORBITAL_>(nalpha, nbeta, {}, {}, 0));
        }

        template <capability_type C, typename=enable_if_t< IS_SUPERSET_OF(C,SPINORBITAL) &&
                                                          !IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        TensorInitializer(const TensorInitializerList<C>& ilist,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          int spin = 0)
        {
            addInitializer(TensorInitializer<IPSYMMETRIC_>(detail::getSymmetry(ncreation,
                                                                               nannihilation)));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getLengths(ilist.template as<SPINORBITAL_>().nalpha,
                                                                          ilist.template as<SPINORBITAL_>().nbeta,
                                                                          ncreation,
                                                                          nannihilation)));
            addInitializer(TensorInitializer<SPINORBITAL_>(ilist.template as<SPINORBITAL_>().nalpha,
                                                           ilist.template as<SPINORBITAL_>().nbeta,
                                                           ncreation, nannihilation, spin));
        }
};

typedef TensorInitializerList<SPINORBITAL> SOInit;

template <> class TensorInitializer<SPINORBITAL|PGSYMMETRIC>
: public TensorInitializerList<SPINORBITAL|PGSYMMETRIC>
{
    public:
        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& nalpha,
                          const vector<vector<int>>& nbeta,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          int spin = 0)
        {
            addInitializer(TensorInitializer<PGSYMMETRIC_>(group));
            addInitializer(TensorInitializer<IPSYMMETRIC_>(detail::getSymmetry(ncreation,
                                                                               nannihilation)));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getLengths(detail::getNumTotal(nalpha),
                                                                          detail::getNumTotal(nbeta),
                                                                          ncreation, nannihilation)));
            addInitializer(TensorInitializer<SPINORBITAL_>(nalpha, nbeta, ncreation,
                                                           nannihilation, spin));
        }

        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& nalpha,
                          const vector<vector<int>>& nbeta)
        {
            addInitializer(TensorInitializer<PGSYMMETRIC_>(group));
            addInitializer(TensorInitializer<IPSYMMETRIC_>({}));
            addInitializer(TensorInitializer<BOUNDED_>({}));
            addInitializer(TensorInitializer<SPINORBITAL_>(nalpha, nbeta, {}, {}, 0));
        }

        TensorInitializer(const symmetry::PointGroup& group,
                          const symmetry::Representation& rep,
                          const vector<vector<int>>& nalpha,
                          const vector<vector<int>>& nbeta,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          int spin = 0)
        {
            addInitializer(TensorInitializer<PGSYMMETRIC_>(group, rep));
            addInitializer(TensorInitializer<IPSYMMETRIC_>(detail::getSymmetry(ncreation,
                                                                               nannihilation)));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getLengths(detail::getNumTotal(nalpha),
                                                                          detail::getNumTotal(nbeta),
                                                                          ncreation, nannihilation)));
            addInitializer(TensorInitializer<SPINORBITAL_>(nalpha, nbeta, ncreation,
                                                           nannihilation, spin));
        }

        TensorInitializer(const symmetry::PointGroup& group,
                          const symmetry::Representation& rep,
                          const vector<vector<int>>& nalpha,
                          const vector<vector<int>>& nbeta)
        {
            addInitializer(TensorInitializer<PGSYMMETRIC_>(group, rep));
            addInitializer(TensorInitializer<IPSYMMETRIC_>({}));
            addInitializer(TensorInitializer<BOUNDED_>({}));
            addInitializer(TensorInitializer<SPINORBITAL_>(nalpha, nbeta, {}, {}, 0));
        }

        template <capability_type C, typename=enable_if_t<IS_SUPERSET_OF(C,SPINORBITAL) &&
                                                          IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        TensorInitializer(const TensorInitializerList<C>& ilist,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          int spin = 0)
        {
                addInitializer(TensorInitializer<IPSYMMETRIC_>(detail::getSymmetry(ncreation,
                                                                                   nannihilation)));
                addInitializer(TensorInitializer<BOUNDED_>(detail::getLengths(ilist.template as<SPINORBITAL_>().nalpha,
                                                                              ilist.template as<SPINORBITAL_>().nbeta,
                                                                              ncreation,
                                                                              nannihilation)));
                addInitializer(TensorInitializer<SPINORBITAL_>(ilist.template as<SPINORBITAL_>().nalpha_per_irrep,
                                                               ilist.template as<SPINORBITAL_>().nbeta_per_irrep,
                                                               ncreation, nannihilation, spin));
        }
};

typedef TensorInitializerList<SPINORBITAL|PGSYMMETRIC> SOPGInit;

TENSOR_INTERFACE(SPINORBITAL_)
{
    public:
        virtual int getNumClasses() const = 0;

        virtual int getSpin() const = 0;

        virtual vector<vector<int>> getLengthsPerSpin() const = 0;

        virtual vector<vector<vector<int>>> getLengthsPerSpinAndIrrep() const = 0;

        virtual const vector<int>& getNumAlpha() const = 0;

        virtual const vector<int>& getNumBeta() const = 0;

        virtual const vector<int>& getNumCreation() const = 0;

        virtual const vector<int>& getNumAnnihilation() const = 0;

        virtual const vector<vector<int>>& getNumAlphaPerIrrep() const = 0;

        virtual const vector<vector<int>>& getNumBetaPerIrrep() const = 0;

        virtual void getAllKeysBySpin(const vector<int>& spins, KeyVector& keys) const = 0;

        virtual void getAllDataBySpin(const vector<int>& spins, KeyValueVector& kv) const = 0;

        virtual void getDataBySpin(const vector<int>& spins, key_type n, key_type* keys, void* values) const = 0;

        virtual void setDataBySpin(const vector<int>& spins, key_type n, const key_type* keys, const void* values) = 0;

        virtual void addDataBySpin(const vector<int>& spins, key_type n, const Scalar& alpha, const key_type* keys,
                                    const void* values, const Scalar& beta) = 0;

        virtual void getLocalKeysBySpin(const vector<int>& spins, KeyVector& keys) const { assert(0); }

        virtual void getLocalDataBySpin(const vector<int>& spins, KeyValueVector& kv) const { assert(0); }

        virtual void setLocalDataBySpin(const vector<int>& spins, key_type n, const key_type* keys, const void* values) { assert(0); }

        virtual void addLocalDataBySpin(const vector<int>& spins, key_type n, const Scalar& alpha, const key_type* keys,
                                  const void* values, const Scalar& beta) { assert(0); }

        virtual void getRemoteDataBySpin(const vector<int>& spins, key_type n, key_type* keys, void* values) const { assert(0); }

        virtual void setRemoteDataBySpin(const vector<int>& spins, key_type n, const key_type* keys, const void* values) { assert(0); }

        virtual void addRemoteDataBySpin(const vector<int>& spins, key_type n, const Scalar& alpha, const key_type* keys,
                                   const void* values, const Scalar& beta) { assert(0); }

        virtual void getRemoteDataBySpin(const vector<int>& spins) const { assert(0); }

        virtual void setRemoteDataBySpin(const vector<int>& spins) { assert(0); }

        virtual void addRemoteDataBySpin(const vector<int>& spins, const Scalar& alpha, const Scalar& beta) { assert(0); }

        virtual void getAllKeysBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyVector& keys) const = 0;

        virtual void getAllDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyValueVector& kv) const = 0;

        virtual void getDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, key_type* keys, void* values) const = 0;

        virtual void setDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const key_type* keys, const void* values) = 0;

        virtual void addDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                                    const void* values, const Scalar& beta) = 0;

        virtual void getLocalKeysBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyVector& keys) const { assert(0); }

        virtual void getLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyValueVector& kv) const { assert(0); }

        virtual void setLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const key_type* keys, const void* values) { assert(0); }

        virtual void addLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                                  const void* values, const Scalar& beta) { assert(0); }

        virtual void getRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, key_type* keys, void* values) const { assert(0); }

        virtual void setRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const key_type* keys, const void* values) { assert(0); }

        virtual void addRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                                   const void* values, const Scalar& beta) { assert(0); }

        virtual void getRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps) const { assert(0); }

        virtual void setRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps) { assert(0); }

        virtual void addRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const Scalar& alpha, const Scalar& beta) { assert(0); }

        virtual void sliceBySpin(const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, const TensorImplementation<>& A,
                                 const Scalar&  beta,             const vector<vector<int>>& start_B, const vector<vector<int>>& length) = 0;

        virtual void sliceBySpinAndIrrep(const Scalar& alpha, bool conja, const vector<vector<vector<int>>>& start_A, const TensorImplementation<>& A,
                                         const Scalar&  beta,             const vector<vector<vector<int>>>& start_B, const vector<vector<vector<int>>>& length) = 0;
};

TENSOR_WRAPPER(SPINORBITAL_)
{
    public:
        int getNumClasses() const
        {
            return this->template impl<SPINORBITAL_>().getNumClasses();
        }

        int getSpin() const
        {
            return this->template impl<SPINORBITAL_>().getSpin();
        }

        vector<vector<int>> getLengthsPerSpin() const
        {
            return this->template impl<SPINORBITAL_>().getLengthsPerSpin();
        }

        vector<vector<vector<int>>> getLengthsPerSpinAndIrrep() const
        {
            return this->template impl<SPINORBITAL_>().getLengthsPerSpinAndIrrep();
        }

        const vector<int>& getNumAlpha() const
        {
            return this->template impl<SPINORBITAL_>().getNumAlpha();
        }

        const vector<int>& getNumBeta() const
        {
            return this->template impl<SPINORBITAL_>().getNumBeta();
        }

        const vector<int>& getNumCreation() const
        {
            return this->template impl<SPINORBITAL_>().getNumCreation();
        }

        const vector<int>& getNumAnnihilation() const
        {
            return this->template impl<SPINORBITAL_>().getNumAnnihilation();
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        const vector<vector<int>>& getNumAlphaPerIrrep() const
        {
            return this->template impl<SPINORBITAL_>().getNumAlphaPerIrrep();
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        const vector<vector<int>>& getNumBetaPerIrrep() const
        {
            return this->template impl<SPINORBITAL_>().getNumBetaPerIrrep();
        }

        KeyVector getAllKeysBySpin(const vector<int>& spins) const
        {
            KeyVector keys;
            this->template impl<SPINORBITAL_>().getAllKeysBySpin(spins, keys);
            return keys;
        }

        void getAllKeysBySpin(const vector<int>& spins, KeyVector& keys) const
        {
            this->template impl<SPINORBITAL_>().getAllKeysBySpin(spins, keys);
        }

        KeyValueVector getAllDataBySpin(const vector<int>& spins) const
        {
            KeyValueVector kv(this->impl().F);
            this->template impl<SPINORBITAL_>().getAllDataBySpin(spins, kv);
            return kv;
        }

        void getAllDataBySpin(const vector<int>& spins, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().getAllDataBySpin(spins, kv);
        }

        void getDataBySpin(const vector<int>& spins, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().getDataBySpin(spins, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T> typename enable_if<is_field<T>::value>::type
        getDataBySpin(const vector<int>& spins, KeyVector& keys, vector<T>& values) const
        {
            assert(this->impl().F == Field(T()));
            values.resize(keys.size());
            this->template impl<SPINORBITAL_>().getDataBySpin(spins, keys.size(), keys.data(), static_cast<void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void setDataBySpin(const vector<int>& spins, const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().setDataBySpin(spins, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        enable_if_field_t<T>
        setDataBySpin(const vector<int>& spins, const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().setDataBySpin(spins, keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void addDataBySpin(const vector<int>& spins, const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().addDataBySpin(spins, kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        enable_if_field_t<T>
        addDataBySpin(const vector<int>& spins, const Scalar& alpha, const KeyVector& keys,
                const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().addDataBySpin(spins, keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        KeyVector getLocalKeysBySpin(const vector<int>& spins) const
        {
            KeyVector keys;
            getLocalKeysBySpin(spins, keys);
            return keys;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void getLocalKeysBySpin(const vector<int>& spins, KeyVector& keys) const
        {
            this->template impl<SPINORBITAL_>().getLocalKeysBySpin(spins, keys);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        KeyValueVector getLocalDataBySpin(const vector<int>& spins) const
        {
            KeyValueVector kv(this->impl().F);
            getLocalDataBySpin(spins, kv);
            return kv;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void getLocalDataBySpin(const vector<int>& spins, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().getLocalDataBySpin(spins, kv);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void setLocalDataBySpin(const vector<int>& spins, const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().setLocalDataBySpin(spins, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        enable_if_field_t<T>
        setLocalDataBySpin(const vector<int>& spins, const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().setLocalDataBySpin(spins, keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void addLocalDataBySpin(const vector<int>& spins, const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().addLocalDataBySpin(spins, kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        enable_if_field_t<T>
        addLocalDataBySpin(const vector<int>& spins, const Scalar& alpha, const KeyVector& keys,
                     const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().addLocalDataBySpin(spins, keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void getRemoteDataBySpin(const vector<int>& spins) const
        {
            this->template impl<SPINORBITAL_>().getRemoteDataBySpin(spins);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void setRemoteDataBySpin(const vector<int>& spins)
        {
            this->template impl<SPINORBITAL_>().setRemoteDataBySpin(spins);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void addRemoteDataBySpin(const vector<int>& spins, const Scalar& alpha, const Scalar& beta)
        {
            this->template impl<SPINORBITAL_>().addRemoteDataBySpin(spins, alpha, beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void getRemoteDataBySpin(const vector<int>& spins, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().getRemoteDataBySpin(spins, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        enable_if_field_t<T>
        getRemoteDataBySpin(const vector<int>& spins, KeyVector& keys, vector<T>& values) const
        {
            assert(this->impl().F == Field(T()));
            values.resize(keys.size());
            this->template impl<SPINORBITAL_>().getRemoteDataBySpin(spins, keys.size(), keys.data(), static_cast<void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void setRemoteDataBySpin(const vector<int>& spins, const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().setRemoteDataBySpin(spins, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        enable_if_field_t<T>
        setRemoteDataBySpin(const vector<int>& spins, const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().setRemoteDataBySpin(spins, keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        void addRemoteDataBySpin(const vector<int>& spins, const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().addRemoteDataBySpin(spins, kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,DISTRIBUTED)>>
        enable_if_field_t<T>
        addRemoteDataBySpin(const vector<int>& spins, const Scalar& alpha, const KeyVector& keys,
                      const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().addRemoteDataBySpin(spins, keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        KeyVector getAllKeysBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps) const
        {
            KeyVector keys;
            this->template impl<SPINORBITAL_>().getAllKeysBySpinAndIrrep(spins, irreps, keys);
            return keys;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        void getAllKeysBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyVector& keys) const
        {
            this->template impl<SPINORBITAL_>().getAllKeysBySpinAndIrrep(spins, irreps, keys);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        KeyValueVector getAllDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps) const
        {
            KeyValueVector kv(this->impl().F);
            this->template impl<SPINORBITAL_>().getAllDataBySpinAndIrrep(spins, irreps, kv);
            return kv;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        void getAllDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().getAllDataBySpinAndIrrep(spins, irreps, kv);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        void getDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().getDataBySpinAndIrrep(spins, irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        enable_if_field_t<T>
        getDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyVector& keys, vector<T>& values) const
        {
            assert(this->impl().F == Field(T()));
            values.resize(keys.size());
            this->template impl<SPINORBITAL_>().getDataBySpinAndIrrep(spins, irreps, keys.size(), keys.data(), static_cast<void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        void setDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().setDataBySpinAndIrrep(spins, irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        enable_if_field_t<T>
        setDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().setDataBySpinAndIrrep(spins, irreps, keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        void addDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().addDataBySpinAndIrrep(spins, irreps, kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC)>>
        enable_if_field_t<T>
        addDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const Scalar& alpha, const KeyVector& keys,
                const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().addDataBySpinAndIrrep(spins, irreps, keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        KeyVector getLocalKeysBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps) const
        {
            KeyVector keys;
            getLocalKeysBySpinAndIrrep(spins, irreps, keys);
            return keys;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void getLocalKeysBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyVector& keys) const
        {
            this->template impl<SPINORBITAL_>().getLocalKeysBySpinAndIrrep(spins, irreps, keys);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        KeyValueVector getLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps) const
        {
            KeyValueVector kv(this->impl().F);
            getLocalDataBySpinAndIrrep(spins, irreps, kv);
            return kv;
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void getLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().getLocalDataBySpinAndIrrep(spins, irreps, kv);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void setLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().setLocalDataBySpinAndIrrep(spins, irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        enable_if_field_t<T>
        setLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().setLocalDataBySpinAndIrrep(spins, irreps, keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void addLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().addLocalDataBySpinAndIrrep(spins, irreps, kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        enable_if_field_t<T>
        addLocalDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const Scalar& alpha, const KeyVector& keys,
                     const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().addLocalDataBySpinAndIrrep(spins, irreps, keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void getRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps) const
        {
            this->template impl<SPINORBITAL_>().getRemoteDataBySpinAndIrrep(spins, irreps);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void setRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps)
        {
            this->template impl<SPINORBITAL_>().setRemoteDataBySpinAndIrrep(spins, irreps);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void addRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const Scalar& alpha, const Scalar& beta)
        {
            this->template impl<SPINORBITAL_>().addRemoteDataBySpinAndIrrep(spins, irreps, alpha, beta);
        }

        template <capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void getRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyValueVector& kv) const
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().getRemoteDataBySpinAndIrrep(spins, irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        enable_if_field_t<T>
        getRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, KeyVector& keys, vector<T>& values) const
        {
            assert(this->impl().F == Field(T()));
            values.resize(keys.size());
            this->template impl<SPINORBITAL_>().getRemoteDataBySpinAndIrrep(spins, irreps, keys.size(), keys.data(), static_cast<void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void setRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const KeyValueVector& kv)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().setRemoteDataBySpinAndIrrep(spins, irreps, kv.size(), kv.keys().data(), kv.data<void>());
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        enable_if_field_t<T>
        setRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const KeyVector& keys, const vector<T>& values)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().setRemoteDataBySpinAndIrrep(spins, irreps, keys.size(), keys.data(), static_cast<const void*>(values.data()));
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        void addRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const Scalar& alpha, const KeyValueVector& kv, const Scalar& beta)
        {
            assert(this->impl().F == kv.field());
            this->template impl<SPINORBITAL_>().addRemoteDataBySpinAndIrrep(spins, irreps, kv.size(), alpha, kv.keys().data(), kv.data<void>(), beta);
        }

        template <typename T, capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C_,PGSYMMETRIC|DISTRIBUTED)>>
        enable_if_field_t<T>
        addRemoteDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, const Scalar& alpha, const KeyVector& keys,
                      const vector<T>& values, const Scalar& beta)
        {
            assert(this->impl().F == Field(T()));
            assert(keys.size() == values.size());
            this->template impl<SPINORBITAL_>().addRemoteDataBySpinAndIrrep(spins, irreps, keys.size(), alpha, keys.data(), static_cast<const void*>(values.data()), beta);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, ConstTensor<SPINORBITAL> A,
                         const Scalar&  beta,             const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            this->template impl<SPINORBITAL_>().sliceBySpin(alpha, conja, start_A, A.impl(),
                                                             beta,        start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(const Scalar& alpha, const vector<vector<int>>& start_A, ConstTensor<SPINORBITAL> A,
                         const Scalar&  beta, const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            sliceBySpin(alpha, false, start_A, A, beta, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(bool conja, const vector<vector<int>>& start_A, ConstTensor<SPINORBITAL> A,
                                     const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            sliceBySpin(1, conja, start_A, A, 0, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(const vector<vector<int>>& start_A, ConstTensor<SPINORBITAL> A,
                         const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            sliceBySpin(1, false, start_A, A, 0, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(const Scalar& alpha, bool conja, ConstTensor<SPINORBITAL> A,
                         const Scalar&  beta,             const vector<vector<int>>& start_B)
        {
            sliceBySpin(alpha, conja, vector<vector<int>>(this->getDimension(), vector<int>(2)), A, beta, start_B, A.getLengthsPerSpin());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(const Scalar& alpha, ConstTensor<SPINORBITAL> A,
                         const Scalar&  beta, const vector<vector<int>>& start_B)
        {
            sliceBySpin(alpha, false, vector<vector<int>>(this->getDimension(), vector<int>(2)), A, beta, start_B, A.getLengthsPerSpin());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(bool conja, ConstTensor<SPINORBITAL> A, const vector<vector<int>>& start_B)
        {
            sliceBySpin(1, conja, vector<vector<int>>(this->getDimension(), vector<int>(2)), A, 0, start_B, A.getLengthsPerSpin());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(ConstTensor<SPINORBITAL> A, const vector<vector<int>>& start_B)
        {
            sliceBySpin(1, false, vector<vector<int>>(this->getDimension(), vector<int>(2)), A, 0, start_B, A.getLengthsPerSpin());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, ConstTensor<SPINORBITAL> A, const Scalar&  beta)
        {
            sliceBySpin(alpha, conja, start_A, A, beta, vector<vector<int>>(this->getDimension(), vector<int>(2)), this->getLengthsPerSpin());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(const Scalar& alpha, const vector<vector<int>>& start_A, ConstTensor<SPINORBITAL> A, const Scalar&  beta)
        {
            sliceBySpin(alpha, false, start_A, A, beta, vector<vector<int>>(this->getDimension(), vector<int>(2)), this->getLengthsPerSpin());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(bool conja, const vector<vector<int>>& start_A, ConstTensor<SPINORBITAL> A)
        {
            sliceBySpin(1, conja, start_A, A, 0, vector<vector<int>>(this->getDimension(), vector<int>(2)), this->getLengthsPerSpin());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)>>
        void sliceBySpin(const vector<vector<int>>& start_A, ConstTensor<SPINORBITAL> A)
        {
            sliceBySpin(1, false, start_A, A, 0, vector<vector<int>>(this->getDimension(), vector<int>(2)), this->getLengthsPerSpin());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep (const Scalar& alpha, bool conja, const vector<vector<vector<int>>>& start_A, ConstTensor<SPINORBITAL|PGSYMMETRIC> A,
                                  const Scalar&  beta,             const vector<vector<vector<int>>>& start_B, const vector<vector<vector<int>>>& length)
        {
            this->template impl<SPINORBITAL_>().sliceBySpinAndIrrep(alpha, conja, start_A, A.impl(),
                                                                     beta,        start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(const Scalar& alpha, const vector<vector<vector<int>>>& start_A, ConstTensor<SPINORBITAL|PGSYMMETRIC> A,
                                 const Scalar&  beta, const vector<vector<vector<int>>>& start_B, const vector<vector<vector<int>>>& length)
        {
            sliceBySpinAndIrrep(alpha, false, start_A, A, beta, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(bool conja, const vector<vector<vector<int>>>& start_A, ConstTensor<SPINORBITAL|PGSYMMETRIC> A,
                                             const vector<vector<vector<int>>>& start_B, const vector<vector<vector<int>>>& length)
        {
            sliceBySpinAndIrrep(1, conja, start_A, A, 0, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(const vector<vector<vector<int>>>& start_A, ConstTensor<SPINORBITAL|PGSYMMETRIC> A,
                                 const vector<vector<vector<int>>>& start_B, const vector<vector<vector<int>>>& length)
        {
            sliceBySpinAndIrrep(1, false, start_A, A, 0, start_B, length);
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(const Scalar& alpha, bool conja, ConstTensor<SPINORBITAL|PGSYMMETRIC> A,
                                 const Scalar&  beta,             const vector<vector<vector<int>>>& start_B)
        {
            sliceBySpinAndIrrep(alpha, conja, vector<vector<vector<int>>>(this->getDimension(), vector<vector<int>>(2, vector<int>(this->getPointGroup().getNumIrreps()))), A, beta, start_B, A.getLengthsPerSpinAndIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(const Scalar& alpha, ConstTensor<SPINORBITAL|PGSYMMETRIC> A,
                                 const Scalar&  beta, const vector<vector<vector<int>>>& start_B)
        {
            sliceBySpinAndIrrep(alpha, false, vector<vector<vector<int>>>(this->getDimension(), vector<vector<int>>(2, vector<int>(this->getPointGroup().getNumIrreps()))), A, beta, start_B, A.getLengthsPerSpinAndIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(bool conja, ConstTensor<SPINORBITAL|PGSYMMETRIC> A, const vector<vector<vector<int>>>& start_B)
        {
            sliceBySpinAndIrrep(1, conja, vector<vector<vector<int>>>(this->getDimension(), vector<vector<int>>(2, vector<int>(this->getPointGroup().getNumIrreps()))), A, 0, start_B, A.getLengthsPerSpinAndIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(ConstTensor<SPINORBITAL|PGSYMMETRIC> A, const vector<vector<vector<int>>>& start_B)
        {
            sliceBySpinAndIrrep(1, false, vector<vector<vector<int>>>(this->getDimension(), vector<vector<int>>(2, vector<int>(this->getPointGroup().getNumIrreps()))), A, 0, start_B, A.getLengthsPerSpinAndIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(const Scalar& alpha, bool conja, const vector<vector<vector<int>>>& start_A, ConstTensor<SPINORBITAL|PGSYMMETRIC> A, const Scalar&  beta)
        {
            sliceBySpinAndIrrep(alpha, conja, start_A, A, beta, vector<vector<vector<int>>>(this->getDimension(), vector<vector<int>>(2, vector<int>(this->getPointGroup().getNumIrreps()))), this->getLengthsPerSpinAndIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(const Scalar& alpha, const vector<vector<vector<int>>>& start_A, ConstTensor<SPINORBITAL|PGSYMMETRIC> A, const Scalar&  beta)
        {
            sliceBySpinAndIrrep(alpha, false, start_A, A, beta, vector<vector<vector<int>>>(this->getDimension(), vector<vector<int>>(2, vector<int>(this->getPointGroup().getNumIrreps()))), this->getLengthsPerSpinAndIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(bool conja, const vector<vector<vector<int>>>& start_A, ConstTensor<SPINORBITAL|PGSYMMETRIC> A)
        {
            sliceBySpinAndIrrep(1, conja, start_A, A, 0, vector<vector<vector<int>>>(this->getDimension(), vector<vector<int>>(2, vector<int>(this->getPointGroup().getNumIrreps()))), this->getLengthsPerSpinAndIrrep());
        }

        template <capability_type C_=C, typename=enable_if_t<!(C_&CONST_)&&IS_SUPERSET_OF(C,PGSYMMETRIC)>>
        void sliceBySpinAndIrrep(const vector<vector<vector<int>>>& start_A, ConstTensor<SPINORBITAL|PGSYMMETRIC> A)
        {
            sliceBySpinAndIrrep(1, false, start_A, A, 0, vector<vector<vector<int>>>(this->getDimension(), vector<vector<int>>(2, vector<int>(this->getPointGroup().getNumIrreps()))), this->getLengthsPerSpinAndIrrep());
        }
};

TENSOR_DEFINITION(SPINORBITAL_)
{
    public:
        TensorDefinition(const INITIALIZER_TYPE(C)& ilist) : Base(ilist) {}

        int getNumClasses() const
        {
            return this->nalpha.size();
        }

        int getSpin() const
        {
            return this->spin;
        }

        vector<vector<int>> getLengthsPerSpin() const
        {
            vector<vector<int>> lens(this->ndim, vector<int>(2));
            int k = 0;
            for (int i = 0;i < this->ncreation.size();i++)
            {
                for (int j = 0;j < this->ncreation[i];j++, k++)
                {
                    lens[k][0] = this->nalpha[i];
                    lens[k][1] = this->nbeta[i];
                }
            }
            for (int i = 0;i < this->nannihilation.size();i++)
            {
                for (int j = 0;j < this->nannihilation[i];j++, k++)
                {
                    lens[k][0] = this->nalpha[i];
                    lens[k][1] = this->nbeta[i];
                }
            }
            return lens;
        }

        vector<vector<vector<int>>> getLengthsPerSpinAndIrrep() const
        {
            vector<vector<vector<int>>> lens(this->ndim, vector<vector<int>>(2));
            int k = 0;
            for (int i = 0;i < this->ncreation.size();i++)
            {
                for (int j = 0;j < this->ncreation[i];j++, k++)
                {
                    lens[k][0] = this->nalpha_per_irrep[i];
                    lens[k][1] = this->nbeta_per_irrep[i];
                }
            }
            for (int i = 0;i < this->nannihilation.size();i++)
            {
                for (int j = 0;j < this->nannihilation[i];j++, k++)
                {
                    lens[k][0] = this->nalpha_per_irrep[i];
                    lens[k][1] = this->nbeta_per_irrep[i];
                }
            }
            return lens;
        }

        const vector<int>& getNumAlpha() const
        {
            return this->nalpha;
        }

        const vector<int>& getNumBeta() const
        {
            return this->nbeta;
        }

        const vector<int>& getNumCreation() const
        {
            return this->ncreation;
        }

        const vector<int>& getNumAnnihilation() const
        {
            return this->nannihilation;
        }

        const vector<vector<int>>& getNumAlphaPerIrrep() const
        {
            return this->nalpha_per_irrep;
        }

        const vector<vector<int>>& getNumBetaPerIrrep() const
        {
            return this->nbeta_per_irrep;
        }

        void slice(const Scalar& alpha, bool conja, const vector<int>& start_A, const TensorImplementation<>& A,
                   const Scalar&  beta,             const vector<int>& start_B, const vector<int>& length)
        {
            int na = A.as<INDEXABLE_>().getDimension();
            int nb = this->ndim;
            vector<vector<int>> la = {A.as<SPINORBITAL_>().getNumAlpha(),
                                      A.as<SPINORBITAL_>().getNumBeta()};
            vector<vector<int>> lb = {this->nalpha, this->nbeta};

            vector<vector<int>> start_A_spin(na, vector<int>(2));
            vector<vector<int>> start_B_spin(nb, vector<int>(2));
            vector<vector<int>> length_spin(nb, vector<int>(2));

            for (int i = 0;i < na;i++)
            {
                int s = start_A[i];
                for (int j = 0;j < 2;j++)
                {
                    if (s < la[i][j] && s >= 0) start_A_spin[i][j] = s;
                    s -= la[i][j];
                }
            }

            for (int i = 0;i < nb;i++)
            {
                int s = start_B[i];
                int num = 0;
                for (int j = 0;j < 2;j++)
                {
                    if (s < lb[i][j] && s >= 0)
                    {
                        start_B_spin[i][j] = s;
                        num = length[i];
                    }
                    s -= lb[i][j];
                    length_spin[i][j] = min(num,lb[i][j]);
                    num -= length_spin[i][j];
                }
            }

            this->sliceBySpin(alpha, conja, start_A_spin, A, beta, start_B_spin, length_spin);
        }

        void sliceByIrrep(const Scalar& alpha, bool conja, const vector<vector<int>>& start_A, const TensorImplementation<>& A,
                          const Scalar&  beta,             const vector<vector<int>>& start_B, const vector<vector<int>>& length)
        {
            int n = this->getPointGroup().getNumIrreps();
            int na = A.as<INDEXABLE_>().getDimension();
            int nb = this->ndim;
            vector<vector<vector<int>>> la = {A.as<SPINORBITAL_>().getNumAlphaPerIrrep(),
                                              A.as<SPINORBITAL_>().getNumBetaPerIrrep()};
            vector<vector<vector<int>>> lb = {this->nalpha_per_irrep, this->nbeta_per_irrep};

            vector<vector<vector<int>>> start_A_spin(na, vector<vector<int>>(2, vector<int>(n)));
            vector<vector<vector<int>>> start_B_spin(nb, vector<vector<int>>(2, vector<int>(n)));
            vector<vector<vector<int>>> length_spin(nb, vector<vector<int>>(2, vector<int>(n)));

            for (int i = 0;i < na;i++)
            {
                for (int j = 0;j < n;j++)
                {
                    if (start_A[i][j] >= la[i][0][j])
                    {
                        start_A_spin[i][1][j] = start_A[i][j]-la[i][0][j];
                    }
                    else
                    {
                        start_A_spin[i][0][j] = start_A[i][j];
                    }
                }
            }

            for (int i = 0;i < nb;i++)
            {
                for (int j = 0;j < n;j++)
                {
                    if (start_B[i][j] >= lb[i][0][j])
                    {
                        start_B_spin[i][1][j] = start_B[i][j]-lb[i][0][j];
                        length_spin[i][1][j] = length[i][j];
                    }
                    else
                    {
                        start_B_spin[i][0][j] = start_B[i][j];
                        length_spin[i][0][j] = min(length[i][j], lb[i][0][j]-start_B[i][j]);
                        length_spin[i][1][j] = length[i][j]-length_spin[i][0][j];
                    }
                }
            }

            this->sliceBySpinAndIrrep(alpha, conja, start_A_spin, A, beta, start_B_spin, length_spin);
        }

        void getDataBySpin(const vector<int>& spins, key_type n, key_type* keys, void* values) const
        {
            const Arena& arena = this->template as<DISTRIBUTED_>().getArena();

            if (arena.rank == 0)
            {
                this->getRemoteDataBySpin(spins, n, keys, values);
            }
            else
            {
                this->getRemoteDataBySpin(spins);
            }
            arena.comm().Bcast(values, n, 0, Field::MPI_TYPE(this->F));
        }

        void setDataBySpin(const vector<int>& spins, key_type n, const key_type* keys, const void* values)
        {
            auto& dist = this->template as<DISTRIBUTED_>();

            if (dist.getArena().rank == 0)
            {
                this->setRemoteDataBySpin(spins, n, keys, values);
            }
            else
            {
                this->setRemoteDataBySpin(spins);
            }
        }

        void addDataBySpin(const vector<int>& spins, key_type n, const Scalar& alpha, const key_type* keys,
                     const void* values, const Scalar& beta)
        {
            auto& dist = this->template as<DISTRIBUTED_>();

            if (dist.getArena().rank == 0)
            {
                this->addRemoteDataBySpin(spins, n, alpha, keys, values, beta);
            }
            else
            {
                this->addRemoteDataBySpin(spins, alpha, beta);
            }
        }

        void getDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, key_type* keys, void* values) const
        {
            const Arena& arena = this->template as<DISTRIBUTED_>().getArena();

            if (arena.rank == 0)
            {
                this->getRemoteDataBySpinAndIrrep(spins, irreps, n, keys, values);
            }
            else
            {
                this->getRemoteDataBySpinAndIrrep(spins, irreps);
            }
            arena.comm().Bcast(values, n, 0, Field::MPI_TYPE(this->F));
        }

        void setDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const key_type* keys, const void* values)
        {
            auto& dist = this->template as<DISTRIBUTED_>();

            if (dist.getArena().rank == 0)
            {
                this->setRemoteDataBySpinAndIrrep(spins, irreps, n, keys, values);
            }
            else
            {
                this->setRemoteDataBySpinAndIrrep(spins, irreps);
            }
        }

        void addDataBySpinAndIrrep(const vector<int>& spins, const vector<int>& irreps, key_type n, const Scalar& alpha, const key_type* keys,
                     const void* values, const Scalar& beta)
        {
            auto& dist = this->template as<DISTRIBUTED_>();

            if (dist.getArena().rank == 0)
            {
                this->addRemoteDataBySpinAndIrrep(spins, irreps, n, alpha, keys, values, beta);
            }
            else
            {
                this->addRemoteDataBySpinAndIrrep(spins, irreps, alpha, beta);
            }
        }
};

}
}

#endif
