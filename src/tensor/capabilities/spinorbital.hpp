#ifndef _AQUARIUS_TENSOR_SPINORBITAL_TENSOR_HPP_
#define _AQUARIUS_TENSOR_SPINORBITAL_TENSOR_HPP_

#include "util/global.hpp"

#include "tensor/tensor.hpp"

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

template <> class TensorInitializer<SPINORBITAL_> : public Destructible
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
: TensorInitializerList<SPINORBITAL>
{
    public:
        TensorInitializer(const vector<int>& nalpha,
                          const vector<int>& nbeta,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          int spin = 0)
        {
                addInitializer(TensorInitializer<IPSYMMETRIC>(detail::getSymmetry(ncreation,
                                                                                  nannihilation)));
                addInitializer(TensorInitializer<BOUNDED_>(detail::getLengths(nalpha, nbeta,
                                                                              ncreation,
                                                                              nannihilation)));
                addInitializer(TensorInitializer<SPINORBITAL_>(nalpha, nbeta, ncreation,
                                                               nannihilation, spin));
        }
};

template <> class TensorInitializer<SPINORBITAL|PGSYMMETRIC>
: TensorInitializerList<SPINORBITAL|PGSYMMETRIC>
{
    public:
        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& nalpha,
                          const vector<vector<int>>& nbeta,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          int spin = 0)
        {
            addInitializer(TensorInitializer<PGSYMMETRIC>(group));
            addInitializer(TensorInitializer<IPSYMMETRIC>(detail::getSymmetry(ncreation,
                                                                              nannihilation)));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getLengths(detail::getNumTotal(nalpha),
                                                                          detail::getNumTotal(nbeta),
                                                                          ncreation, nannihilation)));
            addInitializer(TensorInitializer<SPINORBITAL_>(nalpha, nbeta, ncreation,
                                                           nannihilation, spin));
        }

        TensorInitializer(const symmetry::PointGroup& group,
                          const vector<vector<int>>& nalpha,
                          const vector<vector<int>>& nbeta,
                          const vector<int>& ncreation,
                          const vector<int>& nannihilation,
                          const symmetry::Representation& rep,
                          int spin = 0)
        {
            addInitializer(TensorInitializer<PGSYMMETRIC>(group, rep));
            addInitializer(TensorInitializer<IPSYMMETRIC>(detail::getSymmetry(ncreation,
                                                                              nannihilation)));
            addInitializer(TensorInitializer<BOUNDED_>(detail::getLengths(detail::getNumTotal(nalpha),
                                                                          detail::getNumTotal(nbeta),
                                                                          ncreation, nannihilation)));
            addInitializer(TensorInitializer<SPINORBITAL_>(nalpha, nbeta, ncreation,
                                                           nannihilation, spin));
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
};

TENSOR_INTERFACE(SPINORBITAL_)
{
    public:
        virtual int getNumClasses() const = 0;

        virtual int getSpin() const = 0;

        virtual const vector<int>& getNumAlpha() const = 0;

        virtual const vector<int>& getNumBeta() const = 0;

        virtual const vector<int>& getNumCreation() const = 0;

        virtual const vector<int>& getNumAnnihilation() const = 0;

        virtual const vector<vector<int>>& getNumAlphaPerIrrep() const = 0;

        virtual const vector<vector<int>>& getNumBetaPerIrrep() const = 0;
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

        const vector<vector<int>>& getNumAlphaPerIrrep() const
        {
            static_assert(C&PGSYMMETRIC, "The operand must be PGSYMMETRIC.");
            return this->template impl<SPINORBITAL_>().getNumAlphaPerIrrep();
        }

        const vector<vector<int>>& getNumBetaPerIrrep() const
        {
            static_assert(C&PGSYMMETRIC, "The operand must be PGSYMMETRIC.");
            return this->template impl<SPINORBITAL_>().getNumBetaPerIrrep();
        }
};

}
}

#endif
