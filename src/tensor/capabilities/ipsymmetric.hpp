#ifndef _AQUARIUS_TENSOR_IPSYMMETRIC_TENSOR_HPP_
#define _AQUARIUS_TENSOR_IPSYMMETRIC_TENSOR_HPP_

#include "util/global.hpp"

#include "tensor/tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <> class TensorInitializer<IPSYMMETRIC_> : public Destructible
{
    public:
        vector<int> sym;

        TensorInitializer(const vector<int>& sym) : sym(sym) {}
};

template <> class TensorInitializer<IPSYMMETRIC> : public TensorInitializerList<IPSYMMETRIC>
{
    public:
        TensorInitializer(const vector<int>& sym)
    	{
            addInitializer(TensorInitializer<INDEXABLE>(sym.size()));
            addInitializer(TensorInitializer<IPSYMMETRIC_>(sym));
    	}
};

TENSOR_DEFINITION(IPSYMMETRIC_)
{
    public:
        TensorDefinition(const INITIALIZER_TYPE(C)& ilist) : Base(ilist) {}

        const vector<int>& getSymmetry() const { return this->sym; }
};

TENSOR_INTERFACE(IPSYMMETRIC_)
{
    public:
        virtual const vector<int>& getSymmetry() const = 0;
};

TENSOR_WRAPPER(IPSYMMETRIC_)
{
    public:
        const vector<int>& getSymmetry() const
        {
            return this->template impl<IPSYMMETRIC_>().getSymmetry();
        }
};

}
}

#endif
