#ifndef _AQUARIUS_OPERATOR_EXCITATIONOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_EXCITATIONOPERATOR_HPP_

#include "util/global.hpp"

#include "tensor/composite_tensor.hpp"
#include "tensor/spinorbital_tensor.hpp"

#include "mooperator.hpp"
#include "denominator.hpp"

namespace aquarius
{
namespace op
{

template <typename T, int np, int nh=np>
class ExcitationOperator
: public MOOperator,
  public tensor::CompositeTensor< ExcitationOperator<T,np,nh>,
                                  tensor::SpinorbitalTensor<T>, T >
{
    INHERIT_FROM_COMPOSITE_TENSOR(CONCAT(ExcitationOperator<T,np,nh>),
                                  tensor::SpinorbitalTensor<T>, T)

    protected:
        const int spin;

    public:
        ExcitationOperator(const ExcitationOperator& other)
        : MOOperator(other.arena, other.occ, other.vrt),
          tensor::CompositeTensor< ExcitationOperator<T,np,nh>,
           tensor::SpinorbitalTensor<T>, T >(other.name, max(np,nh)+1),
          spin(other.spin)
        {
            for (int ex = 0;ex <= min(np,nh);ex++)
            {
                int nv = ex+(np > nh ? np-nh : 0);
                int no = ex+(nh > np ? nh-np : 0);

                tensors[ex+abs(np-nh)].isAlloced = true;
                tensors[ex+abs(np-nh)].tensor =
                    new tensor::SpinorbitalTensor<T>(*other.tensors[ex+abs(np-nh)].tensor);
            }
        }

        ExcitationOperator(const string& name, const Arena& arena, const Space& occ, const Space& vrt, int spin=0)
        : MOOperator(arena, occ, vrt),
          tensor::CompositeTensor< ExcitationOperator<T,np,nh>,
           tensor::SpinorbitalTensor<T>, T >(name, max(np,nh)+1),
          spin(spin)
        {
            for (int ex = 0;ex <= min(np,nh);ex++)
            {
                int nv = ex+(np > nh ? np-nh : 0);
                int no = ex+(nh > np ? nh-np : 0);

                tensors[ex+abs(np-nh)].isAlloced = true;
                tensors[ex+abs(np-nh)].tensor =
                    new tensor::SpinorbitalTensor<T>(name, arena, occ.group, {vrt,occ}, {nv,0}, {0,no}, spin);
            }
        }

        ExcitationOperator(const string& name, const Arena& arena, const Space& occ, const Space& vrt,
                           const symmetry::Representation& rep, int spin=0)
        : MOOperator(arena, occ, vrt),
          tensor::CompositeTensor< ExcitationOperator<T,np,nh>,
           tensor::SpinorbitalTensor<T>, T >(name, max(np,nh)+1),
          spin(spin)
        {
            for (int ex = 0;ex <= min(np,nh);ex++)
            {
                int nv = ex+(np > nh ? np-nh : 0);
                int no = ex+(nh > np ? nh-np : 0);

                tensors[ex+abs(np-nh)].isAlloced = true;
                tensors[ex+abs(np-nh)].tensor =
                    new tensor::SpinorbitalTensor<T>(name, arena, occ.group, rep, {vrt,occ}, {nv,0}, {0,no}, spin);
            }
        }

        const symmetry::PointGroup& getGroup() const { return tensors[abs(np-nh)].tensor->getGroup(); }

        const symmetry::Representation& getRepresentation() const { return tensors[abs(np-nh)].tensor->getRepresentation(); }

        void weight(const Denominator<T>& d, double shift = 0)
        {
            vector<const vector<vector<T>>*> da{&d.getDA(), &d.getDI()};
            vector<const vector<vector<T>>*> db{&d.getDa(), &d.getDi()};

            for (int ex = 0;ex <= min(np,nh);ex++)
            {
                if (ex== 0 && np == nh) continue;
                tensors[ex+abs(np-nh)].tensor->weight(da, db, shift);
            }
        }

        T dot(bool conja, const op::ExcitationOperator<T,np,nh>& A, bool conjb) const
        {
            T s = (T)0;

            for (int i = abs(np-nh);i <= max(np,nh);i++)
            {
                s += (*this)(i).dot(conja, A(i), conjb)/(T)factorial(i)/(T)factorial(i-abs(np-nh));
            }

            return s;
        }

        /*
         * Return the largest p-norm of the constituent operators
         */
        real_type_t<T> norm(int p) const
        {
            real_type_t<T> nrm = 0;

            for (int i = abs(np-nh);i <= max(np,nh);i++)
            {
                nrm = max(nrm,(*this)(i).norm(p));
            }

            return nrm;
        }

        template <int np_, int nh_, typename=enable_if<(np_ < np) && (np-np_ == nh-nh_)>>
        operator ExcitationOperator<T,np_,nh_>&()
        {
            return reinterpret_cast<ExcitationOperator<T,np_,nh_>&>(*this);
        }

        template <int np_, int nh_, typename=enable_if<(np_ < np) && (np-np_ == nh-nh_)>>
        operator const ExcitationOperator<T,np_,nh_>&() const
        {
            return reinterpret_cast<const ExcitationOperator<T,np_,nh_>&>(*this);
        }
};

}
}

#endif
