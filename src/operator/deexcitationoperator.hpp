#ifndef _AQUARIUS_OPERATOR_DEEXCITATIONOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_DEEXCITATIONOPERATOR_HPP_

#include "util/global.hpp"

#include "tensor/composite_tensor.hpp"
#include "tensor/spinorbital_tensor.hpp"

#include "mooperator.hpp"
#include "denominator.hpp"
#include "excitationoperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T, int np, int nh=np>
class DeexcitationOperator
: public MOOperator,
  public tensor::CompositeTensor< DeexcitationOperator<T,np,nh>,
                                  tensor::SpinorbitalTensor<T>, T >
{
    INHERIT_FROM_COMPOSITE_TENSOR(CONCAT(DeexcitationOperator<T,np,nh>),
                                  tensor::SpinorbitalTensor<T>, T)

    protected:
        const int spin;

    public:
        DeexcitationOperator(const DeexcitationOperator& other)
        : MOOperator(other.arena, other.occ, other.vrt),
          tensor::CompositeTensor< DeexcitationOperator<T,np,nh>,
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

        DeexcitationOperator(const string& name, const Arena& arena, const Space& occ, const Space& vrt, int spin=0)
        : MOOperator(arena, occ, vrt),
          tensor::CompositeTensor< DeexcitationOperator<T,np,nh>,
           tensor::SpinorbitalTensor<T>, T >(name, max(np,nh)+1),
          spin(spin)
        {
            for (int ex = 0;ex <= min(np,nh);ex++)
            {
                int nv = ex+(np > nh ? np-nh : 0);
                int no = ex+(nh > np ? nh-np : 0);

                tensors[ex+abs(np-nh)].isAlloced = true;
                tensors[ex+abs(np-nh)].tensor =
                    new tensor::SpinorbitalTensor<T>(name, arena, occ.group, {vrt,occ}, {0,no}, {nv,0}, spin);
            }
        }

        DeexcitationOperator(const string& name, const Arena& arena, const Space& occ, const Space& vrt,
                             const symmetry::Representation& rep, int spin=0)
        : MOOperator(arena, occ, vrt),
          tensor::CompositeTensor< DeexcitationOperator<T,np,nh>,
           tensor::SpinorbitalTensor<T>, T >(name, max(np,nh)+1),
          spin(spin)
        {
            for (int ex = 0;ex <= min(np,nh);ex++)
            {
                int nv = ex+(np > nh ? np-nh : 0);
                int no = ex+(nh > np ? nh-np : 0);

                tensors[ex+abs(np-nh)].isAlloced = true;
                tensors[ex+abs(np-nh)].tensor =
                    new tensor::SpinorbitalTensor<T>(name, arena, occ.group, rep, {vrt,occ}, {0,no}, {nv,0}, spin);
            }
        }

        const symmetry::PointGroup& getGroup() const { return tensors[abs(np-nh)].tensor->getGroup(); }

        const symmetry::Representation& getRepresentation() const { return tensors[abs(np-nh)].tensor->getGroup(); }

        void weight(const Denominator<T>& d, double shift = 0)
        {
            vector<const vector<vector<T>>*> da{&d.getDA(), &d.getDI()};
            vector<const vector<vector<T>>*> db{&d.getDa(), &d.getDi()};

            for (int ex = 0;ex <= min(np,nh);ex++)
            {
                if (ex == 0 && np == nh) continue;
                tensors[ex+abs(np-nh)].tensor->weight(da, db, shift);
            }
        }

        T dot(bool conja, const op::DeexcitationOperator<T,np,nh>& A, bool conjb) const
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
        operator DeexcitationOperator<T,np_,nh_>&()
        {
            return reinterpret_cast<DeexcitationOperator<T,np_,nh_>&>(*this);
        }

        template <int np_, int nh_, typename=enable_if<(np_ < np) && (np-np_ == nh-nh_)>>
        operator const DeexcitationOperator<T,np_,nh_>&() const
        {
            return reinterpret_cast<const DeexcitationOperator<T,np_,nh_>&>(*this);
        }
};

namespace detail
{

template <typename T, int np, int nh>
constexpr int get_np(const ExcitationOperator<T,np,nh>&)
{
    return np;
}

template <typename T, int np, int nh>
constexpr int get_nh(const ExcitationOperator<T,np,nh>&)
{
    return nh;
}

template <typename T1, typename T2, typename T>
T exdot(const tensor::ScaledTensor<T1,T>& t1,
        const tensor::ScaledTensor<T2,T>& t2)
{
    T s = (T)0;

    int np = get_np(t2.tensor_);
    int nh = get_nh(t2.tensor_);

    for (int i = abs(np-nh);i <= max(np,nh);i++)
    {
        int nv = i - abs(np-nh) + (np > nh ? np-nh : 0);
        int no = i - abs(np-nh) + (nh > np ? nh-np : 0);

        string out, in;
        for (int j = 0;j < nv;j++) out.push_back(j);
        for (int j = 0;j < no;j++) in.push_back(j+nv);

        s += t1.factor_*t2.factor_*t1.tensor_(i).dot(t2.conj_, t2.tensor_(i), out+in, t1.conj_, in+out)/
               (T)factorial(i)/(T)factorial(i-abs(np-nh));
    }

    return s;
}

template <typename T>
struct is_excitationoperator : std::false_type {};

template <typename T, int np, int nh>
struct is_excitationoperator<ExcitationOperator<T,np,nh>> : std::true_type { typedef T type; };

template <typename T, int np, int nh>
struct is_excitationoperator<tensor::ScaledTensor<ExcitationOperator<T,np,nh>,T>> : std::true_type { typedef T type; };

template <typename T, int np, int nh>
struct is_excitationoperator<tensor::ScaledTensor<const ExcitationOperator<T,np,nh>,T>> : std::true_type { typedef T type; };

template <typename T>
struct is_deexcitationoperator : std::false_type {};

template <typename T, int np, int nh>
struct is_deexcitationoperator<DeexcitationOperator<T,np,nh>> : std::true_type { typedef T type; };

template <typename T, int np, int nh>
struct is_deexcitationoperator<tensor::ScaledTensor<DeexcitationOperator<T,np,nh>,T>> : std::true_type { typedef T type; };

template <typename T, int np, int nh>
struct is_deexcitationoperator<tensor::ScaledTensor<const DeexcitationOperator<T,np,nh>,T>> : std::true_type { typedef T type; };

}

template <typename T, typename U>
auto operator*(const T& t1, const U& t2) ->
    enable_if_t<detail::is_deexcitationoperator<T>::value &&
                detail::is_excitationoperator<U>::value,
                typename detail::is_deexcitationoperator<T>::type>
{
    return detail::exdot(t1*1, t2*1);
}

template <typename T, typename U>
auto operator*(const T& t1, const U& t2) ->
    enable_if_t<detail::is_excitationoperator<T>::value &&
                detail::is_deexcitationoperator<U>::value,
                typename detail::is_excitationoperator<T>::type>
{
    return detail::exdot(t2*1, t1*1);
}

}
}

#endif
