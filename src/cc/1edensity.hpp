#ifndef _AQUARIUS_CC_1EDENSITY_HPP_
#define _AQUARIUS_CC_1EDENSITY_HPP_

#include "util/global.hpp"

#include "operator/excitationoperator.hpp"
#include "operator/deexcitationoperator.hpp"
#include "operator/1eoperator.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class OneElectronDensity : public op::OneElectronOperator<U>
{
    public:
        /*
         * Form the SCF density
         */
        OneElectronDensity(const string& name, const op::MOSpace<U>& occ, const op::MOSpace<U>& vrt,
                           const tensor::SymmetryBlockedTensor<U>& Da,
                           const tensor::SymmetryBlockedTensor<U>& Db);

        /*
         * Form the unrelaxed CCSD density
         */
        OneElectronDensity(const string& name, const op::ExcitationOperator<U,2>& T);

        /*
         * Form the partial perturbed CCSD Density
         */
        OneElectronDensity(const string& name, const op::DeexcitationOperator<U,2>& L,
                           const op::ExcitationOperator<U,2>& T,
                           const op::ExcitationOperator<U,2>& TA);

        /*
         * Form the relaxed CCSD density
         */
        OneElectronDensity(const string& name, const op::DeexcitationOperator<U,2>& L,
                           const op::ExcitationOperator<U,2>& T);

        /*
         * Form the relaxed perturbed CCSD Density
         */
        OneElectronDensity(const string& name, const op::DeexcitationOperator<U,2>& L,
                           const op::DeexcitationOperator<U,2>& LA,
                           const op::ExcitationOperator<U,2>& T,
                           const op::ExcitationOperator<U,2>& TA);
};

}
}

#endif
