#ifndef _AQUARIUS_CC_DENSITY_HPP_
#define _AQUARIUS_CC_DENSITY_HPP_

#include "util/global.hpp"

#include "operator/2eoperator.hpp"

#include "1edensity.hpp"

namespace aquarius
{
namespace cc
{

template <typename U>
class TwoElectronDensity : public op::TwoElectronOperator<U>
{
    public:
        /*
         * Form the SCF density
         *
         * P1_ab = D_ab (2D?)
         *
         * P2_abcd = P1_ab P1_cd + 1/2 P1_ac P1_bd
         */
        TwoElectronDensity(const string& name, const op::MOSpace<U>& occ, const op::MOSpace<U>& vrt,
                           const tensor::SymmetryBlockedTensor<U>& Da,
                           const tensor::SymmetryBlockedTensor<U>& Db);

        /*
         * Form the unrelaxed CCSD Density
         */
        TwoElectronDensity(const string& name, const op::ExcitationOperator<U,2>& T);

        /*
         * Form the mixed perturbed CCSD Density
         */
        TwoElectronDensity(const string& name, const op::DeexcitationOperator<U,2>& L,
                           const op::ExcitationOperator<U,2>& T,
                           const op::ExcitationOperator<U,2>& TA);

        /*
         * Form the relaxed CCSD Density
         */
        TwoElectronDensity(const string& name, const op::DeexcitationOperator<U,2>& L,
                           const op::ExcitationOperator<U,2>& T);

        /*
         * Form the relaxed perturbed CCSD Density
         */
        TwoElectronDensity(const string& name, const op::DeexcitationOperator<U,2>& L,
                           const op::DeexcitationOperator<U,2>& LA,
                           const op::ExcitationOperator<U,2>& T,
                           const op::ExcitationOperator<U,2>& TA);

        double getS2() const;

        double getProjectedMultiplicity() const;
};

}
}

#endif
