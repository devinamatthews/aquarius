#ifndef _AQUARIUS_SCF_UHF_COMMON_HPP_
#define _AQUARIUS_SCF_UHF_COMMON_HPP_

#include "../../frameworks/convergence/diis.hpp"
#include "../../frameworks/input/config.hpp"
#include "../../frameworks/input/molecule.hpp"
#include "../../frameworks/integrals/1eints.hpp"
#include "../../frameworks/operator/space.hpp"
#include "../../frameworks/task/task.hpp"
#include "../../frameworks/util/global.hpp"
#include "../../frameworks/util/iterative.hpp"
#include "tensor/symblocked_tensor.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class UHF : public Iterative<T>
{
    protected:
        bool frozen_core;
        T damping;
        vector<int> occ_alpha, occ_beta;
        vector<vector<real_type_t<T>>> E_alpha, E_beta;
        convergence::DIIS<tensor::SymmetryBlockedTensor<T>> diis;

    public:
        UHF(const string& name, input::Config& config);

        void iterate(const Arena& arena);

        bool run(task::TaskDAG& dag, const Arena& arena);

    protected:
        virtual void calcSMinusHalf() = 0;

        void calcS2();

        virtual void diagonalizeFock() = 0;

        virtual void buildFock() = 0;

        void calcEnergy();

        void calcDensity();

        void DIISExtrap();
};

}
}

#endif
