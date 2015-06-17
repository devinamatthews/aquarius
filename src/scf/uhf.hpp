#ifndef _AQUARIUS_SCF_UHF_COMMON_HPP_
#define _AQUARIUS_SCF_UHF_COMMON_HPP_

#include "util/global.hpp"

#include "tensor/tensor.hpp"
#include "integrals/1eints.hpp"
#include "input/molecule.hpp"
#include "input/config.hpp"
#include "util/iterative.hpp"
#include "convergence/diis.hpp"
#include "task/task.hpp"

namespace aquarius
{
namespace scf
{

class UHF : public Iterative
{
    protected:
        bool frozen_core;
        double damping;
        vector<int> occ_alpha, occ_beta;
        vector<vector<tensor::Scalar>> E_alpha, E_beta;
        convergence::DIIS diis;

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
