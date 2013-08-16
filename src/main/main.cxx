/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "input/config.hpp"
#include "input/molecule.hpp"
#include "integrals/1eints.hpp"
#include "integrals/2eints.hpp"
#include "scf/aoscf.hpp"
#include "operator/aomoints.hpp"
#include "cc/ccsd.hpp"
#include "cc/lambdaccsd.hpp"
#include "cc/2edensity.hpp"
#include "operator/st2eoperator.hpp"
#include "time/time.hpp"
#include "task/task.hpp"

#ifdef USE_ELEMENTAL
#include "elemental.hpp"
using namespace elem;
#endif

using namespace std;
using namespace MPI;
using namespace aquarius;
using namespace aquarius::integrals;
using namespace aquarius::input;
using namespace aquarius::scf;
using namespace aquarius::cc;
using namespace aquarius::op;
using namespace aquarius::time;
using namespace aquarius::tensor;
using namespace aquarius::task;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
#ifdef USE_ELEMENTAL
    elem::Initialize(argc, argv);
#endif

    {
        int i;
        double dt;
        Arena world;

        assert(argc > 1);
        Schema schema(TOPDIR "/input_schema");
        Config config(argv[1]);
        schema.apply(config);

        Molecule* mol_ = new Molecule(world, config);
        Molecule& mol = *mol_;
        Product p("molecule", "molecule");
        p.put(mol_);

        Task* t1 = Task::createTask("1eints", "1eints", config);
        Task* t2 = Task::createTask("2eints", "2eints", config);
        t1->getProduct("S").getRequirements()[0].fulfil(p);
        t1->getProduct("T").getRequirements()[0].fulfil(p);
        t1->getProduct("G").getRequirements()[0].fulfil(p);
        t1->getProduct("H").getRequirements()[0].fulfil(p);
        t2->getProduct("I").getRequirements()[0].fulfil(p);
        Product S = t1->getProduct("S");
        Product T = t1->getProduct("T");
        Product G = t1->getProduct("G");
        Product H = t1->getProduct("H");
        Product I = t2->getProduct("I");
        TaskDAG dag;
        dag.addTask(t1);
        dag.addTask(t2);

        PRINT("nA: %d\n", mol.getNumOrbitals()-mol.getNumAlphaElectrons());
        PRINT("na: %d\n", mol.getNumOrbitals()-mol.getNumBetaElectrons());
        PRINT("nI: %d\n", mol.getNumAlphaElectrons());
        PRINT("ni: %d\n", mol.getNumBetaElectrons());

        tic();
        dag.execute(world);
        dt = todouble(toc());
        PRINT("\nAO integrals took: %8.3f s\n", dt);

        tic();
        AOUHF<double> scf(config.get("scf"), mol,
                          I.get<ERI>(),
                          S.get<OVI>(),
                          H.get<OneElectronHamiltonian>());

        PRINT("\nUHF-SCF\n\n");
        PRINT("It.            SCF Energy     Residual Walltime\n");
        tic();
        for (i = 0;scf.iterate();i++)
        {
            dt = todouble(toc());
            PRINT("%3d % 21.15f %12.6e %8.3f\n", i+1, scf.getEnergy(), scf.getConvergence(), dt);
            tic();
        }
        toc();

        PRINT("\nUHF Orbital Energies\n\n");
        for (int i = 0;i < mol.getNumOrbitals();i++)
        {
            PRINT("%4d ", i+1);

            if (i < mol.getNumAlphaElectrons())
            {
                PRINT("%21.15f a ", scf.getAlphaEigenvalues()[i]);
            }
            else
            {
                PRINT("%21.15f   ", scf.getAlphaEigenvalues()[i]);
            }

            if (i < mol.getNumBetaElectrons())
            {
                PRINT("%21.15f b ", scf.getBetaEigenvalues()[i]);
            }
            else
            {
                PRINT("%21.15f   ", scf.getBetaEigenvalues()[i]);
            }

            PRINT("\n");
        }

        dt = todouble(toc());
        PRINT("\nAO SCF took: %8.3f s (%8.3f s/it.)\n", dt, dt/i);

        double s2 = scf.getS2();
        double mult = scf.getMultiplicity();
        double na = scf.getAvgNumAlpha();
        double nb = scf.getAvgNumBeta();

        PRINT("\n");
        PRINT("<0|S^2|0>     = %f\n", s2);
        PRINT("<0|2S+1|0>    = %f\n", mult);
        PRINT("<0|n_alpha|0> = %f\n", na);
        PRINT("<0|n_beta|0>  = %f\n", nb);

        tic();
        AOMOIntegrals<double> moints(scf);
        dt = todouble(toc());
        PRINT("\nAO MO took: %8.3f s\n", dt);

        CCSD<double> ccsd(config.get("cc"), moints);

        PRINT("\nUHF-MP2 Energy: %.15f\n", ccsd.getEnergy());

        s2 = ccsd.getProjectedS2();
        mult = ccsd.getProjectedMultiplicity();

        PRINT("\n");
        PRINT("<0|S^2|MP2>  = %f\n", s2);
        PRINT("<0|2S+1|MP2> = %f\n", mult);

        PRINT("\nUHF-CCSD\n\n");
        PRINT("It.   Correlation Energy     Residual Walltime\n");
        tic();
        tic();
        for (i = 0;ccsd.iterate();i++)
        {
            dt = todouble(toc());
            PRINT("%3d % 20.15f %12.6e %8.3f\n", i+1, ccsd.getEnergy(), ccsd.getConvergence(), dt);
            tic();
        }
        toc();

        dt = todouble(toc());
        PRINT("\nCCSD took: %8.3f s (%8.3f s/it.)\n", dt, dt/i);

        s2 = ccsd.getProjectedS2();
        mult = ccsd.getProjectedMultiplicity();

        PRINT("\n");
        PRINT("<0|S^2|CC>  = %f\n", s2);
        PRINT("<0|2S+1|CC> = %f\n", mult);

        tic();
        STTwoElectronOperator<double,2> Hbar(moints, ccsd);
        dt = todouble(toc());
        PRINT("             _\n");
        PRINT("Formation of H took: %8.3f s\n", dt);

        LambdaCCSD<double> lambda(config.get("cc"), Hbar, ccsd, ccsd.getEnergy());

        PRINT("\nUHF-Lambda-CCSD\n\n");
        PRINT("It.   Correlation Energy     Residual Walltime\n");
        tic();
        tic();
        for (i = 0;lambda.iterate();i++)
        {
            dt = todouble(toc());
            PRINT("%3d % 20.15f %12.6e %8.3f\n", i+1, lambda.getEnergy(), lambda.getConvergence(), dt);
            tic();
        }
        toc();

        dt = todouble(toc());
        PRINT("\nLambda-CCSD took: %8.3f s (%8.3f s/it.)\n", dt, dt/i);

        PRINT("\nFinal Energy: %.15f\n\n", scf.getEnergy()+ccsd.getEnergy());

        print_timers();
    }

#ifdef USE_ELEMENTAL
    elem::Finalize();
#endif
    MPI_Finalize();
}
