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
#include "scf/aoints.hpp"
#include "scf/aoscf.hpp"
#include "scf/aomoints.hpp"
#include "scf/choleskyscf.hpp"
#include "scf/choleskymoints.hpp"
#include "cc/ccsd.hpp"
#include "cc/ccsdt.hpp"
#include "cc/lambdaccsd.hpp"
#include "operator/st2eoperator.hpp"
#include "time/time.hpp"
#include "tensor/dist_tensor.hpp"

using namespace std;
using namespace elem;
using namespace MPI;
using namespace aquarius::slide;
using namespace aquarius::input;
using namespace aquarius::scf;
using namespace aquarius::cc;
using namespace aquarius::op;
using namespace aquarius::time;
using namespace aquarius::tensor;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    SLIDE::init();
    elem::Initialize(argc, argv);

    {
        tCTF_World<double> ctf;

        assert(argc > 1);
        Schema schema(TOPDIR "/input_schema");
        Config config(argv[1]);
        schema.apply(config);

        Molecule mol(config);

        CholeskyIntegrals<double> chol(ctf, config.get("cholesky"), mol);
        CholeskyUHF<double> scf(config.get("scf"), chol);

        //chol.test();

        PRINT("UHF-SCF\n\n");
        PRINT("It.            SCF Energy     Residual\n");
        for (int i = 0;scf.iterate();i++)
        {
            PRINT("%3d % 21.15f %12.6e\n", i+1, scf.getEnergy(), scf.getConvergence());
        }

        double s2 = scf.getS2();
        double mult = scf.getMultiplicity();
        double na = scf.getAvgNumAlpha();
        double nb = scf.getAvgNumBeta();

        PRINT("\n");
        PRINT("<0|S^2|0>     = %f\n", s2);
        PRINT("<0|2S+1|0>    = %f\n", mult);
        PRINT("<0|n_alpha|0> = %f\n", na);
        PRINT("<0|n_beta|0>  = %f\n", nb);
        PRINT("\n");

        AOIntegrals<double> ao(ctf, mol);
        AOUHF<double> aoscf(config.get("scf"), ao);

        PRINT("UHF-SCF\n\n");
        PRINT("It.            SCF Energy     Residual\n");
        for (int i = 0;aoscf.iterate();i++)
        {
            PRINT("%3d % 21.15f %12.6e\n", i+1, aoscf.getEnergy(), aoscf.getConvergence());
        }

        s2 = aoscf.getS2();
        mult = aoscf.getMultiplicity();
        na = aoscf.getAvgNumAlpha();
        nb = aoscf.getAvgNumBeta();

        PRINT("\n");
        PRINT("<0|S^2|0>     = %f\n", s2);
        PRINT("<0|2S+1|0>    = %f\n", mult);
        PRINT("<0|n_alpha|0> = %f\n", na);
        PRINT("<0|n_beta|0>  = %f\n", nb);
        PRINT("\n");

        CholeskyMOIntegrals<double> moints(scf);
        AOMOIntegrals<double> aomo(aoscf);

        CCSDT<double> ccsdt(config.get("cc"), moints);

        PRINT("UHF-MP2 Energy: %.15f\n", ccsdt.getEnergy());

        PRINT("\nUHF-CCSDT\n\n");
        PRINT("It.   Correlation Energy     Residual\n");
        for (int i = 0;ccsdt.iterate();i++)
        {
            PRINT("%3d % 20.15f %12.6e\n", i+1, ccsdt.getEnergy(), ccsdt.getConvergence());
        }

        s2 = ccsdt.getProjectedS2();
        mult = ccsdt.getProjectedMultiplicity();

        PRINT("\n");
        PRINT("<0|S^2|CC>  = %f\n", s2);
        PRINT("<0|2S+1|CC> = %f\n", mult);
        PRINT("\n");

        /*
        STTwoElectronOperator<double,2> H(moints, ccsd);
        LambdaCCSD<double> lambda(config.get("cc"), H, ccsd, ccsd.getEnergy());

        PRINT("UHF-Lambda-CCSD\n\n");
        PRINT("It.   Correlation Energy     Residual\n");
        for (int i = 0;lambda.iterate();i++)
        {
            PRINT("%3d % 20.15f %12.6e\n", i+1, lambda.getEnergy(), lambda.getConvergence());
        }

        PRINT("\nFinal Energy: %.15f\n\n", scf.getEnergy()+ccsd.getEnergy());
        */

        print_timers();
    }

    elem::Finalize();
    SLIDE::finish();
    MPI_Finalize();
}
