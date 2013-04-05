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
#include "cc/2edensity.hpp"
#include "operator/st2eoperator.hpp"
#include "time/time.hpp"
#include "tensor/dist_tensor.hpp"
#include "tensor/spinorbital.hpp"

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

        double diff;

        /*
        SpinorbitalTensor<DistTensor<double> > ij(moints.getIJ());
        ij -= aomo.getIJ();
        diff = ij(0,1,0,1).reduce(CTF_OP_MAXABS);
        PRINT("IJ:   %e\n", diff);
        diff = ij(0,0,0,0).reduce(CTF_OP_MAXABS);
        PRINT("ij:   %e\n", diff);

        SpinorbitalTensor<DistTensor<double> > ai(moints.getAI());
        ai -= aomo.getAI();
        diff = ai(1,0,0,1).reduce(CTF_OP_MAXABS);
        PRINT("AI:   %e\n", diff);
        diff = ai(0,0,0,0).reduce(CTF_OP_MAXABS);
        PRINT("ai:   %e\n", diff);

        SpinorbitalTensor<DistTensor<double> > ab(moints.getAB());
        ab -= aomo.getAB();
        diff = ab(1,0,1,0).reduce(CTF_OP_MAXABS);
        PRINT("AB:   %e\n", diff);
        diff = ab(0,0,0,0).reduce(CTF_OP_MAXABS);
        PRINT("ab:   %e\n", diff);

        SpinorbitalTensor<DistTensor<double> > abcd(moints.getABCD());
        abcd -= aomo.getABCD();
        diff = abcd(2,0,2,0).reduce(CTF_OP_MAXABS);
        PRINT("ABCD:   %e\n", diff);
        diff = abcd(1,0,1,0).reduce(CTF_OP_MAXABS);
        PRINT("AbCd:   %e\n", diff);
        diff = abcd(0,0,0,0).reduce(CTF_OP_MAXABS);
        PRINT("abcd:   %e\n", diff);

        SpinorbitalTensor<DistTensor<double> > abci(moints.getABCI());
        abci -= aomo.getABCI();
        diff = abci(2,0,1,1).reduce(CTF_OP_MAXABS);
        PRINT("ABCI:   %e\n", diff);
        diff = abci(1,0,1,0).reduce(CTF_OP_MAXABS);
        PRINT("AbCi:   %e\n", diff);
        diff = abci(1,0,0,1).reduce(CTF_OP_MAXABS);
        PRINT("aBcI:   %e\n", diff);
        diff = abci(0,0,0,0).reduce(CTF_OP_MAXABS);
        PRINT("abci:   %e\n", diff);

        SpinorbitalTensor<DistTensor<double> > abij(moints.getABIJ());
        abij -= aomo.getABIJ();
        diff = abij(2,0,0,2).reduce(CTF_OP_MAXABS);
        PRINT("ABIJ:   %e\n", diff);
        diff = abij(1,0,0,1).reduce(CTF_OP_MAXABS);
        PRINT("AbIj:   %e\n", diff);
        diff = abij(0,0,0,0).reduce(CTF_OP_MAXABS);
        PRINT("abij:   %e\n", diff);
        */

        //moints.getAIBJ()(1,0,0,1)["AibJ"] = -moints.getABIJ()(1,0,0,1)["AbJi"];
        //moints.getAIBJ()(0,1,1,0)["aIBj"] = -moints.getABIJ()(1,0,0,1)["BaIj"];
        //aomo.getAIBJ()(1,0,0,1)["AibJ"] = -aomo.getABIJ()(1,0,0,1)["AbJi"];
        //aomo.getAIBJ()(0,1,1,0)["aIBj"] = -aomo.getABIJ()(1,0,0,1)["BaIj"];

        /*
        SpinorbitalTensor<DistTensor<double> > aibj(moints.getAIBJ());
        aibj -= aomo.getAIBJ();
        diff = aibj(1,1,1,1).reduce(CTF_OP_MAXABS);
        PRINT("AIBJ:   %e\n", diff);
        diff = aibj(1,0,1,0).reduce(CTF_OP_MAXABS);
        PRINT("AiBj:   %e\n", diff);
        diff = aibj(0,1,0,1).reduce(CTF_OP_MAXABS);
        PRINT("aIbJ:   %e\n", diff);
        diff = aibj(0,0,0,0).reduce(CTF_OP_MAXABS);
        PRINT("aibj:   %e\n", diff);
        diff = aibj(1,0,0,1).reduce(CTF_OP_MAXABS);
        PRINT("AibJ:   %e\n", diff);
        diff = aibj(0,1,1,0).reduce(CTF_OP_MAXABS);
        PRINT("aIBj:   %e\n", diff);

        SpinorbitalTensor<DistTensor<double> > ijka(moints.getIJKA());
        ijka -= aomo.getIJKA();
        diff = ijka(0,2,1,1).reduce(CTF_OP_MAXABS);
        PRINT("IJKA:   %e\n", diff);
        diff = ijka(0,1,0,1).reduce(CTF_OP_MAXABS);
        PRINT("IjKa:   %e\n", diff);
        diff = ijka(0,1,1,0).reduce(CTF_OP_MAXABS);
        PRINT("iJkA:   %e\n", diff);
        diff = ijka(0,0,0,0).reduce(CTF_OP_MAXABS);
        PRINT("ijka:   %e\n", diff);

        SpinorbitalTensor<DistTensor<double> > ijkl(moints.getIJKL());
        ijkl -= aomo.getIJKL();
        diff = ijkl(0,2,0,2).reduce(CTF_OP_MAXABS);
        PRINT("IJKL:   %e\n", diff);
        diff = ijkl(0,1,0,1).reduce(CTF_OP_MAXABS);
        PRINT("IjKl:   %e\n", diff);
        diff = ijkl(0,0,0,0).reduce(CTF_OP_MAXABS);
        PRINT("ijkl:   %e\n", diff);
        */

        DistTensor<double> AibJ1(moints.getAIBJ()(1,0,0,1));
        AibJ1["AibJ"] += moints.getABIJ()(1,0,0,1)["AbJi"];
        diff = AibJ1.reduce(CTF_OP_MAXABS);
        PRINT("AibJ1:   %e\n", diff);

        DistTensor<double> AibJ2(aomo.getAIBJ()(1,0,0,1));
        AibJ2["AibJ"] += aomo.getABIJ()(1,0,0,1)["AbJi"];
        diff = AibJ2.reduce(CTF_OP_MAXABS);
        PRINT("AibJ2:   %e\n", diff);

        DistTensor<double> aIBj1(moints.getAIBJ()(0,1,1,0));
        aIBj1["aIBj"] += moints.getABIJ()(1,0,0,1)["BaIj"];
        diff = aIBj1.reduce(CTF_OP_MAXABS);
        PRINT("aIBj1:   %e\n", diff);

        DistTensor<double> aIBj2(aomo.getAIBJ()(0,1,1,0));
        aIBj2["aIBj"] += aomo.getABIJ()(1,0,0,1)["BaIj"];
        diff = aIBj2.reduce(CTF_OP_MAXABS);
        PRINT("aIBj2:   %e\n", diff);

        /*
        CCSD<double> ccsd(config.get("cc"), moints);

        PRINT("UHF-MP2 Energy: %.15f\n", ccsd.getEnergy());

        s2 = ccsd.getProjectedS2();
        mult = ccsd.getProjectedMultiplicity();

        PRINT("\n");
        PRINT("<0|S^2|MP2>  = %f\n", s2);
        PRINT("<0|2S+1|MP2> = %f\n", mult);
        PRINT("\n");

        PRINT("\nUHF-CCSD\n\n");
        PRINT("It.   Correlation Energy     Residual Walltime\n");
        tic();
        for (int i = 0;ccsd.iterate();i++)
        {
            double dt = todouble(toc());
            PRINT("%3d % 20.15f %12.6e %8.3f\n", i+1, ccsd.getEnergy(), ccsd.getConvergence(), dt);
            tic();
        }

        s2 = ccsd.getProjectedS2();
        mult = ccsd.getProjectedMultiplicity();

        PRINT("\n");
        PRINT("<0|S^2|CC>  = %f\n", s2);
        PRINT("<0|2S+1|CC> = %f\n", mult);
        PRINT("\n");

        STTwoElectronOperator<double,2> H(moints, ccsd);
        LambdaCCSD<double> lambda(config.get("cc"), H, ccsd, ccsd.getEnergy());

        PRINT("UHF-Lambda-CCSD\n\n");
        PRINT("It.   Correlation Energy     Residual Walltime\n");
        tic();
        for (int i = 0;lambda.iterate();i++)
        {
            double dt = todouble(toc());
            PRINT("%3d % 20.15f %12.6e %8.3f\n", i+1, lambda.getEnergy(), lambda.getConvergence(), dt);
            tic();
        }

        PRINT("\nFinal Energy: %.15f\n\n", scf.getEnergy()+ccsd.getEnergy());

        TwoElectronDensity<double> Ds(scf);
        TwoElectronDensity<double> Du(ccsd);
        TwoElectronDensity<double> Dr(lambda, ccsd);

        double es = scalar(Ds*moints);
        double eu = scalar(Du*moints);
        double er = scalar(Dr*moints);

        cout << es << " " << eu << " " << er << endl;
        */

        print_timers();
    }

    elem::Finalize();
    SLIDE::finish();
    MPI_Finalize();
}
