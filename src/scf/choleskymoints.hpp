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

#ifndef _AQUARIUS_SCF_CHOLESKYMOINTS_HPP_
#define _AQUARIUS_SCF_CHOLESKYMOINTS_HPP_

#include "choleskyscf.hpp"
#include "moints.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class CholeskyMOIntegrals : public MOIntegrals<T>
{
    public:
        CholeskyMOIntegrals(const CholeskyUHF<T>& uhf)
        : MOIntegrals<T>(uhf)
        {
            doTransformation(uhf.chol);
        }

    protected:
        void doTransformation(const CholeskyIntegrals<T>& chol)
        {
            const tensor::DistTensor<T>& cA = this->uhf.getCA();
            const tensor::DistTensor<T>& ca = this->uhf.getCa();
            const tensor::DistTensor<T>& cI = this->uhf.getCI();
            const tensor::DistTensor<T>& ci = this->uhf.getCi();
            const tensor::DistTensor<T>& Lpq = chol.getL();
            const tensor::DistTensor<T>& D = chol.getD();

            int N = this->uhf.getMolecule().getNumOrbitals();
            int nI = this->uhf.getMolecule().getNumAlphaElectrons();
            int ni = this->uhf.getMolecule().getNumBetaElectrons();
            int nA = N-nI;
            int na = N-ni;
            int R = chol.getRank();

            int sizeIIR[] = {nI, nI, R};
            int sizeiiR[] = {ni, ni, R};
            int sizeAAR[] = {nA, nA, R};
            int sizeaaR[] = {na, na, R};
            int sizeAIR[] = {nA, nI, R};
            int sizeaiR[] = {na, ni, R};

            int shapeNNN[] = {NS, NS, NS};

            tensor::DistTensor<T> LIJ(this->ctf, 3, sizeIIR, shapeNNN, false);
            tensor::DistTensor<T> Lij(this->ctf, 3, sizeiiR, shapeNNN, false);
            tensor::DistTensor<T> LAB(this->ctf, 3, sizeAAR, shapeNNN, false);
            tensor::DistTensor<T> Lab(this->ctf, 3, sizeaaR, shapeNNN, false);
            tensor::DistTensor<T> LAI(this->ctf, 3, sizeAIR, shapeNNN, false);
            tensor::DistTensor<T> Lai(this->ctf, 3, sizeaiR, shapeNNN, false);

            {
                int sizeNIR[] = { N, nI, R};
                int sizeNiR[] = { N, ni, R};
                int sizeNAR[] = { N, nA, R};
                int sizeNaR[] = { N, na, R};

                tensor::DistTensor<T> LpI(this->ctf, 3, sizeNIR, shapeNNN, false);
                tensor::DistTensor<T> Lpi(this->ctf, 3, sizeNiR, shapeNNN, false);
                tensor::DistTensor<T> LpA(this->ctf, 3, sizeNAR, shapeNNN, false);
                tensor::DistTensor<T> Lpa(this->ctf, 3, sizeNaR, shapeNNN, false);

                LpI["pIR"] = Lpq["pqR"]*cI["qI"];
                Lpi["piR"] = Lpq["pqR"]*ci["qi"];
                LpA["pAR"] = Lpq["pqR"]*cA["qA"];
                Lpa["paR"] = Lpq["pqR"]*ca["qa"];

                LIJ["IJR"] = LpI["pJR"]*cI["pI"];
                Lij["ijR"] = Lpi["pjR"]*ci["pi"];
                LAI["AIR"] = LpI["pIR"]*cA["pA"];
                Lai["aiR"] = Lpi["piR"]*ca["pa"];
                LAB["ABR"] = LpA["pBR"]*cA["pA"];
                Lab["abR"] = Lpa["pbR"]*ca["pa"];
            }

            tensor::DistTensor<T> LDIJ(this->ctf, 3, sizeIIR, shapeNNN, false);
            tensor::DistTensor<T> LDij(this->ctf, 3, sizeiiR, shapeNNN, false);
            tensor::DistTensor<T> LDAB(this->ctf, 3, sizeAAR, shapeNNN, false);
            tensor::DistTensor<T> LDab(this->ctf, 3, sizeaaR, shapeNNN, false);
            tensor::DistTensor<T> LDAI(this->ctf, 3, sizeAIR, shapeNNN, false);
            tensor::DistTensor<T> LDai(this->ctf, 3, sizeaiR, shapeNNN, false);

            LDIJ["IJR"] = D["R"]*LIJ["IJR"];
            LDij["ijR"] = D["R"]*Lij["ijR"];
            LDAI["AIR"] = D["R"]*LAI["AIR"];
            LDai["aiR"] = D["R"]*Lai["aiR"];
            LDAB["ABR"] = D["R"]*LAB["ABR"];
            LDab["abR"] = D["R"]*Lab["abR"];

            (*this->IJKL_)["IJKL"] = 0.5*LDIJ["IKR"]*LIJ["JLR"];
            (*this->IjKl_)["IjKl"] =     LDIJ["IKR"]*Lij["jlR"];
            (*this->ijkl_)["ijkl"] = 0.5*LDij["ikR"]*Lij["jlR"];

            (*this->IJKA_)["IJKA"] = LDIJ["IKR"]*LAI["AJR"];
            (*this->IjKa_)["IjKa"] = LDIJ["IKR"]*Lai["ajR"];
            (*this->iJkA_)["iJkA"] = LDij["ikR"]*LAI["AJR"];
            (*this->ijka_)["ijka"] = LDij["ikR"]*Lai["ajR"];

            (*this->ABIJ_)["ABIJ"] = 0.5*LDAI["AIR"]*LAI["BJR"];
            (*this->AbIj_)["AbIj"] =     LDAI["AIR"]*Lai["bjR"];
            (*this->abij_)["abij"] = 0.5*LDai["aiR"]*Lai["bjR"];

            (*this->AIBJ_)["AIBJ"]  = LDAB["ABR"]*LIJ["IJR"];
            (*this->AIBJ_)["AIBJ"] -= LDAI["AJR"]*LAI["BIR"];
            (*this->AiBj_)["AiBj"]  = LDAB["ABR"]*Lij["ijR"];
            (*this->aIbJ_)["aIbJ"]  = LDab["abR"]*LIJ["IJR"];
            (*this->aibj_)["aibj"]  = LDab["abR"]*Lij["ijR"];
            (*this->aibj_)["aibj"] -= LDai["ajR"]*Lai["biR"];
            (*this->AibJ_)["AibJ"] -= (*this->AbIj_)["AbJi"];

            (*this->ABCI_)["ABCI"] = LDAB["ACR"]*LAI["BIR"];
            (*this->AbCi_)["AbCi"] = LDAB["ACR"]*Lai["biR"];
            (*this->aBcI_)["aBcI"] = LDab["acR"]*LAI["BIR"];
            (*this->abci_)["abci"] = LDab["acR"]*Lai["biR"];

            (*this->ABCD_)["ABCD"] = 0.5*LDAB["ACR"]*LAB["BDR"];
            (*this->AbCd_)["AbCd"] =     LDAB["ACR"]*Lab["bdR"];
            (*this->abcd_)["abcd"] = 0.5*LDab["acR"]*Lab["bdR"];
        }
};

}
}

#endif
