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

#ifndef _AQUARIUS_SCF_MOINTS_HPP_
#define _AQUARIUS_SCF_MOINTS_HPP_

#include "autocc/autocc.hpp"
#include "scf.hpp"

namespace aquarius
{
namespace scf
{

class MOIntegrals
{
    protected:
        const scf::UHF& uhf;
        autocc::SpinorbitalTensor<DistTensor> fab;
        autocc::SpinorbitalTensor<DistTensor> fai;
        autocc::SpinorbitalTensor<DistTensor> fij;
        autocc::SpinorbitalTensor<DistTensor> ijkl;
        autocc::SpinorbitalTensor<DistTensor> ijka;
        autocc::SpinorbitalTensor<DistTensor> abij;
        autocc::SpinorbitalTensor<DistTensor> aibj;
        autocc::SpinorbitalTensor<DistTensor> abci;
        autocc::SpinorbitalTensor<DistTensor> abcd;
        DistTensor *fAB_;
        DistTensor *fab_;
        DistTensor *fAI_;
        DistTensor *fai_;
        DistTensor *fIJ_;
        DistTensor *fij_;
        DistTensor *IJKL_;
        DistTensor *IjKl_;
        DistTensor *ijkl_;
        DistTensor *IJKA_;
        DistTensor *IjKa_;
        DistTensor *iJkA_;
        DistTensor *ijka_;
        DistTensor *ABIJ_;
        DistTensor *AbIj_;
        DistTensor *abij_;
        DistTensor *AIBJ_;
        DistTensor *AiBj_;
        DistTensor *aIbJ_;
        DistTensor *aibj_;
        DistTensor *ABCI_;
        DistTensor *AbCi_;
        DistTensor *aBcI_;
        DistTensor *abci_;
        DistTensor *ABCD_;
        DistTensor *AbCd_;
        DistTensor *abcd_;
        DistWorld *dw;
        MPI::Intracomm comm;

        MOIntegrals(DistWorld *dw, const scf::UHF& uhf);

    public:
        const scf::UHF& getSCF() const { return uhf; }

        autocc::SpinorbitalTensor<DistTensor>& getFAB()   { return  fab; }
        autocc::SpinorbitalTensor<DistTensor>& getFAI()   { return  fai; }
        autocc::SpinorbitalTensor<DistTensor>& getFIJ()   { return  fij; }
        autocc::SpinorbitalTensor<DistTensor>& getVIJKL() { return ijkl; }
        autocc::SpinorbitalTensor<DistTensor>& getVIJKA() { return ijka; }
        autocc::SpinorbitalTensor<DistTensor>& getVABIJ() { return abij; }
        autocc::SpinorbitalTensor<DistTensor>& getVAIBJ() { return aibj; }
        autocc::SpinorbitalTensor<DistTensor>& getVABCI() { return abci; }
        autocc::SpinorbitalTensor<DistTensor>& getVABCD() { return abcd; }

        const autocc::SpinorbitalTensor<DistTensor>& getFAB() const   { return  fab; }
        const autocc::SpinorbitalTensor<DistTensor>& getFAI() const   { return  fai; }
        const autocc::SpinorbitalTensor<DistTensor>& getFIJ() const   { return  fij; }
        const autocc::SpinorbitalTensor<DistTensor>& getVIJKL() const { return ijkl; }
        const autocc::SpinorbitalTensor<DistTensor>& getVIJKA() const { return ijka; }
        const autocc::SpinorbitalTensor<DistTensor>& getVABIJ() const { return abij; }
        const autocc::SpinorbitalTensor<DistTensor>& getVAIBJ() const { return aibj; }
        const autocc::SpinorbitalTensor<DistTensor>& getVABCI() const { return abci; }
        const autocc::SpinorbitalTensor<DistTensor>& getVABCD() const { return abcd; }
};

}
}

#endif
