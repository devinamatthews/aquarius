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

#include "util/distributed.hpp"
#include "tensor/dist_tensor.hpp"
#include "tensor/spinorbital.hpp"
#include "scf.hpp"

namespace aquarius
{
namespace scf
{

class MOIntegrals : public Distributed<double>
{
    protected:
        const scf::UHF& uhf;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > fab;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > fai;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > fij;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > ijkl;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > ijka;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > abij;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > aibj;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > abci;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > abcd;
        tensor::DistTensor<double> *fAB_;
        tensor::DistTensor<double> *fab_;
        tensor::DistTensor<double> *fAI_;
        tensor::DistTensor<double> *fai_;
        tensor::DistTensor<double> *fIJ_;
        tensor::DistTensor<double> *fij_;
        tensor::DistTensor<double> *IJKL_;
        tensor::DistTensor<double> *IjKl_;
        tensor::DistTensor<double> *ijkl_;
        tensor::DistTensor<double> *IJKA_;
        tensor::DistTensor<double> *IjKa_;
        tensor::DistTensor<double> *iJkA_;
        tensor::DistTensor<double> *ijka_;
        tensor::DistTensor<double> *ABIJ_;
        tensor::DistTensor<double> *AbIj_;
        tensor::DistTensor<double> *abij_;
        tensor::DistTensor<double> *AIBJ_;
        tensor::DistTensor<double> *AiBj_;
        tensor::DistTensor<double> *aIbJ_;
        tensor::DistTensor<double> *aibj_;
        tensor::DistTensor<double> *ABCI_;
        tensor::DistTensor<double> *AbCi_;
        tensor::DistTensor<double> *aBcI_;
        tensor::DistTensor<double> *abci_;
        tensor::DistTensor<double> *ABCD_;
        tensor::DistTensor<double> *AbCd_;
        tensor::DistTensor<double> *abcd_;

        MOIntegrals(const scf::UHF& uhf);

    public:
        const scf::UHF& getSCF() const { return uhf; }

        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFAB()   { return  fab; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFAI()   { return  fai; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFIJ()   { return  fij; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVIJKL() { return ijkl; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVIJKA() { return ijka; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVABIJ() { return abij; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVAIBJ() { return aibj; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVABCI() { return abci; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVABCD() { return abcd; }

        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFAB() const   { return  fab; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFAI() const   { return  fai; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFIJ() const   { return  fij; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVIJKL() const { return ijkl; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVIJKA() const { return ijka; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVABIJ() const { return abij; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVAIBJ() const { return aibj; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVABCI() const { return abci; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getVABCD() const { return abcd; }
};

}
}

#endif
