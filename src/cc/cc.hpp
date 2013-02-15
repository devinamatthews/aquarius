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

#ifndef _AQUARIUS_CC_HPP_
#define _AQUARIUS_CC_HPP_

#include "mpi.h"

#include "time/time.hpp"
#include "autocc/spinorbital.hpp"
#include "scf/moints.hpp"
#include "util/iterative.hpp"
#include "hamiltonian.hpp"

namespace aquarius
{
namespace cc
{

class CCD : public Iterative
{
    protected:
        scf::MOIntegrals& moints;
        autocc::SpinorbitalTensor<DistTensor> T2, E2, D2, Z2;

    public:
        CCD(const input::Config& config, scf::MOIntegrals& moints);

        void _iterate();
};

class CCSD : public Iterative
{
    protected:
        scf::MOIntegrals& moints;
        autocc::SpinorbitalTensor<DistTensor> T1, E1, D1, Z1;
        autocc::SpinorbitalTensor<DistTensor> T2, E2, D2, Z2;
        diis::DIIS< autocc::SpinorbitalTensor<DistTensor> > diis;

    public:
        CCSD(const input::Config& config, scf::MOIntegrals& moints);

        void _iterate();
};

class CCSDT : public Iterative
{
    protected:
        scf::MOIntegrals& moints;
        autocc::SpinorbitalTensor<DistTensor> T1, E1, D1, Z1;
        autocc::SpinorbitalTensor<DistTensor> T2, E2, D2, Z2;
        autocc::SpinorbitalTensor<DistTensor> T3, E3, D3, Z3;

    public:
        CCSDT(const input::Config& config, scf::MOIntegrals& moints);

        void _iterate();
};

}
}

#endif
