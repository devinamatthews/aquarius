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

#ifndef _AQUARIUS_CC_HAMILTONIAN_HPP_
#define _AQUARIUS_CC_HAMILTONIAN_HPP_

#include "tensor.hpp"
#include "dist_tensor.hpp"

#include "autocc/spinorbital.hpp"
#include "scf/moints.hpp"

namespace aquarius
{
namespace cc
{

class Hamiltonian
{
    protected:
        scf::MOIntegrals& moints;
        autocc::SpinorbitalTensor<DistTensor> fae;
        autocc::SpinorbitalTensor<DistTensor> fmi;
        autocc::SpinorbitalTensor<DistTensor> fme;
        autocc::SpinorbitalTensor<DistTensor> wmnij;
        autocc::SpinorbitalTensor<DistTensor> wmbij;
        autocc::SpinorbitalTensor<DistTensor> wmnie;
        autocc::SpinorbitalTensor<DistTensor> wmnef;
        autocc::SpinorbitalTensor<DistTensor> wmbej;
        autocc::SpinorbitalTensor<DistTensor> wamef;
        autocc::SpinorbitalTensor<DistTensor> wabej;
        autocc::SpinorbitalTensor<DistTensor> wabef;

    public:
        enum
        {
            NONE  = 0x000,
            FAE   = 0x001,
            FMI   = 0x002,
            FME   = 0x004,
            WMNIJ = 0x008,
            WMBIJ = 0x010,
            WMNIE = 0x020,
            WMNEF = 0x040,
            WMBEJ = 0x080,
            WAMEF = 0x100,
            WABEJ = 0x200,
            WABEF = 0x400
        };

        Hamiltonian(scf::MOIntegrals& moints, int copy = NONE);

        autocc::SpinorbitalTensor<DistTensor>& getFAE()   { return fae; }
        autocc::SpinorbitalTensor<DistTensor>& getFMI()   { return fmi; }
        autocc::SpinorbitalTensor<DistTensor>& getFME()   { return fme; }
        autocc::SpinorbitalTensor<DistTensor>& getWMNIJ() { return wmnij; }
        autocc::SpinorbitalTensor<DistTensor>& getWMBIJ() { return wmbij; }
        autocc::SpinorbitalTensor<DistTensor>& getWMNIE() { return wmnie; }
        autocc::SpinorbitalTensor<DistTensor>& getWMNEF() { return wmnef; }
        autocc::SpinorbitalTensor<DistTensor>& getWMBEJ() { return wmbej; }
        autocc::SpinorbitalTensor<DistTensor>& getWAMEF() { return wamef; }
        autocc::SpinorbitalTensor<DistTensor>& getWABEJ() { return wabej; }
        autocc::SpinorbitalTensor<DistTensor>& getWABEF() { return wabef; }

        const autocc::SpinorbitalTensor<DistTensor>& getFAE()   const { return fae; }
        const autocc::SpinorbitalTensor<DistTensor>& getFMI()   const { return fmi; }
        const autocc::SpinorbitalTensor<DistTensor>& getFME()   const { return fme; }
        const autocc::SpinorbitalTensor<DistTensor>& getWMNIJ() const { return wmnij; }
        const autocc::SpinorbitalTensor<DistTensor>& getWMBIJ() const { return wmbij; }
        const autocc::SpinorbitalTensor<DistTensor>& getWMNIE() const { return wmnie; }
        const autocc::SpinorbitalTensor<DistTensor>& getWMNEF() const { return wmnef; }
        const autocc::SpinorbitalTensor<DistTensor>& getWMBEJ() const { return wmbej; }
        const autocc::SpinorbitalTensor<DistTensor>& getWAMEF() const { return wamef; }
        const autocc::SpinorbitalTensor<DistTensor>& getWABEJ() const { return wabej; }
        const autocc::SpinorbitalTensor<DistTensor>& getWABEF() const { return wabef; }
};

}
}

#endif
