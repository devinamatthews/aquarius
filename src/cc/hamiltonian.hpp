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

#include "tensor/dist_tensor.hpp"
#include "tensor/spinorbital.hpp"
#include "scf/moints.hpp"

namespace aquarius
{
namespace cc
{

class Hamiltonian
{
    protected:
        scf::MOIntegrals& moints;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > fae;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > fmi;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > fme;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > wmnij;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > wmbij;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > wmnie;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > wmnef;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > wmbej;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > wamef;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > wabej;
        tensor::SpinorbitalTensor< tensor::DistTensor<double> > wabef;

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

        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFAE()   { return fae; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFMI()   { return fmi; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFME()   { return fme; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMNIJ() { return wmnij; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMBIJ() { return wmbij; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMNIE() { return wmnie; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMNEF() { return wmnef; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMBEJ() { return wmbej; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWAMEF() { return wamef; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWABEJ() { return wabej; }
        tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWABEF() { return wabef; }

        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFAE()   const { return fae; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFMI()   const { return fmi; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getFME()   const { return fme; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMNIJ() const { return wmnij; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMBIJ() const { return wmbij; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMNIE() const { return wmnie; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMNEF() const { return wmnef; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWMBEJ() const { return wmbej; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWAMEF() const { return wamef; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWABEJ() const { return wabej; }
        const tensor::SpinorbitalTensor< tensor::DistTensor<double> >& getWABEF() const { return wabef; }
};

}
}

#endif
