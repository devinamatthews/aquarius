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
 * ARE DISCLAIMED. IN NO EVENT SHALL EDGAR SOLOMONIK BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "hamiltonian.hpp"

using namespace std;
using namespace aquarius::scf;
using namespace aquarius::autocc;

namespace aquarius
{
namespace cc
{

Hamiltonian::Hamiltonian(MOIntegrals& moints, int copy)
: moints(moints),
  fae("a,e"),
  fmi("m,i"),
  fme("m,e"),
  wmnij("mn,ij"),
  wmbij("mb,ij"),
  wmnie("mn,ie"),
  wmnef("mn,ef"),
  wmbej("mb,ej"),
  wamef("am,ef"),
  wabej("ab,ej"),
  wabef("ab,ef")
{
    if (copy&FAE)
    {
        fae.addSpinCase(new DistTensor(moints.getFAB().getSpinCase(0)), "A,E", "AE");
        fae.addSpinCase(new DistTensor(moints.getFAB().getSpinCase(1)), "a,e", "ae");
    }
    else
    {
        fae.addSpinCase(moints.getFAB().getSpinCase(0), "A,E", "AE");
        fae.addSpinCase(moints.getFAB().getSpinCase(1), "a,e", "ae");
    }

    if (copy&FMI)
    {
        fmi.addSpinCase(new DistTensor(moints.getFIJ().getSpinCase(0)), "M,I", "MI");
        fmi.addSpinCase(new DistTensor(moints.getFIJ().getSpinCase(1)), "m,i", "mi");
    }
    else
    {
        fmi.addSpinCase(moints.getFIJ().getSpinCase(0), "M,I", "MI");
        fmi.addSpinCase(moints.getFIJ().getSpinCase(1), "m,i", "mi");
    }

    if (copy&FME)
    {
        fme.addSpinCase(new DistTensor(moints.getFAI().getSpinCase(0)), "M,E", "EM");
        fme.addSpinCase(new DistTensor(moints.getFAI().getSpinCase(1)), "m,e", "em");
    }
    else
    {
        fme.addSpinCase(moints.getFAI().getSpinCase(0), "M,E", "EM");
        fme.addSpinCase(moints.getFAI().getSpinCase(1), "m,e", "em");
    }

    if (copy&WMNIJ)
    {
        wmnij.addSpinCase(new DistTensor(moints.getVIJKL().getSpinCase(0)), "MN,IJ", "MNIJ");
        wmnij.addSpinCase(new DistTensor(moints.getVIJKL().getSpinCase(1)), "Mn,Ij", "MnIj");
        wmnij.addSpinCase(new DistTensor(moints.getVIJKL().getSpinCase(2)), "mn,ij", "mnij");
    }
    else
    {
        wmnij.addSpinCase(moints.getVIJKL().getSpinCase(0), "MN,IJ", "MNIJ");
        wmnij.addSpinCase(moints.getVIJKL().getSpinCase(1), "Mn,Ij", "MnIj");
        wmnij.addSpinCase(moints.getVIJKL().getSpinCase(2), "mn,ij", "mnij");
    }

    if (copy&WMBIJ)
    {
        wmbij.addSpinCase(new DistTensor(moints.getVIJKA().getSpinCase(0)), "MB,IJ", "IJMB");
        wmbij.addSpinCase(new DistTensor(moints.getVIJKA().getSpinCase(1)), "Mb,Ij", "IjMb");
        wmbij.addSpinCase(new DistTensor(moints.getVIJKA().getSpinCase(2)), "mB,iJ", "iJmB");
        wmbij.addSpinCase(new DistTensor(moints.getVIJKA().getSpinCase(3)), "mb,ij", "ijmb");
    }
    else
    {
        wmbij.addSpinCase(moints.getVIJKA().getSpinCase(0), "MB,IJ", "IJMB");
        wmbij.addSpinCase(moints.getVIJKA().getSpinCase(1), "Mb,Ij", "IjMb");
        wmbij.addSpinCase(moints.getVIJKA().getSpinCase(2), "mB,iJ", "iJmB");
        wmbij.addSpinCase(moints.getVIJKA().getSpinCase(3), "mb,ij", "ijmb");
    }

    if (copy&WMNIE)
    {
        wmnie.addSpinCase(new DistTensor(moints.getVIJKA().getSpinCase(0)), "MN,IE", "MNIE");
        wmnie.addSpinCase(new DistTensor(moints.getVIJKA().getSpinCase(1)), "Mn,Ie", "MnIe");
        wmnie.addSpinCase(new DistTensor(moints.getVIJKA().getSpinCase(2)), "mN,iE", "mNiE");
        wmnie.addSpinCase(new DistTensor(moints.getVIJKA().getSpinCase(3)), "mn,ie", "mnie");
    }
    else
    {
        wmnie.addSpinCase(moints.getVIJKA().getSpinCase(0), "MN,IE", "MNIE");
        wmnie.addSpinCase(moints.getVIJKA().getSpinCase(1), "Mn,Ie", "MnIe");
        wmnie.addSpinCase(moints.getVIJKA().getSpinCase(2), "mN,iE", "mNiE");
        wmnie.addSpinCase(moints.getVIJKA().getSpinCase(3), "mn,ie", "mnie");
    }

    if (copy&WMNEF)
    {
        wmnef.addSpinCase(new DistTensor(moints.getVABIJ().getSpinCase(0)), "MN,EF", "EFMN");
        wmnef.addSpinCase(new DistTensor(moints.getVABIJ().getSpinCase(1)), "Mn,Ef", "EfMn");
        wmnef.addSpinCase(new DistTensor(moints.getVABIJ().getSpinCase(2)), "mn,ef", "efmn");
    }
    else
    {
        wmnef.addSpinCase(moints.getVABIJ().getSpinCase(0), "MN,EF", "EFMN");
        wmnef.addSpinCase(moints.getVABIJ().getSpinCase(1), "Mn,Ef", "EfMn");
        wmnef.addSpinCase(moints.getVABIJ().getSpinCase(2), "mn,ef", "efmn");
    }

    if (copy&WMBEJ)
    {
        wmbej.addSpinCase(new DistTensor(moints.getVAIBJ().getSpinCase(0)), "MB,EJ", "BMEJ", -1.0);
        wmbej.addSpinCase(new DistTensor(moints.getVAIBJ().getSpinCase(1)), "mB,Ej", "BmEj", -1.0);
        wmbej.addSpinCase(new DistTensor(moints.getVAIBJ().getSpinCase(2)), "Mb,eJ", "bMeJ", -1.0);
        wmbej.addSpinCase(new DistTensor(moints.getVAIBJ().getSpinCase(3)), "mb,ej", "bmej", -1.0);
        wmbej.addSpinCase(new DistTensor(moints.getVABIJ().getSpinCase(1)), "Mb,Ej", "EbMj");
        wmbej.addSpinCase(new DistTensor(moints.getVABIJ().getSpinCase(1)), "mB,eJ", "BeJm");
    }
    else
    {
        wmbej.addSpinCase(moints.getVAIBJ().getSpinCase(0), "MB,EJ", "BMEJ", -1.0);
        wmbej.addSpinCase(moints.getVAIBJ().getSpinCase(1), "mB,Ej", "BmEj", -1.0);
        wmbej.addSpinCase(moints.getVAIBJ().getSpinCase(2), "Mb,eJ", "bMeJ", -1.0);
        wmbej.addSpinCase(moints.getVAIBJ().getSpinCase(3), "mb,ej", "bmej", -1.0);
        wmbej.addSpinCase(moints.getVABIJ().getSpinCase(1), "Mb,Ej", "EbMj");
        wmbej.addSpinCase(moints.getVABIJ().getSpinCase(1), "mB,eJ", "BeJm");
    }

    if (copy&WAMEF)
    {
        wamef.addSpinCase(new DistTensor(moints.getVABCI().getSpinCase(0)), "AM,EF", "EFAM");
        wamef.addSpinCase(new DistTensor(moints.getVABCI().getSpinCase(1)), "Am,Ef", "EfAm");
        wamef.addSpinCase(new DistTensor(moints.getVABCI().getSpinCase(2)), "aM,eF", "eFaM");
        wamef.addSpinCase(new DistTensor(moints.getVABCI().getSpinCase(3)), "am,ef", "efam");
    }
    else
    {
        wamef.addSpinCase(moints.getVABCI().getSpinCase(0), "AM,EF", "EFAM");
        wamef.addSpinCase(moints.getVABCI().getSpinCase(1), "Am,Ef", "EfAm");
        wamef.addSpinCase(moints.getVABCI().getSpinCase(2), "aM,eF", "eFaM");
        wamef.addSpinCase(moints.getVABCI().getSpinCase(3), "am,ef", "efam");
    }

    if (copy&WABEJ)
    {
        wabej.addSpinCase(new DistTensor(moints.getVABCI().getSpinCase(0)), "AB,EJ", "ABEJ");
        wabej.addSpinCase(new DistTensor(moints.getVABCI().getSpinCase(1)), "Ab,Ej", "AbEj");
        wabej.addSpinCase(new DistTensor(moints.getVABCI().getSpinCase(2)), "aB,eJ", "aBeJ");
        wabej.addSpinCase(new DistTensor(moints.getVABCI().getSpinCase(3)), "ab,ej", "abej");
    }
    else
    {
        wabej.addSpinCase(moints.getVABCI().getSpinCase(0), "AB,EJ", "ABEJ");
        wabej.addSpinCase(moints.getVABCI().getSpinCase(1), "Ab,Ej", "AbEj");
        wabej.addSpinCase(moints.getVABCI().getSpinCase(2), "aB,eJ", "aBeJ");
        wabej.addSpinCase(moints.getVABCI().getSpinCase(3), "ab,ej", "abej");
    }

    if (copy&WABEF)
    {
        wabef.addSpinCase(new DistTensor(moints.getVABCD().getSpinCase(0)), "AB,EF", "ABEF");
        wabef.addSpinCase(new DistTensor(moints.getVABCD().getSpinCase(1)), "Ab,Ef", "AbEf");
        wabef.addSpinCase(new DistTensor(moints.getVABCD().getSpinCase(2)), "ab,ef", "abef");
    }
    else
    {
        wabef.addSpinCase(moints.getVABCD().getSpinCase(0), "AB,EF", "ABEF");
        wabef.addSpinCase(moints.getVABCD().getSpinCase(1), "Ab,Ef", "AbEf");
        wabef.addSpinCase(moints.getVABCD().getSpinCase(2), "ab,ef", "abef");
    }
}

}
}
