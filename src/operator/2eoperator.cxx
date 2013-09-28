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

#include "2eoperator.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::tensor;

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(const Arena& arena, const Space& occ, const Space& vrt)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(arena, occ, vrt),
  ijkl(this->addTensor(new SpinorbitalTensor<T>(arena, occ.group, std::vec(vrt, occ), std::vec(0,2), std::vec(0,2)))),
  aijk(this->addTensor(new SpinorbitalTensor<T>(arena, occ.group, std::vec(vrt, occ), std::vec(1,1), std::vec(0,2)))),
  ijak(this->addTensor(new SpinorbitalTensor<T>(arena, occ.group, std::vec(vrt, occ), std::vec(0,2), std::vec(1,1)))),
  abij(this->addTensor(new SpinorbitalTensor<T>(arena, occ.group, std::vec(vrt, occ), std::vec(2,0), std::vec(0,2)))),
  ijab(this->addTensor(new SpinorbitalTensor<T>(arena, occ.group, std::vec(vrt, occ), std::vec(0,2), std::vec(2,0)))),
  aibj(this->addTensor(new SpinorbitalTensor<T>(arena, occ.group, std::vec(vrt, occ), std::vec(1,1), std::vec(1,1)))),
  aibc(this->addTensor(new SpinorbitalTensor<T>(arena, occ.group, std::vec(vrt, occ), std::vec(1,1), std::vec(2,0)))),
  abci(this->addTensor(new SpinorbitalTensor<T>(arena, occ.group, std::vec(vrt, occ), std::vec(2,0), std::vec(1,1)))),
  abcd(this->addTensor(new SpinorbitalTensor<T>(arena, occ.group, std::vec(vrt, occ), std::vec(2,0), std::vec(2,0)))) {}

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(OneElectronOperator<T>& other, int copy)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(other, copy),
  ijkl(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(0,2), std::vec(0,2)))),
  aijk(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(1,1), std::vec(0,2)))),
  ijak(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(0,2), std::vec(1,1)))),
  abij(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(2,0), std::vec(0,2)))),
  ijab(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(0,2), std::vec(2,0)))),
  aibj(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(1,1), std::vec(1,1)))),
  aibc(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(1,1), std::vec(2,0)))),
  abci(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(2,0), std::vec(1,1)))),
  abcd(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(2,0), std::vec(2,0)))) {}

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(const OneElectronOperator<T>& other)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(other),
  ijkl(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(0,2), std::vec(0,2)))),
  aijk(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(1,1), std::vec(0,2)))),
  ijak(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(0,2), std::vec(1,1)))),
  abij(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(2,0), std::vec(0,2)))),
  ijab(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(0,2), std::vec(2,0)))),
  aibj(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(1,1), std::vec(1,1)))),
  aibc(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(1,1), std::vec(2,0)))),
  abci(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(2,0), std::vec(1,1)))),
  abcd(this->addTensor(new SpinorbitalTensor<T>(other.arena, other.occ.group, std::vec(other.vrt, other.occ), std::vec(2,0), std::vec(2,0)))) {}

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(TwoElectronOperator<T>& other, int copy)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(other, copy),
  ijkl(copy&IJKL ? this->addTensor(new SpinorbitalTensor<T>(other.getIJKL())) : this->addTensor(other.getIJKL())),
  aijk(copy&AIJK ? this->addTensor(new SpinorbitalTensor<T>(other.getAIJK())) : this->addTensor(other.getAIJK())),
  ijak(copy&IJAK ? this->addTensor(new SpinorbitalTensor<T>(other.getIJAK())) : this->addTensor(other.getIJAK())),
  abij(copy&ABIJ ? this->addTensor(new SpinorbitalTensor<T>(other.getABIJ())) : this->addTensor(other.getABIJ())),
  ijab(copy&IJAB ? this->addTensor(new SpinorbitalTensor<T>(other.getIJAB())) : this->addTensor(other.getIJAB())),
  aibj(copy&AIBJ ? this->addTensor(new SpinorbitalTensor<T>(other.getAIBJ())) : this->addTensor(other.getAIBJ())),
  aibc(copy&AIBC ? this->addTensor(new SpinorbitalTensor<T>(other.getAIBC())) : this->addTensor(other.getAIBC())),
  abci(copy&ABCI ? this->addTensor(new SpinorbitalTensor<T>(other.getABCI())) : this->addTensor(other.getABCI())),
  abcd(copy&ABCD ? this->addTensor(new SpinorbitalTensor<T>(other.getABCD())) : this->addTensor(other.getABCD())) {}

template <typename T>
TwoElectronOperator<T>::TwoElectronOperator(const TwoElectronOperator<T>& other)
: OneElectronOperatorBase<T,TwoElectronOperator<T> >(other),
  ijkl(this->addTensor(new SpinorbitalTensor<T>(other.getIJKL()))),
  aijk(this->addTensor(new SpinorbitalTensor<T>(other.getAIJK()))),
  ijak(this->addTensor(new SpinorbitalTensor<T>(other.getIJAK()))),
  abij(this->addTensor(new SpinorbitalTensor<T>(other.getABIJ()))),
  ijab(this->addTensor(new SpinorbitalTensor<T>(other.getIJAB()))),
  aibj(this->addTensor(new SpinorbitalTensor<T>(other.getAIBJ()))),
  aibc(this->addTensor(new SpinorbitalTensor<T>(other.getAIBC()))),
  abci(this->addTensor(new SpinorbitalTensor<T>(other.getABCI()))),
  abcd(this->addTensor(new SpinorbitalTensor<T>(other.getABCD()))) {}

template <typename T>
T TwoElectronOperator<T>::dot(bool conja, const TwoElectronOperator<T>& A, bool conjb) const
{
    T sum = 0;

    sum += this->ab.dot(conja, A.ab, conjb);
    sum += this->ai.dot(conja, A.ai, conjb);
    sum += this->ia.dot(conja, A.ia, conjb);
    sum += this->ij.dot(conja, A.ij, conjb);

    sum += (T)0.25*ijkl.dot(conja, A.ijkl, conjb);
    sum += (T)0.25*abcd.dot(conja, A.abcd, conjb);
    sum += (T)0.25*abij.dot(conja, A.abij, conjb);
    sum += (T)0.25*ijab.dot(conja, A.ijab, conjb);
    sum +=  (T)0.5*abci.dot(conja, A.abci, conjb);
    sum +=  (T)0.5*aibc.dot(conja, A.aibc, conjb);
    sum +=  (T)0.5*ijak.dot(conja, A.ijak, conjb);
    sum +=  (T)0.5*aijk.dot(conja, A.aijk, conjb);
    sum +=         aibj.dot(conja, A.aibj, conjb);

    return sum;
}

INSTANTIATE_SPECIALIZATIONS(TwoElectronOperator);
