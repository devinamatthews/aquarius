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

#ifndef _AQUARIUS_CC_PERTURBEDLAMBDACCSD_HPP_
#define _AQUARIUS_CC_PERTURBEDLAMBDACCSD_HPP_

#include "operator/2eoperator.hpp"
#include "operator/deexcitationoperator.hpp"
#include "operator/exponentialoperator.hpp"

namespace aquarius
{
namespace cc
{

/*
 * Solve the frequency-dependent left-hand eigenfunction response equations:
 *
 *     A             _                    _A
 * <0|L (w)|Phi><Phi|H  + w|Phi> = - <0|L H    |Phi> = - <0|N|Phi>
 *                    open                 open
 *
 *       _    -T   T       T       _A   _    _  A
 * where X = e  X e  = (X e )  and H  = A + (H T )
 *                           c                    c
 *
 * As in LambdaCCSD, the full left-hand eigenfunction is used instead of lambda.
 */
template <typename U>
class PerturbedLambdaCCSD : public Iterative, public op::DeexcitationOperator<U,2>
{
    protected:
        op::DeexcitationOperator<U,2>& LA;
        op::DeexcitationOperator<U,2> Z, D, N;
        op::STTwoElectronOperator<U,2> H;
        U omega;
        convergence::DIIS< op::DeexcitationOperator<U,2> > diis;

    public:
        PerturbedLambdaCCSD(const input::Config& config, const op::STTwoElectronOperator<U,2>& H,
                            const op::DeexcitationOperator<U,2>& L,
                            const op::PerturbedSTTwoElectronOperator<U,2>& A, const U omega=0)
        : Iterative(config), op::DeexcitationOperator<U,2>(H.getSCF()), LA(*this),
          Z(this->uhf), D(this->uhf), N(this->uhf),
          H(const_cast<op::TwoElectronOperator<U>&>(H),
            op::TwoElectronOperator<U>::AB|op::TwoElectronOperator<U>::IJ),
          omega(omega), diis(config.get("diis"))
        {
            D(0) = 1;
            D(1)["ia"]  = H.getIJ()["ii"];
            D(1)["ia"] -= H.getAB()["aa"];
            D(2)["ijab"]  = H.getIJ()["ii"];
            D(2)["ijab"] += H.getIJ()["jj"];
            D(2)["ijab"] -= H.getAB()["aa"];
            D(2)["ijab"] -= H.getAB()["bb"];

            D -= omega;
            D = 1/D;

            N(1) = 0;
            N(2) = 0;
            A.contract(L, N);
            N(0) = 0;

            LA = N*D;
        }

        void _iterate()
        {
            Z = N;
            H.contract(L, Z);

            energy = Ecc + tensor::scalar(Z*conj(L))/tensor::scalar(L*conj(L));

             Z *= D;
            //Z -= LA;
            LA += Z;

            conv =               Z(1)(0).reduce(CTF_OP_MAXABS);
            conv = std::max(conv,Z(1)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(0).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(1).reduce(CTF_OP_MAXABS));
            conv = std::max(conv,Z(2)(2).reduce(CTF_OP_MAXABS));

            diis.extrapolate(L, Z);
        }
};

}
}

#endif
