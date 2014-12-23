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

#ifndef _AQUARIUS_CC_LAMBDACCSDT_HPP_
#define _AQUARIUS_CC_LAMBDACCSDT_HPP_

#include "operator/st2eoperator.hpp"
#include "operator/deexcitationoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "convergence/diis.hpp"
#include "util/iterative.hpp"
#include "task/task.hpp"

namespace aquarius
{
namespace cc
{

/*
 * Solve the left-hand coupled cluster eigenvalue equation:
 *
 *               _
 * <0|L|Phi><Phi|H    |Phi> = 0
 *                open
 *
 *       _    -T   T       T
 * where X = e  X e  = (X e )
 *                           c
 */
template <typename U>
class LambdaCCSDT : public Iterative<U>
{
    protected:
        convergence::DIIS< op::DeexcitationOperator<U,3> > diis;

    public:
        LambdaCCSDT(const std::string& name, const input::Config& config);

        void run(task::TaskDAG& dag, const Arena& arena);

        void iterate(const Arena& arena);
};

}
}

#endif
