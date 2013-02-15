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

#ifndef _AQUARIUS_UTIL_ITERATIVE_HPP_
#define _AQUARIUS_UTIL_ITERATIVE_HPP_

#include <limits>

#include "input/config.hpp"

namespace aquarius
{

class Iterative
{
    protected:
        double energy;
        double conv;
        double convtol;
        int iter;
        int maxiter;

        virtual void _iterate() = 0;

    public:
        Iterative(const input::Config& config)
        : energy(0),
          conv(std::numeric_limits<double>::infinity()),
          convtol(config.get<double>("convergence")),
          iter(0),
          maxiter(config.get<int>("max_iterations")) {}

        virtual ~Iterative() {}

        bool iterate()
        {
            iter++;
            if (iter > maxiter || isConverged()) return false;
            _iterate();
            return true;
        }

        double getEnergy() const { return energy; }

        double getConvergence() const { return conv; }

        bool isConverged() const { return conv < convtol; }
};

}

#endif
