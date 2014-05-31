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
#include <string>

#include "time/time.hpp"
#include "task/task.hpp"

#include "distributed.hpp"

namespace aquarius
{

class Iterative : public task::Task
{
    public:
        enum ConvergenceType {MAX_ABS, RMSD, MAD};

    protected:
        double energy;
        double conv;
        double convtol;
        ConvergenceType convtype;
        int iter;
        int maxiter;

        virtual void iterate(const Arena& arena) = 0;

    public:
        Iterative(const std::string& type, const std::string& name, const input::Config& config)
        : Task(type, name),
          energy(0),
          conv(std::numeric_limits<double>::infinity()),
          convtol(config.get<double>("convergence")),
          iter(0),
          maxiter(config.get<int>("max_iterations"))
        {
            std::string sconv = config.get<std::string>("conv_type");

            if (sconv == "MAXE")
            {
                convtype = MAX_ABS;
            }
            else if (sconv == "RMSE")
            {
                convtype = RMSD;
            }
            else if (sconv == "MAE")
            {
                convtype = MAD;
            }
        }

        virtual ~Iterative() {}

        void run(task::TaskDAG& dag, const Arena& arena)
        {
            for (iter = 1;iter <= maxiter && !isConverged();iter++)
            {
                time::Timer timer;
                timer.start();
                iterate(arena);
                timer.stop();
                double dt = timer.seconds(arena);

                int ndigit = (int)(ceil(-log10(convtol))+0.5);

                log(arena) << "Iteration " << iter << " took " << std::fixed <<
                              std::setprecision(3) << dt << " s" << std::endl;
                log(arena) << "Iteration " << iter <<
                              " energy = " << std::fixed << std::setprecision(ndigit) << energy <<
                              ", convergence = " << std::scientific << std::setprecision(3) << conv << std::endl;
            }

            if (!isConverged())
            {
                log(arena) << "Did not converge in " << maxiter << " iterations" << std::endl;
            }
        }

        double getEnergy() const { return energy; }

        double getConvergence() const { return conv; }

        bool isConverged() const { return conv < convtol; }
};

class NonIterative : public task::Task
{
    protected:
        double energy;

    public:
        NonIterative(const std::string& type, const std::string& name, const input::Config& config)
        : Task(type, name),
          energy(0)
        { }

        virtual ~NonIterative() { }

        void run(task::TaskDAG& dag, const Arena& arena) { }

        double getEnergy() const { return energy; }

};

}

#endif
