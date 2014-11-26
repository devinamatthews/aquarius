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
#include <vector>

#include "time/time.hpp"
#include "task/task.hpp"

#include "distributed.hpp"
#include "tensor/ctf_tensor.hpp"

namespace aquarius
{

template <typename U>
class Iterative : public task::Task
{
    public:
        enum ConvergenceType {MAX_ABS, RMSD, MAD};

    private:
        std::vector<U> energy_;
        std::vector<double> conv_;
        double convtol;
        int iter_;
        int maxiter;
        int nsolution_;

        static ConvergenceType getConvType(const input::Config& config)
        {
            std::string sconv = config.get<std::string>("conv_type");

            if (sconv == "MAXE")
            {
                return MAX_ABS;
            }
            else if (sconv == "RMSE")
            {
                return RMSD;
            }
            else if (sconv == "MAE")
            {
                return MAD;
            }

            assert(0);
            return MAX_ABS;
        }

    protected:
        const ConvergenceType convtype;

        U& energy()
        {
            assert(energy_.size() == 1);
            return energy_[0];
        }

        U& energy(int i)
        {
            assert(i >= 0 && i < energy_.size());
            return energy_[i];
        }

        double& conv()
        {
            assert(conv_.size() == 1);
            return conv_[0];
        }

        double& conv(int i)
        {
            assert(i >= 0 && i < conv_.size());
            return conv_[i];
        }

        const U& energy() const
        {
            assert(energy_.size() == 1);
            return energy_[0];
        }

        const U& energy(int i) const
        {
            assert(i >= 0 && i < energy_.size());
            return energy_[i];
        }

        const double& conv() const
        {
            assert(conv_.size() == 1);
            return conv_[0];
        }

        const double& conv(int i) const
        {
            assert(i >= 0 && i < conv_.size());
            return conv_[i];
        }

        int nsolution() const
        {
            return nsolution_;
        }

        int iter() const
        {
            return iter_;
        }

        virtual void iterate(const Arena& arena) = 0;

    public:
        Iterative(const std::string& type, const std::string& name, const input::Config& config)
        : Task(type, name),
          convtol(config.get<double>("convergence")),
          maxiter(config.get<int>("max_iterations")),
          nsolution_(0),
          convtype(getConvType(config)) {}

        virtual ~Iterative() {}

        void run(task::TaskDAG& dag, const Arena& arena, int nsolution = 1)
        {
            nsolution_ = nsolution;
            energy_.resize(nsolution);
            conv_.resize(nsolution, std::numeric_limits<U>::max());

            for (iter_ = 1;iter_ <= maxiter && !isConverged();iter_++)
            {
                time::Timer timer;
                timer.start();
                iterate(arena);
                timer.stop();
                double dt = timer.seconds(arena);

                int ndigit = (int)(ceil(-log10(convtol))+0.5);

                log(arena) << "Iteration " << iter_ << " took " << std::fixed <<
                              std::setprecision(3) << dt << " s" << std::endl;

                for (int i = 0;i < nsolution;i++)
                {
                    if (nsolution > 1)
                    {
                        log(arena) << "Iteration " << iter_ << " sol'n " << (i+1) <<
                                      " energy = " << std::fixed << std::setprecision(ndigit) << energy_[i] <<
                                      ", convergence = " << std::scientific << std::setprecision(3) << conv_[i] << std::endl;
                    }
                    else if (!isConverged(i))
                    {
                        log(arena) << "Iteration " << iter_ <<
                                      " energy = " << std::fixed << std::setprecision(ndigit) << energy_[i] <<
                                      ", convergence = " << std::scientific << std::setprecision(3) << conv_[i] << std::endl;
                    }
                }

            }

            if (!isConverged())
            {
                log(arena) << "Did not converge in " << maxiter << " iterations" << std::endl;
            }
        }

        double getConvergence() const
        {
            return conv();
        }

        double getConvergence(int i) const
        {
            return conv(i);
        }

        bool isConverged() const
        {
            bool converged = true;

            for (int i = 0;i < conv_.size();i++)
            {
                if (conv_[i] >= convtol) converged = false;
            }

            return converged;
        }

        bool isConverged(int i) const
        {
            return conv(i) < convtol;
        }
};

}

#endif
