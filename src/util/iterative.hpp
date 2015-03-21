#ifndef _AQUARIUS_UTIL_ITERATIVE_HPP_
#define _AQUARIUS_UTIL_ITERATIVE_HPP_

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
        vector<U> energy_;
        vector<double> conv_;
        double convtol;
        int iter_;
        int maxiter;
        int nsolution_;

        static ConvergenceType getConvType(const input::Config& config)
        {
            string sconv = config.get<string>("conv_type");

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
        Iterative(const string& name, input::Config& config)
        : Task(name, config),
          convtol(config.get<double>("convergence")),
          maxiter(config.get<int>("max_iterations")),
          nsolution_(0),
          convtype(getConvType(config)) {}

        virtual ~Iterative() {}

        bool run(task::TaskDAG& dag, const Arena& arena)
        {
            return run(dag, arena, 1);
        }

        bool run(task::TaskDAG& dag, const Arena& arena, int nsolution)
        {
            nsolution_ = nsolution;
            energy_.resize(nsolution);
            conv_.assign(nsolution, numeric_limits<U>::max());

            for (iter_ = 1;iter_ <= maxiter && !isConverged();iter_++)
            {
                time::Timer timer;
                timer.start();
                iterate(arena);
                timer.stop();
                double dt = timer.seconds(arena);

                int ndigit = (int)(ceil(-log10(convtol))+0.5);

                log(arena) << "Iteration " << iter_ << " took " << fixed <<
                              setprecision(3) << dt << " s" << endl;

                for (int i = 0;i < nsolution;i++)
                {
                    if (nsolution > 1)
                    {
                        log(arena) << "Iteration " << iter_ << " sol'n " << (i+1) <<
                                      " energy = " << fixed << setprecision(ndigit) << energy_[i] <<
                                      ", convergence = " << scientific << setprecision(3) << conv_[i] << endl;
                    }
                    else
                    {
                        log(arena) << "Iteration " << iter_ <<
                                      " energy = " << fixed << setprecision(ndigit) << energy_[i] <<
                                      ", convergence = " << scientific << setprecision(3) << conv_[i] << endl;
                    }
                }

            }

            if (!isConverged())
            {
                log(arena) << "Did not converge in " << maxiter << " iterations" << endl;
            }

            return true;
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
