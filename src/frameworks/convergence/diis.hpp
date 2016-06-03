#ifndef _AQUARIUS_FRAMEWORKS_CONVERGENCE_DIIS_HPP_
#define _AQUARIUS_FRAMEWORKS_CONVERGENCE_DIIS_HPP_

#include "frameworks/util.hpp"
#include "frameworks/task.hpp"
#include "frameworks/tensor.hpp"

namespace aquarius
{
namespace convergence
{

struct InnerProd
{
    virtual ~InnerProd() {}

    virtual double operator()(const vector<tensor::Tensor<>>& a,
                              const vector<tensor::Tensor<>>& b) const
    {
        double p = 0;
        for (int j = 0;j < a.size();j++)
        {
            p += scalar(conj(a[j])*b[j]).to<double>();
        }
        return p;
    }
};

class DIIS
{
    protected:
        vector<vector<tensor::Tensor<>>> old_x;
        vector<vector<tensor::Tensor<>>> old_dx;
        row<double> c;
        matrix<double> e;
        int nextrap, start;
        double damping;
        int nx, ndx;
        unique_ptr<InnerProd> innerProd;

    public:
        DIIS(const task::Config& config, int nx = 1, int ndx = 1, unique_ptr<InnerProd> innerProd = new InnerProd())
        : nx(nx), ndx(ndx), innerProd(move(innerProd))
        {
            nextrap = config.get<int>("order");
            start = config.get<int>("start");
            damping = config.get<double>("damping");

            e.resize(nextrap+1, nextrap+1);
            c.resize(nextrap+1);

            old_x.resize(nextrap);
            old_dx.resize(nextrap);
        }

        void extrapolate(tensor::Tensor<> x, tensor::Tensor<> dx)
        {
            extrapolate(make_vector(x), make_vector(dx));
        }

        void extrapolate(vector<tensor::Tensor<>>&& x, vector<tensor::Tensor<>>&& dx)
        {
            assert(x.size() == nx);
            assert(dx.size() == ndx);

            if (nextrap <= 1) return;

            /*
             * Move things around such that in iteration n, the data from
             * iteration n-k is in slot k
             */
            rotate( old_x.begin(),  old_x.end()-1,  old_x.end());
            rotate(old_dx.begin(), old_dx.end()-1, old_dx.end());
            e.rotate(-1, -1);
            c.rotate(-1);

            for (int i = 0;i < nx;i++)
            {
                if (i >= old_x[0].size()) old_x[0].push_back(x[i].construct());
                old_x[0][i] = x[i];
            }

            for (int i = 0;i < ndx;i++)
            {
                if (i >= old_dx[0].size()) old_dx[0].push_back(dx[i].construct());
                old_dx[0][i] = dx[i];
            }

            e[0][0] = innerProd(dx, dx);

            /*
             * Get the new off-diagonal error matrix elements for all
             * previous vectors which exist. There may be fewer than nextrap of them
             * (e.g. in iterations 1 to nextrap-1), so save this number.
             */
            int nextrap_real = 1;
            for (int i = 1;i < nextrap && !old_dx[i].empty();i++)
            {
                e[i][0] = (*innerProd)(dx, old_dx[i]);
                e[0][i] = e[i][0];
                nextrap_real++;
            }

            /*
             * Set the elements corresponding to the unity
             * constraints and the solution vector
             */
            for (int i = 0;i < nextrap_real;i++)
            {
                e[i][nextrap_real] = -1.0;
                e[nextrap_real][i] = -1.0;
                c[i] = 0.0;
            }

            e[nextrap_real][nextrap_real] = 0.0;
            c[nextrap_real] = -1.0;

            if (nextrap_real == 1) return;

            if (--start > 1)
            {
                if (damping > 0.0)
                {
                    for (int i = 0;i < nx;i++)
                    {
                        (damping-1)*old_x[0][i] += damping*old_x[1][i];
                    }
                }

                return;
            }

            {
                int info;
                matrix<double> tmp = e;
                vector<integer> ipiv(nextrap+1);

                info = hesv('U', nextrap_real+1, 1, tmp.data(), nextrap+1, ipiv.data(),
                            c.data(), nextrap+1);

                /*
                 * Attempt to stave off singularity due to "exact" convergence
                 */
                if (info > 0)
                {
                    double eps = e[nextrap_real-1][nextrap_real-1]*
                                 numeric_limits<double>::epsilon();
                    tmp = e;
                    axpy(nextrap_real, 1.0, &eps, 0, tmp.data(), nextrap+1);
                    info = hesv('U', nextrap_real+1, 1, tmp.data(), nextrap+1, ipiv.data(),
                                c.data(), nextrap+1);
                }

                if (info != 0) throw runtime_error(str("DIIS: Info in hesv: %d", info));
            }

            /*
            for (int i = 0;i < nextrap_real;i++)
            {
                for (int j = 0;j < nextrap_real;j++)
                {
                    printf("%+.3e ", e[i][j]);
                }

                printf(": %+.3e\n", c[i]);
            }
            printf("\n");
            */

            //for (int i = 0;i <= nextrap_real;i++) printf("%+11e ", c[i]); printf("\n");

            for (int i = 0;i < ndx;i++)
            {
                dx[i] = old_dx[0][i]*c[0];
            }

            for (int i = 0;i < nx;i++)
            {
                x[i] = old_x[0][i]*c[0];
            }

            for (int i = 1;i < nextrap_real;i++)
            {
                for (int j = 0;j < ndx;j++)
                {
                    dx[j] += old_dx[i][j]*c[i];
                }

                for (int j = 0;j < nx;j++)
                {
                    x[j] += old_x[i][j]*c[i];
                }
            }
        }
};

}
}

#endif
