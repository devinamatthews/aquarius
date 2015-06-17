#ifndef _AQUARIUS_DIIS_HPP_
#define _AQUARIUS_DIIS_HPP_

#include "util/global.hpp"

#include "input/config.hpp"
#include "task/task.hpp"
#include "tensor/tensor.hpp"

namespace aquarius
{
namespace convergence
{

class DIIS
{
    protected:
        vector<vector<tensor::Tensor<>>> old_x;
        vector<vector<tensor::Tensor<>>> old_dx;
        marray<tensor::Scalar,1> c;
        marray<tensor::Scalar,2> e;
        int nextrap, start, nx, ndx;
        double damping;

    public:
        DIIS(const input::Config& config, int nx = 1, int ndx = 1)
        : nx(nx), ndx(ndx)
        {
            nextrap = config.get<int>("order");
            start = config.get<int>("start");
            damping = config.get<double>("damping");

            e.resize(nextrap+1, nextrap+1);
            c.resize(nextrap+1);

            old_x.resize(nextrap);
            old_dx.resize(nextrap);
        }

        void extrapolate(tensor::Tensor<>& x, tensor::Tensor<>& dx)
        {
            extrapolate({x}, {dx});
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

            /*
             * Lazily allocate elements of old_x etc. so that we can
             * just use the copy ctor and subclasses do not have to
             * worry about allocation/deallocation
             */
            if (old_x[0].empty())
            {
                for (int i = 0;i < nx;i++)
                {
                    old_x[0].push_back(x[i]);
                }
            }
            else
            {
                for (int i = 0;i < nx;i++)
                {
                    old_x[0][i] = x[i];
                }
            }

            if (old_dx[0].empty())
            {
                for (int i = 0;i < ndx;i++)
                {
                    old_dx[0].push_back(dx[i]);
                }
            }
            else
            {
                for (int i = 0;i < ndx;i++)
                {
                    old_dx[0][i] = dx[i];
                }
            }

            /*
             * Shift the previous error matrix (the last row and column will be discarded)
             */
            for (int i = nextrap-1;i > 0;i--)
            {
                for (int j = nextrap-1;j > 0;j--)
                {
                    e[i][j] = e[i-1][j-1];
                }

                c[i] = c[i-1];
            }

            e[0][0] = 0;
            for (int i = 0;i < ndx;i++)
            {
                e[0][0] += scalar(conj(dx[i])*dx[i]);
            }

            /*
             * Get the new off-diagonal error matrix elements for all
             * previous vectors which exist. There may be fewer than nextrap of them
             * (e.g. in iterations 1 to nextrap-1), so save this number.
             */
            int nextrap_real = 1;
            for (int i = 1;i < nextrap && !old_dx[i].empty();i++)
            {
                e[i][0] = 0;
                for (int j = 0;j < ndx;j++)
                {
                    e[i][0] += scalar(conj(dx[j])*old_dx[i][j]);
                }
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

            /*
            for (int i = 0;i < nextrap_real;i++)
            {
                for (int j = 0;j < nextrap_real;j++)
                {
                    printf("%+.3e ", e[i+j*(nextrap+1)]);
                }

                printf(": %+.3e\n", c[i]);
            }
            printf("\n");
            */

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
                marray<double,2> tmp(nextrap+1,nextrap+1);
                marray<double,1> ctmp(nextrap+1);
                vector<integer> ipiv(nextrap+1);

                for (int i = 0;i < nextrap+1;i++)
                    for (int j = 0;j < nextrap+1;j++)
                        tmp[i][j] = e[i][j].to<double>();
                info = hesv('U', nextrap_real+1, 1, tmp.data(), nextrap+1, ipiv.data(),
                            ctmp.data(), nextrap+1);
                for (int i = 0;i < nextrap+1;i++)
                    c[i] = ctmp[i];

                /*
                 * Attempt to stave off singularity due to "exact" convergence
                 */
                if (info > 0)
                {
                    double eps = abs(e[nextrap_real-1][nextrap_real-1]).to<double>()*
                                 numeric_limits<double>::epsilon();

                    for (int i = 0;i < nextrap+1;i++)
                        for (int j = 0;j < nextrap+1;j++)
                            tmp[i][j] = e[i][j].to<double>() + (i == j ? eps : 0);
                    info = hesv('U', nextrap_real+1, 1, tmp.data(), nextrap+1, ipiv.data(),
                                ctmp.data(), nextrap+1);
                    for (int i = 0;i < nextrap+1;i++)
                        c[i] = ctmp[i];
                }

                if (info != 0) throw runtime_error(strprintf("DIIS: Info in hesv: %d", info));
            }

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
