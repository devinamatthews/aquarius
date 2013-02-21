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

#ifndef _AQUARIUS_DIIS_HPP_
#define _AQUARIUS_DIIS_HPP_

#include "tensor/tensor.hpp"
#include "input/config.hpp"
#include "util/lapack.h"

#include <vector>
#include <cassert>
#include <algorithm>

namespace aquarius
{
namespace diis
{

template<typename T, typename U = T>
class DIIS
{
    protected:
        std::vector< std::vector<T*> > old_x;
        std::vector< std::vector<U*> > old_dx;
        std::vector<double> c, e;
        int nextrap, start, nx, ndx;

    public:
        DIIS(const input::Config& config, const int nx = 1, const int ndx = 1)
        : nx(nx), ndx(ndx)
        {
            nextrap = config.get<int>("order");
            start = config.get<int>("start");

            e.resize((nextrap+1)*(nextrap+1));
            c.resize(nextrap+1);

            old_x.resize(nextrap, std::vector<T*>(nx, (T*)NULL));
            old_dx.resize(nextrap, std::vector<U*>(ndx, (U*)NULL));
        }

        ~DIIS()
        {
            for (typename std::vector< std::vector<T*> >::iterator j = old_x.begin();j != old_x.end();++j)
            {
                for (typename std::vector<T*>::iterator i = j->begin();i != j->end();i++)
                {
                    if (*i != NULL) delete *i;
                }
            }

            for (typename std::vector< std::vector<U*> >::iterator j = old_dx.begin();j != old_dx.end();++j)
            {
                for (typename std::vector<U*>::iterator i = j->begin();i != j->end();++i)
                {
                    if (*i != NULL) delete *i;
                }
            }
        }

        void extrapolate(T& x, U& dx)
        {
            extrapolate(std::vector<T*>(1, &x), std::vector<U*>(1, &dx));
        }

        void extrapolate(const std::vector<T*>& x, const std::vector<U*>& dx)
        {
            assert(x.size() == nx);
            assert(dx.size() == ndx);

            if (nextrap <= 1) return;

            for (int i = 0;i < nx;i++) assert(x[i] != NULL);
            for (int i = 0;i < ndx;i++) assert(dx[i] != NULL);

            /*
             * Move things around such that in iteration n, the data from
             * iteration n-k is in slot k
             */
            std::rotate( old_x.begin(),  old_x.end()-1,  old_x.end());
            std::rotate(old_dx.begin(), old_dx.end()-1, old_dx.end());

            /*
             * Lazily allocate elements of old_x etc. so that we can
             * just use the copy ctor and subclasses do not have to
             * worry about allocation/deallocation
             */
            if (old_x[0][0] == NULL)
            {
                for (int i = 0;i < nx;i++)
                {
                    old_x[0][i] = new T(*x[i]);
                }
            }
            else
            {
                for (int i = 0;i < nx;i++)
                {
                    (*old_x[0][i]) = *x[i];
                }
            }

            if (old_dx[0][0] == NULL)
            {
                for (int i = 0;i < ndx;i++)
                {
                    old_dx[0][i] = new U(*dx[i]);
                }
            }
            else
            {
                for (int i = 0;i < ndx;i++)
                {
                    (*old_dx[0][i]) = *dx[i];
                }
            }

            /*
             * Shift the previous error matrix (the last row and column will be discarded)
             */
            for (int i = nextrap-1;i > 0;i--)
            {
                for (int j = nextrap-1;j > 0;j--)
                {
                    e[i+j*(nextrap+1)] = e[(i-1)+(j-1)*(nextrap+1)];
                }

                c[i] = c[i-1];
            }

            e[0] = 0;
            for (int i = 0;i < ndx;i++)
            {
                e[0] += tensor::scalar((*dx[i])*(*dx[i]));
            }

            /*
             * Get the new off-diagonal error matrix elements for all
             * previous vectors which exist. There may be fewer than nextrap of them
             * (e.g. in iterations 1 to nextrap-1), so save this number.
             */
            int nextrap_real = 1;
            for (int i = 1;i < nextrap && old_dx[i][0] != NULL;i++)
            {
                e[i] = 0;
                for (int j = 0;j < ndx;j++)
                {
                    e[i] += tensor::scalar((*dx[j])*(*old_dx[i][j]));
                }
                e[i*(nextrap+1)] = e[i];
                nextrap_real++;
            }

            /*
             * Set the elements corresponding to the unity
             * constraints and the solution vector
             */
            for (int i = 0;i < nextrap_real;i++)
            {
                e[i+nextrap_real*(nextrap+1)] = -1.0;
                e[nextrap_real+i*(nextrap+1)] = -1.0;
                c[i] = 0.0;
            }

            e[nextrap_real+nextrap_real*(nextrap+1)] = 0.0;
            c[nextrap_real] = -1.0;

            /*
            for (int i = 0;i <= nextrap_real;i++)
            {
                for (int j = 0;j <= nextrap_real;j++)
                {
                    printf("%+11e ", e[i+j*(nextrap+1)]);
                }

                printf(": %+11e\n", c[i]);
            }
            printf("\n");
            */

            if (nextrap_real == 1 || --start > 1) return;

            {
                int info;
                std::vector<double> tmp((nextrap+1)*(nextrap+1));
                std::vector<double> work(nextrap+1);
                std::vector<int> ipiv(nextrap+1);

                std::copy(e.begin(), e.end(), tmp.begin());
                info = dsysv('U', nextrap_real+1, 1, tmp.data(), nextrap+1, ipiv.data(),
                             c.data(), nextrap+1, work.data(), nextrap+1);
                ASSERT(info == 0, "failure in dsysv, info = %d", info);
            }

            for (int i = 0;i < ndx;i++)
            {
                *dx[i] = (*old_dx[0][i])*c[0];
            }

            for (int i = 0;i < nx;i++)
            {
                *x[i] = (*old_x[0][i])*c[0];
            }

            for (int i = 1;i < nextrap_real;i++)
            {
                for (int j = 0;j < ndx;j++)
                {
                    *dx[j] += (*old_dx[i][j])*c[i];
                }

                for (int j = 0;j < nx;j++)
                {
                    *x[j] += (*old_x[i][j])*c[i];
                }
            }
        }
};

}
}

#endif
