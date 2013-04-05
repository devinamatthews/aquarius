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

#ifndef _AQUARIUS_DAVIDSON_HPP_
#define _AQUARIUS_DAVIDSON_HPP_

#include "tensor/tensor.hpp"
#include "input/config.hpp"
#include "util/lapack.h"
#include "util/util.h"

#include <vector>
#include <cassert>
#include <algorithm>
#include <cfloat>

namespace aquarius
{
namespace convergence
{

template<typename T>
class Davidson
{
    protected:
        std::vector< std::vector<T*> > old_c;
        std::vector< std::vector<T*> > old_hc;
        std::vector<double> c, e;
        int nvec, nextrap;
        int mode;
        bool lock;
        double target;

        enum {GUESS_OVERLAP, LOWEST_ENERGY, CLOSEST_ENERGY};

    public:
        Davidson(const input::Config& config, const int nvec=1)
        : nvec(nvec), mode(GUESS_OVERLAP), lock(false)
        {
            nextrap = config.get<int>("order");

            e.resize(nextrap*nextrap);
            c.resize(nextrap);

            old_c.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
            old_hc.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
        }

        ~Davidson()
        {
            for (typename std::vector< std::vector<T*> >::iterator j = old_c.begin();j != old_c.end();++j)
            {
                for (typename std::vector<T*>::iterator i = j->begin();i != j->end();i++)
                {
                    if (*i != NULL) delete *i;
                }
            }

            for (typename std::vector< std::vector<T*> >::iterator j = old_hc.begin();j != old_hc.end();++j)
            {
                for (typename std::vector<T*>::iterator i = j->begin();i != j->end();++i)
                {
                    if (*i != NULL) delete *i;
                }
            }
        }

        double extrapolate(T& c, T& hc, T& D)
        {
            return extrapolate(std::vector<T*>(1, &c), std::vector<T*>(1, &hc), std::vector<T*>(1, &D));
        }

        double extrapolate(std::vector<T*>& c, std::vector<T*>& hc, std::vector<T*>& D)
        {
            assert(c.size() == nvec);
            assert(hc.size() == nvec);

            for (int i = 0;i < nvec;i++) assert(c[i] != NULL);
            for (int i = 0;i < nvec;i++) assert(hc[i] != NULL);

            int nextrap_real;
            for (nextrap_real = 0;nextrap_real < nextrap && old_c[nextrap_real][0] != NULL;nextrap_real++);

            if (nextrap_real == nextrap)
            {
                //TODO: compact
                assert(0);
            }

            nextrap_real++;

            /*
             * Lazily allocate elements of old_c etc. so that we can
             * just use the copy ctor and subclasses do not have to
             * worry about allocation/deallocation
             */
            if (old_c[nextrap_real][0] == NULL)
            {
                for (int i = 0;i < nvec;i++)
                    old_c[nextrap_real][i] = new T(*c[i]);
            }
            else
            {
                for (int i = 0;i < nvec;i++)
                    *old_c[nextrap_real][i] = *c[i];
            }

            if (old_hc[nextrap_real][0] == NULL)
            {
                for (int i = 0;i < nvec;i++)
                    old_hc[nextrap_real][i] = new T(*hc[i]);
            }
            else
            {
                for (int i = 0;i < nvec;i++)
                    *old_hc[nextrap_real][i] = *hc[i];
            }

            e[nextrap_real+nextrap_real*nextrap] = 0;
            for (int i = 0;i < nvec;i++)
                e[nextrap_real+nextrap_real*nextrap] += scalar((*c[i])*(*hc[i]));

            /*
             * Get the new off-diagonal subspace matrix elements for all
             * previous vectors.
             */
            for (int i = 0;i < nextrap_real;i++)
            {
                e[i+nextrap_real*nextrap] = 0;
                e[nextrap_real+i*nextrap] = 0;
                for (int j = 0;j < nvec;i++)
                {
                    e[i+nextrap_real*nextrap] += scalar((*old_c[i][j])*(*hc[j]));
                    e[nextrap_real+i*nextrap] += scalar((*c[j])*(*old_hc[i][j]));
                }
            }

            /*
            for (int i = 0;i < nextrap_real;i++)
            {
                for (int j = 0;j < nextrap_real;j++)
                {
                    printf("%+11e ", e[i+j*nextrap]);
                }

                printf(": %+11e\n", c[i]);
            }
            printf("\n");
            */

            int info;
            std::vector<double> tmp(nextrap*nextrap);
            std::vector<double> vr(nextrap*nextrap);
            std::vector<double> lr(nextrap);
            std::vector<double> lc(nextrap);
            std::vector<double> work(3*nextrap);

            std::copy(e.begin(), e.end(), tmp.begin());
            info = dgeev('N', 'V', nextrap_real, tmp.data(), nextrap, lr.data(), lc.data(),
                         NULL, 1, vr.data(), nextrap, work.data(), 4*nextrap);
            ASSERT(info == 0, "failure in dgeev, info = %d", info);

            int bestev = 0;
            double mincrit = DBL_MAX;
            for (int i = 0;i < nextrap_real;i++)
            {
                double crit = 0;

                switch (mode)
                {
                    case GUESS_OVERLAP:
                        crit = 1-vr[i*nextrap_real];
                        break;
                    case LOWEST_ENERGY:
                        crit = abs(lr[i]-lr[0]);
                        break;
                    case CLOSEST_ENERGY:
                        crit = abs(lr[i]-target);
                        break;
                }

                if (crit < mincrit)
                {
                    mincrit = crit;
                    bestev = i;
                }
            }

            ASSERT(abs(lc[bestev]) < 1e-10, "Complex eigenvalue");

            for (int j = 0;j < nvec;j++)
            {
                *hc[j] = (*old_c[0][j])*vr[0+bestev*nextrap_real];
                *c[j] = (*old_hc[0][j])*vr[0+bestev*nextrap_real];

                for (int i = 1;i < nextrap_real;i++)
                {
                    *hc[j] += (*old_c[i][j])*vr[i+bestev*nextrap_real];
                    *c[j] += (*old_hc[i][j])*vr[i+bestev*nextrap_real];
                }

                *c[j] -= lr[bestev]*(*hc[j]);

                *D[j] -= lr[bestev];
                *c[j] /= *D[j];
                *D[j] += lr[bestev];
            }

            return lr[bestev];
        }
};

}
}

#endif
