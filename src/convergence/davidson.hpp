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

#include "util/stl_ext.hpp"
#include "input/config.hpp"
#include "util/lapack.h"
#include "util/util.h"
#include "task/task.hpp"
#include "operator/denominator.hpp"

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
        typedef typename T::dtype dtype;
        std::vector< std::vector<T*> > old_c; // hold all R[k][i] where i is over nvec and k is over nextrap
        std::vector< std::vector<T*> > old_hc; // hold all H*R[i] = Z[i] at every iteration aka Z[k][i]
        std::vector<dtype> c, e; // e will hold chc[k]
        int nvec, nextrap; // number of energies, number of iterations
        int mode;
        bool lock;
        double target0;
        double target1;

        enum {GUESS_OVERLAP, LOWEST_ENERGY, CLOSEST_ENERGY};

    public:
        Davidson(const input::Config& config, int nvec=1)
        : nvec(nvec), mode(CLOSEST_ENERGY), lock(false)
        {
            nextrap = config.get<int>("order"); // max number of iterations

            e.resize(nextrap*nextrap*nvec*nvec);
            c.resize(nvec);

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

        std::vector<double> extrapolate(T& c, T& hc, const op::Denominator<typename T::dtype>& D, double target0=0.0)
        {
            return extrapolate(std::vector<T*>(1, &c), std::vector<T*>(1, &hc), D, std::vector<double>(1,target0));
        }

        std::vector<double> extrapolate(const std::vector<T*>& c, const std::vector<T*>& hc, const op::Denominator<typename T::dtype>& D, std::vector<double> targets)
        {
            using namespace std;

            if (targets[0] == 0.0)
            {
                assert(mode != CLOSEST_ENERGY);
            }

            // std::cout << "nvec = " << nvec << std::endl;

            assert(c.size() == nvec);
            assert(hc.size() == nvec);

            // std::cout << "c = " << c << std::endl;
            // std::cout << "hc = " << hc << std::endl;

            for (int i = 0;i < nvec;i++)
            {
                assert(c[i] != NULL);
                assert(hc[i] != NULL);

                double norm = sqrt(std::abs(scalar(conj(*c[i])*(*c[i]))));

                // std::cout << "norm = " << norm << std::endl;

                *c[i] /= norm;
                *hc[i] /= norm;
            }

            int nextrap_real; // What iteration are we on
            for (nextrap_real = 0;nextrap_real < nextrap && old_c[nextrap_real][0] != NULL;nextrap_real++);

            if (nextrap_real == nextrap) // aka we've reached our maximum iteration
            {
                //TODO: compact
                assert(0);
            }

            /*
             * Lazily allocate elements of old_c etc. so that we can
             * just use the copy ctor and subclasses do not have to
             * worry about allocation/deallocation
             */

            // std::cout << "nextrap_real = " << nextrap_real << std::endl;
            // std::cout << "nextrap = " << nextrap << std::endl;

            if (old_c[nextrap_real][0] == NULL)
            {
                for (int i = 0;i < nvec;i++)
                    old_c[nextrap_real][i] = new T(*c[i]);
            }
            else // Not sure how we'd get to this else... let's put a print statement to catch if we do
            { // Maybe we'd get here if we were compacting?
                std::cout << "Got to a place I never thought we would get!" << std::endl;
                for (int i = 0;i < nvec;i++)
                    *old_c[nextrap_real][i] = *c[i];
            }

            // std::cout << "*old_c[nextrap_real] = " << old_c[nextrap_real] << std::endl;

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

            // e[nextrap_real+nextrap_real*nextrap] = 0;
            for (int i = 0;i < nvec;i++) // Set diagonals of e
            {
                // e[i+i*nextrap*nvec+nextrap_real*nvec+nextrap_real*nextrap*nvec*nvec] = 0;
                e[i+i*nextrap*nvec+nextrap_real*nvec+nextrap_real*nextrap*nvec*nvec] = scalar(conj(*c[i])*(*hc[i])); 
            }  

            // Get the new off-diagonal subspace matrix elements for the current vectors

            int min_index = nvec * nextrap_real;
            int max_index = nvec * (nextrap_real + 1);

            for (int vert = min_index;vert < max_index;vert++)
            {
                for (int horiz = vert + 1;horiz < max_index;horiz++) 
                {
                    int upper_tri_ind = vert*nextrap*nvec + horiz;
                    int lower_tri_ind = horiz*nextrap*nvec + vert;
                    e[upper_tri_ind] = scalar(conj(*c[horiz%nvec])*(*hc[vert%nvec]));
                    e[lower_tri_ind] = scalar(conj(*c[vert%nvec])*(*hc[horiz%nvec]));
                }
            }

            // Get the new off-diagonal subspace matrix elements for the past vectors

            if (nextrap_real > 0)
            {
                for (int vert = min_index;vert < max_index;vert++)
                {
                    for (int horiz = 0;horiz < min_index;horiz++) 
                    {
                        int upper_tri_ind = horiz*nextrap*nvec + vert;
                        int lower_tri_ind = vert*nextrap*nvec + horiz;
                        e[upper_tri_ind] = scalar(conj(*c[vert%nvec])*(*old_hc[horiz/nvec][horiz%nvec]));
                        e[lower_tri_ind] = scalar(conj(*old_c[horiz/nvec][horiz%nvec])*(*hc[vert%nvec]));
                    }
                }
            }

            

            // for (int i = 0;i < nextrap_real;i++)
            // {
            //     // std::cout << "Got here" << std::endl;
            //     e[i+nextrap_real*nextrap] = 0;
            //     e[nextrap_real+i*nextrap] = 0;
            //     //for (int j = 0;j < nvec;i++)
            //     for (int j = 0;j < nvec;j++) // I don't think we should be adding these together
            //     {
            //         //std::cout << "top" << std::endl;
            //         ////std::cout << "old_c " << old_c[i][j] << std::endl;
            //         //std::cout << "hc " << hc[j] << std::endl;
            //         e[i+nextrap_real*nextrap] += scalar(conj(*old_c[i][j])*(*hc[j])); // index wrong
            //         //std::cout << "middle" << std::endl;
            //         e[nextrap_real+i*nextrap] += scalar(conj(*c[j])*(*old_hc[i][j])); // index wrong
            //         //std::cout << "bot" << std::endl;
            //     }
            // }

            // std::cout << "e2 = " << e << std::endl;

            nextrap_real++; // DONT FORGET THIS GOT PLUSED.

            // DID YOU FORGET? DON'T FORGET.

            // This is because now we're constructing the next state

            /*
            for (int i = 0;i < nextrap_real;i++)
            {
                for (int j = 0;j < nextrap_real;j++)
                {
                    printf("%+11e ", e[i+j*nextrap]);
                }
                printf("\n");
            }
            printf("\n");
            */

            //std::cout << "Check 2" << std::endl;

            int info;
            std::vector<dtype> tmp(nextrap*nextrap*nvec*nvec);
            std::vector<dtype> vr(nextrap*nextrap*nvec*nvec); // Eigenvectors
            std::vector<typename std::complex_type<dtype>::type> l(nextrap*nvec); // Eigenvalues

            std::copy(e.begin(), e.end(), tmp.begin());
            info = geev('N', 'V', nextrap_real*nvec, tmp.data(), nextrap*nvec, l.data(),
                        NULL, 1, vr.data(), nextrap*nvec);
            if (info != 0) throw std::runtime_error(std::strprintf("davidson: Info in geev: %d", info));

            std::cout << "l.size() = " << l.size() << std::endl;
            for (int i = 0; i < nextrap_real*nvec; i++)
            {
                std::cout << "l[" << i << "] = " << l[i] << std::endl;
            }

            std::vector<int> bestevs;
            std::vector<double> mincrit;
            std::vector<double> crit;
            for (int i = 0; i < nvec; i++) {
                bestevs.push_back(-1);
                mincrit.push_back(DBL_MAX);
                crit.push_back(0.0);
            }
            for (int i = 0;i < nextrap_real*nvec;i++)
            {
                
                for (int j = 0; j < nvec; j++)
                    crit[j] = 0.0;

                switch (mode)
                {
                    case GUESS_OVERLAP:
                        for (int j = 0; j < nvec; j++)
                        {
                            crit[j] = 1-std::abs(vr[i*nextrap_real]);
                            if (crit[j] < mincrit[j] && std::abs(std::imag(l[i])) < 1e-9)
                            {
                                mincrit[j] = crit[j];
                                bestevs[j] = i;
                            }
                        }
                        break;
                    case LOWEST_ENERGY:
                        crit[0] = std::real(l[i]);
                        if (crit[0] < mincrit[0] && std::abs(std::imag(l[i])) < 1e-9)
                        {
                            mincrit[0] = crit[0];
                            bestevs[0] = i;
                        }
                        for (int j = 1; j < nvec; j++)
                        {
                            crit[j] = std::real(l[i]);
                            if (crit[j] < mincrit[j] && std::abs(std::imag(l[i])) < 1e-9 && crit[j] > mincrit[j-1])
                            {
                                mincrit[j] = crit[j];
                                bestevs[j] = i;
                            }
                        }
                        break;
                    case CLOSEST_ENERGY:
                        for (int j = 0; j < nvec; j++)
                        {
                            crit[j] = std::abs(std::real(l[i])-targets[j]);
                            std::cout << "crit[" << i << "," << j << "] = " << crit[j] << std::endl;
                            if (crit[j] < mincrit[j] && std::abs(std::imag(l[i])) < 1e-9)
                            {
                                std::cout << "changed!" << std::endl;
                                mincrit[j] = crit[j];
                                bestevs[j] = i;
                            }
                        }
                        break;
                }

                // std::cout << "crit = " << crit << std::endl;
                // std::cout << "mincrit = " << mincrit << std::endl;
                // std::cout << "bestev = " << bestev << std::endl;
                // std::cout << "real = " << std::real(l[i]) << std::endl;
                // std::cout << "imag = " << std::imag(l[i]) << std::endl;

            }

            assert(bestevs[0] != -1);
            assert(bestevs[0] != -1);

            std::cout << "bestevs[0] = " << bestevs[0] << std::endl;
            std::cout << "bestevs[1] = " << bestevs[1] << std::endl;
            // if (bestevs[0] == bestevs[1])
            // {
            //     bestevs[1] = bestevs[0] + 1;
            // }
            //std::cout << "Check 9" << std::endl;

            //if (std::abs(std::imag(l[bestev])) > 1e-5)
            //    throw std::runtime_error("davidson: complex eigenvalue");

            for (int j = 0;j < nvec;j++)
            {
                *hc[j] = (*old_c[0][j])*vr[0+bestevs[j]*nextrap*nvec]; // Not sure why using hc for c and vice versa
                *c[j] = (*old_hc[0][j])*vr[0+bestevs[j]*nextrap*nvec]; // original times the evec element for it

                for (int i = 1;i < nextrap_real;i++)
                {
                    *hc[j] += (*old_c[i][j])*vr[i+bestevs[j]*nextrap*nvec]; // weight each old c by its evec value
                    *c[j] += (*old_hc[i][j])*vr[i+bestevs[j]*nextrap*nvec]; // same for old_hc. making the new state state
                }

                // so now hc is V*y = x and c is A*V*y = A*x

                *c[j] -= std::real(l[bestevs[j]])*(*hc[j]); 

                // now c = A*x - mu*x = -r

                *hc[j] = *c[j]; // This is what we norm to determine convergence, which is r, makes sense.

                c[j]->weight(D, std::real(l[bestevs[j]])); // Look into weight function

                for (int i = 0;i < nextrap_real;i++)
                {
                    for (int k = 0; k < nvec; k++)
                    {
                        double olap = scalar(conj(*c[j])*(*old_c[i][k])); // orthogonalize against all previous cs
                        *c[j] -= (*old_c[i][k])*olap;
                    }
                }
                double norm = sqrt(std::abs(scalar(conj(*c[j])*(*c[j]))));
                *c[j] /= norm;
            }
            std::vector<double> myreturn(nvec);
            for (int i = 0; i < nvec; i++)
                myreturn[i] = std::real(l[bestevs[i]]);
            // if (nvec == 2)
            // {
            //     myreturn[1] = std::real(l[bestevs[1]]);
            // }
            return myreturn;
        }
};

}
}

#endif
