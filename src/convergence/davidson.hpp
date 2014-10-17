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
#include <limits>

namespace aquarius
{
namespace convergence
{

template<typename T>
class Davidson : public task::Destructible
{
    private:
        Davidson(const Davidson& other);

        Davidson& operator=(const Davidson& other);

    protected:
        typedef typename T::dtype dtype;
        std::vector< std::vector<T*> > old_c; // hold all R[k][i] where i is over nvec and k is over nextrap
        std::vector< std::vector<T*> > old_hc; // hold all H*R[i] = Z[i] at every iteration aka Z[k][i]
        std::vector<T*> guess;
        std::matrix<dtype> overlap;
        std::tensor<dtype,4> e; // e will hold chc[k]
        std::vector<dtype> c;
        int nvec, nextrap; // number of energies, number of iterations
        int mode;
        std::vector<dtype> target;

        enum {GUESS_OVERLAP, LOWEST_ENERGY, CLOSEST_ENERGY};

    public:
        Davidson(const input::Config& config, int nvec=1)
        : nvec(nvec), mode(LOWEST_ENERGY), target(nvec)
        {
            nextrap = config.get<int>("order"); // max number of iterations

            assert(nvec > 0);
            assert(nextrap > 0);

            overlap.resize(nvec, nextrap);
            e.resize(std::vec(nvec, nextrap, nvec, nextrap));
            c.resize(nextrap*nvec);

            old_c.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
            old_hc.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
        }

        Davidson(const input::Config& config, dtype target)
        : nvec(1), mode(CLOSEST_ENERGY), target(1,target)
        {
            nextrap = config.get<int>("order"); // max number of iterations

            assert(nvec > 0);
            assert(nextrap > 0);

            overlap.resize(nvec, nextrap);
            e.resize(std::vec(nvec, nextrap, nvec, nextrap));
            c.resize(nextrap*nvec);

            old_c.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
            old_hc.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
        }

        Davidson(const input::Config& config, const std::vector<dtype>& target)
        : nvec(target.size()), mode(CLOSEST_ENERGY), target(target)
        {
            nextrap = config.get<int>("order"); // max number of iterations

            assert(nvec > 0);
            assert(nextrap > 0);

            overlap.resize(nvec, nextrap);
            e.resize(std::vec(nvec, nextrap, nvec, nextrap));
            c.resize(nextrap*nvec);

            old_c.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
            old_hc.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
        }

        Davidson(const input::Config& config, const T& guess)
        : nvec(1), mode(GUESS_OVERLAP), target(1)
        {
            nextrap = config.get<int>("order"); // max number of iterations

            assert(nvec > 0);
            assert(nextrap > 0);

            overlap.resize(nvec, nextrap);
            e.resize(std::vec(nvec, nextrap, nvec, nextrap));
            c.resize(nextrap*nvec);

            old_c.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
            old_hc.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));

            guess.resize(1, new T(guess));
        }

        Davidson(const input::Config& config, const std::vector<T*>& guess)
        : nvec(guess.size()), mode(GUESS_OVERLAP), target(guess.size())
        {
            nextrap = config.get<int>("order"); // max number of iterations

            overlap.resize(nvec, nextrap);
            assert(nvec > 0);
            assert(nextrap > 0);

            e.resize(std::vec(nvec, nextrap, nvec, nextrap));
            c.resize(nextrap*nvec);

            old_c.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));
            old_hc.resize(nextrap, std::vector<T*>(nvec, (T*)NULL));

            guess.resize(nvec);
            for (int i = 0;i < nvec;i++) guess[i] = new T(*guess[i]);
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

            for (typename std::vector<T*>::iterator i = guess.begin();i != guess.end();++i)
            {
                if (*i != NULL) delete *i;
            }
        }

        dtype extrapolate(T& c, T& hc, const op::Denominator<dtype>& D)
        {
            assert(nvec == 1);
            return extrapolate(std::vec(&c), std::vec(&hc), D)[0];
        }

        std::vector<dtype> extrapolate(const std::vector<T*>& c, const std::vector<T*>& hc, const op::Denominator<dtype>& D)
        {
            using namespace std;

            assert(nvec == c.size() && nvec == hc.size());

            // std::cout << "nvec = " << nvec << std::endl;

            assert(c.size() == nvec);
            assert(hc.size() == nvec);

            // std::cout << "c = " << c << std::endl;
            // std::cout << "hc = " << hc << std::endl;

            /*
             * Check and normalize incoming c and H*c vectors
             */
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
            else
            {
                // should only happen after compaction
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
                // should only happen after compaction
                for (int i = 0;i < nvec;i++)
                    *old_hc[nextrap_real][i] = *hc[i];
            }

            /*
             * Augment the subspace matrix with the new vectors
             */
            for (int i = 0;i < nvec;i++)
            {
                /*
                 * Compute and save overlap with guess for later
                 */
                if (mode == GUESS_OVERLAP)
                    overlap[i][nextrap_real] = std::abs(scalar(conj(*c[i])*(*guess[i])));

                for (int j = 0;j < nvec;j++)
                {
                    // "Diagonal"
                    e[i][nextrap_real][j][nextrap_real] = scalar(conj(*c[i])*(*hc[j]));

                    // "Off-diagonal"
                    for (int k = 0;k < nextrap_real;k++)
                    {
                        e[i][nextrap_real][j][k] = scalar(conj(*c[i])*(*old_hc[k][j]));
                        e[i][k][j][nextrap_real] = scalar(conj(*old_c[k][i])*(*hc[j]));
                    }
                }
            }

            /*
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
            */



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

            /*
             * Diagonalize the subspace matrix to obtain approximate solutions
             */
            int info;
            std::tensor<dtype,4> tmp(std::vec(nvec,nextrap_real,nvec,nextrap_real));
            std::tensor<dtype,3> vr(std::vec(nvec,nextrap_real,nvec*nextrap_real)); // Eigenvectors
            std::vector<typename std::complex_type<dtype>::type> l(nextrap_real*nvec); // Eigenvalues

            for (int m = 0;m < nextrap_real;m++)
                for (int k = 0;k < nvec;k++)
                    for (int j = 0;j < nextrap_real;j++)
                        for (int i = 0;i < nvec;i++)
                            tmp[i][j][k][m] = e[i][j][k][m];

            info = geev('N', 'V', nextrap_real*nvec, tmp.data(), nextrap_real*nvec, l.data(),
                        NULL, 1, vr.data(), nextrap_real*nvec);
            if (info != 0) throw std::runtime_error(std::strprintf("davidson: Info in geev: %d", info));

            std::cout << "l.size() = " << l.size() << std::endl;
            for (int i = 0; i < nextrap_real*nvec; i++)
            {
                std::cout << "l[" << i << "] = " << l[i] << std::endl;
            }

            /*
             * Assign eigenvalues (exclusively) to states by the selected criterion
             */
            std::vector<int> bestevs(nvec, -1);
            for (int j = 0; j < nvec; j++)
            {
                dtype crit;
                dtype mincrit = std::numeric_limits<dtype>::max();

                for (int i = 0;i < nextrap_real*nvec;i++)
                {
                    bool found = false;
                    for (int k = 0;k < j;k++)
                    {
                        if (bestevs[k] == i) found = true;
                    }
                    if (found) continue;

                    switch (mode)
                    {
                        case GUESS_OVERLAP:
                            crit = 0.0;
                            for (int m = 0;m < nextrap_real;m++)
                                for (int k = 0;k < nvec;k++)
                                    crit -= vr[k][m][i]*overlap[k][m];
                            break;
                        case LOWEST_ENERGY:
                            crit = std::real(l[i]);
                            break;
                        case CLOSEST_ENERGY:
                            crit = std::abs(std::real(l[i])-target[j]);
                            break;
                    }

                    if (crit < mincrit && std::abs(std::imag(l[i])) < 1e-9)
                    {
                        mincrit = crit;
                        bestevs[j] = i;
                    }
                }

                assert(bestevs[j] != -1);
            }

            std::cout << "bestevs[0] = " << bestevs[0] << std::endl;
            std::cout << "bestevs[1] = " << bestevs[1] << std::endl;
            // if (bestevs[0] == bestevs[1])
            // {
            //     bestevs[1] = bestevs[0] + 1;
            // }
            //std::cout << "Check 9" << std::endl;

            //if (std::abs(std::imag(l[bestev])) > 1e-5)
            //    throw std::runtime_error("davidson: complex eigenvalue");

            /*
             * Calculate residuals and apply Davidson correction
             */
            for (int j = 0;j < nvec;j++)
            {
                /*
                 * Form current trial vector y and H*y = x
                 */
                *c [j] = (*old_c [0][j])*vr[j][0][bestevs[j]];
                *hc[j] = (*old_hc[0][j])*vr[j][0][bestevs[j]]; // original times the evec element for it

                for (int i = 1;i < nextrap_real;i++)
                {
                    *c [j] += (*old_c [i][j])*vr[j][i][bestevs[j]]; // weight each old c by its evec value
                    *hc[j] += (*old_hc[i][j])*vr[j][i][bestevs[j]]; // same for old_hc. making the new state state
                }

                // so now hc is V*y = x and c is A*V*y = A*x

                *hc[j] -= std::real(l[bestevs[j]])*(*c[j]);

                /*
                 * If we were to determine convergence at this point, then c
                 * would be the solution vector
                 *
                 * This needs to be fixed, since right now the caller has no way
                 * of obtaining the solution!
                 */

                // now c = A*x - mu*x = -r

                *c[j] = *hc[j]; // This is what we norm to determine convergence, which is r, makes sense.
                c[j]->weight(D, std::real(l[bestevs[j]])); // Look into weight function

                /*
                 * Orthogonalize new guess vector against existing subspace and against
                 * other guess vectors
                 */
                for (int k = 0; k < j; k++)
                {
                    dtype olap = scalar(conj(*c[j])*(*c[k])); // orthogonalize against new cs which have just been constructed
                    *c[j] -= (*c[k])*olap;
                }

                for (int i = 0;i < nextrap_real;i++)
                {
                    for (int k = 0; k < nvec; k++)
                    {
                        dtype olap = scalar(conj(*c[j])*(*old_c[i][k])); // orthogonalize against all previous cs
                        *c[j] -= (*old_c[i][k])*olap;
                    }
                }

                /*
                 * Normalize
                 */
                double norm = sqrt(std::abs(scalar(conj(*c[j])*(*c[j]))));
                *c[j] /= norm;
            }

            std::vector<dtype> myreturn(nvec);
            for (int i = 0; i < nvec; i++)
                myreturn[i] = std::real(l[bestevs[i]]);

            return myreturn;
        }
};

}
}

#endif
