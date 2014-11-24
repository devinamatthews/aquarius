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

namespace
{

template <typename T> bool absGreaterThan(const T& a, const T& b)
{
    return std::abs(a) > std::abs(b);
}

}

template<typename T>
class Davidson : public task::Destructible
{
    private:
        Davidson(const Davidson& other);

        Davidson& operator=(const Davidson& other);

    protected:
        typedef typename T::dtype dtype;
        std::vector< std::vector<T*> > old_c; // hold all R[k][i] where i is over nvec and k is over maxextrap
        std::vector< std::vector<T*> > old_hc; // hold all H*R[i] = Z[i] at every iteration aka Z[k][i]
        std::vector<T*> guess;
        std::tensor<dtype,3> guess_overlap;
        std::tensor<dtype,4> s, e; // e will hold chc[k]
        std::vector<dtype> c;
        int nvec, maxextrap, nextrap; // number of energies, number of iterations
        std::vector<int> mode;
        std::vector<dtype> target;
        std::vector<dtype> previous;

        enum {GUESS_OVERLAP, LOWEST_ENERGY, CLOSEST_ENERGY};

        void init(const input::Config& config)
        {
            nextrap = 0;
            maxextrap = config.get<int>("order"); // max number of iterations

            assert(nvec > 0);
            assert(maxextrap > 0);

            guess_overlap.resize(std::vec(nvec, maxextrap, nvec));
            s.resize(std::vec(nvec, maxextrap, nvec, maxextrap));
            e.resize(std::vec(nvec, maxextrap, nvec, maxextrap));
            c.resize(maxextrap*nvec);

            old_c.resize(maxextrap, std::vector<T*>(nvec, (T*)NULL));
            old_hc.resize(maxextrap, std::vector<T*>(nvec, (T*)NULL));

            for (int i = 0;i < guess.size();i++)
            {
                guess[i] = new T(*guess[i]);
            }
        }

    public:
        Davidson(const input::Config& config, int nvec=1)
        : nvec(nvec), mode(nvec,LOWEST_ENERGY), target(nvec), previous(nvec)
        { init(config); }

        Davidson(const input::Config& config, dtype target)
        : nvec(1), mode(1,CLOSEST_ENERGY), target(1,target), previous(1)
        { init(config); }

        Davidson(const input::Config& config, const std::vector<dtype>& target)
        : nvec(target.size()), mode(target.size(),CLOSEST_ENERGY), target(target), previous(target.size())
        { init(config); }

        Davidson(const input::Config& config, const T& guess)
        : guess(1,&guess), nvec(1), mode(1,GUESS_OVERLAP), target(1), previous(1)
        { init(config); }

        Davidson(const input::Config& config, const std::vector<T*>& guess)
        : guess(guess), nvec(guess.size()), mode(guess.size(),GUESS_OVERLAP), target(guess.size()), previous(guess.size())
        { init(config); }

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

            // std::cout << setprecision(10) <<"Inf Norm hc = " << hc[0]->norm(00) << std::endl;
            // std::cout << setprecision(10) <<"Inf Norm c = " << c[0]->norm(00) << std::endl;

            assert(nvec == c.size() && nvec == hc.size());

            assert(c.size() == nvec);
            assert(hc.size() == nvec);

            /*
             * Check and normalize incoming c and H*c vectors
             */
            for (int i = 0;i < nvec;i++)
            {
                assert(c[i] != NULL);
                assert(hc[i] != NULL);

                double norm = sqrt(std::abs(scalar(conj(*c[i])*(*c[i]))));

                *c[i] /= norm;
                *hc[i] /= norm;

                // printf("Norm R %d %18.15f\n", i+1, (*c[i])(1)(vec(1,0),vec(0,1)).norm(2));

                // {
                //     vector<dtype> values;
                //     (*c[i])(1)(vec(1,0),vec(0,1))(vec(0,0)).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int j = 0;j < 30;j++)
                //     {
                //         int ia = keys[j]%19+6;
                //         int ii = keys[j]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[j]);
                //     }
                // }

                // printf("Norm H*R %d %18.15f\n", i+1, (*hc[i])(1)(vec(1,0),vec(0,1)).norm(2));

                // {
                //     vector<dtype> values;
                //     (*hc[i])(1)(vec(1,0),vec(0,1))(vec(0,0)).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int j = 0;j < 10;j++)
                //     {
                //         int ia = keys[j]%19+6;
                //         int ii = keys[j]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[j]);
                //     }
                // }
            }

            if (nextrap == maxextrap) // aka we've reached our maximum iteration
            {
                //TODO: compact
                assert(0);
            }

            /*
             * Lazily allocate elements of old_c etc. so that we can
             * just use the copy ctor and subclasses do not have to
             * worry about allocation/deallocation
             */

            if (old_c[nextrap][0] == NULL)
            {
                for (int i = 0;i < nvec;i++)
                {
                    old_c[nextrap][i] = new T(*c[i]);
                    {
                        vector<dtype> values;
                        (*old_c[nextrap][i])(1)(vec(1,0),vec(0,1))(vec(0,0)).getAllData(values);
                        vector<int64_t> keys = range<int64_t>(5*19);
                        cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                               absGreaterThan<dtype>);
                        //printf("Badness old c %d %d %d %15.12g\n", i+1, i+1, nextrap+1, std::abs(values[19-4*i]));
                    }
                }
            }
            else
            {
                // should only happen after compaction
                for (int i = 0;i < nvec;i++)
                    *old_c[nextrap][i] = *c[i];
            }

            if (old_hc[nextrap][0] == NULL)
            {
                for (int i = 0;i < nvec;i++)
                {
                    old_hc[nextrap][i] = new T(*hc[i]);
                    {
                        vector<dtype> values;
                        (*old_hc[nextrap][i])(1)(vec(1,0),vec(0,1))(vec(0,0)).getAllData(values);
                        vector<int64_t> keys = range<int64_t>(5*19);
                        cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                               absGreaterThan<dtype>);
                        //printf("Badness old hc %d %d %d %15.12g\n", i+1, i+1, nextrap+1, std::abs(values[19-4*i]));
                    }
                }
            }
            else
            {
                // should only happen after compaction
                for (int i = 0;i < nvec;i++)
                    *old_hc[nextrap][i] = *hc[i];
            }

            /*
             * Augment the subspace matrix with the new vectors
             */
            for (int i = 0;i < nvec;i++)
            {
                /*
                 * Compute and save overlap with guess for later
                 */
                if (!guess.empty())
                {
                    for (int j = 0;j < nvec;j++)
                    {
                        guess_overlap[i][nextrap][j] = std::abs(scalar(conj(*c[i])*(*guess[j])));
                    }
                }

                for (int j = 0;j < nvec;j++)
                {
                    // "Diagonal"
                    e[i][nextrap][j][nextrap] = scalar(conj(*hc[i])*(*c[j]));
                    s[i][nextrap][j][nextrap] = scalar(conj(*c[i])*(*c[j]));
                    //printf("s %d %d %d %18.15f\n", nextrap, i, j, s[i][nextrap][j][nextrap]);

                    // "Off-diagonal"
                    for (int k = 0;k < nextrap;k++)
                    {
                        e[i][nextrap][j][k] = scalar(conj(*old_hc[k][i])*(*c[j]));
                        e[i][k][j][nextrap] = scalar(conj(*hc[i])*(*old_c[k][j]));
                        s[i][nextrap][j][k] = scalar(conj(*old_c[k][i])*(*c[j]));
                        s[i][k][j][nextrap] = scalar(conj(*c[i])*(*old_c[k][j]));
                        //printf("s %d %d %d %18.15f\n", k, i, j, s[i][nextrap][j][nextrap]);
                    }
                }
            }

            nextrap++;

            /*
             * Diagonalize the subspace matrix to obtain approximate solutions
             */
            int info;
            std::vector<dtype> beta(nvec*nextrap);
            std::tensor<dtype,4> tmp1(std::vec(nvec,nextrap,nvec,nextrap));
            std::tensor<dtype,4> tmp2(std::vec(nvec,nextrap,nvec,nextrap));
            std::tensor<dtype,3> vr(std::vec(nvec,nextrap,nvec*nextrap));
            std::vector<typename std::complex_type<dtype>::type> l(nextrap*nvec); // Eigenvalues

            for (int m = 0;m < nextrap;m++)
            {
                for (int k = 0;k < nvec;k++)
                {
                    for (int j = 0;j < nextrap;j++)
                    {
                        for (int i = 0;i < nvec;i++)
                        {
                            // printf("%15.12f ", e[i][j][k][m]);
                            tmp1[i][j][k][m] = e[i][j][k][m];
                            tmp2[i][j][k][m] = s[i][j][k][m];
                        }
                    }
                    // printf("\n");
                }
            }

            //info = ggev('N', 'V', nextrap*nvec, tmp1.data(), nextrap*nvec,
            //            tmp2.data(), nextrap*nvec, l.data(), beta.data(), NULL, 1,
            //            vr.data(), nextrap*nvec);
            info = geev('N', 'V', nextrap*nvec, tmp1.data(), nextrap*nvec,
                        l.data(), NULL, 1,
                        vr.data(), nextrap*nvec);
            if (info != 0) throw std::runtime_error(std::strprintf("davidson: Info in ggev: %d", info));

            /*
             * Fix sign of eigenvectors
             */
            for (int i = 0;i < nextrap*nvec;i++)
            {
                //l[i] /= beta[i];
                // printf("%15.12f\n", std::real(l[i]));

                for (int m = 0;m < nextrap;m++)
                {
                    for (int k = 0;k < nvec;k++)
                    {
                        if (std::abs(vr[k][m][i]) > 1e-10)
                        {
                            if (std::real(vr[k][m][i]) < 0)
                            {
                                scal(nextrap*nvec, -1, &vr[0][0][i], 1);
                            }
                            m = nextrap;
                            break;
                        }
                    }
                }
            }

            // std::cout << "evec check:" << std::endl;
            //  for (int i = 0;i < nextrap*nvec;i++)
            // {
            //     for (int m = 0;m < nextrap;m++)
            //     {
            //         for (int k = 0;k < nvec;k++) 
            //         {
            //             std::cout << vr[k][m][i] << std::endl;
            //         }
            //     }
            // }

            /*
             * Assign eigenvalues (exclusively) to states by the selected criterion
             */
            std::vector<int> bestevs(nvec, -1);
            for (int j = 0; j < nvec; j++)
            {
                dtype crit;
                dtype mincrit = std::numeric_limits<dtype>::max();

                for (int i = 0;i < nextrap*nvec;i++)
                {
                    //if (j == 0) printf("%15.12f\n", std::real(l[i]));

                    bool found = false;
                    for (int k = 0;k < j;k++)
                    {
                        if (bestevs[k] == i) found = true;
                    }
                    if (found) continue;

                    switch (mode[j])
                    {
                        case GUESS_OVERLAP:
                            crit = 0.0;
                            for (int m = 0;m < nextrap;m++)
                                for (int k = 0;k < nvec;k++)
                                    crit -= vr[k][m][i]*guess_overlap[k][m][j];
                            break;
                        case LOWEST_ENERGY:
                            crit = std::real(l[i]);
                            break;
                        case CLOSEST_ENERGY:
                            crit = std::abs(std::real(l[i])-target[j]);
                            break;
                    }

                    if (crit < mincrit && std::abs(std::imag(l[i])) < 1e-12)
                    {
                        mincrit = crit;
                        bestevs[j] = i;
                    }
                }

                assert(bestevs[j] != -1);

                if (nextrap > 1 && mode[j] != CLOSEST_ENERGY &&
                    std::abs(previous[j]-std::real(l[bestevs[j]])) < 1e-6)
                {
                    cout << "Locking root " << (j+1) << endl;
                    mode[j] = CLOSEST_ENERGY;
                    target[j] = std::real(l[bestevs[j]]);
                }
                previous[j] = std::real(l[bestevs[j]]);
            }

            /*
             * Calculate residuals and apply Davidson correction
             */
            for (int j = 0;j < nvec;j++)
            {
                /*
                 * Form current trial vector y and H*y = x
                 */
                *c [j] = 0;
                *hc[j] = 0;

                for (int i = nextrap-1;i >= 0;i--)
                {
                    for (int k = nvec-1;k >= 0;k--)
                    {
                        *c [j] += (*old_c [i][k])*vr[k][i][bestevs[j]]; // weight each old c by its evec value
                        *hc[j] += (*old_hc[i][k])*vr[k][i][bestevs[j]]; // same for old_hc. making the new state state
                    }
                }

                // std::cout << setprecision(10) <<"Inf Norm hc = " << hc[j]->norm(00) << std::endl;
                // std::cout << setprecision(10) <<"Inf Norm c = " << c[j]->norm(00) << std::endl;

                // so now c is V*y = x and hc is A*V*y = A*x

                // printf("Norm c %d %15.12f\n", j+1, (*c[j])(1)(vec(1,0),vec(0,1)).norm(2));
                // {
                //     vector<dtype> values;
                //     (*c[j])(1)(vec(1,0),vec(0,1))(vec(0,0)).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 30;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }
                //     printf("Badness     c %d %15.12g\n", j+1, std::abs(values[19-4*j]));
                // }
                // printf("Norm H*c %d %15.12f\n", j+1, (*hc[j])(1)(vec(1,0),vec(0,1)).norm(2));
                // {
                //     vector<dtype> values;
                //     (*hc[j])(1)(vec(1,0),vec(0,1))(vec(0,0)).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 30;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }
                //     printf("Badness   H*c %d %15.12g\n", j+1, std::abs(values[19-4*j]));
                // }
                *hc[j] -= std::real(l[bestevs[j]])*(*c[j]);
                // std::cout << setprecision(10) <<"Inf Norm hc = " << hc[j]->norm(00) << std::endl;

                /*
                 * If we were to determine convergence at this point, then c
                 * would be the solution vector
                 *
                 * This needs to be fixed, since right now the caller has no way
                 * of obtaining the solution!
                 */

                // now hc = A*x - mu*x = -r

                *c[j] = -*hc[j]; // This is what we norm to determine convergence, which is r, makes sense.
                // std::cout << setprecision(10) <<"Inf Norm c = " << c[j]->norm(00) << std::endl;
                // printf("Norm r %d %15.12f\n", j+1, (*c[j])(1)(vec(1,0),vec(0,1)).norm(2));
                // {
                //     vector<dtype> values;
                //     (*c[j])(1)(vec(1,0),vec(0,1))(vec(0,0)).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 10;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }
                //     printf("Badness     r %d %15.12g\n", j+1, std::abs(values[19-4*j]));
                // }
                c[j]->weight(D, std::real(l[bestevs[j]])); // Look into weight function
                // std::cout << setprecision(10) <<"Inf Norm c = " << c[j]->norm(00) << std::endl;
                // printf("Norm d %d %15.12f\n", j+1, (*c[j])(1)(vec(1,0),vec(0,1)).norm(2));
                // {
                //     vector<dtype> values;
                //     (*c[j])(1)(vec(1,0),vec(0,1))(vec(0,0)).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 10;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }
                //     printf("Badness     d %d %15.12g\n", j+1, std::abs(values[19-4*j]));
                // }

                double norm = sqrt(std::abs(scalar(conj(*c[j])*(*c[j]))));
                *c[j] /= norm;
                

                //printf("%d %15.12g\n", j+1, norm);
            }

            for (int j = 0;j < nvec;j++)
            {
                for (int k = 0;k < nvec;k++)
                {
                    for (int i = 0;i < nextrap;i++)
                    {
                        dtype olap = scalar(conj(*c[j])*(*old_c[i][k]));
                        *c[j] -= olap*(*old_c[i][k]);
                    }

                    if (k == j) continue;
                    dtype olap = scalar(conj(*c[j])*(*c[k]));
                    *c[j] -= olap*(*c[k]);
                }

                for (int k = nvec-1;k >= 0;k--)
                {
                    for (int i = nextrap-1;i >= 0;i--)
                    {
                        dtype olap = scalar(conj(*c[j])*(*old_c[i][k]));
                        *c[j] -= olap*(*old_c[i][k]);
                    }

                    if (k == j) continue;
                    dtype olap = scalar(conj(*c[j])*(*c[k]));
                    *c[j] -= olap*(*c[k]);
                }

                double norm = sqrt(std::abs(scalar(conj(*c[j])*(*c[j]))));
                *c[j] /= norm;

                // std::cout << setprecision(10) <<"Inf Norm c = " << c[j]->norm(00) << std::endl;

                // printf("Norm new c %d %15.12f\n", j+1, (*c[j])(1)(vec(1,0),vec(0,1)).norm(2));
                // {
                //     vector<dtype> values;
                //     (*c[j])(1)(vec(1,0),vec(0,1))(vec(0,0)).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 10;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }

                //     printf("Badness new c %d %15.12g\n", j+1, std::abs(values[19-4*j]));
                // }
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
