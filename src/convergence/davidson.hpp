#ifndef _AQUARIUS_DAVIDSON_HPP_
#define _AQUARIUS_DAVIDSON_HPP_

#include "util/global.hpp"

#include "input/config.hpp"
#include "task/task.hpp"
#include "operator/denominator.hpp"

namespace aquarius
{
namespace convergence
{

namespace
{

template <typename T> bool absGreaterThan(const T& a, const T& b)
{
    return abs(a) > abs(b);
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
        vector<unique_vector<T>> old_c; // hold all R[k][i] where i is over nvec and k is over maxextrap
        vector<unique_vector<T>> old_hc; // hold all H*R[i] = Z[i] at every iteration aka Z[k][i]
        unique_vector<T> guess;
        marray<dtype,3> guess_overlap;
        marray<dtype,4> s, e; // e will hold chc[k]
        vector<dtype> c;
        int nvec, maxextrap, nextrap; // number of energies, number of iterations
        vector<int> mode;
        vector<dtype> target;
        vector<dtype> previous;
        vector<int> root;
        vector<typename complex_type<dtype>::type> l;
        marray<dtype,3> vr;

        enum {GUESS_OVERLAP, LOWEST_ENERGY, CLOSEST_ENERGY};

        void init(const input::Config& config)
        {
            nextrap = 0;
            maxextrap = config.get<int>("order"); // max number of iterations

            assert(nvec > 0);
            assert(maxextrap > 0);

            guess_overlap.resize({nvec, maxextrap, nvec});
            s.resize({nvec, maxextrap, nvec, maxextrap});
            e.resize({nvec, maxextrap, nvec, maxextrap});
            c.resize(maxextrap*nvec);

            old_c.resize(maxextrap);
            old_hc.resize(maxextrap);

            root.resize(nvec);
            l.resize(nvec*maxextrap);
            vr.resize({nvec,maxextrap,nvec*maxextrap});
        }

    public:
        Davidson(const input::Config& config, int nvec=1)
        : nvec(nvec), mode(nvec,LOWEST_ENERGY), target(nvec), previous(nvec)
        { init(config); }

        Davidson(const input::Config& config, dtype target)
        : nvec(1), mode(1,CLOSEST_ENERGY), target(1,target), previous(1)
        { init(config); }

        Davidson(const input::Config& config, vector<dtype>&& target)
        : nvec(target.size()), mode(target.size(),CLOSEST_ENERGY),
          target(forward<vector<dtype>>(target)), previous(target.size())
        { init(config); }

        Davidson(const input::Config& config, T&& guess)
        : guess(1,forward<T>(guess)), nvec(1), mode(1,GUESS_OVERLAP), target(1), previous(1)
        { init(config); }

        Davidson(const input::Config& config, unique_vector<T>&& guess)
        : guess(move(guess)), nvec(guess.size()), mode(guess.size(),GUESS_OVERLAP),
          target(guess.size()), previous(guess.size())
        { init(config); }

        template <typename guess_container>
        Davidson(const input::Config& config, const guess_container& guess)
        : nvec(guess.size()), mode(guess.size(),GUESS_OVERLAP),
          target(guess.size()), previous(guess.size())
        {
            init(config);
            for (auto& g : guess)
            {
                this->guess.push_back(g);
            }
        }

        template <typename guess_container>
        Davidson(const input::Config& config, guess_container&& guess)
        : nvec(guess.size()), mode(guess.size(),GUESS_OVERLAP),
          target(guess.size()), previous(guess.size())
        {
            init(config);
            for (auto& g : guess)
            {
                this->guess.push_back(move(g));
            }
        }

        dtype extrapolate(T& c, T& hc, const op::Denominator<dtype>& D)
        {
            assert(nvec == 1);
            return extrapolate(ptr_vector<T>{&c}, ptr_vector<T>{&hc}, D)[0];
        }

        template <typename c_container, typename hc_container>
        vector<dtype> extrapolate(c_container&& c, hc_container&& hc, const op::Denominator<dtype>& D)
        {
            // cout << setprecision(10) <<"Inf Norm hc = " << hc[0]->norm(00) << endl;
            // cout << setprecision(10) <<"Inf Norm c = " << c[0]->norm(00) << endl;

            assert(nvec == c.size() && nvec == hc.size());

            assert(c.size() == nvec);
            assert(hc.size() == nvec);

            /*
             * Check and normalize incoming c and H*c vectors
             */
            for (int i = 0;i < nvec;i++)
            {
                double norm = sqrt(abs(scalar(conj(c[i])*c[i])));

                c[i] /= norm;
                hc[i] /= norm;

                // printf("Norm R %d %18.15f\n", i+1, (*c[i])(1)({1,0},{0,1}).norm(2));

                // {
                //     vector<dtype> values;
                //     (*c[i])(1)({1,0},{0,1})({0,0}).getAllData(values);
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

                // printf("Norm H*R %d %18.15f\n", i+1, (*hc[i])(1)({1,0},{0,1}).norm(2));

                // {
                //     vector<dtype> values;
                //     (*hc[i])(1)({1,0},{0,1})({0,0}).getAllData(values);
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

            if (old_c[nextrap].empty())
            {
                for (int i = 0;i < nvec;i++)
                {
                    old_c[nextrap].emplace_back(c[i]);
                    //{
                    //    vector<dtype> values;
                    //    (*old_c[nextrap][i])(1)({1,0},{0,1})({0,0}).getAllData(values);
                    //    vector<int64_t> keys = range<int64_t>(5*19);
                    //    cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                    //           absGreaterThan<dtype>);
                    //    //printf("Badness old c %d %d %d %15.12g\n", i+1, i+1, nextrap+1, abs(values[19-4*i]));
                    //}
                }
            }
            else
            {
                // should only happen after compaction
                for (int i = 0;i < nvec;i++)
                {
                    old_c[nextrap][i] = c[i];
                }
            }

            if (old_hc[nextrap].empty())
            {
                for (int i = 0;i < nvec;i++)
                {
                    old_hc[nextrap].emplace_back(hc[i]);
                    //{
                    //    vector<dtype> values;
                    //    (*old_hc[nextrap][i])(1)({1,0},{0,1})({0,0}).getAllData(values);
                    //    vector<int64_t> keys = range<int64_t>(5*19);
                    //    cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                    //           absGreaterThan<dtype>);
                    //    //printf("Badness old hc %d %d %d %15.12g\n", i+1, i+1, nextrap+1, abs(values[19-4*i]));
                    //}
                }
            }
            else
            {
                // should only happen after compaction
                for (int i = 0;i < nvec;i++)
                {
                    old_hc[nextrap][i] = hc[i];
                }
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
                        guess_overlap[i][nextrap][j] = abs(scalar(conj(c[i])*guess[j]));
                    }
                }

                for (int j = 0;j < nvec;j++)
                {
                    // "Diagonal"
                    e[i][nextrap][j][nextrap] = scalar(conj(hc[i])*c[j]);
                    s[i][nextrap][j][nextrap] = scalar(conj( c[i])*c[j]);
                    //printf("s %d %d %d %18.15f\n", nextrap, i, j, s[i][nextrap][j][nextrap]);

                    // "Off-diagonal"
                    for (int k = 0;k < nextrap;k++)
                    {
                        e[i][nextrap][j][k] = scalar(conj(old_hc[k][i])*    c   [j]);
                        e[i][k][j][nextrap] = scalar(conj(    hc   [i])*old_c[k][j]);
                        s[i][nextrap][j][k] = scalar(conj( old_c[k][i])*    c   [j]);
                        s[i][k][j][nextrap] = scalar(conj(     c   [i])*old_c[k][j]);
                        //printf("s %d %d %d %18.15f\n", k, i, j, s[i][nextrap][j][nextrap]);
                    }
                }
            }

            nextrap++;

            /*
             * Diagonalize the subspace matrix to obtain approximate solutions
             */
            int info;
            vector<dtype> beta(nvec*nextrap);
            marray<dtype,4> tmp1({nvec,nextrap,nvec,nextrap});
            marray<dtype,4> tmp2({nvec,nextrap,nvec,nextrap});
            marray<dtype,3> tmp3({nvec,nextrap,nvec*nextrap});

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
                        tmp3.data(), nextrap*nvec);
            if (info != 0) throw runtime_error(strprintf("davidson: Info in ggev: %d", info));

            for (int k = 0;k < nvec*nextrap;k++)
            {
                for (int j = 0;j < nextrap;j++)
                {
                    for (int i = 0;i < nvec;i++)
                    {
                        vr[i][j][k] = tmp3[i][j][k];
                    }
                }
            }

            /*
             * Fix sign of eigenvectors
             */
            for (int i = 0;i < nextrap*nvec;i++)
            {
                //l[i] /= beta[i];
                // printf("%15.12f\n", real(l[i]));

                for (int m = 0;m < nextrap;m++)
                {
                    for (int k = 0;k < nvec;k++)
                    {
                        if (abs(vr[k][m][i]) > 1e-10)
                        {
                            if (real(vr[k][m][i]) < 0)
                            {
                                scal(nextrap*nvec, -1, &vr[0][0][i], 1);
                            }
                            m = nextrap;
                            break;
                        }
                    }
                }
            }

            // cout << "evec check:" << endl;
            //  for (int i = 0;i < nextrap*nvec;i++)
            // {
            //     for (int m = 0;m < nextrap;m++)
            //     {
            //         for (int k = 0;k < nvec;k++)
            //         {
            //             cout << vr[k][m][i] << endl;
            //         }
            //     }
            // }

            /*
             * Assign eigenvalues (exclusively) to states by the selected criterion
             */
            for (int j = 0; j < nvec; j++)
            {
                root[j] = -1;

                dtype crit = numeric_limits<dtype>::max();
                dtype mincrit = numeric_limits<dtype>::max();

                for (int i = 0;i < nextrap*nvec;i++)
                {
                    //if (j == 0) printf("%15.12f\n", real(l[i]));

                    bool found = false;
                    for (int k = 0;k < j;k++)
                    {
                        if (root[k] == i) found = true;
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
                            crit = real(l[i]);
                            break;
                        case CLOSEST_ENERGY:
                            crit = abs(real(l[i])-target[j]);
                            break;
                    }

                    if (crit < mincrit && abs(imag(l[i])) < 1e-12)
                    {
                        mincrit = crit;
                        root[j] = i;
                    }
                }

                assert(root[j] != -1);

                if (nextrap > 1 && mode[j] != CLOSEST_ENERGY &&
                    abs(previous[j]-real(l[root[j]])) < 1e-4)
                {
                    cout << "Locking root " << (j+1) << endl;
                    mode[j] = CLOSEST_ENERGY;
                    target[j] = real(l[root[j]]);
                }
                previous[j] = real(l[root[j]]);
            }

            /*
             * Calculate residuals and apply Davidson correction
             */
            for (int j = 0;j < nvec;j++)
            {
                /*
                 * Form current trial vector y and H*y = x
                 */
                c [j] = 0;
                hc[j] = 0;

                for (int i = nextrap-1;i >= 0;i--)
                {
                    for (int k = nvec-1;k >= 0;k--)
                    {
                        c [j] +=  old_c[i][k]*vr[k][i][root[j]]; // weight each old c by its evec value
                        hc[j] += old_hc[i][k]*vr[k][i][root[j]]; // same for old_hc. making the new state state
                    }
                }

                // cout << setprecision(10) <<"Inf Norm hc = " << hc[j]->norm(00) << endl;
                // cout << setprecision(10) <<"Inf Norm c = " << c[j]->norm(00) << endl;

                // so now c is V*y = x and hc is A*V*y = A*x

                // printf("Norm c %d %15.12f\n", j+1, (*c[j])(1)({1,0},{0,1}).norm(2));
                // {
                //     vector<dtype> values;
                //     (*c[j])(1)({1,0},{0,1})({0,0}).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 30;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }
                //     printf("Badness     c %d %15.12g\n", j+1, abs(values[19-4*j]));
                // }
                // printf("Norm H*c %d %15.12f\n", j+1, (*hc[j])(1)({1,0},{0,1}).norm(2));
                // {
                //     vector<dtype> values;
                //     (*hc[j])(1)({1,0},{0,1})({0,0}).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 30;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }
                //     printf("Badness   H*c %d %15.12g\n", j+1, abs(values[19-4*j]));
                // }
                hc[j] -= real(l[root[j]])*c[j];
                // cout << setprecision(10) <<"Inf Norm hc = " << hc[j]->norm(00) << endl;

                // now hc = A*x - mu*x = -r

                c[j] = -hc[j]; // This is what we norm to determine convergence, which is r, makes sense.
                // cout << setprecision(10) <<"Inf Norm c = " << c[j]->norm(00) << endl;
                // printf("Norm r %d %15.12f\n", j+1, (*c[j])(1)({1,0},{0,1}).norm(2));
                // {
                //     vector<dtype> values;
                //     (*c[j])(1)({1,0},{0,1})({0,0}).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 10;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }
                //     printf("Badness     r %d %15.12g\n", j+1, abs(values[19-4*j]));
                // }
                c[j].weight(D, real(l[root[j]])); // Look into weight function
                // cout << setprecision(10) <<"Inf Norm c = " << c[j]->norm(00) << endl;
                // printf("Norm d %d %15.12f\n", j+1, (*c[j])(1)({1,0},{0,1}).norm(2));
                // {
                //     vector<dtype> values;
                //     (*c[j])(1)({1,0},{0,1})({0,0}).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 10;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }
                //     printf("Badness     d %d %15.12g\n", j+1, abs(values[19-4*j]));
                // }

                c[j] /= sqrt(abs(scalar(conj(c[j])*c[j])));

                //if (scalar((*c[j])(1)({1,0},{0,1})*(*c[j])(1)({0,0},{0,0})) < 0)
                //{
                //    cout << "NOOOOO!" << endl;
                //}

                //0.5*(*c[j])(1)({0,0},{0,0})[  "ai"] += 0.5*(*c[j])(1)({1,0},{0,1})[  "ai"];
                //    (*c[j])(1)({1,0},{0,1})[  "ai"]  =     (*c[j])(1)({0,0},{0,0})[  "ai"];
                //0.5*(*c[j])(2)({1,0},{0,1})["abij"] += 0.5*(*c[j])(2)({1,0},{0,1})["baji"];
                //0.5*(*c[j])(2)({0,0},{0,0})["abij"] += 0.5*(*c[j])(2)({2,0},{0,2})["abij"];
                //    (*c[j])(2)({2,0},{0,2})["abij"]  =     (*c[j])(2)({0,0},{0,0})["abij"];

                //printf("%d %15.12g\n", j+1, norm);
            }

            for (int j = 0;j < nvec;j++)
            {
                for (int k = 0;k < nvec;k++)
                {
                    for (int i = 0;i < nextrap;i++)
                    {
                        c[j] -= scalar(conj(c[j])*old_c[i][k])*old_c[i][k];
                    }

                    if (k != j)
                    {
                        c[j] -= scalar(conj(c[j])*c[k])*c[k];
                    }
                }

                for (int k = nvec-1;k >= 0;k--)
                {
                    for (int i = nextrap-1;i >= 0;i--)
                    {
                        c[j] -= scalar(conj(c[j])*old_c[i][k])*old_c[i][k];
                    }

                    if (k != j)
                    {
                        c[j] -= scalar(conj(c[j])*c[k])*c[k];
                    }
                }

                c[j] /= sqrt(abs(scalar(conj(c[j])*c[j])));

                // cout << setprecision(10) <<"Inf Norm c = " << c[j]->norm(00) << endl;

                // printf("Norm new c %d %15.12f\n", j+1, (*c[j])(1)({1,0},{0,1}).norm(2));
                // {
                //     vector<dtype> values;
                //     (*c[j])(1)({1,0},{0,1})({0,0}).getAllData(values);
                //     vector<int64_t> keys = range<int64_t>(5*19);
                //     cosort(values.begin(), values.end(), keys.begin(), keys.end(),
                //            absGreaterThan<dtype>);

                //     for (int k = 0;k < 10;k++)
                //     {
                //         int ia = keys[k]%19+6;
                //         int ii = keys[k]/19+1;
                //         //printf("%2d %2d %18.15f\n", ii, ia, values[k]);
                //     }

                //     printf("Badness new c %d %15.12g\n", j+1, abs(values[19-4*j]));
                // }
            }

            vector<dtype> myreturn(nvec);
            for (int i = 0; i < nvec; i++)
                myreturn[i] = real(l[root[i]]);

            return myreturn;
        }

        void getSolution(int j, T& s)
        {
            s = 0;
            for (int i = nextrap-1;i >= 0;i--)
            {
                for (int k = nvec-1;k >= 0;k--)
                {
                    s += old_c[i][k]*vr[k][i][root[j]];
                }
            }
        }

        void clear()
        {
            nextrap = 0;
        }
};

}
}

#endif
