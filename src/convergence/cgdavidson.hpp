#ifndef _AQUARIUS_CGDAVIDSON_HPP_
#define _AQUARIUS_CGDAVIDSON_HPP_

#include "util/global.hpp"

#include "input/config.hpp"
#include "task/task.hpp"
#include "operator/denominator.hpp"

#include "davidson.hpp"

namespace aquarius
{
namespace convergence
{

template<typename T>
class CGDavidson : public task::Destructible
{
    private:
        CGDavidson(const CGDavidson& other);

        CGDavidson& operator=(const CGDavidson& other);

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
        vector<complex_type_t<dtype>> l;
        marray<dtype,3> vr;
        vector<bool> lock;
        vector<dtype> lock_e;

        enum {GUESS_OVERLAP, LOWEST_ENERGY, CLOSEST_ENERGY};

        void init(const input::Config& config)
        {
            nextrap = 0;
            maxextrap = config.get<int>("order"); // max number of iterations

            assert(nvec > 0);
            assert(maxextrap > 0);

            guess_overlap.resize({nvec, nvec, maxextrap});
            s.resize({nvec, maxextrap, nvec, maxextrap});
            e.resize({nvec, maxextrap, nvec, maxextrap});
            c.resize(maxextrap*nvec);

            old_c.resize(maxextrap);
            old_hc.resize(maxextrap);

            root.resize(nvec);
            l.resize(nvec*maxextrap);
            vr.resize({nvec*maxextrap,nvec,maxextrap});

            lock.resize(nvec, false);
            lock_e.resize(nvec);
        }

    public:
        CGDavidson(const input::Config& config, int nvec=1)
        : nvec(nvec), mode(nvec,LOWEST_ENERGY), target(nvec), previous(nvec)
        { init(config); }

        CGDavidson(const input::Config& config, dtype target)
        : nvec(1), mode(1,CLOSEST_ENERGY), target(1,target), previous(1)
        { init(config); }

        CGDavidson(const input::Config& config, vector<dtype>&& target)
        : nvec(target.size()), mode(target.size(),CLOSEST_ENERGY),
          target(forward<vector<dtype>>(target)), previous(target.size())
        { init(config); }

        CGDavidson(const input::Config& config, T&& guess)
        : guess(1,forward<T>(guess)), nvec(1), mode(1,GUESS_OVERLAP), target(1), previous(1)
        { init(config); }

        CGDavidson(const input::Config& config, unique_vector<T>&& guess)
        : guess(move(guess)), nvec(guess.size()), mode(guess.size(),GUESS_OVERLAP),
          target(guess.size()), previous(guess.size())
        { init(config); }

        template <typename guess_container>
        CGDavidson(const input::Config& config, const guess_container& guess)
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
        CGDavidson(const input::Config& config, guess_container&& guess)
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
            assert(nvec == c.size() && nvec == hc.size());
	    
            assert(c.size() == nvec);
            assert(hc.size() == nvec);

            /*
             * Check and normalize incoming c and H*c vectors
             */
            for (int i = 0;i < nvec;i++)
            {
                double norm = sqrt(aquarius::abs(scalar(conj(c[i])*c[i])));
                c[i] /= norm;
                hc[i] /= norm;
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
                        guess_overlap[j][i][nextrap] = aquarius::abs(scalar(conj(c[i])*guess[j]));
                    }
                }

                for (int j = 0;j < nvec;j++)
                {
                    // "Diagonal"
                    e[i][nextrap][j][nextrap] = scalar(conj(c[i])*hc[j]);
                    s[i][nextrap][j][nextrap] = scalar(conj(c[i])* c[j]);
                    //printf("s %d %d %d %18.15f\n", nextrap, i, j, s[i][nextrap][j][nextrap]);

                    // "Off-diagonal"
                    for (int k = 0;k < nextrap;k++)
                    {
                        e[i][k][j][nextrap] = scalar(conj(old_c[k][i])*    hc   [j]);
                        e[i][nextrap][j][k] = scalar(conj(    c   [i])*old_hc[k][j]);
                        s[i][k][j][nextrap] = scalar(conj(old_c[k][i])*     c   [j]);
                        s[i][nextrap][j][k] = scalar(conj(    c   [i])* old_c[k][j]);
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
            marray<dtype,3> tmp3({nvec*nextrap,nvec,nextrap});

            for (int m = 0;m < nextrap;m++)
            {
                for (int k = 0;k < nvec;k++)
                {
                    for (int j = 0;j < nextrap;j++)
                    {
                        for (int i = 0;i < nvec;i++)
                        {
                            //printf("%15.12f ", e[i][j][k][m]);
                            tmp1[i][j][k][m] = e[k][m][i][j];
                            tmp2[i][j][k][m] = s[k][m][i][j];
                        }
                    }
                    //printf("\n");
                }
            }

            info = ggev('N', 'V', nextrap*nvec, tmp1.data(), nextrap*nvec,
                        tmp2.data(), nextrap*nvec, l.data(), beta.data(), NULL, 1,
                        tmp3.data(), nextrap*nvec);
            //info = geev('N', 'V', nextrap*nvec, tmp1.data(), nextrap*nvec,
            //            l.data(), NULL, 1, tmp3.data(), nextrap*nvec);
            if (info != 0) throw runtime_error(strprintf("davidson: Info in ggev: %d", info));

            for (int k = 0;k < nvec*nextrap;k++)
            {
                for (int j = 0;j < nextrap;j++)
                {
                    for (int i = 0;i < nvec;i++)
                    {
                        vr[k][i][j] = tmp3[k][i][j];
                    }
                }
            }

            /*
             * Fix sign of eigenvectors
             */
            for (int i = 0;i < nextrap*nvec;i++)
            {
                l[i] /= beta[i];
                //printf("%15.12f\n", real(l[i]));

                int k = i/nextrap;
                int m = nextrap-1;
                if (real(vr[i][k][m]) < 0)
                {
                    scal(nextrap*nvec, -1, vr[i].data(), 1);
                }

                /*
                for (int m = 0;m < nextrap;m++)
                {
                    for (int k = 0;k < nvec;k++)
                    {
                        if (aquarius::abs(vr[i][k][m]) > 1e-10)
                        {
                            if (real(vr[i][k][m]) < 0)
                            {
                                scal(nextrap*nvec, -1, vr[i].data(), 1);
                            }
                            m = nextrap;
                            break;
                        }
                    }
                }
                */
            }

            // cout << "evec check:" << endl;
            //  for (int i = 0;i < nextrap*nvec;i++)
            // {
            //     for (int m = 0;m < nextrap;m++)
            //     {
            //         for (int k = 0;k < nvec;k++)
            //         {
            //             cout << vr[i][k][m] << endl;
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
                    bool found = false;
                    for (int k = 0;k < j;k++)
                    {
                        if (root[k] == i) found = true;
                    }
                    if (found) continue;

                    if (lock[j])
                    {
                        crit = aquarius::abs(real(l[i])-lock_e[j]);
                    }
                    else if (mode[j] == GUESS_OVERLAP)
                    {
                        crit = 0.0;
                        for (int m = 0;m < nextrap;m++)
                            for (int k = 0;k < nvec;k++)
                                crit -= vr[i][k][m]*guess_overlap[j][k][m];
                    }
                    else if (mode[j] == LOWEST_ENERGY)
                    {
                        crit = real(l[i]);
                    }
                    else if (mode[j] == CLOSEST_ENERGY)
                    {
                        crit = aquarius::abs(real(l[i])-target[j]);
                    }

                    if (crit < mincrit && aquarius::abs(imag(l[i])) < 1e-12)
                    {
                        mincrit = crit;
                        root[j] = i;
                    }
                }

                assert(root[j] != -1);

                if (nextrap > 1 && mode[j] != CLOSEST_ENERGY && !lock[j] &&
                    aquarius::abs(previous[j]-real(l[root[j]])) < 1e-4)
                {
                    if (c[0].arena.rank == 0)
                        cout << "Locking root " << (j+1) << endl;
                    lock[j] = true;
                    lock_e[j] = real(l[root[j]]);
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

                for (int i = 0;i < nextrap;i++)
                {
                    for (int k = 0;k < nvec;k++)
                    {
                        c [j] +=  old_c[i][k]*vr[root[j]][k][i]; // weight each old c by its evec value
                        hc[j] += old_hc[i][k]*vr[root[j]][k][i]; // same for old_hc. making the new state state
                    }
                }


                /*
                 * Save current solution as new vector in the Krylov subspace
                 */
                old_c[nextrap-1][j] = c[j];
                old_hc[nextrap-1][j] = hc[j];

                /*
                 * Form residual and apply Davidson correction
                 */
                hc[j] -= real(l[root[j]])*c[j];
                c[j] = -hc[j];
                c[j].weight(D, real(l[root[j]]));
                c[j] /= sqrt(aquarius::abs(scalar(conj(c[j])*c[j])));
            }

            /*
             * Recompute the error and overlap elements for the last Krylov
             * vector (which was replace with the current solution).
             */
            for (int j = 0;j < nvec;j++)
            {
                for (int jj = 0;jj < nvec;jj++)
                {
                    dtype snew = 0;
                    dtype enew = 0;
                    for (int i = 0;i < nextrap;i++)
                    {
                        for (int k = 0;k < nvec;k++)
                        {
                            for (int ii = 0;ii < nextrap;ii++)
                            {
                                for (int kk = 0;kk < nvec;kk++)
                                {
                                    snew += s[k][i][kk][ii]*conj(vr[root[j]][k][i])*vr[root[jj]][kk][ii];
                                    enew += e[k][i][kk][ii]*conj(vr[root[j]][k][i])*vr[root[jj]][kk][ii];
                                }
                            }
                        }
                    }
                    s[j][nextrap-1][jj][nextrap-1] = snew;
                    e[j][nextrap-1][jj][nextrap-1] = enew;
                }
            }

            /*
             * ...and the cross elements with the other vectors.
             */
            for (int j = 0;j < nvec;j++)
            {
                for (int i = 0;i < nextrap-1;i++)
                {
                    for (int k = 0;k < nvec;k++)
                    {
                        dtype snew = 0;
                        dtype enew = 0;
                        for (int ii = 0;ii < nextrap;ii++)
                        {
                            for (int kk = 0;kk < nvec;kk++)
                            {
                                snew += s[k][i][kk][ii]*vr[root[j]][kk][ii];
                                enew += e[k][i][kk][ii]*vr[root[j]][kk][ii];
                            }
                        }
                        s[k][i][j][nextrap-1] = snew;
                        e[k][i][j][nextrap-1] = enew;

                        snew = 0;
                        enew = 0;
                        for (int ii = 0;ii < nextrap;ii++)
                        {
                            for (int kk = 0;kk < nvec;kk++)
                            {
                                snew += s[kk][ii][k][i]*conj(vr[root[j]][kk][ii]);
                                enew += e[kk][ii][k][i]*conj(vr[root[j]][kk][ii]);
                            }
                        }
                        s[j][nextrap-1][k][i] = snew;
                        e[j][nextrap-1][k][i] = enew;
                    }
                }
            }

            /*
             * In subsequent iterations, again replace the past solution vector
             * with the update vector (u_i = c_i+1 - c_i). Adjust error and
             * overlap matrices to match.
             */
            if (nextrap > 1)
            {
                for (int k = 0;k < nvec;k++)
                {
                    old_c[nextrap-2][k] -= old_c[nextrap-1][k];
                    old_hc[nextrap-2][k] -= old_hc[nextrap-1][k];

                    for (int kk = 0;kk < nvec;kk++)
                    {
                        for (int i = 0;i < nextrap-2;i++)
                        {
                            s[k][i][kk][nextrap-2] = s[k][i][kk][nextrap-2] -
                                                     s[k][i][kk][nextrap-1];
                            e[k][i][kk][nextrap-2] = e[k][i][kk][nextrap-2] -
                                                     e[k][i][kk][nextrap-1];

                            s[k][nextrap-2][kk][i] = s[k][nextrap-2][kk][i] -
                                                     s[k][nextrap-1][kk][i];
                            e[k][nextrap-2][kk][i] = e[k][nextrap-2][kk][i] -
                                                     e[k][nextrap-1][kk][i];
                        }

                        s[k][nextrap-2][kk][nextrap-2] = s[k][nextrap-2][kk][nextrap-2] -
                                                         s[k][nextrap-1][kk][nextrap-2] -
                                                         s[k][nextrap-2][kk][nextrap-1] +
                                                         s[k][nextrap-1][kk][nextrap-1];
                        e[k][nextrap-2][kk][nextrap-2] = e[k][nextrap-2][kk][nextrap-2] -
                                                         e[k][nextrap-1][kk][nextrap-2] -
                                                         e[k][nextrap-2][kk][nextrap-1] +
                                                         e[k][nextrap-1][kk][nextrap-1];

                        s[k][nextrap-1][kk][nextrap-2] = s[k][nextrap-1][kk][nextrap-2] -
                                                         s[k][nextrap-1][kk][nextrap-1];
                        e[k][nextrap-1][kk][nextrap-2] = e[k][nextrap-1][kk][nextrap-2] -
                                                         e[k][nextrap-1][kk][nextrap-1];

                        s[k][nextrap-2][kk][nextrap-1] = s[k][nextrap-2][kk][nextrap-1] -
                                                         s[k][nextrap-1][kk][nextrap-1];
                        e[k][nextrap-2][kk][nextrap-1] = e[k][nextrap-2][kk][nextrap-1] -
                                                         e[k][nextrap-1][kk][nextrap-1];
                    }
                }
            }

            /*
            for (int i = 0;i < nextrap;i++)
            {
                for (int k = 0;k < nvec;k++)
                {
                    for (int ii = 0;ii < nextrap;ii++)
                    {
                        for (int kk = 0;kk < nvec;kk++)
                        {
                            dtype tmps = scalar(conj(old_c[i][k])* old_c[ii][kk]);
                            dtype tmpe = scalar(conj(old_c[i][k])*old_hc[ii][kk]);
                            tmps -= s[k][i][kk][ii];
                            tmpe -= e[k][i][kk][ii];
                            printf("%d %d %d %d %18.15f\n", i, k, ii, kk, tmps);
                            printf("%d %d %d %d %18.15f\n", i, k, ii, kk, tmpe);
                        }
                    }
                }
            }
            */

            /*
             * If the subspace is full, eject the oldest vector.
             */
            if (nextrap == maxextrap)
            {
                rotate(old_c.begin(), old_c.begin()+1, old_c.end());
                rotate(old_hc.begin(), old_hc.begin()+1, old_hc.end());

                for (int i = 0;i < nextrap-1;i++)
                {
                    for (int k = 0;k < nvec;k++)
                    {
                        for (int ii = 0;ii < nextrap-1;ii++)
                        {
                            for (int kk = 0;kk < nvec;kk++)
                            {
                                s[k][i][kk][ii] = s[k][i+1][kk][ii+1];
                                e[k][i][kk][ii] = e[k][i+1][kk][ii+1];
                            }
                        }
                    }
                }

                nextrap--;
            }

            vector<dtype> myreturn(nvec);
            for (int i = 0; i < nvec; i++)
                myreturn[i] = real(l[root[i]]);

            return myreturn;
        }

        void getSolution(int j, T& s)
        {
            s = old_c[nextrap-1][j];
        }

        void clear()
        {
            nextrap = 0;
            for (int k = 0;k < nvec;k++)
            {
                lock[k] = false;
            }
        }

        void orthogonalizeToSolution()
        {
            nextrap--;

            for (int k = 0;k < nvec;k++)
            {
                lock[k] = false;
                for (int i = 0;i < nextrap;i++)
                {
                    dtype olap = scalar(conj(old_c[nextrap][k])*old_c[i][k]);
                    old_c[i][k] -= olap*old_c[nextrap][k];
                    old_hc[i][k] -= olap*old_hc[nextrap][k];
                    dtype nrm = sqrt(aquarius::abs(scalar(conj(old_c[i][k])*old_c[i][k])));
                    old_c[i][k] /= nrm;
                    old_hc[i][k] /= nrm;
                }
            }
        }
};

}
}

#endif
