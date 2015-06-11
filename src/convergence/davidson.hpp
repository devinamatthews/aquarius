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
    return aquarius::abs(a) > aquarius::abs(b);
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
        vector<dtype> soln_e;
        vector<unique_vector<T>> old_c; // hold all R[k][i] where i is over nvec and k is over maxextrap
        vector<unique_vector<T>> old_hc; // hold all H*R[i] = Z[i] at every iteration aka Z[k][i]
        unique_vector<T> guess;
        marray<dtype,3> guess_overlap;
        marray<dtype,4> s, e;
        int nvec, maxextrap, nreduce, nextrap, nsoln; // number of energies, number of iterations
        vector<int> mode;
        vector<dtype> target;
        vector<dtype> previous;
        vector<int> root;
        vector<complex_type_t<dtype>> l;
        marray<dtype,3> vr;
        vector<bool> lock;
        vector<dtype> lock_e;
        bool continuous;

        enum {GUESS_OVERLAP, LOWEST_ENERGY, CLOSEST_ENERGY};

        void parse(const input::Config& config)
        {
            nextrap = 0;
            nsoln = 0;
            maxextrap = config.get<int>("order");
            nreduce = config.get<int>("num_reduce");
            continuous = config.get<string>("compaction") == "continuous";
        }

        void init(int nvec_, int mode_)
        {
            nextrap = nsoln;
            nvec = nvec_;
            mode.assign(nvec, mode_);
            target.resize(nvec);
            previous.assign(nvec, 1e9);
            root.resize(nvec);
            lock.assign(nvec, false);
            lock_e.resize(nvec);

            guess_overlap.resize(nvec, nvec, nextrap);
            s.resize(nvec, nextrap, nvec, nextrap);
            e.resize(nvec, nextrap, nvec, nextrap);
            l.resize(nvec*nextrap);
            vr.resize(nvec*nextrap, nvec, nextrap);
        }

        vector<vector<int>> getBestRoots(int n)
        {
            vector<vector<int>> roots(n, vector<int>(nvec, -1));
            vector<int> solns(nsoln*nvec, -1);

            for (int idx = 0;idx < nsoln*nvec;idx++)
            {
                dtype crit = numeric_limits<dtype>::max();
                dtype mincrit = numeric_limits<dtype>::max();

                for (int rt = 0;rt < nextrap*nvec;rt++)
                {
                    bool found = false;

                    /*
                     * Make sure that this root hasn't already been selected
                     */
                    for (int idx2 = 0;idx2 < nvec*nsoln;idx2++)
                    {
                        if (solns[idx2] == rt)
                        {
                            //if (n == 1) printf("Root %d (%.12f) already picked\n", rt+1, real(l[rt]));
                            found = true;
                        }
                    }

                    if (found) continue;

                    crit = aquarius::abs(real(l[rt])-soln_e[idx]);
                    if (crit < mincrit && aquarius::abs(imag(l[rt])) < 1e-12)
                    {
                        //if (n == 1) printf("Solution %d (%.12f) matches root %d (%.12f)\n", idx+1, soln_e[idx], rt+1, real(l[rt]));
                        mincrit = crit;
                        solns[idx] = rt;
                    }
                }

                assert(solns[idx] != -1);
            }

            for (int idx = 0;idx < n;idx++)
            {
                for (int vec = 0;vec < nvec;vec++)
                {
                    dtype crit = numeric_limits<dtype>::max();
                    dtype mincrit = numeric_limits<dtype>::max();

                    for (int rt = 0;rt < nextrap*nvec;rt++)
                    {
                        bool found = false;

                        /*
                         * Make sure that this root hasn't already been selected
                         */
                        for (int idx2 = 0;idx2 < n;idx2++)
                        {
                            for (int vec2 = 0;vec2 < nvec;vec2++)
                            {
                                if (roots[idx2][vec2] == rt)
                                {
                                    //if (n == 1) printf("Root %d (%.12f) already picked\n", rt+1, real(l[rt]));
                                    found = true;
                                }
                            }
                        }

                        /*
                         * And that it is not a previously converged root
                         */
                        for (int idx2 = 0;idx2 < nvec*nsoln;idx2++)
                        {
                            if (solns[idx2] == rt)
                            {
                                //if (n == 1) printf("Root %d (%.12f) matches solution %d (%.12f)\n", rt+1, real(l[rt]), idx2, soln_e[idx2]);
                                found = true;
                            }
                        }

                        if (found) continue;

                        if (lock[vec])
                        {
                            crit = aquarius::abs(real(l[rt])-lock_e[vec]);
                        }
                        else if (mode[vec] == GUESS_OVERLAP)
                        {
                            crit = 0.0;
                            for (int m = 0;m < nextrap;m++)
                                for (int k = 0;k < nvec;k++)
                                    crit -= vr[rt][k][m]*guess_overlap[vec][k][m];
                        }
                        else if (mode[vec] == LOWEST_ENERGY)
                        {
                            crit = real(l[rt]);
                        }
                        else if (mode[vec] == CLOSEST_ENERGY)
                        {
                            crit = aquarius::abs(real(l[rt])-target[vec]);
                        }

                        if (crit < mincrit && aquarius::abs(imag(l[rt])) < 1e-12)
                        {
                            //if (n == 1) printf("Root %d (%.12f) picked\n", rt+1, real(l[rt]));
                            mincrit = crit;
                            roots[idx][vec] = rt;
                        }
                        else
                        {
                            //if (n == 1) printf("Root %d (%.12f) not picked\n", rt+1, real(l[rt]));
                        }
                    }

                    assert(roots[idx][vec] != -1);
                }
            }

            return roots;
        }

        vector<int> getBestRoot()
        {
            return getBestRoots(1)[0];
        }

        template <typename c_container, typename hc_container>
        void addVectors(c_container&& c, hc_container&& hc)
        {
            nextrap++;

            old_c.resize(max((int)old_c.size(),nextrap));
            old_hc.resize(max((int)old_hc.size(),nextrap));

            guess_overlap.resize(nvec, nvec, nextrap);
            s.resize(nvec, nextrap, nvec, nextrap);
            e.resize(nvec, nextrap, nvec, nextrap);
            l.resize(nvec*nextrap);
            vr.resize(nvec*nextrap, nvec, nextrap);

            for (int vec = 0;vec < nvec;vec++)
            {
                if (vec >= old_c[nextrap-1].size())
                    old_c[nextrap-1].emplace_back(c[vec]);
                else
                    old_c[nextrap-1][vec] = c[vec];

                if (vec >= old_hc[nextrap-1].size())
                    old_hc[nextrap-1].emplace_back(hc[vec]);
                else
                    old_hc[nextrap-1][vec] = hc[vec];
            }

            /*
             * Compute the overlap with the guess vectors
             */
            if (!guess.empty())
            {
                for (int cvec = 0;cvec < nvec;cvec++)
                {
                    for (int gvec = 0;gvec < nvec;gvec++)
                    {
                        guess_overlap[gvec][cvec][nextrap-1] = aquarius::abs(scalar(conj(c[cvec])*guess[gvec]));
                    }
                }
            }

            /*
             * Augment the subspace matrix with the new vectors
             */
            for (int lvec = 0;lvec < nvec;lvec++)
            {
                for (int rvec = 0;rvec < nvec;rvec++)
                {
                    e[lvec][nextrap-1][rvec][nextrap-1] = scalar(conj(c[lvec])*hc[rvec]);
                    s[lvec][nextrap-1][rvec][nextrap-1] = scalar(conj(c[lvec])* c[rvec]);

                    for (int extrap = 0;extrap < nextrap;extrap++)
                    {
                        e[lvec][   extrap][rvec][nextrap-1] = scalar(conj(old_c[extrap][lvec])*    hc        [rvec]);
                        e[lvec][nextrap-1][rvec][   extrap] = scalar(conj(    c        [lvec])*old_hc[extrap][rvec]);
                        s[lvec][   extrap][rvec][nextrap-1] = scalar(conj(old_c[extrap][lvec])*     c        [rvec]);
                        s[lvec][nextrap-1][rvec][   extrap] = scalar(conj(    c        [lvec])* old_c[extrap][rvec]);
                    }
                }
            }
        }

        void getRoot(int rt, T& c, T& hc)
        {
            getRoot(rt, c);

            hc = 0;
            for (int extrap = nextrap-1;extrap >= 0;extrap--)
            {
                for (int vec = nvec-1;vec >= 0;vec--)
                {
                    hc += old_hc[extrap][vec]*vr[rt][vec][extrap];
                }
            }
        }

        void getRoot(int rt, T& c)
        {
            c = 0;
            for (int extrap = nextrap-1;extrap >= 0;extrap--)
            {
                for (int vec = nvec-1;vec >= 0;vec--)
                {
                    c += old_c[extrap][vec]*vr[rt][vec][extrap];
                }
            }
        }

    public:
        template <typename... U>
        Davidson(const input::Config& config, U&&... args)
        {
            parse(config);
            reset(forward<U>(args)...);
        }

        dtype extrapolate(T& c, T& hc, const op::Denominator<dtype>& D)
        {
            assert(nvec == 1);
            return extrapolate(ptr_vector<T>{&c}, ptr_vector<T>{&hc}, D)[0];
        }

        template <typename c_container, typename hc_container>
        vector<dtype> extrapolate(c_container&& c, hc_container&& hc, const op::Denominator<dtype>& D)
        {
            using slice::all;

            assert(nvec == c.size() && nvec == hc.size());
	    
            assert(c.size() == nvec);
            assert(hc.size() == nvec);

            /*
             * Check and normalize incoming c and H*c vectors
             */
            for (int vec = 0;vec < nvec;vec++)
            {
                for (int soln = 0;soln < nsoln;soln++)
                {
                    for (int svec = 0;svec < nvec;svec++)
                    {
                        dtype olap = aquarius::abs(scalar(conj(old_c[soln][svec])*c[vec]));
                        c[vec] -= olap*old_c[soln][svec];
                        hc[vec] -= olap*old_hc[soln][svec];
                    }
                }

                dtype norm = sqrt(aquarius::abs(scalar(conj(c[vec])*c[vec])));
                c[vec] /= norm;
                hc[vec] /= norm;
            }

            /*
             * If the maximum size of the subspace has been reached, a smaller
             * subspace must be constructed which contains the best approximation
             * to the solution.
             */
            if (nextrap == maxextrap+nsoln)
            {
                task::Logger::log(c[0].arena) << "Compacting..." << endl;
                int new_nextrap = max(nsoln+1, nextrap-nreduce);

                vector<vector<int>> roots = getBestRoots(new_nextrap-nsoln);

                vector<unique_vector<T>> new_c(new_nextrap-nsoln);
                vector<unique_vector<T>> new_hc(new_nextrap-nsoln);

                for (int extrap = 0; extrap < new_nextrap-nsoln; extrap++)
                {
                    new_c[extrap].assign(c.begin(), c.end());
                    new_hc[extrap].assign(hc.begin(), hc.end());

                    for (int vec = 0;vec < nvec;vec++)
                    {
                        getRoot(roots[extrap][vec], new_c[extrap][vec], new_hc[extrap][vec]);
                    }
                }

                nextrap = nsoln;
                for (int extrap = 0; extrap < new_nextrap-nsoln; extrap++)
                {
                    addVectors(new_c[extrap], new_hc[extrap]);
                }
            }

            addVectors(c, hc);

            /*
             * Diagonalize the subspace matrix to obtain approximate solutions
             */
            vector<dtype> beta(nvec*nextrap);
            marray<dtype,4> e_tmp = copy(e);
            marray<dtype,4> s_tmp = copy(s);
            marray<dtype,3> vr_tmp = copy(vr);

            int info = ggev('V', 'N', nextrap*nvec, e_tmp.data(), nextrap*nvec,
                        s_tmp.data(), nextrap*nvec, l.data(), beta.data(),
                        vr_tmp.data(), nextrap*nvec, NULL, 1);
            if (info != 0) throw runtime_error(strprintf("davidson: Info in ggev: %d", info));

            vr = vr_tmp;

            /*
             * Fix sign of eigenvectors
             */
            for (int rt = 0;rt < nextrap*nvec;rt++)
            {
                l[rt] /= beta[rt];

                int vec = rt/nextrap;
                if (real(vr[rt][vec][nextrap-1]) < 0)
                {
                    scal(nextrap*nvec, -1, vr[rt].data(), 1);
                }
            }

            /*
             * Assign eigenvalues (exclusively) to states by the selected criterion
             */
            root = getBestRoot();

            /*
             * Check proximity to previous root and lock on if within tolerance
             */
            for (int vec = 0; vec < nvec; vec++)
            {
                if (nextrap > 1 &&
                    mode[vec] != CLOSEST_ENERGY &&
                    !lock[vec] &&
                    aquarius::abs(previous[vec]-real(l[root[vec]])) < 1e-4)
                {
                    task::Logger::log(c[0].arena) << "Locking root " << (vec+1) << endl;
                    lock[vec] = true;
                    lock_e[vec] = real(l[root[vec]]);
                }
                previous[vec] = real(l[root[vec]]);
            }

            /*
             * Calculate residuals and apply Davidson correction
             */
            for (int vec = 0;vec < nvec;vec++)
            {
                /*
                 * Form current trial vector y and H*y = x
                 */
                getRoot(root[vec], c[vec], hc[vec]);

                if (continuous)
                {
                    /*
                     * Save current solution as new vector in the Krylov subspace
                     */
                    old_c[nextrap-1][vec] = c[vec];
                    old_hc[nextrap-1][vec] = hc[vec];
                }

                /*
                 * Form residual and apply Davidson correction
                 */
                hc[vec] -= real(l[root[vec]])*c[vec];
                c[vec] = -hc[vec];
                c[vec].weight(D, real(l[root[vec]]));
                c[vec] /= sqrt(aquarius::abs(scalar(conj(c[vec])*c[vec])));
            }

            if (continuous)
            {
                /*
                 * Recompute the error and overlap elements for the last Krylov
                 * vector (which was replaced with the current solution).
                 */
                for (int lvec = 0;lvec < nvec;lvec++)
                {
                    //TODO: guess_overlap

                    for (int rvec = 0;rvec < nvec;rvec++)
                    {
                        dtype snew = 0;
                        dtype enew = 0;
                        for (int i = 0;i < nvec;i++)
                        {
                            for (int j = 0;j < nextrap;j++)
                            {
                                for (int k = 0;k < nvec;k++)
                                {
                                    for (int l = 0;l < nextrap;l++)
                                    {
                                        snew += conj(vr[root[lvec]][i][j])*s[i][j][k][l]*vr[root[rvec]][k][l];
                                        enew += conj(vr[root[lvec]][i][j])*e[i][j][k][l]*vr[root[rvec]][k][l];
                                    }
                                }
                            }
                        }
                        s[lvec][nextrap-1][rvec][nextrap-1] = snew;
                        e[lvec][nextrap-1][rvec][nextrap-1] = enew;
                    }
                }

                /*
                 * ...and the cross elements with the other vectors.
                 */
                for (int vec = 0;vec < nvec;vec++)
                {
                    //TODO: guess_overlap

                    for (int i = 0;i < nvec;i++)
                    {
                        for (int j = 0;j < nextrap-1;j++)
                        {
                            dtype snew = 0;
                            dtype enew = 0;
                            for (int k = 0;k < nvec;k++)
                            {
                                for (int l = 0;l < nextrap;l++)
                                {
                                    snew += s[i][j][k][l]*vr[root[vec]][k][l];
                                    enew += e[i][j][k][l]*vr[root[vec]][k][l];
                                }
                            }
                            s[i][j][vec][nextrap-1] = snew;
                            e[i][j][vec][nextrap-1] = enew;

                            snew = 0;
                            enew = 0;
                            for (int k = 0;k < nvec;k++)
                            {
                                for (int l = 0;l < nextrap;l++)
                                {
                                    snew += conj(vr[root[vec]][k][l])*s[k][l][i][j];
                                    enew += conj(vr[root[vec]][k][l])*e[k][l][i][j];
                                }
                            }
                            s[vec][nextrap-1][i][j] = snew;
                            e[vec][nextrap-1][i][j] = enew;
                        }
                    }
                }

                /*
                 * In subsequent iterations, again replace the past solution vector
                 * with the update vector (u_i = c_i+1 - c_i). Adjust error and
                 * overlap matrices to match.
                 */
                if (nextrap > nsoln+1)
                {
                    //TODO: guess_overlap

                    for (int lvec = 0;lvec < nvec;lvec++)
                    {
                        old_c[nextrap-2][lvec] -= old_c[nextrap-1][lvec];
                        old_hc[nextrap-2][lvec] -= old_hc[nextrap-1][lvec];

                        for (int rvec = 0;rvec < nvec;rvec++)
                        {
                            for (int extrap = 0;extrap < nextrap-2;extrap++)
                            {
                                s[lvec][extrap][rvec][nextrap-2] = s[lvec][extrap][rvec][nextrap-2] -
                                                                   s[lvec][extrap][rvec][nextrap-1];
                                e[lvec][extrap][rvec][nextrap-2] = e[lvec][extrap][rvec][nextrap-2] -
                                                                   e[lvec][extrap][rvec][nextrap-1];

                                s[lvec][nextrap-2][rvec][extrap] = s[lvec][nextrap-2][rvec][extrap] -
                                                                   s[lvec][nextrap-1][rvec][extrap];
                                e[lvec][nextrap-2][rvec][extrap] = e[lvec][nextrap-2][rvec][extrap] -
                                                                   e[lvec][nextrap-1][rvec][extrap];
                            }

                            s[lvec][nextrap-2][rvec][nextrap-2] = s[lvec][nextrap-2][rvec][nextrap-2] -
                                                                  s[lvec][nextrap-1][rvec][nextrap-2] -
                                                                  s[lvec][nextrap-2][rvec][nextrap-1] +
                                                                  s[lvec][nextrap-1][rvec][nextrap-1];
                            e[lvec][nextrap-2][rvec][nextrap-2] = e[lvec][nextrap-2][rvec][nextrap-2] -
                                                                  e[lvec][nextrap-1][rvec][nextrap-2] -
                                                                  e[lvec][nextrap-2][rvec][nextrap-1] +
                                                                  e[lvec][nextrap-1][rvec][nextrap-1];

                            s[lvec][nextrap-1][rvec][nextrap-2] = s[lvec][nextrap-1][rvec][nextrap-2] -
                                                                  s[lvec][nextrap-1][rvec][nextrap-1];
                            e[lvec][nextrap-1][rvec][nextrap-2] = e[lvec][nextrap-1][rvec][nextrap-2] -
                                                                  e[lvec][nextrap-1][rvec][nextrap-1];

                            s[lvec][nextrap-2][rvec][nextrap-1] = s[lvec][nextrap-2][rvec][nextrap-1] -
                                                                  s[lvec][nextrap-1][rvec][nextrap-1];
                            e[lvec][nextrap-2][rvec][nextrap-1] = e[lvec][nextrap-2][rvec][nextrap-1] -
                                                                  e[lvec][nextrap-1][rvec][nextrap-1];
                        }
                    }
                }

                /*
                 * If the subspace is full, eject the oldest vector.
                 */
                if (nextrap == maxextrap+nsoln)
                {
                    rotate(old_c.begin()+nsoln, old_c.begin()+nsoln+1, old_c.end());
                    rotate(old_hc.begin()+nsoln, old_hc.begin()+nsoln+1, old_hc.end());
                    e[all][range(nsoln,nextrap)][all][all].rotate(0,1,0,0);
                    s[all][range(nsoln,nextrap)][all][all].rotate(0,1,0,0);
                    e[all][all][all][range(nsoln,nextrap)].rotate(0,0,0,1);
                    s[all][all][all][range(nsoln,nextrap)].rotate(0,0,0,1);

                    //TODO: guess_overlap

                    nextrap--;
                }
            }

            vector<dtype> myreturn(nvec);
            for (int i = 0; i < nvec; i++)
                myreturn[i] = real(l[root[i]]);

            return myreturn;
        }

        void getSolution(int j, T& c)
        {
            if (continuous)
            {
                c = old_c[nextrap-1][j];
            }
            else
            {
                getRoot(root[j], c);
            }
        }

        void reset(int nvec = 1)
        {
            init(nvec, LOWEST_ENERGY);
        }

        void reset(dtype t)
        {
            init(1, CLOSEST_ENERGY);
            target[0] = t;
        }

        void reset(const vector<dtype>& t)
        {
            init(target.size(), CLOSEST_ENERGY);
            target = t;
        }

        void reset(const T& g)
        {
            init(1, GUESS_OVERLAP);
            guess.push_back(g);
        }

        template <typename guess_container>
        enable_if_t<is_same<typename guess_container::value_type,T>::value>
        reset(const guess_container& gs)
        {
            init(guess.size(), GUESS_OVERLAP);
            for (auto& g : gs) guess.push_back(g);
        }

        template <typename... Args>
        void nextRoot(Args&&... args)
        {
            //if (!continuous)
            //{
                unique_vector<T> c(old_c[0].begin(), old_c[0].end());
                unique_vector<T> hc(old_hc[0].begin(), old_hc[0].end());

                for (int vec = 0;vec < nvec;vec++)
                {
                    getRoot(root[vec], c[vec], hc[vec]);
                }

                nextrap = nsoln;
                addVectors(c, hc);
            //}

            /*
            using slice::all;

            rotate(old_c.begin()+nsoln, old_c.begin()+nextrap-1, old_c.end());
            rotate(old_hc.begin()+nsoln, old_hc.begin()+nextrap-1, old_hc.end());
            e[all][range(nsoln,nextrap)][all][all].rotate(0,-1,0,0);
            s[all][range(nsoln,nextrap)][all][all].rotate(0,-1,0,0);
            e[all][all][all][range(nsoln,nextrap)].rotate(0,0,0,-1);
            s[all][all][all][range(nsoln,nextrap)].rotate(0,0,0,-1);
            */

            nsoln++;
            for (int vec = 0;vec < nvec;vec++)
            {
                soln_e.push_back(real(l[root[vec]]));
            }

            reset(forward<Args>(args)...);
        }
};

}
}

#endif
