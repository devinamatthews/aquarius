#ifndef _AQUARIUS_DAVIDSON_HPP_
#define _AQUARIUS_DAVIDSON_HPP_

#include "util/global.hpp"

#include "input/config.hpp"
#include "task/task.hpp"
#include "operator/denominator.hpp"

#include "diis.hpp"

namespace aquarius
{
namespace convergence
{

namespace detail
{

template <typename T>
struct DefaultWeight
{
    typedef typename T::dtype dtype;

    template <typename a_container>
    void operator()(a_container& a, const op::Denominator<dtype>& D, dtype omega) const
    {
        for (int j = 0;j < a.size();j++)
        {
            a[j].weight(D, omega);
        }
    }
};

}

template <typename T, typename InnerProd = detail::DefaultInnerProd<T>, typename Weight = detail::DefaultWeight<T>>
class Davidson : public task::Destructible
{
    private:
        Davidson(const Davidson& other);

        Davidson& operator=(const Davidson& other);

    protected:
        typedef typename T::dtype dtype;
        vector<dtype> soln_e;
        vector<vector<unique_vector<T>>> old_c; // hold all R[k][i] where i is over nvec and k is over maxextrap
        vector<vector<unique_vector<T>>> old_hc; // hold all H*R[i] = Z[i] at every iteration aka Z[k][i]
        vector<unique_vector<T>> guess;
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
        int nc;
        InnerProd innerProd;
        Weight weight;

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

        vector<vector<int>> getBestRoots(int n, const Arena& arena)
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
                            //if (n == 1) //printf("Root %d (%.12f) already picked\n", rt+1, real(l[rt]));
                            //    task::Logger::log(arena) << "Root " << rt+1 << " (" << real(l[rt]) << ") already picked" << endl;
                            found = true;
                        }
                    }

                    if (found) continue;

                    crit = aquarius::abs(real(l[rt])-soln_e[idx]);
                    if (crit < mincrit)// && aquarius::abs(imag(l[rt])) < 1e-12)
                    {
                        if (aquarius::abs(imag(l[rt])) > 1e-12)
                        {
                            task::Logger::log(arena) << "WARNING: Root is imaginary! (1)" << endl;
                        }
                        //if (n == 1) //printf("Solution %d (%.12f) matches root %d (%.12f)\n", idx+1, soln_e[idx], rt+1, real(l[rt]));
                        //    task::Logger::log(arena) << "Solution " << idx+1 << " (" << soln_e[idx] << ") matches root " << rt+1 << " (" << real(l[rt]) << ")" << endl;
                        mincrit = crit;
                        solns[idx] = rt;
                    }
                }

                assert(solns[idx] != -1);
                if (solns[idx] == -1)
                    task::Logger::log(arena) << "WARNING: No root selected! (1)" << endl;
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
                                    //if (n == 1) //printf("Root %d (%.12f) already picked\n", rt+1, real(l[rt]));
                                    //    task::Logger::log(arena) << "Root " << rt+1 << " (" << real(l[rt])<< ") already picked" << endl;
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
                                //if (n == 1) //printf("Root %d (%.12f) matches solution %d (%.12f)\n", rt+1, real(l[rt]), idx2, soln_e[idx2]);
                                //    task::Logger::log(arena) << "Root " << rt+1 << " (" << real(l[rt]) << ") matches solution " << idx2+1 << " (" << soln_e[idx2] << ")" << endl;
                                found = true;
                            }
                        }

                        /*
                         * Check against energies of converged roots
                         */
                        for (int idx2 = 0;idx2 < soln_e.size();idx2++)
                        {
                            if (aquarius::abs(real(l[rt])-soln_e[idx2]) < 1e-6)
                            {
                                //if (n == 1) //printf("Root %d (%.12f) matches solution %d (%.12f)\n", rt+1, real(l[rt]), idx2, soln_e[idx2]);
                                //    task::Logger::log(arena) << "Root " << rt+1 << " (" << real(l[rt]) << ") matches solution " << idx2+1 << " (" << soln_e[idx2] << ")" << endl;
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

                        if (crit < mincrit)// && aquarius::abs(imag(l[rt])) < 1e-12)
                        {
                            if (aquarius::abs(imag(l[rt])) > 1e-12)
                            {
                                task::Logger::log(arena) << "WARNING: Root is imaginary! (2)" << endl;
                            }
                            //if (n == 1) //printf("Root %d (%.12f) picked\n", rt+1, real(l[rt]));
                            //    task::Logger::log(arena) << "Root " << rt+1 << " (" << real(l[rt]) << ") picked" << endl;
                            mincrit = crit;
                            roots[idx][vec] = rt;
                        }
                        else
                        {
                            //if (n == 1) //printf("Root %d (%.12f) not picked\n", rt+1, real(l[rt]));
                            //    task::Logger::log(arena) << "Root " << rt+1 << " (" << real(l[rt]) << ") not picked" << endl;
                        }
                    }

                    assert(roots[idx][vec] != -1);
                    if (roots[idx][vec] == -1)
                        task::Logger::log(arena) << "WARNING: No root selected! (2)" << endl;
                }
            }

            return roots;
        }

        vector<int> getBestRoot(const Arena& arena)
        {
            return getBestRoots(1, arena)[0];
        }

        template <typename c_container, typename hc_container>
        void addVectors(c_container&& c, hc_container&& hc)
        {
            nextrap++;

            old_c.resize(max<int>(old_c.size(),nextrap));
            old_hc.resize(max<int>(old_hc.size(),nextrap));

            guess_overlap.resize(nvec, nvec, nextrap);
            s.resize(nvec, nextrap, nvec, nextrap);
            e.resize(nvec, nextrap, nvec, nextrap);
            l.resize(nvec*nextrap);
            vr.resize(nvec*nextrap, nvec, nextrap);

            for (int vec = 0;vec < nvec;vec++)
            {
                if (vec >= old_c[nextrap-1].size()) old_c[nextrap-1].emplace_back();
                if (vec >= old_hc[nextrap-1].size()) old_hc[nextrap-1].emplace_back();

                for (int idx = 0;idx < nc;idx++)
                {
                    if (idx >= old_c[nextrap-1][vec].size())
                        old_c[nextrap-1][vec].emplace_back(c[vec][idx]);
                    else
                        old_c[nextrap-1][vec][idx] = c[vec][idx];

                    if (idx >= old_hc[nextrap-1][vec].size())
                        old_hc[nextrap-1][vec].emplace_back(hc[vec][idx]);
                    else
                        old_hc[nextrap-1][vec][idx] = hc[vec][idx];
                }
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
                        guess_overlap[gvec][cvec][nextrap-1] = innerProd(old_c[nextrap-1][cvec], guess[gvec]);
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
                    e[lvec][nextrap-1][rvec][nextrap-1] = innerProd(old_c[nextrap-1][lvec], old_hc[nextrap-1][rvec]);
                    s[lvec][nextrap-1][rvec][nextrap-1] = innerProd(old_c[nextrap-1][lvec],  old_c[nextrap-1][rvec]);

                    for (int extrap = 0;extrap < nextrap;extrap++)
                    {
                        e[lvec][   extrap][rvec][nextrap-1] = innerProd(old_c[   extrap][lvec], old_hc[nextrap-1][rvec]);
                        e[lvec][nextrap-1][rvec][   extrap] = innerProd(old_c[nextrap-1][lvec], old_hc[   extrap][rvec]);
                        s[lvec][   extrap][rvec][nextrap-1] = innerProd(old_c[   extrap][lvec],  old_c[nextrap-1][rvec]);
                        s[lvec][nextrap-1][rvec][   extrap] = innerProd(old_c[nextrap-1][lvec],  old_c[   extrap][rvec]);
                    }
                }
            }
        }

        void getRoot(int rt, T& c, T& hc, bool normalize = true)
        {
            assert(nc == 1);
            getRoot(rt, ptr_vector<T>{&c}, ptr_vector<T>{&hc}, normalize);
        }

        void getRoot(int rt, T& c, bool normalize = true)
        {
            assert(nc == 1);
            getRoot(rt, ptr_vector<T>{&c}, normalize);
        }

        template <typename c_container, typename hc_container>
        void getRoot(int rt, c_container&& c, hc_container&& hc, bool normalize = true)
        {
            getRoot(rt, c, false);

            for (int idx = 0;idx < nc;idx++)
            {
                hc[idx] = 0;
                for (int extrap = nextrap-1;extrap >= 0;extrap--)
                {
                    for (int vec = nvec-1;vec >= 0;vec--)
                    {
                        hc[idx] += old_hc[extrap][vec][idx]*vr[rt][vec][extrap];
                    }
                }
            }

            if (normalize)
            {
                dtype nrm = sqrt(aquarius::abs(innerProd(c, c)));

                for (int idx = 0;idx < nc;idx++)
                {
                    c[idx] /= nrm;
                    hc[idx] /= nrm;
                }
            }
        }

        template <typename c_container>
        void getRoot(int rt, c_container&& c, bool normalize = true)
        {
            for (int idx = 0;idx < nc;idx++)
            {
                c[idx] = 0;
                for (int extrap = nextrap-1;extrap >= 0;extrap--)
                {
                    for (int vec = nvec-1;vec >= 0;vec--)
                    {
                        c[idx] += old_c[extrap][vec][idx]*vr[rt][vec][extrap];
                    }
                }
            }

            if (normalize)
            {
                dtype nrm = sqrt(aquarius::abs(innerProd(c, c)));
                for (int idx = 0;idx < nc;idx++) c[idx] /= nrm;
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
            assert(nc == 1);
            assert(nvec == 1);
            return extrapolate(ptr_vector<T>{&c}, ptr_vector<T>{&hc}, D)[0];
        }

        template <typename c_container, typename hc_container>
        enable_if_t<is_same<typename decay_t<c_container>::value_type, T>::value, vector<dtype>>
        extrapolate(c_container&& c, hc_container&& hc, const op::Denominator<dtype>& D)
        {
            assert((nc == 1 && nvec  > 1) ||
                   (nc  > 1 && nvec == 1) ||
                   c.size() == 1);

            vector<ptr_vector<T>> new_c, new_hc;

            if (nc == 1)
            {
                for (auto& t : c) new_c.push_back({&t});
                for (auto& t : hc) new_hc.push_back({&t});
            }
            else
            {
                new_c.emplace_back();
                new_hc.emplace_back();
                for (auto& t : c) new_c[0].push_back(&t);
                for (auto& t : hc) new_hc[0].push_back(&t);
            }

            return extrapolate(new_c, new_hc, D);
        }

        template <typename c_container, typename hc_container>
        enable_if_t<is_same<typename decay_t<c_container>::value_type::value_type, T>::value, vector<dtype>>
        extrapolate(c_container&& c, hc_container&& hc, const op::Denominator<dtype>& D)
        {
            using slice::all;

            assert(c.size() == nvec);
            assert(hc.size() == nvec);

            for (int vec = 0;vec < nvec;vec++)
            {
                assert(c[vec].size() == nc);
                assert(hc[vec].size() == nc);
            }

            /*
             * If the maximum size of the subspace has been reached, a smaller
             * subspace must be constructed which contains the best approximation
             * to the solution.
             */
            if (nextrap == maxextrap+nsoln)
            {
                task::Logger::log(c[0][0].arena) << "Compacting..." << endl;
                int new_nextrap = nsoln+nreduce;

                vector<vector<int>> roots = getBestRoots(nreduce,c[0][0].arena);

                vector<vector<unique_vector<T>>> new_c(nreduce);
                vector<vector<unique_vector<T>>> new_hc(nreduce);

                for (int extrap = 0; extrap < nreduce; extrap++)
                {
                    new_c[extrap].resize(nvec);
                    new_hc[extrap].resize(nvec);

                    for (int vec = 0;vec < nvec;vec++)
                    {
                        new_c[extrap][vec].assign(c[vec].begin(), c[vec].end());
                        new_hc[extrap][vec].assign(hc[vec].begin(), hc[vec].end());
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
            marray<dtype,4> e_tmp(e);
            marray<dtype,4> s_tmp(s);
            marray<dtype,3> vr_tmp(vr);

            int info = ggev('V', 'N', nextrap*nvec, e_tmp.data(), nextrap*nvec,
                        s_tmp.data(), nextrap*nvec, l.data(), beta.data(),
                        vr_tmp.data(), nextrap*nvec, NULL, 1);
            if (info != 0) throw runtime_error(str("davidson: Info in ggev: %d", info));

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
            root = getBestRoot(c[0][0].arena);

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
                    task::Logger::log(c[0][0].arena) << "Locking root " << (vec+1) << endl;
                    lock[vec] = true;
                    lock_e[vec] = real(l[root[vec]]);
                }
                else if (lock[vec] &&
                         aquarius::abs(lock_e[vec]-real(l[root[vec]])) > 1e-4)
                {
                    task::Logger::log(c[0][0].arena) << "Re-locking root " << (vec+1) << endl;
                    lock_e[vec] = real(l[root[vec]]);
                }
                previous[vec] = real(l[root[vec]]);
            }

            /*
             * Form the current solution
             */
            for (int vec = 0;vec < nvec;vec++)
            {
                getRoot(root[vec], c[vec], hc[vec]);
            }

            if (continuous)
            {
                /*
                 * Save current solution as new vector in the Krylov subspace
                 */
                nextrap--;
                addVectors(c, hc);
            }

            /*
             * Calculate residuals and apply Davidson correction
             */
            for (int vec = 0;vec < nvec;vec++)
            {
                /*
                 * Form residual and apply Davidson correction
                 */
                for (int idx = 0;idx < nc;idx++)
                {
                    hc[vec][idx] -= real(l[root[vec]])*c[vec][idx];
                    c[vec][idx] = -hc[vec][idx];
                }
                weight(c[vec], D, real(l[root[vec]]));

                /*
                 * Orthogonalize and normalize
                 */
                for (int soln = 0;soln < nsoln;soln++)
                {
                    for (int svec = 0;svec < nvec;svec++)
                    {
                        dtype olap = innerProd(old_c[soln][svec], c[vec]);
                        for (int idx = 0;idx < nc;idx++) c[vec][idx] -= old_c[soln][svec][idx]*olap;
                    }
                }
                dtype nrm = innerProd(c[vec], c[vec]);
                for (int idx = 0;idx < nc;idx++) c[vec][idx] /= sqrt(nrm);
            }

            if (continuous)
            {
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
                        for (int idx = 0;idx < nc;idx++)
                        {
                            old_c[nextrap-2][lvec][idx] -= old_c[nextrap-1][lvec][idx];
                            old_hc[nextrap-2][lvec][idx] -= old_hc[nextrap-1][lvec][idx];
                        }

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
            assert(nc == 1);

            if (continuous)
            {
                c = old_c[nextrap-1][j][0];
            }
            else
            {
                getRoot(root[j], c);
            }
        }

        template <typename c_container>
        void getSolution(int j, c_container&& c)
        {
            assert(nc == c.size());

            if (continuous)
            {
                for (int idx = 0;idx < nc;idx++)
                    c[idx] = old_c[nextrap-1][j][idx];
            }
            else
            {
                getRoot(root[j], c);
            }
        }

        void getSolution(int j, T& c, T& hc)
        {
            assert(nc == 1);

            if (continuous)
            {
                c = old_c[nextrap-1][j][0];
                hc = old_hc[nextrap-1][j][0];
            }
            else
            {
                getRoot(root[j], c, hc);
            }
        }

        template <typename c_container, typename hc_container>
        void getSolution(int j, c_container&& c, hc_container&& hc)
        {
            assert(nc == c.size());

            if (continuous)
            {
                for (int idx = 0;idx < nc;idx++)
                {
                    c[idx] = old_c[nextrap-1][j][idx];
                    hc[idx] = old_hc[nextrap-1][j][idx];
                }
            }
            else
            {
                getRoot(root[j], c, hc);
            }
        }

        void reset(int nvec = 1, int nc = 1, InnerProd innerProd = InnerProd(), Weight weight = Weight())
        {
            this->nc = nc;
            this->innerProd = innerProd;
            this->weight = weight;
            init(nvec, LOWEST_ENERGY);
        }

        void reset(dtype t, int nc = 1, InnerProd innerProd = InnerProd(), Weight weight = Weight())
        {
            this->nc = nc;
            this->innerProd = innerProd;
            this->weight = weight;
            init(1, CLOSEST_ENERGY);
            target[0] = t;
        }

        void reset(const vector<dtype>& t, int nc = 1, InnerProd innerProd = InnerProd(), Weight weight = Weight())
        {
            this->nc = nc;
            this->innerProd = innerProd;
            this->weight = weight;
            init(target.size(), CLOSEST_ENERGY);
            target = t;
        }

        void reset(const T& g, int nc = 1, InnerProd innerProd = InnerProd(), Weight weight = Weight())
        {
            this->nc = nc;
            this->innerProd = innerProd;
            this->weight = weight;
            assert(nc == 1);
            init(1, GUESS_OVERLAP);
            guess.push_back({g});
        }

        template <typename guess_container>
        enable_if_t<is_same<typename decay_t<guess_container>::value_type,T>::value>
        reset(const guess_container& gs, int nc = 1, InnerProd innerProd = InnerProd(), Weight weight = Weight())
        {
            this->nc = nc;
            this->innerProd = innerProd;
            this->weight = weight;

            assert((nc == 1 && nvec  > 1) ||
                   (nc  > 1 && nvec == 1));

            if (nc == 1)
            {
                init(gs.size(), GUESS_OVERLAP);
                guess.clear();
                guess.resize(gs.size());
                for (int vec = 0;vec < nvec;vec++) guess[vec].emplace_back(gs[vec]);
            }
            else
            {
                assert(nc == gs.size());
                init(1, GUESS_OVERLAP);
                guess.resize(1);
                guess[0].assign(gs.begin(), gs.end());
            }
        }

        template <typename guess_container>
        enable_if_t<is_same<typename decay_t<guess_container>::value_type::value_type,T>::value>
        reset(const guess_container& gs, int nc = 1, InnerProd innerProd = InnerProd(), Weight weight = Weight())
        {
            this->nc = nc;
            this->innerProd = innerProd;
            this->weight = weight;

            assert(gs.size() == nvec);
            for (int vec = 0;vec < nvec;vec++)
            {
                assert(gs[vec].size() == nc);
            }

            init(gs.size(), GUESS_OVERLAP);
            guess.resize(gs.size());
            for (int vec = 0;vec < nvec;vec++) guess[vec].assign(gs[vec].begin(), gs[vec].end());
        }

        template <typename... Args>
        void nextRoot(Args&&... args)
        {
            vector<unique_vector<T>> c(nvec);
            vector<unique_vector<T>> hc(nvec);

            for (int vec = 0;vec < nvec;vec++)
            {
                c[vec].assign(old_c[0][vec].begin(), old_c[0][vec].end());
                hc[vec].assign(old_hc[0][vec].begin(), old_hc[0][vec].end());
                getRoot(root[vec], c[vec], hc[vec]);
            }

            nextrap = nsoln;
            addVectors(c, hc);

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
