#ifndef _AQUARIUS_FRAMEWORKS_CONVERGENCE_DAVIDSON_HPP_
#define _AQUARIUS_FRAMEWORKS_CONVERGENCE_DAVIDSON_HPP_

#include "frameworks/util.hpp"
#include "frameworks/task.hpp"
#include "frameworks/logging.hpp"

#include "diis.hpp"

namespace aquarius
{
namespace convergence
{

struct Weight
{
    virtual ~Weight() {}

    virtual void operator()(vector<tensor::Tensor<>>& a, double omega) const = 0;
};

class Davidson
{
    protected:
        vector<double> soln_e;
        vector<vector<vector<tensor::Tensor<>>>> old_c; // hold all R[k][i] where i is over nvec and k is over maxextrap
        vector<vector<vector<tensor::Tensor<>>>> old_hc; // hold all H*R[i] = Z[i] at every iteration aka Z[k][i]
        vector<vector<tensor::ConstTensor<>>> guess;
        marray<double,3> guess_overlap;
        marray<double,4> s, e;
        int nvec, maxextrap, nreduce, nextrap, nsoln; // number of energies, number of iterations
        vector<int> mode;
        vector<double> target;
        vector<double> previous;
        vector<int> root;
        vector<dcomplex> l;
        marray<double,3> vr;
        vector<bool> lock;
        vector<double> lock_e;
        bool continuous;
        int nc;
        unique_ptr<InnerProd> innerProd;
        unique_ptr<Weight> weight;

        enum {GUESS_OVERLAP, LOWEST_ENERGY, CLOSEST_ENERGY};

        void parse(const task::Config& config)
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
                double crit = numeric_limits<double>::max();
                double mincrit = numeric_limits<double>::max();

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

                    crit = std::abs(real(l[rt])-soln_e[idx]);
                    if (crit < mincrit)// && aquarius::abs(imag(l[rt])) < 1e-12)
                    {
                        if (std::abs(imag(l[rt])) > 1e-12)
                        {
                            logging::Logger::log(arena()) << "WARNING: Root is imaginary! (1)" << endl;
                        }
                        //if (n == 1) //printf("Solution %d (%.12f) matches root %d (%.12f)\n", idx+1, soln_e[idx], rt+1, real(l[rt]));
                        //    task::Logger::log(arena) << "Solution " << idx+1 << " (" << soln_e[idx] << ") matches root " << rt+1 << " (" << real(l[rt]) << ")" << endl;
                        mincrit = crit;
                        solns[idx] = rt;
                    }
                }

                assert(solns[idx] != -1);
                if (solns[idx] == -1)
                    logging::Logger::log(arena()) << "WARNING: No root selected! (1)" << endl;
            }

            for (int idx = 0;idx < n;idx++)
            {
                for (int vec = 0;vec < nvec;vec++)
                {
                    double crit = numeric_limits<double>::max();
                    double mincrit = numeric_limits<double>::max();

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
                            if (std::abs(real(l[rt])-soln_e[idx2]) < 1e-6)
                            {
                                //if (n == 1) //printf("Root %d (%.12f) matches solution %d (%.12f)\n", rt+1, real(l[rt]), idx2, soln_e[idx2]);
                                //    task::Logger::log(arena) << "Root " << rt+1 << " (" << real(l[rt]) << ") matches solution " << idx2+1 << " (" << soln_e[idx2] << ")" << endl;
                                found = true;
                            }
                        }

                        if (found) continue;

                        if (lock[vec])
                        {
                            crit = std::abs(real(l[rt])-lock_e[vec]);
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
                            crit = std::abs(real(l[rt])-target[vec]);
                        }

                        if (crit < mincrit)// && aquarius::abs(imag(l[rt])) < 1e-12)
                        {
                            if (std::abs(imag(l[rt])) > 1e-12)
                            {
                                logging::Logger::log(arena()) << "WARNING: Root is imaginary! (2)" << endl;
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
                        logging::Logger::log(arena()) << "WARNING: No root selected! (2)" << endl;
                }
            }

            return roots;
        }

        vector<int> getBestRoot()
        {
            return getBestRoots(1)[0];
        }

        void addVectors(const vector<vector<tensor::Tensor<>>>& c,
                        const vector<vector<tensor::Tensor<>>>& hc)
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
                        old_c[nextrap-1][vec].push_back(c[vec][idx].construct());
                    old_c[nextrap-1][vec][idx] = c[vec][idx];

                    if (idx >= old_hc[nextrap-1][vec].size())
                        old_hc[nextrap-1][vec].push_back(hc[vec][idx].construct());
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
                        guess_overlap[gvec][cvec][nextrap-1] = (*innerProd)(old_c[nextrap-1][cvec], guess[gvec]);
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
                    e[lvec][nextrap-1][rvec][nextrap-1] = (*innerProd)(old_c[nextrap-1][lvec], old_hc[nextrap-1][rvec]);
                    s[lvec][nextrap-1][rvec][nextrap-1] = (*innerProd)(old_c[nextrap-1][lvec],  old_c[nextrap-1][rvec]);

                    for (int extrap = 0;extrap < nextrap;extrap++)
                    {
                        e[lvec][   extrap][rvec][nextrap-1] = (*innerProd)(old_c[   extrap][lvec], old_hc[nextrap-1][rvec]);
                        e[lvec][nextrap-1][rvec][   extrap] = (*innerProd)(old_c[nextrap-1][lvec], old_hc[   extrap][rvec]);
                        s[lvec][   extrap][rvec][nextrap-1] = (*innerProd)(old_c[   extrap][lvec],  old_c[nextrap-1][rvec]);
                        s[lvec][nextrap-1][rvec][   extrap] = (*innerProd)(old_c[nextrap-1][lvec],  old_c[   extrap][rvec]);
                    }
                }
            }
        }

        void getRoot(int rt, tensor::Tensor<> c, tensor::Tensor<> hc, bool normalize = true)
        {
            assert(nc == 1);
            getRoot(rt, make_vector(c), make_vector(hc), normalize);
        }

        void getRoot(int rt, tensor::Tensor<> c, bool normalize = true)
        {
            assert(nc == 1);
            getRoot(rt, make_vector(c), normalize);
        }

        void getRoot(int rt, vector<tensor::Tensor<>>&& c, vector<tensor::Tensor<>>&& hc, bool normalize = true)
        {
            getRoot(rt, c, hc, normalize);
        }

        void getRoot(int rt, vector<tensor::Tensor<>>& c, vector<tensor::Tensor<>>& hc, bool normalize = true)
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
                double nrm = sqrt(std::abs((*innerProd)(c, c)));

                for (int idx = 0;idx < nc;idx++)
                {
                    c[idx] /= nrm;
                    hc[idx] /= nrm;
                }
            }
        }

        void getRoot(int rt, vector<tensor::Tensor<>>&& c, bool normalize = true)
        {
            getRoot(rt, c, normalize);
        }

        void getRoot(int rt, vector<tensor::Tensor<>>& c, bool normalize = true)
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
                double nrm = sqrt(std::abs((*innerProd)(c, c)));
                for (int idx = 0;idx < nc;idx++) c[idx] /= nrm;
            }
        }

    public:
        template <typename... U>
        Davidson(const task::Config& config, U&&... args)
        {
            parse(config);
            reset(forward<U>(args)...);
        }

        Davidson(const Davidson& other) = delete;

        Davidson& operator=(const Davidson& other) = delete;

        double extrapolate(tensor::Tensor<> c, tensor::Tensor<> hc)
        {
            assert(nc == 1);
            assert(nvec == 1);
            return extrapolate(make_vector(c), make_vector(hc))[0];
        }

        vector<double> extrapolate(vector<tensor::Tensor<>>&& c, vector<tensor::Tensor<>>&& hc)
        {
            return extrapolate(c, hc);
        }

        vector<double> extrapolate(vector<tensor::Tensor<>>& c, vector<tensor::Tensor<>>& hc)
        {
            assert((nc == 1 && nvec  > 1) ||
                   (nc  > 1 && nvec == 1));

            vector<vector<tensor::Tensor<>>> new_c, new_hc;

            if (nc == 1)
            {
                assert(c.size() == nvec);
                assert(hc.size() == nvec);
                for (auto& t : c) new_c.push_back({t});
                for (auto& t : hc) new_hc.push_back({t});
            }
            else
            {
                assert(c.size() == nc);
                assert(hc.size() == nc);
                new_c.emplace_back();
                new_hc.emplace_back();
                for (auto& t : c) new_c[0].push_back(t);
                for (auto& t : hc) new_hc[0].push_back(t);
            }

            return extrapolate(new_c, new_hc);
        }

        vector<double> extrapolate(vector<vector<tensor::Tensor<>>>&& c, vector<vector<tensor::Tensor<>>>&& hc)
        {
            return extrapolate(c, hc);
        }

        vector<double> extrapolate(vector<vector<tensor::Tensor<>>>& c, vector<vector<tensor::Tensor<>>>& hc)
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
                logging::Logger::log(arena()) << "Compacting..." << endl;
                int new_nextrap = nsoln+nreduce;

                vector<vector<int>> roots = getBestRoots(nreduce);

                vector<vector<vector<tensor::Tensor<>>>> new_c(nreduce);
                vector<vector<vector<tensor::Tensor<>>>> new_hc(nreduce);

                for (int extrap = 0; extrap < nreduce; extrap++)
                {
                    new_c[extrap].resize(nvec);
                    new_hc[extrap].resize(nvec);

                    for (int vec = 0;vec < nvec;vec++)
                    {
                        for (int i = 0;i < nc;i++)
                        {
                            new_c[extrap][vec].push_back(c[vec][i].construct());
                            new_hc[extrap][vec].push_back(hc[vec][i].construct());
                        }

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
            vector<double> beta(nvec*nextrap);
            marray<double,4> e_tmp = e;
            marray<double,4> s_tmp = s;
            marray<double,3> vr_tmp = vr;

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
            root = getBestRoot();

            /*
             * Check proximity to previous root and lock on if within tolerance
             */
            for (int vec = 0; vec < nvec; vec++)
            {
                if (nextrap > 1 &&
                    mode[vec] != CLOSEST_ENERGY &&
                    !lock[vec] &&
                    std::abs(previous[vec]-real(l[root[vec]])) < 1e-4)
                {
                    logging::Logger::log(arena()) << "Locking root " << (vec+1) << endl;
                    lock[vec] = true;
                    lock_e[vec] = real(l[root[vec]]);
                }
                else if (lock[vec] &&
                         std::abs(lock_e[vec]-real(l[root[vec]])) > 1e-4)
                {
                    logging::Logger::log(arena()) << "Re-locking root " << (vec+1) << endl;
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
                (*weight)(c[vec], real(l[root[vec]]));

                /*
                 * Orthogonalize and normalize
                 */
                for (int soln = 0;soln < nsoln;soln++)
                {
                    for (int svec = 0;svec < nvec;svec++)
                    {
                        double olap = (*innerProd)(old_c[soln][svec], c[vec]);
                        for (int idx = 0;idx < nc;idx++) c[vec][idx] -= old_c[soln][svec][idx]*olap;
                    }
                }
                double nrm = (*innerProd)(c[vec], c[vec]);
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

            vector<double> myreturn(nvec);
            for (int i = 0; i < nvec; i++)
                myreturn[i] = real(l[root[i]]);

            return myreturn;
        }

        void getSolution(int j, tensor::Tensor<> c)
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

        void getSolution(int j, vector<tensor::Tensor<>>&& c)
        {
            getSolution(j, c);
        }

        void getSolution(int j, vector<tensor::Tensor<>>& c)
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

        void getSolution(int j, tensor::Tensor<> c, tensor::Tensor<> hc)
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


        void getSolution(int j, vector<tensor::Tensor<>>&& c, vector<tensor::Tensor<>>&& hc)
        {
            getSolution(j, c, hc);
        }

        void getSolution(int j, vector<tensor::Tensor<>>& c, vector<tensor::Tensor<>>& hc)
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

        void reset(unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            this->nc = 1;
            this->innerProd = move(innerProd);
            this->weight = move(weight);
            init(1, LOWEST_ENERGY);
        }

        void reset(int nvec, unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            this->nc = 1;
            this->innerProd = move(innerProd);
            this->weight = move(weight);
            init(nvec, LOWEST_ENERGY);
        }

        void reset(int nvec, int nc, unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            this->nc = nc;
            this->innerProd = move(innerProd);
            this->weight = move(weight);
            init(nvec, LOWEST_ENERGY);
        }

        void reset(double t, unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            this->nc = 1;
            this->innerProd = move(innerProd);
            this->weight = move(weight);
            init(1, CLOSEST_ENERGY);
            target[0] = t;
        }

        void reset(double t, int nc, unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            this->nc = nc;
            this->innerProd = move(innerProd);
            this->weight = move(weight);
            init(1, CLOSEST_ENERGY);
            target[0] = t;
        }

        void reset(const vector<double>& t, unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            this->nc = 1;
            this->innerProd = move(innerProd);
            this->weight = move(weight);
            init(target.size(), CLOSEST_ENERGY);
            target = t;
        }

        void reset(const vector<double>& t, int nc, unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            this->nc = nc;
            this->innerProd = move(innerProd);
            this->weight = move(weight);
            init(target.size(), CLOSEST_ENERGY);
            target = t;
        }

        void reset(tensor::ConstTensor<> g, unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            this->nc = 1;
            this->innerProd = move(innerProd);
            this->weight = move(weight);
            init(1, GUESS_OVERLAP);
            guess = {{g}};
        }

        void reset(const vector<vector<tensor::ConstTensor<>>>& gs,
                   unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            this->nc = gs[0].size();
            this->innerProd = move(innerProd);
            this->weight = move(weight);

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
            vector<vector<tensor::Tensor<>>> c(nvec);
            vector<vector<tensor::Tensor<>>> hc(nvec);

            for (int vec = 0;vec < nvec;vec++)
            {
                for (int i = 0;i < nc;i++)
                {
                    c[vec].push_back(old_c[0][vec][i].construct());
                    hc[vec].push_back(old_hc[0][vec][i].construct());
                }

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
