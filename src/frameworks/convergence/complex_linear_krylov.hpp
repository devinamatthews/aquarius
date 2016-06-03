#ifndef _AQUARIUS_FRAMEWORKS_CONVERGENCE_COMPLEX_LINEAR_KRYLOV_HPP_
#define _AQUARIUS_FRAMEWORKS_CONVERGENCE_COMPLEX_LINEAR_KRYLOV_HPP_

#include "frameworks/util.hpp"
#include "frameworks/task.hpp"
#include "frameworks/tensor.hpp"
#include "frameworks/logging.hpp"

namespace aquarius
{
namespace convergence
{

struct ComplexWeight
{
    virtual ~ComplexWeight() {}

    virtual void operator()(vector<tensor::Tensor<>>& a, dcomplex omega) const = 0;
};

class ComplexLinearKrylov
{
    protected:
        vector<tensor::Tensor<>> old_c_r;
        vector<tensor::Tensor<>> old_c_i;
        vector<tensor::Tensor<>> old_hc_r;
        vector<tensor::Tensor<>> old_hc_i;
        tensor::ConstTensor<> rhs;
        vector<dcomplex> b;
        matrix<dcomplex> e;
        int maxextrap, nextrap;
        vector<dcomplex> v;
        bool continuous;
        unique_ptr<InnerProd> innerProd;
        unique_ptr<Weight> weight;

        void parse(const task::Config& config)
        {
            nextrap = 0;
            maxextrap = config.get<int>("order");
            continuous = config.get<string>("compaction") == "continuous";
        }

        void init()
        {
            nextrap = 0;
            b.resize(nextrap);
            e.resize(nextrap, nextrap);
            v.resize(nextrap);
        }

        void addVectors(tensor::ConstTensor<> c_r, tensor::ConstTensor<> c_i,
                        tensor::ConstTensor<> hc_r, tensor::ConstTensor<> hc_i)
        {
            nextrap++;

            b.resize(nextrap);
            e.resize(nextrap, nextrap);
            v.resize(nextrap);

            if (old_c_r.size() < nextrap)
                old_c_r.emplace_back(c_r);
            else
                old_c_r[nextrap-1] = c_r;

            if (old_c_i.size() < nextrap)
                old_c_i.emplace_back(c_i);
            else
                old_c_i[nextrap-1] = c_i;

            if (old_hc_r.size() < nextrap)
                old_hc_r.emplace_back(hc_r);
            else
                old_hc_r[nextrap-1] = hc_r;

            if (old_hc_i.size() < nextrap)
                old_hc_i.emplace_back(hc_i);
            else
                old_hc_i[nextrap-1] = hc_i;

            /*
             * Compute the overlap with the rhs vector
             */
            b[nextrap-1] = dcomplex((*innerProd)(c_r, rhs), -(*innerProd)(c_i, rhs));

            /*
             * Augment the subspace matrix with the new vectors
             */
            e[nextrap-1][nextrap-1] =
                dcomplex((*innerProd)(hc_r, c_r) + (*innerProd)(hc_i, c_i),
                         (*innerProd)(hc_i, c_r) - (*innerProd)(hc_r, c_i));

            for (int extrap = 0;extrap < nextrap-1;extrap++)
            {
                e[   extrap][nextrap-1] =
                    dcomplex((*innerProd)(old_hc_r[extrap], c_r) + (*innerProd)(old_hc_i[extrap], c_i),
                             (*innerProd)(old_hc_i[extrap], c_r) - (*innerProd)(old_hc_r[extrap], c_i));
                e[nextrap-1][   extrap] =
                    dcomplex((*innerProd)(hc_r, old_c_r[extrap]) + (*innerProd)(hc_i, old_c_i[extrap]),
                             (*innerProd)(hc_i, old_c_r[extrap]) - (*innerProd)(hc_r, old_c_i[extrap]));
            }
        }

        void getRoot(tensor::Tensor<> c_r, tensor::Tensor<> c_i,
                     tensor::Tensor<> hc_r, tensor::Tensor<> hc_i)
        {
            getRoot(c_r, c_i);

            hc_r = 0;
            hc_i = 0;
            for (int extrap = 0;extrap < nextrap;extrap++)
            {
                hc_r += old_hc_r[extrap]*real(v[extrap]);
                hc_r -= old_hc_i[extrap]*imag(v[extrap]);
                hc_i += old_hc_r[extrap]*imag(v[extrap]);
                hc_i += old_hc_i[extrap]*real(v[extrap]);
            }
        }

        void getRoot(tensor::Tensor<> c_r, tensor::Tensor<> c_i)
        {
            c_r = 0;
            c_i = 0;
            for (int extrap = 0;extrap < nextrap;extrap++)
            {
                c_r += old_c_r[extrap]*real(v[extrap]);
                c_r -= old_c_i[extrap]*imag(v[extrap]);
                c_i += old_c_r[extrap]*imag(v[extrap]);
                c_i += old_c_i[extrap]*real(v[extrap]);
            }
        }

    public:
        template <typename... Args>
        ComplexLinearKrylov(const task::Config& config, Args&&... args)
        {
            parse(config);
            reset(forward<Args>(args)...);
        }

        ComplexLinearKrylov(const ComplexLinearKrylov& other) = delete;

        ComplexLinearKrylov& operator=(const ComplexLinearKrylov& other) = delete;

        void extrapolate(tensor::Tensor<> c_r, tensor::Tensor<> c_i,
                         tensor::Tensor<> hc_r, tensor::Tensor<> hc_i,
                         dcomplex omega)
        {
            using slice::all;

            addVectors(c_r, c_i, hc_r, hc_i);

            /*
            for (int i = 0;i < nextrap;i++)
            {
                for (int j = 0;j < nextrap;j++)
                {
                    printf(" %23.15e", e[i][j].real());
                }
                printf("\n");
            }
            printf("\n");

            for (int i = 0;i < nextrap;i++)
            {
                for (int j = 0;j < nextrap;j++)
                {
                    printf(" %23.15e", e[i][j].imag());
                }
                printf("\n");
            }
            printf("\n");
            */

            /*
             * Solve the linear problem in the subspace
             */
            vector<integer> ipiv(nextrap);
            matrix<dcomplex> tmp = e;
            v = b;

            /*
            for (int j = 0;j < nextrap;j++)
            {
                printf(" %23.15e", v[j].real());
            }
            printf("\n\n");

            for (int j = 0;j < nextrap;j++)
            {
                printf(" %23.15e", v[j].imag());
            }
            printf("\n\n");
            */

            int info = gesv(nextrap, 1, tmp.data(), nextrap, ipiv.data(), v.data(), nextrap);
            if (info != 0) throw runtime_error(str("krylov: Info in gesv: %d", info));

            /*
            for (int j = 0;j < nextrap;j++)
            {
                printf(" %23.15e", v[j].real());
            }
            printf("\n\n");

            for (int j = 0;j < nextrap;j++)
            {
                printf(" %23.15e", v[j].imag());
            }
            printf("\n\n");
            */

            /*
             * Generate solution vector
             */
            getRoot(c_r, c_i, hc_r, hc_i);

            if (continuous)
            {
                /*
                 * Orthogonalize to previous vectors
                 */
                for (int extrap = 0;extrap < nextrap-1;extrap++)
                {
                    Scalar olap_r = (*innerProd)(old_c_r[extrap], c_r) + (*innerProd)(old_c_i[extrap], c_i);
                    Scalar olap_i = (*innerProd)(old_c_r[extrap], c_i) - (*innerProd)(old_c_i[extrap], c_r);
                    c_r -= olap_r*old_c_r[extrap];
                    c_r += olap_i*old_c_i[extrap];
                    c_i -= olap_r*old_c_i[extrap];
                    c_i -= olap_i*old_c_r[extrap];
                }

                /*
                 * Save current solution as new vector in the Krylov subspace
                 */
                nextrap--;
                addVectors(c_r, c_i, hc_r, hc_i);

                /*
                 * If the subspace is full, eject the oldest vector.
                 */
                if (nextrap == maxextrap)
                {
                    rotate(old_c_r.begin(), old_c_r.begin()+1, old_c_r.end());
                    rotate(old_c_i.begin(), old_c_i.begin()+1, old_c_i.end());
                    rotate(old_hc_r.begin(), old_hc_r.begin()+1, old_hc_r.end());
                    rotate(old_hc_i.begin(), old_hc_i.begin()+1, old_hc_i.end());
                    e.rotate(1,1);
                    rotate(b.begin(), b.begin()+1, b.end());
                    nextrap--;
                }
            }
            else if (nextrap == maxextrap)
            {
                /*
                 * If the maximum size of the subspace has been reached, a smaller
                 * subspace must be constructed which contains the best approximation
                 * to the solution.
                 */
                logging::Logger::log(arena()) << "Compacting..." << endl;
                nextrap = 0;
                addVectors(c_r, c_i, hc_r, hc_i);
                v[0] = 1.0;
            }

            /*
             * Form residual and apply Davidson-like correction
             */
            hc_r -= rhs;
            c_r = hc_r;
            c_i = hc_i;
            (*weight)(c_r, c_i, omega);

            //printf("<B|r>: %.15f\n", (*innerProd)((*rhs)*c_r));
            //printf("<B|i>: %.15f\n", (*innerProd)((*rhs)*c_i));

            //printf("<rr|rr>: %.15f\n", (*innerProd)(hc_r*hc_r));
            //printf("<ri|ri>: %.15f\n", (*innerProd)(hc_i*hc_i));

            //printf("<Rr|Rr>: %.15f\n", (*innerProd)(c_r*c_r));
            //printf("<Ri|Ri>: %.15f\n", (*innerProd)(c_i*c_i));

            //printf("<B|r>: %.15f\n", (*innerProd)((*rhs)*c_r));
            //printf("<B|i>: %.15f\n", (*innerProd)((*rhs)*c_i));

            /*
             * Orthogonalize to previous vectors
             */
            for (int extrap = 0;extrap < nextrap;extrap++)
            {
                Scalar olap_r = (*innerProd)(old_c_r[extrap], c_r) + (*innerProd)(old_c_i[extrap], c_i);
                Scalar olap_i = (*innerProd)(old_c_r[extrap], c_i) - (*innerProd)(old_c_i[extrap], c_r);
                c_r -= olap_r*old_c_r[extrap];
                c_r += olap_i*old_c_i[extrap];
                c_i -= olap_r*old_c_i[extrap];
                c_i -= olap_i*old_c_r[extrap];
            }

            Scalar norm = sqrt(aquarius::abs((*innerProd)(c_r, c_r)) +
                               aquarius::abs((*innerProd)(c_i, c_i)));
            c_r /= norm;
            c_i /= norm;

            //printf("<B|r>: %.15f\n", (*innerProd)(rhs, c_r));
            //printf("<B|i>: %.15f\n", (*innerProd)(rhs, c_i));
        }

        void getSolution(tensor::Tensor<> c_r, tensor::Tensor<> c_i)
        {
            if (continuous)
            {
                c_r = old_c_r[nextrap-1];
                c_i = old_c_i[nextrap-1];
            }
            else
            {
                getRoot(c_r, c_i);
            }
        }

        void reset(tensor::ConstTensor<> rhs,
                   unique_ptr<Weight> weight, unique_ptr<InnerProd> innerProd = new InnerProd())
        {
            init();
            this->rhs = move(rhs);
            this->weight = move(weight);
            this->innerProd = move(innerProd);
        }
};

}
}

#endif
