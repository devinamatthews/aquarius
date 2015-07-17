#ifndef _AQUARIUS_COMPLEX_LINEAR_KRYLOV_HPP_
#define _AQUARIUS_COMPLEX_LINEAR_KRYLOV_HPP_

#include "util/global.hpp"

#include "input/config.hpp"
#include "task/task.hpp"
#include "cc/complex_denominator.hpp"

#include "davidson.hpp"

namespace aquarius
{
namespace convergence
{

template<typename T>
class ComplexLinearKrylov : public task::Destructible
{
    private:
        ComplexLinearKrylov(const ComplexLinearKrylov& other);

        ComplexLinearKrylov& operator=(const ComplexLinearKrylov& other);

    protected:
        typedef typename T::dtype U;
        typedef complex_type_t<U> CU;
        unique_vector<T> old_c_r;
        unique_vector<T> old_c_i;
        unique_vector<T> old_hc_r;
        unique_vector<T> old_hc_i;
        unique_ptr<T> rhs;
        vector<CU> b;
        marray<CU,2> e;
        int maxextrap, nextrap;
        vector<CU> v;
        bool continuous;

        void parse(const input::Config& config)
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
            v.resize(nextrap, nextrap);
        }

        void addVectors(const T& c_r, const T& c_i, const T& hc_r, const T& hc_i)
        {
            nextrap++;

            b.resize(nextrap);
            e.resize(nextrap, nextrap);
            v.resize(nextrap, nextrap);

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
            b[nextrap-1].real( scalar(c_r*(*rhs)));
            b[nextrap-1].imag(-scalar(c_i*(*rhs)));

            /*
             * Augment the subspace matrix with the new vectors
             */
            e[nextrap-1][nextrap-1].real(scalar(hc_r*c_r) + scalar(hc_i*c_i));
            e[nextrap-1][nextrap-1].imag(scalar(hc_i*c_r) - scalar(hc_r*c_i));

            for (int extrap = 0;extrap < nextrap-1;extrap++)
            {
                e[   extrap][nextrap-1].real(scalar(old_hc_r[extrap]*c_r) + scalar(old_hc_i[extrap]*c_i));
                e[   extrap][nextrap-1].imag(scalar(old_hc_i[extrap]*c_r) - scalar(old_hc_r[extrap]*c_i));
                e[nextrap-1][   extrap].real(scalar(hc_r*old_c_r[extrap]) + scalar(hc_i*old_c_i[extrap]));
                e[nextrap-1][   extrap].imag(scalar(hc_i*old_c_r[extrap]) - scalar(hc_r*old_c_i[extrap]));
            }
        }

        void getRoot(T& c_r, T& c_i, T& hc_r, T& hc_i)
        {
            getRoot(c_r, c_i);

            hc_r = 0;
            hc_i = 0;
            for (int extrap = 0;extrap < nextrap;extrap++)
            {
                hc_r += old_hc_r[extrap]*v[extrap].real();
                hc_r -= old_hc_i[extrap]*v[extrap].imag();
                hc_i += old_hc_r[extrap]*v[extrap].imag();
                hc_i += old_hc_i[extrap]*v[extrap].real();
            }
        }

        void getRoot(T& c_r, T& c_i)
        {
            c_r = 0;
            c_i = 0;
            for (int extrap = 0;extrap < nextrap;extrap++)
            {
                c_r += old_c_r[extrap]*v[extrap].real();
                c_r -= old_c_i[extrap]*v[extrap].imag();
                c_i += old_c_r[extrap]*v[extrap].imag();
                c_i += old_c_i[extrap]*v[extrap].real();
            }
        }

    public:
        ComplexLinearKrylov(const input::Config& config, const T& g)
        {
            parse(config);
            reset(g);
        }

        void extrapolate(T& c_r, T& c_i, T& hc_r, T& hc_i, const cc::ComplexDenominator<U>& D, const CU& omega)
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
            marray<CU,2> tmp(e, construct_copy);
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
            if (info != 0) throw runtime_error(strprintf("krylov: Info in gesv: %d", info));

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
                    U olap_r = scalar(old_c_r[extrap]*c_r) + scalar(old_c_i[extrap]*c_i);
                    U olap_i = scalar(old_c_r[extrap]*c_i) - scalar(old_c_i[extrap]*c_r);
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
                    rotate(old_c_r.pbegin(), old_c_r.pbegin()+1, old_c_r.pend());
                    rotate(old_c_i.pbegin(), old_c_i.pbegin()+1, old_c_i.pend());
                    rotate(old_hc_r.pbegin(), old_hc_r.pbegin()+1, old_hc_r.pend());
                    rotate(old_hc_i.pbegin(), old_hc_i.pbegin()+1, old_hc_i.pend());
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
                task::Logger::log(c_r.arena) << "Compacting..." << endl;
                nextrap = 0;
                addVectors(c_r, c_i, hc_r, hc_i);
                v[0] = 1.0;
            }

            /*
             * Form residual and apply Davidson-like correction
             */
            hc_r -= *rhs;
            c_r = hc_r;
            c_i = hc_i;
            D.weight(c_r, c_i, omega);

            //printf("<B|r>: %.15f\n", scalar((*rhs)*c_r));
            //printf("<B|i>: %.15f\n", scalar((*rhs)*c_i));

            //printf("<rr|rr>: %.15f\n", scalar(hc_r*hc_r));
            //printf("<ri|ri>: %.15f\n", scalar(hc_i*hc_i));

            //printf("<Rr|Rr>: %.15f\n", scalar(c_r*c_r));
            //printf("<Ri|Ri>: %.15f\n", scalar(c_i*c_i));

            //printf("<B|r>: %.15f\n", scalar((*rhs)*c_r));
            //printf("<B|i>: %.15f\n", scalar((*rhs)*c_i));

            /*
             * Orthogonalize to previous vectors
             */
            for (int extrap = 0;extrap < nextrap;extrap++)
            {
                U olap_r = scalar(old_c_r[extrap]*c_r) + scalar(old_c_i[extrap]*c_i);
                U olap_i = scalar(old_c_r[extrap]*c_i) - scalar(old_c_i[extrap]*c_r);
                c_r -= olap_r*old_c_r[extrap];
                c_r += olap_i*old_c_i[extrap];
                c_i -= olap_r*old_c_i[extrap];
                c_i -= olap_i*old_c_r[extrap];
            }

            U norm = sqrt(aquarius::abs(scalar(c_r*c_r)) + aquarius::abs(scalar(c_i*c_i)));
            c_r /= norm;
            c_i /= norm;

            //printf("<B|r>: %.15f\n", scalar((*rhs)*c_r));
            //printf("<B|i>: %.15f\n", scalar((*rhs)*c_i));
        }

        void getSolution(T& c_r, T& c_i)
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

        void reset(const T& g)
        {
            init();
            rhs.reset(new T(g));
        }
};

}
}

#endif
