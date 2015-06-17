#ifndef _AQUARIUS_INTEGRALS_FMGAMMA_HPP_
#define _AQUARIUS_INTEGRALS_FMGAMMA_HPP_

#include "util/global.hpp"

namespace aquarius
{
namespace integrals
{

class Fm
{
    protected:
        constexpr static int TAYLOR_N = 8;

        const static int TMAX[40];

        static double FMTABLE[2421][40];

        static bool inited;

        static double direct(double T, int m);

        static double taylor(double T, int n);

        static double asymptotic(double T, int m);

    public:
        Fm()
        {
            while (!inited)
            {
                #pragma omp critical
                {
                    if (!inited)
                    {
                        calcTable();
                        inited = true;
                    }
                }
            }
        }

        static void calcTable();

        template <typename Container>
        void operator()(double A0, double T, Container&& array) const
        {
            operator()(A0, T, array.size()-1, array);
        }

        template <typename Container>
        void operator()(double A0, double T, int n, Container&& array) const
        {
            if (T > TMAX[n])
            {
                array[n] = asymptotic(T, n);
            }
            else
            {
                array[n] = taylor(T, n);
            }

            if (n == 0)
            {
                array[0] *= A0;
                return;
            }

            double emt = exp(-T);
            double twoT = 2*T;

            for (int i = n;i > 0;i--)
            {
                array[i-1] = (twoT*array[i] + emt) / (2*i-1);
                array[i] *= A0;
            }
            array[0] *= A0;
        }
};

}
}

#endif
