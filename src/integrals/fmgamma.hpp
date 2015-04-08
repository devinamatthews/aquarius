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

        void operator()(double T, int n, double* array);

        void operator()(double T, vector<double>& array)
        {
            operator()(T, array.size()-1, array.data());
        }
};

}
}

#endif
