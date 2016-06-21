#ifndef _AQUARIUS_FRAMEWORKS_INTEGRALS_FMGAMMA_HPP_
#define _AQUARIUS_FRAMEWORKS_INTEGRALS_FMGAMMA_HPP_

#include "frameworks/util.hpp"

namespace aquarius
{
namespace integrals
{

class FmGamma
{
    protected:
        struct TaylorTable
        {
            static double FMTABLE[2421][40];

            TaylorTable();
        };

        constexpr static int TAYLOR_N = 8;
        const static int TMAX[40];

        static double direct(double T, int m);

        static double asymptotic(double T, int m);

        double taylor(double T, int n);

        TaylorTable& table();

    public:
        void operator()(double T, int n, double* array);

        void operator()(double T, vector<double>& array)
        {
            (*this)(T, array.size()-1, array.data());
        }
};

}
}

#endif
