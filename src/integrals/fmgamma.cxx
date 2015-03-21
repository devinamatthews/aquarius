#include "fmgamma.hpp"

namespace aquarius
{
namespace integrals
{

const int Fm::TMAX[40] = { 33, 37, 40, 43, 46, 49, 51, 53, 56, 58,
                              60, 62, 64, 66, 68, 70, 72, 74, 76, 78,
                              80, 82, 83, 85, 87, 89, 90, 92, 94, 96,
                              97, 99, 101, 102, 104, 106, 108, 110, 112, 114 };

double Fm::FMTABLE[2421][40];

bool Fm::inited = false;

double Fm::direct(double T, int m)
{
    constexpr double epsilon = numeric_limits<double>::epsilon();

    double ap = m + 0.5;

    assert(T >= 0.0);
    if (T < 1e-13)
    {
        return 0.5/ap;
    }

    if (T > TMAX[m]) return asymptotic(T, m);

    double sum = 1.0/(2*ap*exp(T));
    double delt = sum;
    ap += 1.0;

    while (true)
    {
        delt *= T/ap;
        ap += 1.0;
        if (abs(delt) < epsilon*sum) break;
        sum += delt;
    }

    return sum;
}

double Fm::taylor(double T, int n)
{
    int tidx = (int)round(T*20.0);
    double tr = (double)tidx/20.0;
    double trmt = tr-T;

    double ans = 0.0;

    array<double, TAYLOR_N> fac;
    fac[0] = 1.0;
    for (int i = 1;i < TAYLOR_N;i++)
    {
        fac[i] = fac[i-1]*trmt/i;
    }

    for (int i = 0;i < TAYLOR_N;i++)
    {
        ans += FMTABLE[tidx][n+i]*fac[i];
    }

    return ans;
}

double Fm::asymptotic(double T, int m)
{
    double sum = -m*log(2*T);
    for (int i = 3;i < 2*m;i += 2)
    {
        sum += log(i);
    }
    return exp(sum) * sqrt(M_PI/T)/2;
}

void Fm::calcTable()
{
    // store tabulated points for T=0,TMAX[m],0.05, m=0,39
    #pragma omp parallel for
    for (int m = 0;m < 40;m++)
    {
        for (int Tidx = 0;Tidx < (TMAX[m]+TAYLOR_N-1)*20+1;Tidx++)
        {
            FMTABLE[Tidx][m] = direct(((double)Tidx)/20.0, m);
        }
    }
    inited = true;
}

void Fm::operator()(double T, int n, double* array)
{
    if (T > TMAX[n])
    {
        array[n] = asymptotic(T, n);
    }
    else
    {
        array[n] = taylor(T, n);
    }

    if (n == 0) return;

    double emt = exp(-T);
    double twoT = 2*T;

    for (int i = n;i > 0;i--)
    {
        array[i-1] = (twoT*array[i] + emt) / (2*i-1);
    }
}

}
}
