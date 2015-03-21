#include "rys.hpp"

#include "fmgamma.hpp"

namespace aquarius
{
namespace integrals
{

void Rys::operator()(double T, int n, double* restrict rt, double* restrict wt)
{
    Fm fm;
    row<double> ssssm(2*n+1);
    fm(T, 2*n, ssssm.data());

    matrix<double> R(n+1, n+1);
    for (int i = 0;i <= n;i++)
    {
        for (int j = 0;j <= n;j++)
        {
            R[i][j] = ssssm[i+j];
        }
    }
    potrf('U', n+1, R.data(), n+1);

    row<double> a(n), b(n);
    a[0] = R[0][1]/R[0][0];
    for (int i = 0;i < n-1;i++)
    {
        a[i+1] = R[i+1][i+2]/R[i+1][i+1] - R[i][i+1]/R[i][i];
        b[i] = R[i+1][i+1]/R[i][i];
    }

    matrix<double> Z(n,n);
    int info = stev('V', n, a.data(), b.data(), Z.data(), n);
    assert(info == 0);

    for (int i = 0;i < n;i++)
    {
        rt[i] = a[i];
        wt[i] = Z[i][0]*Z[i][0]*ssssm[0];
    }
}

}
}
