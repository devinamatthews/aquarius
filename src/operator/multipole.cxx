#include "multipole.hpp"

using namespace aquarius::tensor;
using namespace aquarius::integrals;
using namespace aquarius::scf;
using namespace aquarius::input;

namespace aquarius
{
namespace op
{

template <typename T>
Multipole<T>::Multipole(const string& name, const UHF<T>& uhf, int Lmin_, int Lmax_)
: MOOperator<T>(uhf), CompositeTensor<Multipole<T>,
  OneElectronOperator<T>,T>(name, Lmax_ == -1 ? (Lmin_+1)*(Lmin_+2)/2 :
          (Lmax_+1)*(Lmax_+2)*(Lmax_+3)/6-Lmin_*(Lmin_+1)*(Lmin_+2)/6),
  Lmin(Lmin_), Lmax(Lmax_ == -1 ? Lmin_ : Lmax_)
{
    Context context;
    const Molecule& m = uhf.getMolecule();
    vector<int> N = m.getNumOrbitals();

    SymmetryBlockedTensor<T> ao(uhf.arena, m.getGroup(), 2, {N,N}, {NS,NS}, false);

    int xyztot = 0;
    for (int L = Lmin;L <= Lmax;L++)
    {
        vector<vector<tkv_pair<T>>> pairs;

        int ij = 0;
        for (Molecule::const_shell_iterator i = m.getShellsBegin();i != m.getShellsEnd();++i)
        {
            for (Molecule::const_shell_iterator j = i;j != m.getShellsEnd();++j)
            {
                if (i < j) continue;

                if (ij%this->arena.nproc == this->arena.rank)
                {
                    //calc integrals

                    //put into pairs
                }

                ij++;
            }
        }

        for (int xyz = 0;xyz < (L+1)*(L+2)/2;xyz++)
        {
            ao = (T)0;
            //ao.writeRemoteData(pairs[xyz]);
            this->tensors[xyztot++].tensor = new OneElectronOperator<T>(name, uhf, ao);
        }
    }
}

template <typename T>
const OneElectronOperator<T>& Multipole<T>::operator()(int L, int xyz) const
{
    assert(L >= Lmin && L <= Lmax);
    assert(xyz >= 0 && xyz < (L+1)*(L+2)/2);
    return (*this)(L*(L+1)*(L+2)/6-Lmin*(Lmin+1)*(Lmin+2)/6+xyz);
}

template <typename T>
const OneElectronOperator<T>& Multipole<T>::operator()(int x, int y, int z) const
{
    int L = x+y+z;
    assert(L >= Lmin && L <= Lmax);
    int xyz = z+(L+1-x)*(L+2-x)/2;
    return (*this)(L, xyz);
}
