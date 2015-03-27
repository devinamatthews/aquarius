#include "uhf_local.hpp"

using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;
using namespace aquarius::time;
using namespace aquarius::op;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace scf
{

template <typename T>
LocalUHF<T>::LocalUHF(const string& name, Config& config)
: UHF<T>(name, config) {}

template <typename T>
void LocalUHF<T>::calcSMinusHalf()
{
    const Molecule& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();

    SymmetryBlockedTensor<T>& S = this->template get<SymmetryBlockedTensor<T> >("S");
    SymmetryBlockedTensor<T>& Smhalf = this->template gettmp<SymmetryBlockedTensor<T> >("S^-1/2");

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        //cout << "S " << (i+1) << endl;
        //vector<T> vals;
        //S({i,i}).getAllData(vals);
        //printmatrix(norb[i], norb[i], vals.data(), 6, 3, 108);
    }

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        if (norb[i] == 0) continue;

        vector<int> irreps(2,i);
        vector<typename real_type<T>::type> E(norb[i]);

        if (S.arena.rank == 0)
        {
            vector<T> s;
            vector<T> smhalf(norb[i]*norb[i]);

            S.getAllData(irreps, s, 0);
            assert(s.size() == norb[i]*norb[i]);

            //PROFILE_FLOPS(26*norb[i]*norb[i]*norb[i]);
            int info = heev('V', 'U', norb[i], s.data(), norb[i], E.data());
            assert(info == 0);

            fill(smhalf.begin(), smhalf.end(), (T)0);
            //PROFILE_FLOPS(2*norb[i]*norb[i]*norb[i]);
            for (int j = 0;j < norb[i];j++)
            {
                ger(norb[i], norb[i], 1/sqrt(E[j]), &s[j*norb[i]], 1, &s[j*norb[i]], 1, smhalf.data(), norb[i]);
            }

            vector< tkv_pair<T> > pairs(norb[i]*norb[i]);

            for (int j = 0;j < norb[i]*norb[i];j++)
            {
                pairs[j].k = j;
                pairs[j].d = smhalf[j];
            }

            Smhalf.writeRemoteData(irreps, pairs);
        }
        else
        {
            S.getAllData(irreps, 0);
            Smhalf.writeRemoteData(irreps);
        }
    }
}

template <typename T>
void LocalUHF<T>::diagonalizeFock()
{
    const Molecule& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();

    SymmetryBlockedTensor<T>& S = this->template get<SymmetryBlockedTensor<T> >("S");
    SymmetryBlockedTensor<T>& Fa = this->template get<SymmetryBlockedTensor<T> >("Fa");
    SymmetryBlockedTensor<T>& Fb = this->template get<SymmetryBlockedTensor<T> >("Fb");
    SymmetryBlockedTensor<T>& Ca = this->template gettmp<SymmetryBlockedTensor<T> >("Ca");
    SymmetryBlockedTensor<T>& Cb = this->template gettmp<SymmetryBlockedTensor<T> >("Cb");

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        //cout << "F " << (i+1) << endl;
        //vector<T> vals;
        //Fa({i,i}).getAllData(vals);
        //printmatrix(norb[i], norb[i], vals.data(), 6, 3, 108);
    }

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        if (norb[i] == 0) continue;

        vector<int> irreps(2,i);

        if (S.arena.rank == 0)
        {
            vector<T> s;
            S.getAllData(irreps, s, 0);
            assert(s.size() == norb[i]*norb[i]);

            for (int spin : {0,1})
            {
                SymmetryBlockedTensor<T>& F = (spin == 0 ? Fa : Fb);
                SymmetryBlockedTensor<T>& C = (spin == 0 ? Ca : Cb);

                int info;
                vector<T> fock, ctsc(norb[i]*norb[i]);
                vector<tkv_pair<T>> pairs(norb[i]*norb[i]);
                vector<T> tmp(s);

                F.getAllData(irreps, fock, 0);
                assert(fock.size() == norb[i]*norb[i]);
                //PROFILE_FLOPS(9*norb[i]*norb[i]*norb[i]);
                info = hegv(AXBX, 'V', 'U', norb[i], fock.data(), norb[i], tmp.data(), norb[i], E_alpha[i].data());
                assert(info == 0);
                S.arena.Bcast(E_alpha[i], 0);

                for (int j = 0;j < norb[i];j++)
                {
                    T sign = 0;
                    for (int k = 0;k < norb[i];k++)
                    {
                        if (aquarius::abs(fock[k+j*norb[i]]) > 1e-10)
                        {
                            sign = (fock[k+j*norb[i]] < 0 ? -1 : 1);
                            break;
                        }
                    }
                    //PROFILE_FLOPS(norb[i]);
                    scal(norb[i], sign, &fock[j*norb[i]], 1);
                }

                for (int j = 0;j < norb[i]*norb[i];j++)
                {
                    pairs[j].k = j;
                    pairs[j].d = fock[j];
                }

                C.writeRemoteData(irreps, pairs);
            }
        }
        else
        {
            S.getAllData(irreps, 0);
            Fa.getAllData(irreps, 0);
            S.arena.Bcast(E_alpha[i], 0);
            Ca.writeRemoteData(irreps);
            Fb.getAllData(irreps, 0);
            S.arena.Bcast(E_beta[i], 0);
            Cb.writeRemoteData(irreps);
        }
    }
}

INSTANTIATE_SPECIALIZATIONS(LocalUHF);

}
}
