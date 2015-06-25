#include "uhf_local.hpp"

using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;
using namespace aquarius::time;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace scf
{

LocalUHF::LocalUHF(const string& name, Config& config)
: UHF(name, config) {}

void LocalUHF::calcSMinusHalf()
{
    const Molecule& molecule = this->template get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();

    const vector<int>& norb = molecule.getNumOrbitals();

    Tensor<BOUNDED|PGSYMMETRIC> S      = this->template get   <Tensor<>>("S");
    Tensor<BOUNDED|PGSYMMETRIC> Smhalf = this->template gettmp<Tensor<>>("S^-1/2");

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        //cout << "S " << (i+1) << endl;
        //vector<T> vals;
        //S({i,i}).getAllData(vals);
        //printmatrix(norb[i], norb[i], vals.data(), 6, 3, 108);
    }

    for (int i = 0;i < group.getNumIrreps();i++)
    {
        if (norb[i] == 0) continue;

        vector<int> irreps = {i,i};

        vector<double> s(norb[i]*norb[i]);
        vector<double> e(norb[i]);
        vector<double> smhalf(norb[i]*norb[i]);

        KeyValueVector kv = S.getAllDataByIrrep(irreps);
        assert(kv.size() == norb[i]*norb[i]);

        for (int j = 0;j < norb[i]*norb[i];j++)
        {
            s[kv.key(j)] = kv.value<double>(j);
        }

        int info = heev('V', 'U', norb[i], s.data(), norb[i], e.data());
        assert(info == 0);

        fill(smhalf.begin(), smhalf.end(), 0);
        for (int j = 0;j < norb[i];j++)
        {
            ger(norb[i], norb[i], 1/sqrt(e[j]), &s[j*norb[i]], 1, &s[j*norb[i]], 1, smhalf.data(), norb[i]);
        }

        for (int j = 0;j < norb[i]*norb[i];j++)
        {
            kv.key(j) = j;
            kv.value(j, smhalf[j]);
        }

        Smhalf.setDataByIrrep(irreps, kv);
    }
}

void LocalUHF::diagonalizeFock()
{
    const Molecule& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();

    Tensor<BOUNDED|PGSYMMETRIC> S  = this->template get   <Tensor<>>("S");
    Tensor<BOUNDED|PGSYMMETRIC> Fa = this->template get   <Tensor<>>("Fa");
    Tensor<BOUNDED|PGSYMMETRIC> Fb = this->template get   <Tensor<>>("Fb");
    Tensor<BOUNDED|PGSYMMETRIC> Ca = this->template gettmp<Tensor<>>("Ca");
    Tensor<BOUNDED|PGSYMMETRIC> Cb = this->template gettmp<Tensor<>>("Cb");

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        if (norb[i] == 0) continue;

        vector<int> irreps = {i,i};
        vector<double> e(norb[i]);

        vector<double> s(norb[i]*norb[i]);
        KeyValueVector kv = S.getAllDataByIrrep(irreps);
        assert(kv.size() == norb[i]*norb[i]);
        for (int j = 0;j < norb[i]*norb[i];j++)
        {
            s[kv.key(j)] = kv.value<double>(j);
        }

        for (int spin : {0,1})
        {
            auto& F = (spin == 0 ? Fa : Fb);
            auto& C = (spin == 0 ? Ca : Cb);
            auto& E = (spin == 0 ? E_alpha[i] : E_beta[i]);

            int info;
            vector<double> fock(norb[i]*norb[i]);
            vector<double> tmp(s);

            F.getAllDataByIrrep(irreps, kv);
            assert(kv.size() == norb[i]*norb[i]);
            double nrm = 0;
            for (int j = 0;j < norb[i]*norb[i];j++)
            {
                fock[kv.key(j)] = kv.value<double>(j);
            }

            info = hegv(AXBX, 'V', 'U', norb[i], fock.data(), norb[i], tmp.data(), norb[i], e.data());
            assert(info == 0);
            arena().comm().Bcast(e);
            for (int j = 0;j < norb[i];j++) E[j] = e[j];

            for (int j = 0;j < norb[i];j++)
            {
                double sign = 0;
                for (int k = 0;k < norb[i];k++)
                {
                    if (aquarius::abs(fock[k+j*norb[i]]) > 1e-10)
                    {
                        sign = (fock[k+j*norb[i]] < 0 ? -1 : 1);
                        break;
                    }
                }
                scal(norb[i], sign, &fock[j*norb[i]], 1);
            }

            for (int j = 0;j < norb[i]*norb[i];j++)
            {
                kv.key(j) = j;
                kv.value(j, fock[j]);
            }

            C.setDataByIrrep(irreps, kv);
        }
    }
}

}
}
