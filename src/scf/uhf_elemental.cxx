#include "uhf_elemental.hpp"

using namespace El;

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

ElementalUHF::ElementalUHF(const string& name, Config& config)
: UHF(name, config) {}

void ElementalUHF::calcSMinusHalf()
{
    const auto& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();

    Tensor<BOUNDED|PGSYMMETRIC|DISTRIBUTED> S      = this->template get   <Tensor<>>("S");
    Tensor<BOUNDED|PGSYMMETRIC|DISTRIBUTED> Smhalf = this->template gettmp<Tensor<>>("S^-1/2");

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        if (norb[i] == 0) continue;

        vector<int> irreps = {i,i};
        vector<double> e(norb[i]);

        DistMatrix<double> S_elem(norb[i], norb[i]);

        KeyValueVector kv(Field::DOUBLE);
        KeyValueVector kv2(Field::DOUBLE);

        int cshift = S_elem.ColShift();
        int rshift = S_elem.RowShift();
        int cstride = S_elem.ColStride();
        int rstride = S_elem.RowStride();

        for (int k = 0;k < S_elem.LocalHeight();k++)
        {
            for (int j = 0;j < S_elem.LocalWidth();j++)
            {
                int c = cshift+k*cstride;
                int r = rshift+j*rstride;
                if (r > c) continue;
                kv.push_back(r*norb[i]+c);
                if (r == c) continue;
                kv2.push_back(c*norb[i]+r);
            }
        }

        S.getRemoteDataByIrrep(irreps, kv);

        for (int p = 0;p < kv.size();p++)
        {
            int r = kv.key(p)/norb[i];
            int c = kv.key(p)%norb[i];
            int k = (c-cshift)/cstride;
            int j = (r-rshift)/rstride;
            S_elem.SetLocal(k, j, kv.value<double>(p));
        }

        DistMatrix<double> Smhalf_elem(S_elem);
        HPSDSquareRoot(LOWER, Smhalf_elem);
        HPDInverse(LOWER, Smhalf_elem);

        for (int p = 0;p < kv.size();p++)
        {
            int r = kv.key(p)/norb[i];
            int c = kv.key(p)%norb[i];
            int k = (c-cshift)/cstride;
            int j = (r-rshift)/rstride;
            kv.value(p, Smhalf_elem.GetLocal(k, j));
        }

        Smhalf.setRemoteDataByIrrep(irreps, kv);

        for (int p = 0;p < kv2.size();p++)
        {
            int c = kv2.key(p)/norb[i];
            int r = kv2.key(p)%norb[i];
            int k = (c-cshift)/cstride;
            int j = (r-rshift)/rstride;
            kv2.value(p, Smhalf_elem.GetLocal(k, j));
        }

        Smhalf.setRemoteDataByIrrep(irreps, kv2);
    }
}

void ElementalUHF::diagonalizeFock()
{
    const auto& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();

    Tensor<BOUNDED|PGSYMMETRIC|DISTRIBUTED> S  = this->template get   <Tensor<>>("S");
    Tensor<BOUNDED|PGSYMMETRIC|DISTRIBUTED> Fa = this->template get   <Tensor<>>("Fa");
    Tensor<BOUNDED|PGSYMMETRIC|DISTRIBUTED> Fb = this->template get   <Tensor<>>("Fb");
    Tensor<BOUNDED|PGSYMMETRIC|DISTRIBUTED> Ca = this->template gettmp<Tensor<>>("Ca");
    Tensor<BOUNDED|PGSYMMETRIC|DISTRIBUTED> Cb = this->template gettmp<Tensor<>>("Cb");

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        if (norb[i] == 0) continue;

        vector<int> irreps = {i,i};

        KeyValueVector kv(Field::DOUBLE);

        DistMatrix<double> S_elem(norb[i], norb[i]);
        DistMatrix<double> F_elem(norb[i], norb[i]);
        DistMatrix<double> C_elem(norb[i], norb[i]);
        DistMatrix<double> S_tmp(norb[i], norb[i]);
        DistMatrix<double> E_elem;
        DistMatrix<double,STAR,STAR> E_local;

        int cshift = S_elem.ColShift();
        int rshift = S_elem.RowShift();
        int cstride = S_elem.ColStride();
        int rstride = S_elem.RowStride();
        for (int k = 0;k < S_elem.LocalHeight();k++)
        {
            for (int j = 0;j < S_elem.LocalWidth();j++)
            {
                int c = cshift+k*cstride;
                int r = rshift+j*rstride;
                kv.push_back(r*norb[i]+c);
            }
        }

        S.getRemoteDataByIrrep(irreps, kv);
        for (int p = 0;p < kv.size();p++)
        {
            int r = kv.key(p)/norb[i];
            int c = kv.key(p)%norb[i];
            int k = (c-cshift)/cstride;
            int j = (r-rshift)/rstride;
            S_elem.SetLocal(k, j, kv.value<double>(p));
        }

        for (int spin : {0,1})
        {
            auto& F = (spin == 0 ? Fa : Fb);
            auto& C = (spin == 0 ? Ca : Cb);
            auto& E = (spin == 0 ? E_alpha : E_beta);

            F.getRemoteDataByIrrep(irreps, kv);
            for (int p = 0;p < kv.size();p++)
            {
                int r = kv.key(p)/norb[i];
                int c = kv.key(p)%norb[i];
                int k = (c-cshift)/cstride;
                int j = (r-rshift)/rstride;
                F_elem.SetLocal(k, j, kv.value<double>(p));
            }

            S_tmp = S_elem;
            HermitianGenDefEig(El::AXBX, LOWER, F_elem, S_tmp, E_elem, C_elem);

            E_local = E_elem;
            for (int j = 0;j < norb[i];j++) E[i][j] = E_local.GetLocal(j,0);

            for (int p = 0;p < kv.size();p++)
            {
                int r = kv.key(p)/norb[i];
                int c = kv.key(p)%norb[i];
                int k = (c-cshift)/cstride;
                int j = (r-rshift)/rstride;
                kv.value(p, C_elem.GetLocal(k, j));
            }

            C.setRemoteDataByIrrep(irreps, kv);
        }
    }
}

}
}
