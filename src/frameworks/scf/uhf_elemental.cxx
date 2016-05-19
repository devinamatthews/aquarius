#include "uhf_elemental.hpp"

using namespace El;

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
ElementalUHF<T>::ElementalUHF(const string& name, Config& config)
: UHF<T>(name, config) {}

template <typename T>
void ElementalUHF<T>::calcSMinusHalf()
{
    const auto& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();

    auto& S = this->template get<SymmetryBlockedTensor<T>>("S");
    auto& Smhalf = this->template gettmp<SymmetryBlockedTensor<T>>("S^-1/2");

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        //cout << "S " << (i+1) << endl;
        //vector<T> vals;
        //S({i,i}).getAllData(vals);
        //printmatrix(norb[i], norb[i], vals.data(), 6, 3, 108);
    }

    SymmetryBlockedTensor<T> Smhalf2(Smhalf);
    SymmetryBlockedTensor<T> tmp(Smhalf);
    SymmetryBlockedTensor<T> tmp2(Smhalf);

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        if (norb[i] == 0) continue;

        vector<int> irreps = {i,i};
        vector<real_type_t<T>> E(norb[i]);

        DistMatrix<T> S_elem(norb[i], norb[i]);
        vector<tkv_pair<T>> pairs, pairs2;

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
                pairs.emplace_back(r*norb[i]+c, 0);
                if (r == c) continue;
                pairs2.emplace_back(c*norb[i]+r, 0);
            }
        }

        S.getRemoteData(irreps, pairs);

        for (auto& p : pairs)
        {
            int r = p.k/norb[i];
            int c = p.k%norb[i];
            int k = (c-cshift)/cstride;
            int j = (r-rshift)/rstride;
            S_elem.SetLocal(k, j, p.d);
        }

        DistMatrix<T> Smhalf_elem(S_elem);
        HPSDSquareRoot(LOWER, Smhalf_elem);
        HPDInverse(LOWER, Smhalf_elem);

        for (auto& p : pairs)
        {
            int r = p.k/norb[i];
            int c = p.k%norb[i];
            int k = (c-cshift)/cstride;
            int j = (r-rshift)/rstride;
            p.d = Smhalf_elem.GetLocal(k, j);
        }

        Smhalf.writeRemoteData(irreps, pairs);

        for (auto& p : pairs2)
        {
            int c = p.k/norb[i];
            int r = p.k%norb[i];
            int k = (c-cshift)/cstride;
            int j = (r-rshift)/rstride;
            p.d = Smhalf_elem.GetLocal(k, j);
        }

        Smhalf.writeRemoteData(irreps, pairs2);

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

            pairs.resize(norb[i]*norb[i]);
            for (int j = 0;j < norb[i]*norb[i];j++)
            {
                pairs[j].k = j;
                pairs[j].d = smhalf[j];
            }

            Smhalf2.writeRemoteData(irreps, pairs);
        }
        else
        {
            S.getAllData(irreps, 0);
            Smhalf2.writeRemoteData(irreps);
        }
    }

    tmp["il"] = S["ik"]*Smhalf["kl"];
    tmp2["ij"] = tmp["il"]*Smhalf2["lj"];
    double nrm = tmp2.norm(1)-sum(norb);
    if (S.arena.rank == 0)
    printf("S^-1/2 diff: %g\n", nrm);
}

template <typename T>
void ElementalUHF<T>::diagonalizeFock()
{
    const auto& molecule = this->template get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();

    auto& S  = this->template get   <SymmetryBlockedTensor<T>>("S");
    auto& Fa = this->template get   <SymmetryBlockedTensor<T>>("Fa");
    auto& Fb = this->template get   <SymmetryBlockedTensor<T>>("Fb");
    auto& Ca = this->template gettmp<SymmetryBlockedTensor<T>>("Ca");
    auto& Cb = this->template gettmp<SymmetryBlockedTensor<T>>("Cb");

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

        vector<tkv_pair<T>> pairs;

        DistMatrix<T> S_elem(norb[i], norb[i]);
        DistMatrix<T> F_elem(norb[i], norb[i]);
        DistMatrix<T> C_elem(norb[i], norb[i]);
        DistMatrix<T> S_tmp(norb[i], norb[i]);
        DistMatrix<T> E_elem;
        DistMatrix<T,STAR,STAR> E_local;

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
                pairs.emplace_back(r*norb[i]+c, 0);
            }
        }

        S.getRemoteData(irreps, pairs);
        for (int p = 0;p < pairs.size();p++)
        {
            int r = pairs[p].k/norb[i];
            int c = pairs[p].k%norb[i];
            int k = (c-cshift)/cstride;
            int j = (r-rshift)/rstride;
            S_elem.SetLocal(k, j, pairs[p].d);
        }

        vector<T> s;
        if (S.arena.rank == 0)
        {
            S.getAllData(irreps, s, 0);
            assert(s.size() == norb[i]*norb[i]);
        }
        else
        {
            S.getAllData(irreps, 0);
        }

        SymmetryBlockedTensor<T> Ca2(Ca);
        SymmetryBlockedTensor<T> Cb2(Cb);
        vector<vector<real_type_t<T>>> E_alpha2(E_alpha);
        vector<vector<real_type_t<T>>> E_beta2(E_beta);

        for (int spin : {0,1})
        {
            auto& F = (spin == 0 ? Fa : Fb);
            auto& C = (spin == 0 ? Ca : Cb);
            auto& C2 = (spin == 0 ? Ca2 : Cb2);
            auto& E = (spin == 0 ? E_alpha : E_beta);
            auto& E2 = (spin == 0 ? E_alpha2 : E_beta2);

            F.getRemoteData(irreps, pairs);
            for (int p = 0;p < pairs.size();p++)
            {
                int r = pairs[p].k/norb[i];
                int c = pairs[p].k-norb[i]*r;
                int k = (c-cshift)/cstride;
                int j = (r-rshift)/rstride;
                F_elem.SetLocal(k, j, pairs[p].d);
            }

            S_tmp = S_elem;
            HermitianGenDefEig(El::AXBX, LOWER, F_elem, S_tmp, E_elem, C_elem);

            E_local = E_elem;
            for (int j = 0;j < norb[i];j++) E[i][j] = E_local.GetLocal(j,0);

            for (int p = 0;p < pairs.size();p++)
            {
                int r = pairs[p].k/norb[i];
                int c = pairs[p].k-norb[i]*r;
                int k = (c-cshift)/cstride;
                int j = (r-rshift)/rstride;
                pairs[p].d = C_elem.GetLocal(k, j);
            }

            C.writeRemoteData(irreps, pairs);

            if (S.arena.rank == 0)
            {
                vector<T> fock;
                vector<T> tmp(s);

                F.getAllData(irreps, fock, 0);
                assert(fock.size() == norb[i]*norb[i]);
                //PROFILE_FLOPS(9*norb[i]*norb[i]*norb[i]);
                int info = hegv(LAWrap::AXBX, 'V', 'U', norb[i], fock.data(), norb[i], tmp.data(), norb[i], E2[i].data());
                assert(info == 0);
                S.arena.comm().Bcast(E2[i], 0);

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

                vector<tkv_pair<T>> pairs2(norb[i]*norb[i]);
                for (int j = 0;j < norb[i]*norb[i];j++)
                {
                    pairs2[j].k = j;
                    pairs2[j].d = fock[j];
                }

                C2.writeRemoteData(irreps, pairs2);
            }
            else
            {
                F.getAllData(irreps, 0);
                S.arena.comm().Bcast(E2[i], 0);
                C2.writeRemoteData(irreps);
            }
        }

        SymmetryBlockedTensor<T> C2C(Ca);
        SymmetryBlockedTensor<T> tmp(Ca);

        tmp["kj"] = S["kl"]*Ca2["lj"];
        C2C["ij"] = Ca["ki"]*tmp["kj"];

        double nrm = C2C.norm(1)-sum(norb);
        if (S.arena.rank == 0)
        {
            printf("C2Ca: %g\n", nrm);

            double dE = 0.0;
            for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
            {
                for (int j = 0;j < norb[i];j++)
                {
                    dE += aquarius::abs(E_alpha[i][j]-E_alpha2[i][j]);
                    dE += aquarius::abs(E_beta[i][j]-E_beta2[i][j]);
                }
            }
            printf("E diff: %g\n", sqrt(dE));
        }
    }
}

INSTANTIATE_SPECIALIZATIONS(ElementalUHF);

}
}
