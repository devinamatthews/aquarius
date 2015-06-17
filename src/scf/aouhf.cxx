#include "aouhf.hpp"

using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;

namespace aquarius
{
namespace scf
{

template <class WhichUHF>
AOUHF<WhichUHF>::AOUHF(const string& name, Config& config)
: WhichUHF(name, config)
{
    for (vector<Product>::iterator i = this->products.begin();i != this->products.end();++i)
    {
        i->addRequirement(Requirement("eri", "I"));
    }
}

template <class WhichUHF>
void AOUHF<WhichUHF>::buildFock()
{
    const Molecule& molecule =this->template get<Molecule>("molecule");
    const ERI& ints = this->template get<ERI>("I");

    const vector<int>& norb = molecule.getNumOrbitals();
    int nirrep = molecule.getGroup().getNumIrreps();

    vector<int> irrep;
    for (int i = 0;i < nirrep;i++) irrep += vector<int>(norb[i],i);

    vector<int> start(nirrep,0);
    for (int i = 1;i < nirrep;i++) start[i] = start[i-1]+norb[i-1];

    Tensor<BOUNDED|PGSYMMETRIC> H  = this->template get<Tensor<>>("H");
    Tensor<BOUNDED|PGSYMMETRIC> Da = this->template get<Tensor<>>("Da");
    Tensor<BOUNDED|PGSYMMETRIC> Db = this->template get<Tensor<>>("Db");
    Tensor<BOUNDED|PGSYMMETRIC> Fa = this->template get<Tensor<>>("Fa");
    Tensor<BOUNDED|PGSYMMETRIC> Fb = this->template get<Tensor<>>("Fb");

    vector<vector<double>> focka(nirrep), fockb(nirrep);
    vector<vector<double>> densa(nirrep), densb(nirrep);
    vector<vector<double>> densab(nirrep);

    for (int i = 0;i < nirrep;i++)
    {
        vector<int> irreps = {i,i};

        focka[i].resize(norb[i]*norb[i]);
        KeyValueVector kv = H.getAllDataByIrrep(irreps);
        for (int j = 0;j < kv.size();j++)
        {
            focka[i][kv.key(j)] = kv.value<double>(j);
        }
        if (this->arena().rank != 0)
        {
            fill(focka[i].begin(), focka[i].end(), 0);
        }
        fockb[i] = focka[i];

        Da.getAllDataByIrrep(irreps, kv);
        assert(kv.size() == norb[i]*norb[i]);
        for (int j = 0;j < kv.size();j++)
        {
            densa[i][kv.key(j)] = kv.value<double>(j);
        }

        Db.getAllDataByIrrep(irreps, kv);
        assert(kv.size() == norb[i]*norb[i]);
        for (int j = 0;j < kv.size();j++)
        {
            densb[i][kv.key(j)] = kv.value<double>(j);
        }

        densab[i] = densa[i];
        axpy(norb[i]*norb[i], 1.0, densb[i].data(), 1, densab[i].data(), 1);
    }

    auto& eris = ints.ints;
    auto& idxs = ints.idxs;
    size_t neris = eris.size();
    assert(eris.size() == idxs.size());

    int64_t flops = 0;
    #pragma omp parallel reduction(+:flops)
    {
        int nt = omp_get_num_threads();
        int tid = omp_get_thread_num();
        size_t n0 = (neris*tid)/nt;
        size_t n1 = (neris*(tid+1))/nt;

        vector<vector<double>> focka_local(nirrep);
        vector<vector<double>> fockb_local(nirrep);

        for (int i = 0;i < nirrep;i++)
        {
            focka_local[i].resize(norb[i]*norb[i]);
            fockb_local[i].resize(norb[i]*norb[i]);
        }

        auto iidx = idxs.begin()+n0;
        auto iend = idxs.begin()+n1;
        auto iint = eris.begin()+n0;
        for (;iidx != iend;++iidx, ++iint)
        {
            int irri = irrep[iidx->i];
            int irrj = irrep[iidx->j];
            int irrk = irrep[iidx->k];
            int irrl = irrep[iidx->l];

            if (irri != irrj && irri != irrk && irri != irrl) continue;

            int i = iidx->i-start[irri];
            int j = iidx->j-start[irrj];
            int k = iidx->k-start[irrk];
            int l = iidx->l-start[irrl];

            /*
            if (i < j)
            {
                swap(i, j);
            }
            if (k < l)
            {
                swap(k, l);
            }
            if (i < k || (i == k && j < l))
            {
                swap(i, k);
                swap(j, l);
            }
            printf("%d %d %d %d %25.15e\n", i+1, j+1, k+1, l+1, eris[n].value);
            */

            bool ieqj = i == j && irri == irrj;
            bool keql = k == l && irrk == irrl;
            bool ijeqkl = i == k && irri == irrk && j == l && irrj == irrl;

            //cout << irri << " " << irrj << " " << irrk << " " << irrl << " "
            //        << i << " " << j << " " << k << " " << l << endl;

            /*
             * Exchange contribution: Fa(ac) -= Da(bd)*(ab|cd)
             */

            double e = 2.0*(*iint)*(ijeqkl ? 0.5 : 1.0);

            if (irri == irrk && irrj == irrl)
            {
                flops += 4;;
                focka_local[irri][i+k*norb[irri]] -= densa[irrj][j+l*norb[irrj]]*e;
                fockb_local[irri][i+k*norb[irri]] -= densb[irrj][j+l*norb[irrj]]*e;
            }
            if (!keql && irri == irrl && irrj == irrk)
            {
                flops += 4;;
                focka_local[irri][i+l*norb[irri]] -= densa[irrj][j+k*norb[irrj]]*e;
                fockb_local[irri][i+l*norb[irri]] -= densb[irrj][j+k*norb[irrj]]*e;
            }
            if (!ieqj)
            {
                if (irri == irrl && irrj == irrk)
                {
                    flops += 4;;
                    focka_local[irrj][j+k*norb[irrj]] -= densa[irri][i+l*norb[irri]]*e;
                    fockb_local[irrj][j+k*norb[irrj]] -= densb[irri][i+l*norb[irri]]*e;
                }
                if (!keql && irri == irrk && irrj == irrl)
                {
                    flops += 4;;
                    focka_local[irrj][j+l*norb[irrj]] -= densa[irri][i+k*norb[irri]]*e;
                    fockb_local[irrj][j+l*norb[irrj]] -= densb[irri][i+k*norb[irri]]*e;
                }
            }

            /*
             * Coulomb contribution: Fa(ab) += [Da(cd)+Db(cd)]*(ab|cd)
             */

            e = 2.0*e*(keql ? 0.5 : 1.0)*(ieqj ? 0.5 : 1.0);

            if (irri == irrj && irrk == irrl)
            {
                flops += 6;;
                focka_local[irri][i+j*norb[irri]] += densab[irrk][k+l*norb[irrk]]*e;
                fockb_local[irri][i+j*norb[irri]] += densab[irrk][k+l*norb[irrk]]*e;
                focka_local[irrk][k+l*norb[irrk]] += densab[irri][i+j*norb[irri]]*e;
                fockb_local[irrk][k+l*norb[irrk]] += densab[irri][i+j*norb[irri]]*e;
            }
        }

        #pragma omp critical
        {
            for (int irr = 0;irr < nirrep;irr++)
            {
                flops += 2*norb[irr]*norb[irr];
                axpy(norb[irr]*norb[irr], 1.0, focka_local[irr].data(), 1, focka[irr].data(), 1);
                axpy(norb[irr]*norb[irr], 1.0, fockb_local[irr].data(), 1, fockb[irr].data(), 1);
            }
        }
    }

    for (int irr = 0;irr < nirrep;irr++)
    {
        for (int i = 0;i < norb[irr];i++)
        {
            for (int j = 0;j < i;j++)
            {
                focka[irr][i+j*norb[irr]] = 0.5*(focka[irr][i+j*norb[irr]]+focka[irr][j+i*norb[irr]]);
                focka[irr][j+i*norb[irr]] = focka[irr][i+j*norb[irr]];
                fockb[irr][i+j*norb[irr]] = 0.5*(fockb[irr][i+j*norb[irr]]+fockb[irr][j+i*norb[irr]]);
                fockb[irr][j+i*norb[irr]] = fockb[irr][i+j*norb[irr]];
            }
        }
    }

    for (int i = 0;i < nirrep;i++)
    {
        vector<int> irreps = {i,i};

        this->arena().comm().Reduce(focka[i], MPI_SUM, 0);
        this->arena().comm().Reduce(fockb[i], MPI_SUM, 0);

        KeyValueVector kv(Field::DOUBLE, norb[i]*norb[i]);

        for (int p = 0;p < norb[i]*norb[i];p++)
        {
            kv.key(p) = p;
            kv.value(p, focka[i][p]);
        }

        Fa.setDataByIrrep(irreps, kv);

        for (int p = 0;p < norb[i]*norb[i];p++)
        {
            kv.key(p) = p;
            kv.value(p, fockb[i][p]);
        }

        Fb.setDataByIrrep(irreps, kv);
    }
}

}
}

static const char* spec = R"(

    frozen_core?
        bool false,
    convergence?
        double 1e-12,
    max_iterations?
        int 150,
    conv_type?
        enum { MAXE, RMSE, MAE },
    diis?
    {
        damping?
            double 0.0,
        start?
            int 8,
        order?
            int 6,
        jacobi?
            bool false
    }

)";

template class aquarius::scf::AOUHF<aquarius::scf::LocalUHF>;
REGISTER_TASK(CONCAT(aquarius::scf::AOUHF<aquarius::scf::LocalUHF>), "localaoscf",spec);

#if HAVE_ELEMENTAL
template class aquarius::scf::AOUHF<aquarius::scf::ElementalUHF>;
REGISTER_TASK(CONCAT(aquarius::scf::AOUHF<aquarius::scf::ElementalUHF>), "elementalaoscf",spec);
#endif
