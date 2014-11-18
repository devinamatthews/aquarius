/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "aoscf.hpp"

#include "util/blas.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::scf;
using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;

template <typename T>
AOUHF<T>::AOUHF(const string& name, const Config& config)
: UHF<T>("aoscf", name, config)
{
    for (vector<Product>::iterator i = this->products.begin();i != this->products.end();++i)
    {
        i->addRequirement(Requirement("eri", "I"));
    }
}

template <typename T>
void AOUHF<T>::buildFock()
{
    const Molecule& molecule =this->template get<Molecule>("molecule");
    const ERI& ints = this->template get<ERI>("I");

    const vector<int>& norb = molecule.getNumOrbitals();
    int nirrep = molecule.getGroup().getNumIrreps();

    vector<int> irrep;
    for (int i = 0;i < nirrep;i++) irrep += vector<int>(norb[i],i);

    vector<int> start(nirrep,0);
    for (int i = 1;i < nirrep;i++) start[i] = start[i-1]+norb[i-1];

    SymmetryBlockedTensor<T>& H = this->template get<SymmetryBlockedTensor<T> >("H");
    SymmetryBlockedTensor<T>& Da = this->template get<SymmetryBlockedTensor<T> >("Da");
    SymmetryBlockedTensor<T>& Db = this->template get<SymmetryBlockedTensor<T> >("Db");
    SymmetryBlockedTensor<T>& Fa = this->template get<SymmetryBlockedTensor<T> >("Fa");
    SymmetryBlockedTensor<T>& Fb = this->template get<SymmetryBlockedTensor<T> >("Fb");

    Arena& arena = H.arena;

    vector<vector<T> > focka(nirrep), fockb(nirrep);
    vector<vector<T> > densa(nirrep), densb(nirrep);
    vector<vector<T> > densab(nirrep);

    for (int i = 0;i < nirrep;i++)
    {
        vector<int> irreps(2,i);

        if (arena.rank == 0)
        {
            H.getAllData(irreps, focka[i], 0);
            assert(focka[i].size() == norb[i]*norb[i]);
            fockb[i] = focka[i];
        }
        else
        {
            H.getAllData(irreps, 0);
            focka[i].resize(norb[i]*norb[i], (T)0);
            fockb[i].resize(norb[i]*norb[i], (T)0);
        }

        Da.getAllData(irreps, densa[i]);
        assert(densa[i].size() == norb[i]*norb[i]);
        Db.getAllData(irreps, densb[i]);
        assert(densa[i].size() == norb[i]*norb[i]);

        densab[i] = densa[i];
        PROFILE_FLOPS(norb[i]*norb[i]);
        axpy(norb[i]*norb[i], 1.0, densb[i].data(), 1, densab[i].data(), 1);

        if (Da.norm(2) > 1e-10)
        {
            //fill(focka[i].begin(), focka[i].end(), 0.0);
            //fill(fockb[i].begin(), fockb[i].end(), 0.0);
        }
    }

    const vector<T>& eris = ints.ints;
    const vector<idx4_t>& idxs = ints.idxs;
    size_t neris = eris.size();
    assert(eris.size() == idxs.size());

    int64_t flops = 0;
    #pragma omp parallel reduction(+:flops)
    {
        int nt = omp_get_num_threads();
        int tid = omp_get_thread_num();
        size_t n0 = (neris*tid)/nt;
        size_t n1 = (neris*(tid+1))/nt;

        vector<vector<T> > focka_local(nirrep);
        vector<vector<T> > fockb_local(nirrep);

        for (int i = 0;i < nirrep;i++)
        {
            focka_local[i].resize(norb[i]*norb[i], (T)0);
            fockb_local[i].resize(norb[i]*norb[i], (T)0);
        }

        for (size_t n = n0;n < n1;n++)
        {
            int irri = irrep[idxs[n].i];
            int irrj = irrep[idxs[n].j];
            int irrk = irrep[idxs[n].k];
            int irrl = irrep[idxs[n].l];

            if (irri != irrj && irri != irrk && irri != irrl) continue;

            int i = idxs[n].i-start[irri];
            int j = idxs[n].j-start[irrj];
            int k = idxs[n].k-start[irrk];
            int l = idxs[n].l-start[irrl];

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

            T e = 2.0*eris[n]*(ijeqkl ? 0.5 : 1.0);

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
                axpy(norb[irr]*norb[irr], (T)1, focka_local[irr].data(), 1, focka[irr].data(), 1);
                axpy(norb[irr]*norb[irr], (T)1, fockb_local[irr].data(), 1, fockb[irr].data(), 1);
            }
        }
    }
    PROFILE_FLOPS(flops);

    for (int irr = 0;irr < nirrep;irr++)
    {
        PROFILE_FLOPS(2*norb[irr]*(norb[irr]-1));
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
        vector<int> irreps(2,i);

        if (arena.rank == 0)
        {
            PROFILE_FLOPS(2*norb[i]*norb[i]);
            arena.Reduce(focka[i], MPI::SUM);
            arena.Reduce(fockb[i], MPI::SUM);

            vector<tkv_pair<T> > pairs(norb[i]*norb[i]);

            for (int p = 0;p < norb[i]*norb[i];p++)
            {
                pairs[p].d = focka[i][p];
                pairs[p].k = p;
            }

            Fa.writeRemoteData(irreps, pairs);

            for (int p = 0;p < norb[i]*norb[i];p++)
            {
                pairs[p].d = fockb[i][p];
                pairs[p].k = p;
            }

            Fb.writeRemoteData(irreps, pairs);
        }
        else
        {
            PROFILE_FLOPS(2*norb[i]*norb[i]);
            arena.Reduce(focka[i], MPI::SUM, 0);
            arena.Reduce(fockb[i], MPI::SUM, 0);

            Fa.writeRemoteData(irreps);
            Fb.writeRemoteData(irreps);
        }
    }
}

INSTANTIATE_SPECIALIZATIONS(AOUHF);
REGISTER_TASK(AOUHF<double>, "aoscf");
