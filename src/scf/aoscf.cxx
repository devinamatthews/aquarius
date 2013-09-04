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

    int norb = molecule.getNumOrbitals();

    DistTensor<T>& H = this->template get<DistTensor<T> >("H");
    DistTensor<T>& Da = this->template get<DistTensor<T> >("Da");
    DistTensor<T>& Db = this->template get<DistTensor<T> >("Db");
    DistTensor<T>& Fa = this->template get<DistTensor<T> >("Fa");
    DistTensor<T>& Fb = this->template get<DistTensor<T> >("Fb");

    Arena& arena = H.arena;

    vector<T> focka, fockb;
    vector<T> densa, densb;

    H.getAllData(focka, 0);
    assert(arena.rank != 0 || focka.size() == norb*norb);
    H.getAllData(fockb, 0);
    assert(arena.rank != 0 || focka.size() == norb*norb);

    if (arena.rank != 0)
    {
        focka.resize(norb*norb, (T)0);
        fockb.resize(norb*norb, (T)0);
    }

    Da.getAllData(densa);
    assert(densa.size() == norb*norb);
    Db.getAllData(densb);
    assert(densa.size() == norb*norb);

    vector<T> densab(densa);
    axpy(norb*norb, 1.0, densb.data(), 1, densab.data(), 1);

    const vector<T>& eris = ints.ints;
    const vector<idx4_t>& idxs = ints.idxs;
    size_t neris = eris.size();

    for (size_t n = 0;n < neris;n++)
    {
        int i = idxs[n].i;
        int j = idxs[n].j;
        int k = idxs[n].k;
        int l = idxs[n].l;

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

        bool ieqj = i == j;
        bool keql = k == l;
        bool ijeqkl = min(i,j) == min(k,l) && max(i,j) == max(k,l);

        /*
         * Exchange contribution: Fa(ac) -= Da(bd)*(ab|cd)
         */

        T e = 2.0*eris[n]*(ijeqkl ? 0.5 : 1.0);

        focka[i+k*norb] -= densa[j+l*norb]*e;
        fockb[i+k*norb] -= densb[j+l*norb]*e;
        if (!keql)
        {
            focka[i+l*norb] -= densa[j+k*norb]*e;
            fockb[i+l*norb] -= densb[j+k*norb]*e;
        }
        if (!ieqj)
        {
            focka[j+k*norb] -= densa[i+l*norb]*e;
            fockb[j+k*norb] -= densb[i+l*norb]*e;
            if (!keql)
            {
                focka[j+l*norb] -= densa[i+k*norb]*e;
                fockb[j+l*norb] -= densb[i+k*norb]*e;
            }
        }

        /*
         * Coulomb contribution: Fa(ab) += [Da(cd)+Db(cd)]*(ab|cd)
         */

        e = 2.0*e*(keql ? 0.5 : 1.0)*(ieqj ? 0.5 : 1.0);

        focka[i+j*norb] += densab[k+l*norb]*e;
        fockb[i+j*norb] += densab[k+l*norb]*e;
        focka[k+l*norb] += densab[i+j*norb]*e;
        fockb[k+l*norb] += densab[i+j*norb]*e;
    }

    for (int i = 0;i < norb;i++)
    {
        for (int j = 0;j < i;j++)
        {
            focka[i+j*norb] = 0.5*(focka[i+j*norb]+focka[j+i*norb]);
            focka[j+i*norb] = focka[i+j*norb];
            fockb[i+j*norb] = 0.5*(fockb[i+j*norb]+fockb[j+i*norb]);
            fockb[j+i*norb] = fockb[i+j*norb];
        }
    }

    if (arena.rank == 0)
    {
        arena.Reduce(focka, MPI::SUM);
        arena.Reduce(fockb, MPI::SUM);

        vector<tkv_pair<T> > pairs(norb*norb);

        for (int p = 0;p < norb*norb;p++)
        {
            pairs[p].d = focka[p];
            pairs[p].k = p;
        }

        Fa.writeRemoteData(pairs);

        for (int p = 0;p < norb*norb;p++)
        {
            pairs[p].d = fockb[p];
            pairs[p].k = p;
        }

        Fb.writeRemoteData(pairs);
    }
    else
    {
        arena.Reduce(focka, MPI::SUM, 0);
        arena.Reduce(fockb, MPI::SUM, 0);

        Fa.writeRemoteData();
        Fb.writeRemoteData();
    }
}

INSTANTIATE_SPECIALIZATIONS(AOUHF);
REGISTER_TASK(AOUHF<double>, "aoscf");
