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

#include "scf.hpp"

#include <algorithm>
#include <vector>
#include <limits>

using namespace std;
using namespace elem;
using namespace MPI;
using namespace aquarius::slide;
using namespace aquarius::input;
using namespace aquarius::diis;
using namespace aquarius::tensor;

namespace aquarius
{
namespace scf
{

UHF::UHF(tCTF_World<double>& ctf, const Molecule& molecule, const Config& config)
: Iterative(config),
  Distributed<double>(ctf),
  molecule(molecule),
  norb(molecule.getNumOrbitals()),
  nalpha(molecule.getNumAlphaElectrons()),
  nbeta(molecule.getNumBetaElectrons()),
  damping(config.get<double>("damping")),
  Ea(norb),
  Eb(norb),
  diis(config.get("diis"), 2),
  grid(comm),
  C_elem(norb, norb, grid),
  S_elem(norb, norb, grid),
  F_elem(norb, norb, grid),
  E_elem(norb, 1   , grid)
{
    energy = molecule.getNuclearRepulsion();

    int shapeNN[] = {NS,NS};

    int sizenn[] = {norb,norb};
    Fa = new DistTensor<double>(ctf, 2, sizenn, shapeNN, false);
    Fb = new DistTensor<double>(ctf, 2, sizenn, shapeNN, false);
    dF = new DistTensor<double>(ctf, 2, sizenn, shapeNN, false);
    int sizenO[] = {norb,nalpha};
    int sizeno[] = {norb,nbeta};
    Ca_occ = new DistTensor<double>(ctf, 2, sizenO, shapeNN, false);
    Cb_occ = new DistTensor<double>(ctf, 2, sizeno, shapeNN, false);
    int sizenV[] = {norb,norb-nalpha};
    int sizenv[] = {norb,norb-nbeta};
    Ca_vrt = new DistTensor<double>(ctf, 2, sizenV, shapeNN, false);
    Cb_vrt = new DistTensor<double>(ctf, 2, sizenv, shapeNN, false);
    Da = new DistTensor<double>(ctf, 2, sizenn, shapeNN, true);
    Db = new DistTensor<double>(ctf, 2, sizenn, shapeNN, true);
    dDa = new DistTensor<double>(ctf, 2, sizenn, shapeNN, false);
    dDb = new DistTensor<double>(ctf, 2, sizenn, shapeNN, false);
    S = new DistTensor<double>(ctf, 2, sizenn, shapeNN, false);
    Smhalf = new DistTensor<double>(ctf, 2, sizenn, shapeNN, false);
    H = new DistTensor<double>(ctf, 2, sizenn, shapeNN, false);

    getOverlap();
    get1eHamiltonian();
}

UHF::~UHF()
{
    delete Fa;
    delete Fb;
    delete dF;
    delete Ca_occ;
    delete Cb_occ;
    delete Ca_vrt;
    delete Cb_vrt;
    delete Da;
    delete Db;
    delete dDa;
    delete dDb;
    delete S;
    delete Smhalf;
    delete H;
}

void UHF::getOverlap()
{
    Context context;

    vector<kv_pair> pairs;

    int pid = comm.Get_rank();
    int nproc = comm.Get_size();
    int block = 0;
    for (Molecule::const_shell_iterator i = molecule.getShellsBegin();i != molecule.getShellsEnd();++i)
    {
        for (Molecule::const_shell_iterator j = molecule.getShellsBegin();j != molecule.getShellsEnd();++j)
        {
            if (i < j) continue;

            if (block%nproc == pid)
            {
                context.calcOVI(1.0, 0.0, *i, *j);

                size_t nint = context.getNumIntegrals();
                vector<double> ints(nint);
                vector<idx2_t> idxs(nint);

                size_t nproc = context.process1eInts(nint, ints.data(), idxs.data(), -1.0);
                for (int i = 0;i < nproc;i++)
                {
                    //printf("%d %d %25.15f\n", idxs[i].i+1, idxs[i].j+1, ints[i]);

                    pairs.push_back(kv_pair(idxs[i].i*norb+idxs[i].j, ints[i]));
                    if (idxs[i].i != idxs[i].j)
                    {
                        pairs.push_back(kv_pair(idxs[i].j*norb+idxs[i].i, ints[i]));
                    }
                }
            }

            block++;
        }
    }

    S->writeRemoteData(pairs.size(), pairs.data());

    pairs.clear();

    int cshift = S_elem.ColShift();
    int rshift = S_elem.RowShift();
    int cstride = S_elem.ColStride();
    int rstride = S_elem.RowStride();
    for (int i = 0;i < S_elem.LocalHeight();i++)
    {
        for (int j = 0;j < S_elem.LocalWidth();j++)
        {
            int c = cshift+i*cstride;
            int r = rshift+j*rstride;

            pairs.push_back(kv_pair(r*norb+c,0));
        }
    }

    S->getRemoteData(pairs.size(), pairs.data());

    for (int p = 0;p < pairs.size();p++)
    {
        int r = pairs[p].k/norb;
        int c = pairs[p].k-norb*r;
        int i = (c-cshift)/cstride;
        int j = (r-rshift)/rstride;
        S_elem.SetLocal(i, j, pairs[p].d);
    }

    DistMatrix<double> Smhalf_elem(S_elem);
    HPSDSquareRoot(UPPER, Smhalf_elem);
    HPDInverse(UPPER, Smhalf_elem);

    for (int p = 0;p < pairs.size();p++)
    {
        int r = pairs[p].k/norb;
        int c = pairs[p].k-norb*r;
        if (r > c) swap(r,c);
        int i = (c-cshift)/cstride;
        int j = (r-rshift)/rstride;
        pairs[p].d = Smhalf_elem.GetLocal(i, j);
    }

    Smhalf->writeRemoteData(pairs.size(), pairs.data());
}

void UHF::get1eHamiltonian()
{
    Context context;

    vector<Center> centers;
    for (vector<Atom>::const_iterator a = molecule.getAtomsBegin();a != molecule.getAtomsEnd();++a)
    {
        centers.push_back(a->getCenter());
    }

    vector<kv_pair> pairs;

    int pid = comm.Get_rank();
    int nproc = comm.Get_size();
    int block = 0;
    for (Molecule::const_shell_iterator i = molecule.getShellsBegin();i != molecule.getShellsEnd();++i)
    {
        for (Molecule::const_shell_iterator j = molecule.getShellsBegin();j != molecule.getShellsEnd();++j)
        {
            if (i < j) continue;

            if (block%nproc == pid)
            {
                context.calcKEI(1.0, 0.0, *i, *j);
                context.calcNAI(1.0, 1.0, *i, *j, centers.data(), centers.size());

                size_t nint = context.getNumIntegrals();
                vector<double> ints(nint);
                vector<idx2_t> idxs(nint);

                size_t nproc = context.process1eInts(nint, ints.data(), idxs.data(), -1.0);
                for (int i = 0;i < nproc;i++)
                {
                    //printf("%d %d %25.15f\n", idxs[i].i+1, idxs[i].j+1, ints[i]);

                    pairs.push_back(kv_pair(idxs[i].i*norb+idxs[i].j, ints[i]));
                    if (idxs[i].i != idxs[i].j)
                    {
                        pairs.push_back(kv_pair(idxs[i].j*norb+idxs[i].i, ints[i]));
                    }
                }
            }

            block++;
        }
    }

    H->writeRemoteData(pairs.size(), pairs.data());
}

void UHF::diagonalizeFock()
{
    vector<kv_pair> pairs;

    int cshift = S_elem.ColShift();
    int rshift = S_elem.RowShift();
    int cstride = S_elem.ColStride();
    int rstride = S_elem.RowStride();
    for (int i = 0;i < F_elem.LocalHeight();i++)
    {
        for (int j = 0;j < F_elem.LocalWidth();j++)
        {
            int c = cshift+i*cstride;
            int r = rshift+j*rstride;

            pairs.push_back(kv_pair(r*norb+c,0));
        }
    }

    Fa->getRemoteData(pairs.size(), pairs.data());

    for (int p = 0;p < pairs.size();p++)
    {
        int r = pairs[p].k/norb;
        int c = pairs[p].k-norb*r;
        int i = (c-cshift)/cstride;
        int j = (r-rshift)/rstride;
        F_elem.SetLocal(i, j, pairs[p].d);
    }

    {
        DistMatrix<double> S_tmp(S_elem);
        HermitianGenDefiniteEig(AXBX, LOWER, F_elem, S_tmp, E_elem, C_elem);
        SortEig(E_elem, C_elem);

        DistMatrix<double,STAR,STAR> E_local(E_elem);
        for (int i = 0;i < norb;i++) Ea[i] = E_local.GetLocal(i,0);

        vector<kv_pair> pairs_occ;
        vector<kv_pair> pairs_vrt;

        for (int i = 0;i < C_elem.LocalHeight();i++)
        {
            for (int j = 0;j < C_elem.LocalWidth();j++)
            {
                int c = cshift+i*cstride;
                int r = rshift+j*rstride;

                if (r < nalpha)
                {
                    pairs_occ.push_back(kv_pair(r*norb+c, C_elem.GetLocal(i,j)));
                }
                else
                {
                    pairs_vrt.push_back(kv_pair((r-nalpha)*norb+c, C_elem.GetLocal(i,j)));
                }
            }
        }

        Ca_occ->writeRemoteData(pairs_occ.size(), pairs_occ.data());
        Ca_vrt->writeRemoteData(pairs_vrt.size(), pairs_vrt.data());
    }

    Fb->getRemoteData(pairs.size(), pairs.data());

    for (int p = 0;p < pairs.size();p++)
    {
        int r = pairs[p].k/norb;
        int c = pairs[p].k-norb*r;
        int i = (c-cshift)/cstride;
        int j = (r-rshift)/rstride;
        F_elem.SetLocal(i, j, pairs[p].d);
    }

    {
        DistMatrix<double> S_tmp(S_elem);
        HermitianGenDefiniteEig(AXBX, LOWER, F_elem, S_tmp, E_elem, C_elem);
        SortEig(E_elem, C_elem);

        DistMatrix<double,STAR,STAR> E_local(E_elem);
        for (int i = 0;i < norb;i++) Eb[i] = E_local.GetLocal(i,0);

        vector<kv_pair> pairs_occ;
        vector<kv_pair> pairs_vrt;

        for (int i = 0;i < C_elem.LocalHeight();i++)
        {
            for (int j = 0;j < C_elem.LocalWidth();j++)
            {
                int c = cshift+i*cstride;
                int r = rshift+j*rstride;

                if (r < nbeta)
                {
                    pairs_occ.push_back(kv_pair(r*norb+c, C_elem.GetLocal(i,j)));
                }
                else
                {
                    pairs_vrt.push_back(kv_pair((r-nbeta)*norb+c, C_elem.GetLocal(i,j)));
                }
            }
        }

        Cb_occ->writeRemoteData(pairs_occ.size(), pairs_occ.data());
        Cb_vrt->writeRemoteData(pairs_vrt.size(), pairs_vrt.data());
    }

    fixPhase(*Ca_occ);
    fixPhase(*Cb_occ);
    fixPhase(*Ca_vrt);
    fixPhase(*Cb_vrt);

    //cout << "Eigenvalues" << endl;
    //for (int i = 0;i < norb;i++)
    //    cout << i << ' ' << Ea[i] << ' ' << Eb[i] << endl;
    //cout << endl;
}

void UHF::fixPhase(DistTensor<double>& C)
{
    int rank = comm.Get_rank();
    int np = comm.Get_size();
    int nr = C.getLengths()[1];

    vector<kv_pair> pairs(norb);

    for (int b = 0;b*np < nr;b++)
    {
        int r = b*np+rank;

        if (r < nr)
        {
            for (int i = 0;i < norb;i++) pairs[i].k = i+r*norb;

            C.getRemoteData(norb, pairs.data());

            sort(pairs.begin(), pairs.end());
            int sign = 0;
            for (int i = 0;i < norb;i++)
            {
                if (sign == 0 && abs(pairs[i].d) > 1e-12)
                {
                    sign = 1;
                    if (pairs[i].d < 0)
                    {
                        pairs[i].d = -pairs[i].d;
                        sign = -1;
                    }
                }
                else
                {
                    pairs[i].d = sign*pairs[i].d;
                }
            }

            C.writeRemoteData(norb, pairs.data());
        }
        else
        {
            C.getRemoteData(0, NULL);
            C.writeRemoteData(0, NULL);
        }
    }
}

void UHF::_iterate()
{
    buildFock();
    DIISExtrap();
    calcEnergy();
    diagonalizeFock();
    calcDensity();

    switch (convtype)
    {
        case MAX_ABS:
            conv = max(dDa->reduce(CTF_OP_MAXABS), dDb->reduce(CTF_OP_MAXABS));
            break;
        case RMSD:
            conv = sqrt((dDa->reduce(CTF_OP_SQNRM2)+dDb->reduce(CTF_OP_SQNRM2))/(2*norb*norb));
            break;
        case MAD:
            conv = (dDa->reduce(CTF_OP_SUMABS)+dDb->reduce(CTF_OP_SUMABS))/(2*norb*norb);
            break;
    }
}

void UHF::DIISExtrap()
{
    /*
     * Generate the residual:
     *
     *  R = FDS - SDF
     *
     * Since F and D commute at convergence, we should have [F,D] = FD - DF = 0,
     * although I'm not totally sure why the S is thrown in there.
     *
     * Then, convert to the orthonormal basis:
     *
     *  ~    -1/2    -1/2
     *  R = S     R S
     *                                    ~   ~
     * This is so that the inner product <R_i|R_j> is what we expect it to be.
     */
    {
        DistTensor<double> tmp1(*Fa);
        DistTensor<double> tmp2(*Fa);

         tmp1["ab"]  =     (*Fa)["ac"]*    (*Da)["cb"];
         tmp2["ab"]  =      tmp1["ac"]*     (*S)["cb"];
         tmp1["ab"]  =      (*S)["ac"]*    (*Da)["cb"];
         tmp2["ab"] -=      tmp1["ac"]*    (*Fa)["cb"];
         tmp1["ab"]  = (*Smhalf)["ac"]*     tmp2["cb"];
        (*dF)["ab"]  =      tmp1["ac"]*(*Smhalf)["cb"];

         tmp1["ab"]  =     (*Fb)["ac"]*    (*Db)["cb"];
         tmp2["ab"]  =      tmp1["ac"]*     (*S)["cb"];
         tmp1["ab"]  =      (*S)["ac"]*    (*Db)["cb"];
         tmp2["ab"] -=      tmp1["ac"]*    (*Fb)["cb"];
         tmp1["ab"]  = (*Smhalf)["ac"]*     tmp2["cb"];
        (*dF)["ab"] +=      tmp1["ac"]*(*Smhalf)["cb"];
    }

    vector<DistTensor<double>*> Fab(2);
    Fab[0] = Fa;
    Fab[1] = Fb;
    diis.extrapolate(Fab, vector<DistTensor<double>*>(1, dF));
}

void UHF::calcEnergy()
{
    /*
     * E = (1/2)Tr[D(F+H)]
     *
     *   = (1/2)Tr[Da*(Fa+H) + Db*(Fb+H)]
     */
    (*Fa)["ab"] += (*H)["ab"];
    (*Fb)["ab"] += (*H)["ab"];
    energy = molecule.getNuclearRepulsion();
    energy += 0.5*scalar((*Da)["ab"]*(*Fa)["ab"]);
    energy += 0.5*scalar((*Db)["ab"]*(*Fb)["ab"]);
    (*Fa)["ab"] -= (*H)["ab"];
    (*Fb)["ab"] -= (*H)["ab"];
}

void UHF::calcDensity()
{
    /*
     * D[ab] = C[ai]*C[bi]
     */
    (*dDa)["ab"] = (*Da)["ab"];
    (*dDb)["ab"] = (*Db)["ab"];
    (*Da)["ab"] = (*Ca_occ)["ai"]*(*Ca_occ)["bi"];
    (*Db)["ab"] = (*Cb_occ)["ai"]*(*Cb_occ)["bi"];
    (*dDa)["ab"] -= (*Da)["ab"];
    (*dDb)["ab"] -= (*Db)["ab"];

    if (damping > 0.0)
    {
        (*Da)["ab"] += damping*(*dDa)["ab"];
        (*Db)["ab"] += damping*(*dDb)["ab"];
    }
}

double UHF::getS2() const
{
    int shapeNN[] = {NS,NS};
    int sizeab[] = {nalpha,nbeta};
    int sizean[] = {nalpha,norb};

    DistTensor<double> Delta(ctf, 2, sizeab, shapeNN, false);
    DistTensor<double> tmp(ctf, 2, sizean, shapeNN, false);

    double ndiff = abs(nalpha-nbeta);
    int nmin = min(nalpha, nbeta);

    double S2 = (ndiff/2)*(ndiff/2+1) + nmin;

    tmp["ai"] = (*Ca_occ)["ja"]*(*S)["ij"];
    Delta["ab"] = tmp["ai"]*(*Cb_occ)["ib"];

    S2 -= scalar(Delta*Delta);

    return fabs(S2);
}

double UHF::getAvgNumAlpha() const
{
    return scalar((*S)*(*Da));
}

double UHF::getAvgNumBeta() const
{
    return scalar((*S)*(*Db));
}

}
}
