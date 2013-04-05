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

#ifndef _AQUARIUS_SCF_HPP_
#define _AQUARIUS_SCF_HPP_

#include <vector>
#include <cmath>
#include <complex>

#include "mpi.h"

#ifdef USE_ELEMENTAL
#include "elemental.hpp"
#endif

#include "tensor/dist_tensor.hpp"
#include "slide/slide.hpp"
#include "input/molecule.hpp"
#include "input/config.hpp"
#include "util/util.h"
#include "util/blas.h"
#include "util/lapack.h"
#include "util/distributed.hpp"
#include "util/iterative.hpp"
#include "convergence/diis.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class UHF : public Iterative, public Distributed<T>
{
    protected:
        const input::Molecule& molecule;
        int norb;
        int nalpha;
        int nbeta;
        T damping;
        std::vector<T> Ea, Eb;
        tensor::DistTensor<T> *Fa, *Fb;
        tensor::DistTensor<T> *dF;
        tensor::DistTensor<T> *Ca_occ, *Cb_occ;
        tensor::DistTensor<T> *Ca_vrt, *Cb_vrt;
        tensor::DistTensor<T> *Da, *Db;
        tensor::DistTensor<T> *dDa, *dDb;
        tensor::DistTensor<T> *S, *Smhalf;
        tensor::DistTensor<T> *H;
        aquarius::convergence::DIIS< tensor::DistTensor<T> > diis;
        #ifdef USE_ELEMENTAL
        elem::Grid grid;
        elem::DistMatrix<T> C_elem;
        elem::DistMatrix<T> S_elem;
        elem::DistMatrix<T> F_elem;
        elem::DistMatrix<T,elem::VR,elem::STAR> E_elem;
        #endif

    public:
        UHF(tCTF_World<T>& ctf, const input::Config& config, const input::Molecule& molecule)
        : Iterative(config),
          Distributed<T>(ctf),
          molecule(molecule),
          norb(molecule.getNumOrbitals()),
          nalpha(molecule.getNumAlphaElectrons()),
          nbeta(molecule.getNumBetaElectrons()),
          damping(config.get<T>("damping")),
          Ea(norb),
          Eb(norb),
          diis(config.get("diis"), 2)
          #ifdef USE_ELEMENTAL
          ,grid(comm),
          C_elem(norb, norb, grid),
          S_elem(norb, norb, grid),
          F_elem(norb, norb, grid),
          E_elem(norb, 1   , grid)
          #endif
        {
            energy = molecule.getNuclearRepulsion();

            int shapeNN[] = {NS,NS};

            int sizenn[] = {norb,norb};
            Fa = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, false);
            Fb = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, false);
            dF = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, false);
            int sizenO[] = {norb,nalpha};
            int sizeno[] = {norb,nbeta};
            Ca_occ = new tensor::DistTensor<T>(ctf, 2, sizenO, shapeNN, false);
            Cb_occ = new tensor::DistTensor<T>(ctf, 2, sizeno, shapeNN, false);
            int sizenV[] = {norb,norb-nalpha};
            int sizenv[] = {norb,norb-nbeta};
            Ca_vrt = new tensor::DistTensor<T>(ctf, 2, sizenV, shapeNN, false);
            Cb_vrt = new tensor::DistTensor<T>(ctf, 2, sizenv, shapeNN, false);
            Da = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, true);
            Db = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, true);
            dDa = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, false);
            dDb = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, false);
            S = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, false);
            Smhalf = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, false);
            H = new tensor::DistTensor<T>(ctf, 2, sizenn, shapeNN, false);

            calcOverlap();
            calc1eHamiltonian();
        }

        ~UHF()
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

        void _iterate()
        {
            buildFock();
            DIISExtrap();
            calcEnergy();
            diagonalizeFock();
            calcDensity();

            switch (convtype)
            {
                case MAX_ABS:
                    conv = std::max(dDa->reduce(CTF_OP_MAXABS), dDb->reduce(CTF_OP_MAXABS));
                    break;
                case RMSD:
                    conv = sqrt((dDa->reduce(CTF_OP_SQNRM2)+dDb->reduce(CTF_OP_SQNRM2))/(2*norb*norb));
                    break;
                case MAD:
                    conv = (dDa->reduce(CTF_OP_SUMABS)+dDb->reduce(CTF_OP_SUMABS))/(2*norb*norb);
                    break;
            }
        }

        T getMultiplicity() const
        {
            return sqrt(4*getS2()+1);
        }

        T getS2() const
        {
            int shapeNN[] = {NS,NS};
            int sizeab[] = {nalpha,nbeta};
            int sizean[] = {nalpha,norb};

            tensor::DistTensor<T> Delta(this->ctf, 2, sizeab, shapeNN, false);
            tensor::DistTensor<T> tmp(this->ctf, 2, sizean, shapeNN, false);

            T ndiff = abs(nalpha-nbeta);
            int nmin = std::min(nalpha, nbeta);

            T S2 = (ndiff/2)*(ndiff/2+1) + nmin;

            tmp["ai"] = (*Ca_occ)["ja"]*(*S)["ij"];
            Delta["ab"] = tmp["ai"]*(*Cb_occ)["ib"];

            S2 -= scalar(Delta*Delta);

            return std::abs(S2);
        }

        T getAvgNumAlpha() const
        {
            return scalar((*S)*(*Da));
        }

        T getAvgNumBeta() const
        {
            return scalar((*S)*(*Db));
        }

        const T* getAlphaEigenvalues() const { return Ea.data(); }

        const T* getBetaEigenvalues() const { return Eb.data(); }

        const input::Molecule& getMolecule() const { return molecule; }

        const tensor::DistTensor<T>& getOverlap() const { return *S; }
        const tensor::DistTensor<T>& get1eHamiltonian() const { return *H; }

        const tensor::DistTensor<T>& getCA() const { return *Ca_vrt; }
        const tensor::DistTensor<T>& getCI() const { return *Ca_occ; }
        const tensor::DistTensor<T>& getCa() const { return *Cb_vrt; }
        const tensor::DistTensor<T>& getCi() const { return *Cb_occ; }

    protected:
        void calcOverlap()
        {
            slide::Context context;

            std::vector< tkv_pair<T> > pairs;

            int pid = this->comm.Get_rank();
            int nproc = this->comm.Get_size();
            int block = 0;
            for (input::Molecule::const_shell_iterator i = molecule.getShellsBegin();i != molecule.getShellsEnd();++i)
            {
                for (input::Molecule::const_shell_iterator j = molecule.getShellsBegin();j != molecule.getShellsEnd();++j)
                {
                    if (i < j) continue;

                    if (block%nproc == pid)
                    {
                        context.calcOVI(1.0, 0.0, *i, *j);

                        size_t nint = context.getNumIntegrals();
                        std::vector<T> ints(nint);
                        std::vector<idx2_t> idxs(nint);

                        size_t nproc = context.process1eInts(nint, ints.data(), idxs.data(), -1.0);
                        for (int i = 0;i < nproc;i++)
                        {
                            //printf("%d %d %25.15f\n", idxs[i].i+1, idxs[i].j+1, ints[i]);

                            pairs.push_back(tkv_pair<T>(idxs[i].i*norb+idxs[i].j, ints[i]));
                            if (idxs[i].i != idxs[i].j)
                            {
                                pairs.push_back(tkv_pair<T>(idxs[i].j*norb+idxs[i].i, ints[i]));
                            }
                        }
                    }

                    block++;
                }
            }

            S->writeRemoteData(pairs.size(), pairs.data());

            pairs.clear();

            #ifdef USE_ELEMENTAL

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

                    pairs.push_back(tkv_pair<T>(r*norb+c,0));
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

            DistMatrix<T> Smhalf_elem(S_elem);
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

            #else

            int64_t size;
            T *s;
            std::vector<T> work(3*norb);
            std::vector<T> smhalf(norb*norb);

            S->getAllData(size, s);
            assert(size == norb*norb);

            int info = dsyev('V', 'U', norb, s, norb, Ea.data(), work.data(), 3*norb);
            assert(info == 0);

            std::fill(smhalf.begin(), smhalf.end(), 0.0);
            for (int i = 0;i < norb;i++)
            {
                dger(norb, norb, 1/sqrt(Ea[i]), s+i*norb, 1, s+i*norb, 1, smhalf.data(), norb);
            }

            free(s);

            if (this->comm.Get_rank() == 0)
            {
                std::vector< tkv_pair<T> > pairs;

                for (int i = 0;i < norb;i++)
                {
                    for (int j = 0;j < norb;j++)
                    {
                        pairs.push_back(tkv_pair<T>(i+j*norb, smhalf[i+j*norb]));
                    }
                }

                Smhalf->writeRemoteData(norb*norb, pairs.data());
            }
            else
            {
                Smhalf->writeRemoteData(0, NULL);
            }

            #endif
        }

        void calc1eHamiltonian()
        {
            slide::Context context;

            std::vector<slide::Center> centers;
            for (std::vector<input::Atom>::const_iterator a = molecule.getAtomsBegin();a != molecule.getAtomsEnd();++a)
            {
                centers.push_back(a->getCenter());
            }

            std::vector< tkv_pair<T> > pairs;

            int pid = this->comm.Get_rank();
            int nproc = this->comm.Get_size();
            int block = 0;
            for (input::Molecule::const_shell_iterator i = molecule.getShellsBegin();i != molecule.getShellsEnd();++i)
            {
                for (input::Molecule::const_shell_iterator j = molecule.getShellsBegin();j != molecule.getShellsEnd();++j)
                {
                    if (i < j) continue;

                    if (block%nproc == pid)
                    {
                        context.calcKEI(1.0, 0.0, *i, *j);
                        context.calcNAI(1.0, 1.0, *i, *j, centers.data(), centers.size());

                        size_t nint = context.getNumIntegrals();
                        std::vector<T> ints(nint);
                        std::vector<idx2_t> idxs(nint);

                        size_t nproc = context.process1eInts(nint, ints.data(), idxs.data(), -1.0);
                        for (int i = 0;i < nproc;i++)
                        {
                            //printf("%d %d %25.15f\n", idxs[i].i+1, idxs[i].j+1, ints[i]);

                            pairs.push_back(tkv_pair<T>(idxs[i].i*norb+idxs[i].j, ints[i]));
                            if (idxs[i].i != idxs[i].j)
                            {
                                pairs.push_back(tkv_pair<T>(idxs[i].j*norb+idxs[i].i, ints[i]));
                            }
                        }
                    }

                    block++;
                }
            }

            H->writeRemoteData(pairs.size(), pairs.data());
        }

        void diagonalizeFock()
        {
            std::vector< tkv_pair<T> > pairs;

            #ifdef USE_ELEMENTAL

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
                DistMatrix<T> S_tmp(S_elem);
                HermitianGenDefiniteEig(AXBX, LOWER, F_elem, S_tmp, E_elem, C_elem);
                SortEig(E_elem, C_elem);

                DistMatrix<T,STAR,STAR> E_local(E_elem);
                for (int i = 0;i < norb;i++) Ea[i] = E_local.GetLocal(i,0);

                std::vector< tkv_pair<T> > pairs_occ;
                std::vector< tkv_pair<T> > pairs_vrt;

                for (int i = 0;i < C_elem.LocalHeight();i++)
                {
                    for (int j = 0;j < C_elem.LocalWidth();j++)
                    {
                        int c = cshift+i*cstride;
                        int r = rshift+j*rstride;

                        if (r < nalpha)
                        {
                            pairs_occ.push_back(tkv_pair<T>(r*norb+c, C_elem.GetLocal(i,j)));
                        }
                        else
                        {
                            pairs_vrt.push_back(tkv_pair<T>((r-nalpha)*norb+c, C_elem.GetLocal(i,j)));
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
                DistMatrix<T> S_tmp(S_elem);
                HermitianGenDefiniteEig(AXBX, LOWER, F_elem, S_tmp, E_elem, C_elem);
                SortEig(E_elem, C_elem);

                DistMatrix<T,STAR,STAR> E_local(E_elem);
                for (int i = 0;i < norb;i++) Eb[i] = E_local.GetLocal(i,0);

                std::vector< tkv_pair<T> > pairs_occ;
                std::vector< tkv_pair<T> > pairs_vrt;

                for (int i = 0;i < C_elem.LocalHeight();i++)
                {
                    for (int j = 0;j < C_elem.LocalWidth();j++)
                    {
                        int c = cshift+i*cstride;
                        int r = rshift+j*rstride;

                        if (r < nbeta)
                        {
                            pairs_occ.push_back(tkv_pair<T>(r*norb+c, C_elem.GetLocal(i,j)));
                        }
                        else
                        {
                            pairs_vrt.push_back(tkv_pair<T>((r-nbeta)*norb+c, C_elem.GetLocal(i,j)));
                        }
                    }
                }

                Cb_occ->writeRemoteData(pairs_occ.size(), pairs_occ.data());
                Cb_vrt->writeRemoteData(pairs_vrt.size(), pairs_vrt.data());
            }

            #else

            int64_t size;
            int info;
            T *fock, *s;
            std::vector<T> work(3*norb);

            S->getAllData(size, s);
            assert(size == norb*norb);
            Fa->getAllData(size, fock);
            assert(size == norb*norb);
            info = dsygv(1, 'V', 'U', norb, fock, norb, s, norb, Ea.data(), work.data(), 3*norb);
            assert(info == 0);

            if (this->comm.Get_rank() == 0)
            {
                for (int i = 0;i < nalpha;i++)
                {
                    for (int p = 0;p < norb;p++)
                    {
                        pairs.push_back(tkv_pair<T>(p+i*norb, fock[p+i*norb]));
                    }
                }

                Ca_occ->writeRemoteData(norb*nalpha, pairs.data());
                pairs.clear();

                for (int i = 0;i < norb-nalpha;i++)
                {
                    for (int p = 0;p < norb;p++)
                    {
                        pairs.push_back(tkv_pair<T>(p+i*norb, fock[p+(i+nalpha)*norb]));
                    }
                }

                Ca_vrt->writeRemoteData(norb*(norb-nalpha), pairs.data());
                pairs.clear();
            }
            else
            {
                Ca_occ->writeRemoteData(0, NULL);
                Ca_vrt->writeRemoteData(0, NULL);
            }

            free(fock);
            free(s);

            S->getAllData(size, s);
            assert(size == norb*norb);
            Fb->getAllData(size, fock);
            assert(size == norb*norb);
            info = dsygv(1, 'V', 'U', norb, fock, norb, s, norb, Eb.data(), work.data(), 3*norb);
            assert(info == 0);

            if (this->comm.Get_rank() == 0)
            {
                for (int i = 0;i < nbeta;i++)
                {
                    for (int p = 0;p < norb;p++)
                    {
                        pairs.push_back(tkv_pair<T>(p+i*norb, fock[p+i*norb]));
                    }
                }

                Cb_occ->writeRemoteData(norb*nbeta, pairs.data());
                pairs.clear();

                for (int i = 0;i < norb-nbeta;i++)
                {
                    for (int p = 0;p < norb;p++)
                    {
                        pairs.push_back(tkv_pair<T>(p+i*norb, fock[p+(i+nbeta)*norb]));
                    }
                }

                Cb_vrt->writeRemoteData(norb*(norb-nbeta), pairs.data());
            }
            else
            {
                Cb_occ->writeRemoteData(0, NULL);
                Cb_vrt->writeRemoteData(0, NULL);
            }

            free(fock);
            free(s);

            #endif

            fixPhase(*Ca_occ);
            fixPhase(*Cb_occ);
            fixPhase(*Ca_vrt);
            fixPhase(*Cb_vrt);

            /*
            cout << "Eigenvalues" << endl;
            for (int i = 0;i < norb;i++)
                cout << i << ' ' << Ea[i] << ' ' << Eb[i] << endl;
            cout << endl;

            DistTensor<T> test1(*Ca_occ);
            test1 -= *Cb_occ;
            test1.print(stdout);

            DistTensor<T> test2(*Ca_vrt);
            test2 -= *Cb_vrt;
            test2.print(stdout);

            exit(1);
            */
        }

        void fixPhase(tensor::DistTensor<T>& C)
        {
            int rank = this->comm.Get_rank();
            int np = this->comm.Get_size();
            int nr = C.getLengths()[1];

            std::vector< tkv_pair<T> > pairs(norb, tkv_pair<T>(0,0));

            for (int b = 0;b*np < nr;b++)
            {
                int r = b*np+rank;

                if (r < nr)
                {
                    for (int i = 0;i < norb;i++) pairs[i].k = i+r*norb;

                    C.getRemoteData(norb, pairs.data());

                    std::sort(pairs.begin(), pairs.end());
                    int sign = 0;
                    for (int i = 0;i < norb;i++)
                    {
                        if (sign == 0 && std::abs(pairs[i].d) > 1e-12)
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

        virtual void buildFock() = 0;

        void calcEnergy()
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

        void calcDensity()
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

        void DIISExtrap()
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
                tensor::DistTensor<T> tmp1(*Fa);
                tensor::DistTensor<T> tmp2(*Fa);

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

            std::vector< tensor::DistTensor<T>* > Fab(2);
            Fab[0] = Fa;
            Fab[1] = Fb;
            diis.extrapolate(Fab, std::vector< tensor::DistTensor<T>* >(1, dF));
        }
};

}
}

#endif
