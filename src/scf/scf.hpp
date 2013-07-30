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
        UHF(tCTF_World<T>& ctf, const input::Config& config, const input::Molecule& molecule);

        ~UHF();

        void _iterate();

        T getMultiplicity() const;

        T getS2() const;

        T getAvgNumAlpha() const;

        T getAvgNumBeta() const;

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
        void calcOverlap();

        void calc1eHamiltonian();

        void diagonalizeFock();

        void fixPhase(tensor::DistTensor<T>& C);

        virtual void buildFock() = 0;

        void calcEnergy();

        void calcDensity();

        void DIISExtrap();
};

}
}

#endif
