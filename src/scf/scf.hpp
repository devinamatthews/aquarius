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
#include "diis/diis.hpp"

namespace aquarius
{
namespace scf
{

class UHF : public Iterative, public Distributed<double>
{
    protected:
        const input::Molecule& molecule;
        int norb;
        int nalpha;
        int nbeta;
        double damping;
        std::vector<double> Ea, Eb;
        diis::DIIS< tensor::DistTensor<double> > diis;
        tensor::DistTensor<double> *Fa, *Fb;
        tensor::DistTensor<double> *dF;
        tensor::DistTensor<double> *Ca_occ, *Cb_occ;
        tensor::DistTensor<double> *Ca_vrt, *Cb_vrt;
        tensor::DistTensor<double> *Da, *Db;
        tensor::DistTensor<double> *dDa, *dDb;
        tensor::DistTensor<double> *S, *Smhalf;
        tensor::DistTensor<double> *H;
        #ifdef USE_ELEMENTAL
        elem::Grid grid;
        elem::DistMatrix<double> C_elem;
        elem::DistMatrix<double> S_elem;
        elem::DistMatrix<double> F_elem;
        elem::DistMatrix<double,elem::VR,elem::STAR> E_elem;
        #endif

        void getOverlap();

        void get1eHamiltonian();

        void diagonalizeFock();

        void fixPhase(tensor::DistTensor<double>& C);

        virtual void buildFock() = 0;

        void calcEnergy();

        void calcDensity();

        void DIISExtrap();

    public:
        UHF(tCTF_World<double>& ctf, const input::Molecule& molecule, const input::Config& config);

        ~UHF();

        void _iterate();

        double getMultiplicity() const
        {
            return sqrt(4*getS2()+1);
        }

        double getS2() const;

        double getAvgNumAlpha() const;

        double getAvgNumBeta() const;

        const double* getAlphaEigenvalues() const { return Ea.data(); }

        const double* getBetaEigenvalues() const { return Eb.data(); }

        const input::Molecule& getMolecule() const { return molecule; }

        const tensor::DistTensor<double>& getCA() const { return *Ca_vrt; }
        const tensor::DistTensor<double>& getCI() const { return *Ca_occ; }
        const tensor::DistTensor<double>& getCa() const { return *Cb_vrt; }
        const tensor::DistTensor<double>& getCi() const { return *Cb_occ; }
};

}
}

#endif
