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

#ifndef _AQUARIUS_SCF_CHOLESKYSCF_HPP_
#define _AQUARIUS_SCF_CHOLESKYSCF_HPP_

#include "scf.hpp"
#include "cholesky.hpp"

namespace aquarius
{
namespace scf
{

template <typename T> class CholeskyMOIntegrals;

template <typename T>
class CholeskyUHF : public UHF<T>
{
    friend class CholeskyMOIntegrals<T>;

    protected:
        const CholeskyIntegrals<T>& chol;
        tensor::DistTensor<T> *J;
        tensor::DistTensor<T> *JD;
        tensor::DistTensor<T> *La_occ, *Lb_occ;
        tensor::DistTensor<T> *LDa_occ, *LDb_occ;

    public:
        CholeskyUHF(const input::Config& config, const CholeskyIntegrals<T>& chol)
        : UHF<T>(chol.ctf, config, chol.getMolecule()), chol(chol)
        {
            int shapeN[] = {NS};
            int shapeNNN[] = {NS,NS,NS};

            int sizer[] = {chol.getRank()};
            J = new tensor::DistTensor<T>(this->ctf, 1, sizer, shapeN, false);
            JD = new tensor::DistTensor<T>(this->ctf, 1, sizer, shapeN, false);
            int sizenOr[] = {this->norb,this->nalpha,chol.getRank()};
            int sizenor[] = {this->norb,this->nbeta,chol.getRank()};
            La_occ = new tensor::DistTensor<T>(this->ctf, 3, sizenOr, shapeNNN, false);
            Lb_occ = new tensor::DistTensor<T>(this->ctf, 3, sizenor, shapeNNN, false);
            LDa_occ = new tensor::DistTensor<T>(this->ctf, 3, sizenOr, shapeNNN, false);
            LDb_occ = new tensor::DistTensor<T>(this->ctf, 3, sizenor, shapeNNN, false);
        }

        ~CholeskyUHF()
        {
            delete J;
            delete JD;
            delete La_occ;
            delete Lb_occ;
            delete LDa_occ;
            delete LDb_occ;
        }

    protected:
        void buildFock()
        {
            /*
             * Coulomb contribution:
             *
             * F[ab] = (Da[cd]+Db[cd])*(ab|cd)
             *
             *       = (Da[cd]+Db[cd])*L[abJ]*D[J]*L[cdJ]
             *
             *       = L(abJ)*{D[J]*{(Da[cd]+Db[cd])*L[cdJ]}}
             *
             *       = L(abJ)*{D[J]*J[J]}
             */
            (*this->Da) += (*this->Db);
            (*J)["J"] = chol.getL()["cdJ"]*(*this->Da)["cd"];
            (*JD)["J"] = chol.getD()["J"]*(*J)["J"];
            (*this->Da) -= (*this->Db);
            (*this->Fa)["ab"] = (*JD)["J"]*chol.getL()["abJ"];

            /*
             * Core contribution:
             *
             * F += H
             *
             * Up though this point, Fa = Fb
             */
            (*this->Fa) += (*this->H);
            (*this->Fb) = (*this->Fa);

            /*
             * Exchange contribution:
             *
             * Fa[ab] -= Da[cd]*(ac|bd)
             *
             *         = C[ci]*C[di]*L[acJ]*D[J]*L[bdJ]
             *
             *         = {C[ci]*L[acJ]}*{D[J]*{C[di]*L[bdJ]}}
             *
             *         = L[aiJ]*{D[J]*L[biJ]}
             */
            (*La_occ)["aiJ"] = chol.getL()["acJ"]*(*this->Ca_occ)["ci"];
            (*LDa_occ)["aiJ"] = chol.getD()["J"]*(*La_occ)["aiJ"];
            (*this->Fa)["ab"] -= (*LDa_occ)["aiJ"]*(*La_occ)["biJ"];

            (*Lb_occ)["aiJ"] = chol.getL()["acJ"]*(*this->Cb_occ)["ci"];
            (*LDb_occ)["aiJ"] = chol.getD()["J"]*(*Lb_occ)["aiJ"];
            (*this->Fb)["ab"] -= (*LDb_occ)["aiJ"]*(*Lb_occ)["biJ"];
        }
};

}
}

#endif
