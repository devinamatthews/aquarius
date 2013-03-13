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

#ifndef _AQUARIUS_SCF_AOSCF_HPP_
#define _AQUARIUS_SCF_AOSCF_HPP_

#include "scf.hpp"
#include "cholesky.hpp"
#include "aoints.hpp"

namespace aquarius
{
namespace scf
{

template <typename T>
class AOUHF : public UHF<T>
{
    friend class AOMOIntegrals<T>;

    protected:
        const AOIntegrals<T>& ints;

    public:
        AOUHF(tCTF_World<T>& ctf, const input::Config& config, const AOIntegrals& ints)
        : UHF<T>(config, ctf, ints.getMolecule()), ints(ints) {}

    protected:
        void buildFock()
        {
            T *focka, *fockb;
            T *densa, *densb;

            int64_t npair;

            if (this->rank == 0)
            {
                H->getAllData(npair, focka);
                assert(npair == norb*norb);
                H->getAllData(npair, fockb);
                assert(npair == norb*norb);
            }
            else
            {
                focka = SAFE_MALLOC(T, norb*norb);
                fockb = SAFE_MALLOC(T, norb*norb);
                std::fill(focka, focka+norb*norb, 0.0);
                std::fill(fockb, fockb+norb*norb, 0.0);
            }

            Da->getAllData(npair, densa);
            assert(npair == norb*norb);
            Db->getAllData(npair, densb);
            assert(npair == norb*norb);

            size_t neris = ints.getNumInts();
            const integral_t* eris = ints.getInts();

            for (size_t n = 0;n < neris;n++)
            {
                int i = eris[n].idx.i;
                int j = eris[n].idx.j;
                int k = eris[n].idx.k;
                int l = eris[n].idx.l;

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

                T e = eris[n].value;
                if (i == j) e *= 0.5;
                if (k == l) e *= 0.5;
                if (min(i,j) == min(k,l) &&
                    max(i,j) == max(k,l)) e *= 0.5;

                /*
                 * Coulomb contribution: Fa(ab) += [Da(cd)+Db(cd)]*(ab|cd)
                 */

                // (ij|kl)
                focka[i+j*norb] += 2.0*(densa[k+l*norb]+densb[k+l*norb])*e;
                fockb[i+j*norb] += 2.0*(densa[k+l*norb]+densb[k+l*norb])*e;
                // (ji|kl)
                focka[j+i*norb] += 2.0*(densa[k+l*norb]+densb[k+l*norb])*e;
                fockb[j+i*norb] += 2.0*(densa[k+l*norb]+densb[k+l*norb])*e;
                // (kl|ij)
                focka[k+l*norb] += 2.0*(densa[i+j*norb]+densb[i+j*norb])*e;
                fockb[k+l*norb] += 2.0*(densa[i+j*norb]+densb[i+j*norb])*e;
                // (lk|ij)
                focka[l+k*norb] += 2.0*(densa[i+j*norb]+densb[i+j*norb])*e;
                fockb[l+k*norb] += 2.0*(densa[i+j*norb]+densb[i+j*norb])*e;

                /*
                 * Exchange contribution: Fa(ac) -= Da(bd)*(ab|cd)
                 */

                // (ij|kl)
                focka[i+k*norb] -= densa[j+l*norb]*e;
                fockb[i+k*norb] -= densb[j+l*norb]*e;
                // (kl|ij)
                focka[k+i*norb] -= densa[j+l*norb]*e;
                fockb[k+i*norb] -= densb[j+l*norb]*e;
                // (ij|lk)
                focka[i+l*norb] -= densa[j+k*norb]*e;
                fockb[i+l*norb] -= densb[j+k*norb]*e;
                // (lk|ij)
                focka[l+i*norb] -= densa[j+k*norb]*e;
                fockb[l+i*norb] -= densb[j+k*norb]*e;
                // (ji|lk)
                focka[j+l*norb] -= densa[i+k*norb]*e;
                fockb[j+l*norb] -= densb[i+k*norb]*e;
                // (lk|ji)
                focka[l+j*norb] -= densa[i+k*norb]*e;
                fockb[l+j*norb] -= densb[i+k*norb]*e;
                // (ji|kl)
                focka[j+k*norb] -= densa[i+l*norb]*e;
                fockb[j+k*norb] -= densb[i+l*norb]*e;
                // (kl|ji)
                focka[k+j*norb] -= densa[i+l*norb]*e;
                fockb[k+j*norb] -= densb[i+l*norb]*e;
            }

            if (this->rank == 0)
            {
                this->comm.Reduce(MPI::IN_PLACE, focka, norb*norb, this->type, MPI::SUM, 0);
                this->comm.Reduce(MPI::IN_PLACE, fockb, norb*norb, this->type, MPI::SUM, 0);

                tkv_pair<T>* pairs = SAFE_MALLOC(tkv_pair<T>, norb*norb);

                for (int p = 0;p < norb*norb;p++)
                {
                    pairs[p].d = focka[p];
                    pairs[p].k = p;
                }

                Fa->writeRemoteData(norb*norb, pairs);

                for (int p = 0;p < norb*norb;p++)
                {
                    pairs[p].d = fockb[p];
                    pairs[p].k = p;
                }

                Fb->writeRemoteData(norb*norb, pairs);

                FREE(pairs);
            }
            else
            {
                this->comm.Reduce(focka, NULL, norb*norb, this->type, MPI::SUM, 0);
                this->comm.Reduce(fockb, NULL, norb*norb, this->type, MPI::SUM, 0);

                Fa->writeRemoteData(0, NULL);
                Fb->writeRemoteData(0, NULL);
            }

            if (this->rank == 0)
            {
                free(focka);
                free(fockb);
            }
            else
            {
                FREE(focka);
                FREE(fockb);
            }

            free(densa);
            free(densb);
        }
};

}
}

#endif
