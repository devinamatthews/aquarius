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

#include "tda.hpp"
#include "util/lapack.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::symmetry;

template <typename U>
TDA<U>::TDA(const std::string& name, const Config& config)
: Task("tda", name)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("molecule", "molecule"));
    reqs.push_back(Requirement("moints", "H"));
    addProduct(Product("tda.TDAevals", "TDAevals", reqs));
    addProduct(Product("tda.TDAevecs", "TDAevecs", reqs));
}

template <typename U>
void TDA<U>::run(TaskDAG& dag, const Arena& arena)
{
    TwoElectronOperator<U>& W = get<TwoElectronOperator<U> >("H");
    const Space& occ = W.occ;
    const Space& vrt = W.vrt;

    SpinorbitalTensor<U> Hguess("Hguess", W.getAIBJ());
    Hguess = 0;

    SpinorbitalTensor<U>& FAB = W.getAB();
    SpinorbitalTensor<U>& FIJ = W.getIJ();
    SpinorbitalTensor<U>& WAIBJ = W.getAIBJ();

    Hguess["aiaj"] -= FIJ["ij"];
    Hguess["aibi"] += FAB["ab"];
    Hguess["aibj"] -= WAIBJ["aibj"];

    const Molecule& molecule = get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();
    int nirrep = group.getNumIrreps();

    put("TDAevecs", new vector<vector<SpinorbitalTensor<U> > >);
    put("TDAevals", new vector<vector<U> >);
    vector<vector<SpinorbitalTensor<U> > >& TDAevecs = get<vector<vector<SpinorbitalTensor<U> > > >("TDAevecs");
    vector<vector<U> >& TDAevals = get<vector<vector<U> > >("TDAevals");


    vector<int> alphasize(nirrep);
    vector<int> betasize(nirrep);
    vector<int> tempvec1(2);
    vector<vector<int> > subindex(nirrep,tempvec1);
    int count;
    int alphatot;
    int betatot;
    for (int R = 0;R < nirrep;R++) {
        const Representation& irr_R = group.getIrrep(R);
        count = 0;
        alphatot = 0;
        betatot = 0;

        for (int j = 0;j < nirrep;j++) {
            const Representation& irr_j = group.getIrrep(j);
            for (int b = 0;b < nirrep;b++) {
                const Representation& irr_b = group.getIrrep(b);
                if (!(irr_b*irr_j*irr_R).isTotallySymmetric()) continue;
                subindex[count][0] = b;
                subindex[count][1] = j;
                alphasize[count] = vrt.nalpha[b]*occ.nalpha[j];
                betasize[count] = vrt.nbeta[b]*occ.nbeta[j];
                alphatot += vrt.nalpha[b]*occ.nalpha[j];
                betatot += vrt.nbeta[b]*occ.nbeta[j];
                count += 1;
                assert(count <= nirrep);
            }
        }

        assert(count == nirrep);
        int SL = alphatot+betatot;
        vector<U> data(SL*SL);
        int Yoffset = 0;
        int Xoffset = 0;

        for (int Y = 0; Y < nirrep; Y++) {
            int b = subindex[Y][0];
            int j = subindex[Y][1];
            if (Y != 0)
                Yoffset += (alphasize[Y-1] + betasize[Y-1]);
            for (int X = 0; X < nirrep; X++) {
                int a = subindex[X][0];
                int i = subindex[X][1];
                if (X == 0)
                    Xoffset = 0;
                else
                    Xoffset += (alphasize[X-1] + betasize[X-1]);
                for (int spin_bj = 1;spin_bj >= 0;spin_bj--) {
                    for (int spin_ai = 1;spin_ai >= 0;spin_ai--) {
                        CTFTensor<U>& this_tensor = Hguess(vec(spin_ai,spin_bj),vec(spin_bj,spin_ai))(vec(a,j,b,i));
                        // CTFTensor<U> trans_tensor(Hguess(vec(spin_ai,spin_bj),vec(spin_bj,spin_ai))(vec(a,i,b,j))); // Fix for UHF, use norm ctor using vrt, occ
                        CTFTensor<U> trans_tensor("trans_tensor", arena, 4, vec((spin_ai == 1 ? vrt.nalpha[a] : vrt.nbeta[a]),
                                                                                (spin_ai == 1 ? occ.nalpha[i] : occ.nbeta[i]),
                                                                                (spin_bj == 1 ? vrt.nalpha[b] : vrt.nbeta[b]),
                                                                                (spin_bj == 1 ? occ.nalpha[j] : occ.nbeta[j])), vec(NS,NS,NS,NS), true);
                        trans_tensor["ajbi"] = this_tensor["aibj"];
                        vector<U> tempdata;
                        trans_tensor.getAllData(tempdata);
                        int ind1 = (1-spin_bj)*alphasize[Y] + spin_bj*betasize[Y];
                        int ind2 = (1-spin_ai)*alphasize[X] + spin_ai*betasize[X];
                        assert(tempdata.size() == ind1*ind2);
                        for (int k = 0; k < ind1; k++) {
                            for (int l = 0; l < ind2; l++) {
                                data[l+k*SL+Xoffset+Yoffset*SL + spin_ai*alphasize[X] + spin_bj*alphasize[Y]*SL] = tempdata[l+k*ind2];
                            }
                        }
                    }
                }
            }
        }
        vector<U> w(SL);
        heev('V','U',SL,data.data(),SL,w.data());
        TDAevals.push_back(w);
        for (int i = 0; i < w.size(); i++)
            cout << setprecision(10)<< w[i] << endl;
        /* First attempt at evec storage:
        SpinorbitalTensor<U> tempSoT(W.getAI());
        tempSoT = 0;
        vector<SpinorbitalTensor<U> > tempvec2(SL,tempSoT);
        TDAevecs.push_back(tempvec2);
        for (int i = 0; i < SL; i++) {
            int offset = 0;
            int alphadone = 0;
            int betadone = 0;
            for (int j = 0; j < nirrep; j++) {
                for (int k = 0; k < alphasize[j]; k++) {
                    TDAevecs[R][i](vec(0,0),vec(0,0))(vec(subindex[j][0],subindex[j][1])).writeRemoteData(vec(kv_pair(alphadone+k,data[i*SL+offset+k])));
                }
                offset += alphasize[j];
                alphadone += alphasize[j];
                for (int l = 0; l < betasize[j]; l++) {
                    TDAevecs[R][i](vec(1,0),vec(0,1))(vec(subindex[j][0],subindex[j][1])).writeRemoteData(vec(kv_pair(betadone+l,data[i*SL+offset+l])));
                }
                offset += betasize[j];
                betadone += betasize[j];
            }
        }
        */
    }
}

INSTANTIATE_SPECIALIZATIONS(TDA);
REGISTER_TASK(TDA<double>, "tda");
