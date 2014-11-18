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
    cout << "Start run" << endl;
    TwoElectronOperator<U>& W = get<TwoElectronOperator<U> >("H");
    const Space& occ = W.occ;
    const Space& vrt = W.vrt;
    // int mysize = occ.nalpha[0] * vrt.nalpha[0] + occ.nbeta[0] * vrt.nbeta[0];

    SpinorbitalTensor<U> Hguess("Hguess", W.getAIBJ());
    Hguess = 0;

    cout << "t1" << endl;

    SpinorbitalTensor<U>& FAB = W.getAB();
    SpinorbitalTensor<U>& FIJ = W.getIJ();
    SpinorbitalTensor<U>& WAIBJ = W.getAIBJ();

    cout << "t2" << endl;

    Hguess["aiaj"] -= FIJ["ij"];
    Hguess["aibi"] += FAB["ab"];
    Hguess["aibj"] -= WAIBJ["aibj"];

    int count = 0;
    int mysize = 0;

    cout << "t6" << endl;  

    const Molecule& molecule = get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();
    for (int R = 0;R < group.getNumIrreps();R++) {
        mysize = 0;
        cout << "R = " << R << endl;
        const Representation& irr_R = group.getIrrep(R);
        for (int spin_bj = 1;spin_bj >= 0;spin_bj--) {
            for (int j = 0;j < group.getNumIrreps();j++) {
                const Representation& irr_j = group.getIrrep(j);
                for (int b = 0;b < group.getNumIrreps();b++) {
                    const Representation& irr_b = group.getIrrep(b);
                    if (!(irr_b*irr_j*irr_R).isTotallySymmetric()) continue; // Loop over j,b only order N

                    if (spin_bj == 1) {
                        mysize += vrt.nalpha[b]*occ.nalpha[j];
                    }
                    else {
                        mysize += vrt.nbeta[b]*occ.nbeta[j];
                    }
                    for (int spin_ai = 1;spin_ai >= 0;spin_ai--) {
                        for (int i = 0;i < group.getNumIrreps();i++) {
                            const Representation& irr_i = group.getIrrep(i);
                            for (int a = 0;a < group.getNumIrreps();a++) {
                                const Representation& irr_a = group.getIrrep(a);
                                if (!(irr_a*irr_i*irr_R).isTotallySymmetric()) continue; // Loop over i,a only order N
                                CTFTensor<U>& this_tensor = Hguess(vec(spin_ai,spin_bj),vec(spin_bj,spin_ai))(vec(a,j,b,i));
                                CTFTensor<U> trans_this_tensor(this_tensor);
                                trans_this_tensor["ajbi"] = this_tensor["aibj"];
                                count += 1;
                            }
                        }
                    }
                }
            }
        }
        // Do diagonalization!
    }

    cout << "count = " << count << endl;


    put("TDAevals", new CTFTensor<U>("TDAevals", arena, 1, vec(mysize), vec(NS), true));
    put("TDAevecs", new CTFTensor<U>("TDAevecs", arena, 2, vec(mysize,mysize), vec(NS,NS), true));

    CTFTensor<U>& TDAevals = get<CTFTensor<U> >("TDAevals");
    CTFTensor<U>& TDAevecs = get<CTFTensor<U> >("TDAevecs");

    

    vector<U> dataAAAA;
    vector<U> dataAABB;
    vector<U> dataBBBB;
    vector<U> dataBBAA;

    CTFTensor<U>& AAAA = Hguess(vec(0,0),vec(0,0))(vec(0,0,0,0));
    CTFTensor<U>& AABB = Hguess(vec(1,0),vec(0,1))(vec(0,0,0,0));
    CTFTensor<U>& BBBB = Hguess(vec(1,1),vec(1,1))(vec(0,0,0,0));
    CTFTensor<U>& BBAA = Hguess(vec(0,1),vec(1,0))(vec(0,0,0,0));
    // CTFTensor<U> H(AAAA);
    // H["ajbi"]  = AAAA["aibj"];
    // H["ajbi"] += AABB["aibj"];

    CTFTensor<U> TAAAA(AAAA);
    CTFTensor<U> TAABB(AABB);
    CTFTensor<U> TBBBB(BBBB);
    CTFTensor<U> TBBAA(BBAA);
    TAAAA["ajbi"] = AAAA["aibj"];
    TAABB["ajbi"] = AABB["aibj"];
    TBBBB["ajbi"] = BBBB["aibj"];
    TBBAA["ajbi"] = BBAA["aibj"];
    TAAAA.getAllData(dataAAAA);
    TAABB.getAllData(dataAABB);
    TBBBB.getAllData(dataBBBB);
    TBBAA.getAllData(dataBBAA);

    int alphasize = sqrt(dataAAAA.size());
    int betasize = sqrt(dataBBBB.size());
    int rowlength = alphasize+betasize;
    // assert(rowlength == mysize);
    assert(dataAABB.size() == alphasize*betasize);
    assert(dataBBAA.size() == alphasize*betasize);
    vector<U> data(pow(rowlength,2));

    for (int i = 0; i < alphasize; ++i) {
        for (int j = 0; j < alphasize; ++j) {
            data[j+i*rowlength] = dataAAAA[j+i*alphasize];
        }
    }

    for (int i = 0; i < betasize; ++i) {
        for (int j = 0; j < betasize; ++j) {
            data[j+alphasize*rowlength+alphasize+i*rowlength] = dataBBBB[j+i*betasize];
        }
    } 

    for (int i = 0; i < alphasize; ++i) {
        for (int j = 0; j < betasize; ++j) {
            data[j+alphasize+i*rowlength] = dataAABB[j+i*betasize];
        }
    } 

    for (int i = 0; i < betasize; ++i) {
        for (int j = 0; j < alphasize; ++j) {
            data[j+alphasize*rowlength+i*rowlength] = dataBBAA[j+i*alphasize];
        }
    } 

    // cout << "TAAA, size = " << data.size() << endl;
    // for (int i = 0; i < data.size(); ++i) {
    //     cout << data[i] << endl;
    // }

    vector<U> w(mysize);
    heev('V','U',mysize,data.data(),mysize,w.data());

    // cout << " ok " << endl;

    // for (int i = 0; i < data.size(); ++i) {
    //     cout << data[i] << endl;
    // }

    int num_singlets = 0;
    int num_triplets = 0;
    vector<int> checkspin(mysize);
    double thisdot = 0.0;
    for (int i = 0; i < mysize; ++i) {
        thisdot = 0.0;
        for (int j = 0; j < alphasize; ++j) { // What if alphasize != betasize?
            thisdot += data[j+i*rowlength]*data[j+alphasize+i*rowlength];
        }
        if (thisdot > 0.0) {
            checkspin[i] = 1; // aka this is a singlet state
            num_singlets += 1;
        }
        else {
            assert(thisdot < 0.0);
            num_triplets += 1;
            checkspin[i] = 3; // aka this is a triplet state
        }
    }

    // std::cout << num_singlets << " singlet states" << std::endl;
    // std::cout << num_triplets << " triplet states" << std::endl;


    vector<int> evalorder = range(mysize);

    cosort(w.begin(), w.end(),
           evalorder.begin(), evalorder.end());
    printf("%18.15f\n", w[0]);
    printf("%18.15f\n", w[1]);

    vector<tkv_pair<U> > pairs;

    for (int i = 0;i < mysize;i++)
    {
        if (checkspin[evalorder[i]] == 1) {
            pairs.push_back(tkv_pair<U>(i, w[evalorder[i]]));
            // cout << setprecision(10) << w[evalorder[i]] << endl;
        }
    }

    if (arena.rank == 0)
        TDAevals.writeRemoteData(pairs);
    else
        TDAevals.writeRemoteData();

    pairs.clear();

    for (int i = 0;i < mysize;i++)
    {
        for (int j = 0;j < mysize;j++)
        {
            pairs.push_back(kv_pair(j+i*mysize,data[evalorder[i]*mysize + j]));
            // cout << j+i*mysize << " " << data[evalorder[i]*mysize + j] << endl;
        }
    }

    if (arena.rank == 0)
        TDAevecs.writeRemoteData(pairs);
    else
        TDAevecs.writeRemoteData();

    // vector<U> data1;
    // TDAevecs.getAllData(data1);
    // cout << "data1.size() = " << data1.size() << endl;

}

INSTANTIATE_SPECIALIZATIONS(TDA);
REGISTER_TASK(TDA<double>, "tda");
