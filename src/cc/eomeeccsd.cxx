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

#include "eomeeccsd.hpp"
#include "util/lapack.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::cc;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::convergence;
using namespace aquarius::symmetry;

template <typename U>
EOMEECCSD<U>::EOMEECCSD(const std::string& name, const Config& config)
: Iterative<U>("eomeeccsd", name, config), davidson_config(config.get("davidson").clone()),
  nroot(config.get<int>("nroot")), multiroot(config.get<bool>("multiroot"))
{
    vector<Requirement> reqs;
    reqs.emplace_back("ccsd.T", "T");
    reqs.emplace_back("ccsd.Hbar", "Hbar");
    reqs.emplace_back("tda.TDAevals", "TDAevals");
    reqs.emplace_back("tda.TDAevecs", "TDAevecs");
    this->addProduct("eomeeccsd.energy", "energy", reqs);
    this->addProduct("eomeeccsd.convergence", "convergence", reqs);
    this->addProduct("eomeeccsd.R", "R", reqs);
}

template <typename U>
void EOMEECCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    const PointGroup& group = H.getABIJ().getGroup();
    int nirrep = group.getNumIrreps();

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    this->puttmp("D", new Denominator<U>(H));
    auto& TDAevecs = this->template get<vector<unique_vector<SpinorbitalTensor<U>>>>("TDAevecs");
    auto& TDAevals = this->template get<vector<vector<U>>>("TDAevals");

    this->puttmp("XMI", new SpinorbitalTensor<U>("X(mi)", H.getIJ()));
    this->puttmp("XAE", new SpinorbitalTensor<U>("X(ae)", H.getAB()));

    vector<tuple<int,U,int,int> > tda_sorted;

    int nsinglet = 0;
    for (int i = 0;i < nirrep;i++)
    {
        vector<int> spin(TDAevals[i].size());

        for (int j = 0;j < TDAevals[i].size();j++)
        {
            if (scalar(TDAevecs[i][j]({1,0},{0,1})*TDAevecs[i][j]({0,0},{0,0})) < 0)
                spin[j] = 1;
            else
                nsinglet++;
        }

        tda_sorted += zip(spin,
                          TDAevals[i],
                          vector<int>(TDAevals[i].size(), i),
                          range<int>(TDAevals[i].size()));
    }

    // for (int i = 0; i < tda_sorted.size(); i++) {
    //     cout << "tda_sorted[i] = " << tda_sorted[i].first << " " << tda_sorted[i].second.first << " " << tda_sorted[i].second.second << endl;
    // }

    sort(tda_sorted);

    // for (int i = 0; i < tda_sorted.size(); i++) {
    //     cout << "tda_sorted[i] = " << tda_sorted[i].first << " " << tda_sorted[i].second.first << " " << tda_sorted[i].second.second << endl;
    // }

    nsinglet = min(nsinglet, nroot);
    vector<vector<int>> root_idx(nirrep);

    for (int i = 0;i < nsinglet;i++)
    {
        root_idx[get<2>(tda_sorted[i])].push_back(get<3>(tda_sorted[i]));
    }

    auto& Rs = this->puttmp("R", new unique_vector<ExcitationOperator<U,2>>());
    auto& Vs = this->puttmp("V", new unique_vector<ExcitationOperator<U,2>>());
    auto& Zs = this->puttmp("Z", new unique_vector<ExcitationOperator<U,2>>());

    for (int i = 0;i < nirrep;i++)
    {
        Vs.clear();
        Rs.clear();
        Zs.clear();

        if (multiroot)
        {
            for (int j = 0;j < root_idx[i].size();j++)
            {
                Rs.emplace_back("R", arena, occ, vrt, group.getIrrep(i));
                Zs.emplace_back("Z", arena, occ, vrt, group.getIrrep(i));
                ExcitationOperator<U,2>& R = Rs.back();
                R(0) = 0;
                R(1) = TDAevecs[i][root_idx[i][j]];
                R(2) = 0;
            }

            auto& davidson = this->puttmp("Davidson",
                new Davidson<ExcitationOperator<U,2> >(davidson_config, (int)root_idx[i].size()));

            Iterative<U>::run(dag, arena, root_idx[i].size());

            for (int j = 0;j < root_idx[i].size();j++)
            {
                if (this->isConverged(j))
                {
                    Vs.emplace_back("V", arena, occ, vrt, group.getIrrep(i));
                    ExcitationOperator<U,2>& V = Vs.back();
                    davidson.getSolution(j, V);
                    V /= sqrt(std::abs(scalar(conj(V)*V)));
                }
            }
        }
        else
        {
            for (int j = 0;j < root_idx[i].size();j++)
            {
                Rs.clear();
                Zs.clear();

                Rs.emplace_back("R", arena, occ, vrt, group.getIrrep(i));
                Zs.emplace_back("Z", arena, occ, vrt, group.getIrrep(i));
                ExcitationOperator<U,2>& R = Rs.back();
                R(0) = 0;
                R(1) = TDAevecs[i][root_idx[i][j]];
                R(2) = 0;

                auto& davidson = this->puttmp("Davidson",
                    new Davidson<ExcitationOperator<U,2> >(davidson_config));

                Iterative<U>::run(dag, arena);

                if (this->isConverged())
                {
                    Vs.emplace_back("V", arena, occ, vrt, group.getIrrep(i));
                    ExcitationOperator<U,2>& V = Vs.back();
                    davidson.getSolution(0, V);
                    V /= sqrt(std::abs(scalar(conj(V)*V)));
                }
            }
        }
    }

    this->put("energy", new CTFTensor<U>("energy", arena, 1, {nroot}, {NS}, true));
    this->put("convergence", new U(this->conv()));
}

template <typename U>
void EOMEECCSD<U>::iterate(const Arena& arena)
{
    auto& H = this->template get<STTwoElectronOperator<U> >("Hbar");

    const SpinorbitalTensor<U>&   FME =   H.getIA();
    const SpinorbitalTensor<U>&   FAE =   H.getAB();
    const SpinorbitalTensor<U>&   FMI =   H.getIJ();
    const SpinorbitalTensor<U>& WMNEF = H.getIJAB();
    const SpinorbitalTensor<U>& WAMEF = H.getAIBC();
    const SpinorbitalTensor<U>& WABEJ = H.getABCI();
    const SpinorbitalTensor<U>& WABEF = H.getABCD();
    const SpinorbitalTensor<U>& WMNIJ = H.getIJKL();
    const SpinorbitalTensor<U>& WMNEJ = H.getIJAK();
    const SpinorbitalTensor<U>& WAMIJ = H.getAIJK();
    const SpinorbitalTensor<U>& WAMEI = H.getAIBJ();

    auto& T = this->template get<ExcitationOperator<U,2> >("T");

    auto& XMI = this->template gettmp<SpinorbitalTensor<U> >("XMI");
    auto& XAE = this->template gettmp<SpinorbitalTensor<U> >("XAE");

    auto& D = this->template gettmp<Denominator<U> >("D");
    auto& davidson = this->template gettmp<Davidson<ExcitationOperator<U,2>>>("Davidson");

    auto& Rs = this->template gettmp<unique_vector<ExcitationOperator<U,2>>>("R");
    auto& Zs = this->template gettmp<unique_vector<ExcitationOperator<U,2>>>("Z");
    auto& Vs = this->template gettmp<unique_vector<ExcitationOperator<U,2>>>("V");

    for (int root = 0;root < this->nsolution();root++)
    {
        // cout << "am I here?" << endl;
        ExcitationOperator<U,2>& R = Rs[root];
        ExcitationOperator<U,2>& Z = Zs[root];
        Z = 0;

        for (auto& V : Vs)
        {
            R -= scalar(conj(R)*V)*V;
        }
        R /= sqrt(std::abs(scalar(conj(R)*R)));

         XMI[  "mi"]  =     WMNEJ["nmei"]*R(1)[  "en"];
         XMI[  "mi"] += 0.5*WMNEF["mnef"]*R(2)["efin"];
         XAE[  "ae"]  =     WAMEF["amef"]*R(1)[  "fm"];
         XAE[  "ae"] -= 0.5*WMNEF["mnef"]*R(2)["afmn"];

        Z(1)[  "ai"] +=       FAE[  "ae"]*R(1)[  "ei"];
        Z(1)[  "ai"] -=       FMI[  "mi"]*R(1)[  "am"];
        Z(1)[  "ai"] -=     WAMEI["amei"]*R(1)[  "em"];
        Z(1)[  "ai"] +=       FME[  "me"]*R(2)["aeim"];
        Z(1)[  "ai"] += 0.5*WAMEF["amef"]*R(2)["efim"];
        Z(1)[  "ai"] -= 0.5*WMNEJ["mnei"]*R(2)["eamn"];

        Z(2)["abij"] +=     WABEJ["abej"]*R(1)[  "ei"];
        Z(2)["abij"] -=     WAMIJ["amij"]*R(1)[  "bm"];
        Z(2)["abij"] +=       FAE[  "ae"]*R(2)["ebij"];
        Z(2)["abij"] -=       FMI[  "mi"]*R(2)["abmj"];
        Z(2)["abij"] +=       XAE[  "ae"]*T(2)["ebij"];
        Z(2)["abij"] -=       XMI[  "mi"]*T(2)["abmj"];
        Z(2)["abij"] += 0.5*WMNIJ["mnij"]*R(2)["abmn"];
        Z(2)["abij"] += 0.5*WABEF["abef"]*R(2)["efij"];
        Z(2)["abij"] -=     WAMEI["amei"]*R(2)["ebmj"];

        //0.5*Z(1)({0,0},{0,0})[  "ai"] += 0.5*Z(1)({1,0},{0,1})[  "ai"];
        //    Z(1)({1,0},{0,1})[  "ai"]  =     Z(1)({0,0},{0,0})[  "ai"];
        //0.5*Z(2)({1,0},{0,1})["abij"] += 0.5*Z(2)({1,0},{0,1})["baji"];
        //0.5*Z(2)({0,0},{0,0})["abij"] += 0.5*Z(2)({2,0},{0,2})["abij"];
        //    Z(2)({2,0},{0,2})["abij"]  =     Z(2)({0,0},{0,0})["abij"];
    }

    vector<U> energies = davidson.extrapolate(Rs, Zs, D);

    for (int i = 0;i < this->nsolution();i++)
    {
        this->energy(i) = energies[i];
        this->conv(i) = Zs[i].norm(00);
    }
}

INSTANTIATE_SPECIALIZATIONS(EOMEECCSD);
REGISTER_TASK(EOMEECCSD<double>, "eomeeccsd");
