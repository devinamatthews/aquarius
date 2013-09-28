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

using namespace std;
using namespace aquarius;
using namespace aquarius::scf;
using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::integrals;
using namespace aquarius::task;
using namespace aquarius::time;
using namespace aquarius::op;
using namespace aquarius::symmetry;

template <typename T>
UHF<T>::UHF(const std::string& type, const std::string& name, const Config& config)
: Iterative(type, name, config),
	frozen_core(config.get<bool>("frozen_core")),
  damping(config.get<T>("damping")),
  diis(config.get("diis"), 2)
{
    vector<Requirement> reqs;
    reqs += Requirement("molecule", "molecule");
    reqs += Requirement("ovi", "S");
    reqs += Requirement("1ehamiltonian", "H");
    addProduct(Product("double", "energy", reqs));
    addProduct(Product("double", "convergence", reqs));
    addProduct(Product("double", "S2", reqs));
    addProduct(Product("double", "multiplicity", reqs));
    addProduct(Product("occspace", "occ", reqs));
    addProduct(Product("vrtspace", "vrt", reqs));
    addProduct(Product("Fa", "Fa", reqs));
    addProduct(Product("Fb", "Fb", reqs));
    addProduct(Product("Da", "Da", reqs));
    addProduct(Product("Db", "Db", reqs));
}

template <typename T>
void UHF<T>::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();

    energy = molecule.getNuclearRepulsion();

    vector<int> shapeNN = vec(NS,NS);
    vector<vector<int> > sizenn = vec(norb,norb);

    put("Fa", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, false));
    put("Fb", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, false));
    put("Da", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, true));
    put("Db", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, true));

    puttmp("dF", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, false));
    puttmp("Ca", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, false));
    puttmp("Cb", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, false));
    puttmp("dDa", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, false));
    puttmp("dDb", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, false));
    puttmp("S^-1/2", new SymmetryBlockedTensor<T>(arena, group, 2, sizenn, shapeNN, false));

    occ_alpha.resize(group.getNumIrreps());
    occ_beta.resize(group.getNumIrreps());

    E_alpha.resize(group.getNumIrreps());
    E_beta.resize(group.getNumIrreps());

    for (int i = 0;i < group.getNumIrreps();i++)
    {
        E_alpha[i].resize(norb[i]);
        E_beta[i].resize(norb[i]);
    }

    calcSMinusHalf();

    Iterative::run(dag, arena);

    if (isUsed("S2") || isUsed("multiplicity"))
    {
        calcS2();
    }

    put("energy", new Scalar(arena, energy));
    put("convergence", new Scalar(arena, conv));

    int nfrozen = 0;
    if (frozen_core)
    {
        for (vector<Atom>::const_iterator a = molecule.getAtomsBegin();a != molecule.getAtomsEnd();++a)
        {
            int Z = a->getCenter().getElement().getAtomicNumber();
            if      (Z > 86) nfrozen += 31;
            else if (Z > 54) nfrozen += 22;
            else if (Z > 36) nfrozen += 13;
            else if (Z > 18) nfrozen += 9;
            else if (Z > 10) nfrozen += 5;
            else if (Z >  2) nfrozen += 1;
        }
    }

    if (nfrozen > nalpha || nfrozen > nbeta)
        Logger::error(arena) << "There are not enough valence electrons for this multiplicity" << endl;

    vector<pair<typename std::real_type<T>::type,int> > E_alpha_occ;
    vector<pair<typename std::real_type<T>::type,int> > E_beta_occ;
    for (int i = 0;i < group.getNumIrreps();i++)
    {
        for (int j = 0;j < occ_alpha[i];j++)
        {
            E_alpha_occ.push_back(make_pair(E_alpha[i][j],i));
        }
        for (int j = 0;j < occ_beta[i];j++)
        {
            E_beta_occ.push_back(make_pair(E_beta[i][j],i));
        }
    }
    assert(E_alpha_occ.size() == nalpha);
    assert(E_beta_occ.size() == nbeta);

    sort(E_alpha_occ.begin(), E_alpha_occ.end());
    sort(E_beta_occ.begin(), E_beta_occ.end());

    vector<int> nfrozen_alpha(group.getNumIrreps());
    vector<int> nfrozen_beta(group.getNumIrreps());
    for (int i = 0;i < nfrozen;i++)
    {
        nfrozen_alpha[E_alpha_occ[i].second]++;
        nfrozen_beta[E_beta_occ[i].second]++;
    }

    vector<int> vrt_alpha(group.getNumIrreps());
    vector<int> vrt_beta(group.getNumIrreps());
    for (int i = 0;i < group.getNumIrreps();i++)
    {
        vrt_alpha[i] = norb[i]-occ_alpha[i];
        vrt_beta[i] = norb[i]-occ_beta[i];
        occ_alpha[i] -= nfrozen_alpha[i];
        occ_beta[i] -= nfrozen_beta[i];
    }

    vector<int> zero(norb.size(), 0);
    put("occ", new MOSpace<T>(new SymmetryBlockedTensor<T>(gettmp<SymmetryBlockedTensor<T> >("Ca"),
                                                           vec(zero,nfrozen_alpha),
                                                           vec(norb,occ_alpha)),
                              new SymmetryBlockedTensor<T>(gettmp<SymmetryBlockedTensor<T> >("Cb"),
                                                           vec(zero,nfrozen_beta),
                                                           vec(norb,occ_beta))));

   	put("vrt", new MOSpace<T>(new SymmetryBlockedTensor<T>(gettmp<SymmetryBlockedTensor<T> >("Ca"),
                                                           vec(zero,occ_alpha),
                                                           vec(norb,vrt_alpha)),
   	                          new SymmetryBlockedTensor<T>(gettmp<SymmetryBlockedTensor<T> >("Cb"),
                                                           vec(zero,occ_beta),
                                                           vec(norb,vrt_beta))));
}

template <typename T>
void UHF<T>::iterate()
{
    buildFock();
    DIISExtrap();
    calcEnergy();
    diagonalizeFock();
    calcDensity();

    const Molecule& molecule = get<Molecule>("molecule");
    int norb = sum(molecule.getNumOrbitals());

    SymmetryBlockedTensor<T>& dDa = gettmp<SymmetryBlockedTensor<T> >("dDa");
    SymmetryBlockedTensor<T>& dDb = gettmp<SymmetryBlockedTensor<T> >("dDb");

    switch (convtype)
    {
        case MAX_ABS:
            conv = max(dDa.norm(00), dDb.norm(00));
            break;
        case RMSD:
            conv = (dDa.norm(2)+dDb.norm(2))/sqrt(2*norb*norb);
            break;
        case MAD:
            conv = (dDa.norm(1)+dDb.norm(1))/(2*norb*norb);
            break;
    }
}

template <typename T>
void UHF<T>::calcS2()
{
    const Molecule& molecule = get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();

    SymmetryBlockedTensor<T>& S = get<SymmetryBlockedTensor<T> >("S");

    vector<int> zero(norb.size(), 0);
    SymmetryBlockedTensor<T> Ca_occ(gettmp<SymmetryBlockedTensor<T> >("Ca"),
                                    vec(zero,zero), vec(norb,occ_alpha));
    SymmetryBlockedTensor<T> Cb_occ(gettmp<SymmetryBlockedTensor<T> >("Cb"),
                                    vec(zero,zero), vec(norb,occ_beta));

    SymmetryBlockedTensor<T> Delta(S.arena, group, 2, vec(vec(nalpha),vec(nbeta)), vec(NS,NS), false);
    SymmetryBlockedTensor<T> tmp(S.arena, group, 2, vec(vec(nalpha),norb), vec(NS,NS), false);

    int ndiff = abs(nalpha-nbeta);
    int nmin = min(nalpha, nbeta);

    double S2 = ((ndiff/2)*(ndiff/2+1) + nmin);

    tmp["ai"] = Ca_occ["ja"]*S["ij"];
    Delta["ab"] = tmp["ai"]*Cb_occ["ib"];

    S2 -= abs(scalar(Delta*conj(Delta)));

    put("S2", new Scalar(S.arena, S2));
    put("multiplicity", new Scalar(S.arena, sqrt(4*S2+1)));
}

template <typename T>
void UHF<T>::calcSMinusHalf()
{
    const Molecule& molecule = get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();

    SymmetryBlockedTensor<T>& S = get<SymmetryBlockedTensor<T> >("S");
    SymmetryBlockedTensor<T>& Smhalf = gettmp<SymmetryBlockedTensor<T> >("S^-1/2");

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        //cout << "S " << (i+1) << endl;
        //vector<T> vals;
        //S(vec(i,i)).getAllData(vals);
        //printmatrix(norb[i], norb[i], vals.data(), 6, 3, 108);
    }

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        if (norb[i] == 0) continue;

        vector<int> irreps(2,i);
        vector<typename real_type<T>::type> E(norb[i]);

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

        if (S.arena.rank == 0)
        {
            vector<T> s;
            vector<T> smhalf(norb[i]*norb[i]);

            S(irreps).getAllData(s,0);
            assert(s.size() == norb[i]*norb[i]);

            int info = heev('V', 'U', norb[i], s.data(), norb[i], E.data());
            assert(info == 0);

            fill(smhalf.begin(), smhalf.end(), (T)0);
            for (int j = 0;j < norb[i];j++)
            {
                ger(norb[i], norb[i], 1/sqrt(E[j]), &s[j*norb[i]], 1, &s[j*norb[i]], 1, smhalf.data(), norb[i]);
            }

            vector< tkv_pair<T> > pairs(norb[i]*norb[i]);

            for (int j = 0;j < norb[i]*norb[i];j++)
            {
                pairs[j].k = j;
                pairs[j].d = smhalf[j];
            }

            Smhalf(irreps).writeRemoteData(pairs);
        }
        else
        {
            S(irreps).getAllData(0);
            Smhalf(irreps).writeRemoteData();
        }

        #endif
    }
}

template <typename T>
void UHF<T>::diagonalizeFock()
{
    const Molecule& molecule = get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();

    SymmetryBlockedTensor<T>& S = get<SymmetryBlockedTensor<T> >("S");
    SymmetryBlockedTensor<T>& Fa = get<SymmetryBlockedTensor<T> >("Fa");
    SymmetryBlockedTensor<T>& Fb = get<SymmetryBlockedTensor<T> >("Fb");
    SymmetryBlockedTensor<T>& Ca = gettmp<SymmetryBlockedTensor<T> >("Ca");
    SymmetryBlockedTensor<T>& Cb = gettmp<SymmetryBlockedTensor<T> >("Cb");

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        //cout << "F " << (i+1) << endl;
        //vector<T> vals;
        //Fa(vec(i,i)).getAllData(vals);
        //printmatrix(norb[i], norb[i], vals.data(), 6, 3, 108);
    }

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

        vector< tkv_pair<T> > pairs_occ;
        vector< tkv_pair<T> > pairs_vrt;

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

        vector< tkv_pair<T> > pairs_occ;
        vector< tkv_pair<T> > pairs_vrt;

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

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        if (norb[i] == 0) continue;

        vector<int> irreps(2,i);

        if (S.arena.rank == 0)
        {
            int info;
            vector<T> fock, s;
            vector< tkv_pair<T> > pairs(norb[i]*norb[i]);

            S(irreps).getAllData(s,0);
            assert(s.size() == norb[i]*norb[i]);
            vector<T> tmp(s);

            Fa(irreps).getAllData(fock,0);
            assert(fock.size() == norb[i]*norb[i]);
            info = hegv(AXBX, 'V', 'U', norb[i], fock.data(), norb[i], tmp.data(), norb[i], E_alpha[i].data());
            assert(info == 0);

            for (int j = 0;j < norb[i];j++)
            {
                T sign = 0;
                for (int k = 0;k < norb[i];k++)
                {
                    if (abs(fock[k+j*norb[i]]) > 1e-10)
                    {
                        sign = (fock[k+j*norb[i]] < 0 ? -1 : 1);
                        break;
                    }
                }
                scal(norb[i], sign, &fock[j*norb[i]], 1);
            }

            for (int j = 0;j < norb[i]*norb[i];j++)
            {
                pairs[j].k = j;
                pairs[j].d = fock[j];
            }

            Ca(irreps).writeRemoteData(pairs);

            Fb(irreps).getAllData(fock,0);
            assert(fock.size() == norb[i]*norb[i]);
            info = hegv(AXBX, 'V', 'U', norb[i], fock.data(), norb[i], s.data(), norb[i], E_beta[i].data());
            assert(info == 0);

            for (int j = 0;j < norb[i];j++)
            {
                T sign = 0;
                for (int k = 0;k < norb[i];k++)
                {
                    if (abs(fock[k+j*norb[i]]) > 1e-10)
                    {
                        sign = (fock[k+j*norb[i]] < 0 ? -1 : 1);
                        break;
                    }
                }
                scal(norb[i], sign, &fock[j*norb[i]], 1);
            }

            for (int j = 0;j < norb[i]*norb[i];j++)
            {
                pairs[j].k = j;
                pairs[j].d = fock[j];
            }

            Cb(irreps).writeRemoteData(pairs);
        }
        else
        {
            S(irreps).getAllData(0);
            Fa(irreps).getAllData(0);
            Fb(irreps).getAllData(0);
            Ca(irreps).writeRemoteData();
            Cb(irreps).writeRemoteData();
        }
    }

    #endif

    vector<pair<typename real_type<T>::type,int> > E_alpha_sorted;
    vector<pair<typename real_type<T>::type,int> > E_beta_sorted;
    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        occ_alpha[i] = 0;
        occ_beta[i] = 0;
        for (int j = 0;j < norb[i];j++)
        {
            E_alpha_sorted.push_back(make_pair(E_alpha[i][j],i));
        }
        for (int j = 0;j < norb[i];j++)
        {
            E_beta_sorted.push_back(make_pair(E_beta[i][j],i));
        }
    }

    sort(E_alpha_sorted.begin(), E_alpha_sorted.end());
    sort(E_beta_sorted.begin(), E_beta_sorted.end());

    for (int i = 0;i < nalpha;i++)
    {
        occ_alpha[E_alpha_sorted[i].second]++;
    }
    for (int i = 0;i < nbeta;i++)
    {
        occ_beta[E_beta_sorted[i].second]++;
    }

    log(S.arena) << "Iteration " << iter << " occupation = " << occ_alpha << ", " << occ_beta << endl;
}

template <typename T>
void UHF<T>::calcEnergy()
{
    const Molecule& molecule = get<Molecule>("molecule");

    SymmetryBlockedTensor<T>& H = get<SymmetryBlockedTensor<T> >("H");
    SymmetryBlockedTensor<T>& Fa = get<SymmetryBlockedTensor<T> >("Fa");
    SymmetryBlockedTensor<T>& Fb = get<SymmetryBlockedTensor<T> >("Fb");
    SymmetryBlockedTensor<T>& Da = get<SymmetryBlockedTensor<T> >("Da");
    SymmetryBlockedTensor<T>& Db = get<SymmetryBlockedTensor<T> >("Db");

    /*
     * E = (1/2)Tr[D(F+H)]
     *
     *   = (1/2)Tr[Da*(Fa+H) + Db*(Fb+H)]
     */
    Fa["ab"] += H["ab"];
    Fb["ab"] += H["ab"];
    energy = molecule.getNuclearRepulsion();
    energy += 0.5*scalar(Da["ab"]*Fa["ab"]);
    energy += 0.5*scalar(Db["ab"]*Fb["ab"]);
    Fa["ab"] -= H["ab"];
    Fb["ab"] -= H["ab"];
}

template <typename T>
void UHF<T>::calcDensity()
{
    const Molecule& molecule = get<Molecule>("molecule");

    const vector<int>& norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumBetaElectrons();

    SymmetryBlockedTensor<T>& dDa = gettmp<SymmetryBlockedTensor<T> >("dDa");
    SymmetryBlockedTensor<T>& dDb = gettmp<SymmetryBlockedTensor<T> >("dDb");
    SymmetryBlockedTensor<T>& Da = get<SymmetryBlockedTensor<T> >("Da");
    SymmetryBlockedTensor<T>& Db = get<SymmetryBlockedTensor<T> >("Db");

    vector<int> zero(norb.size(), 0);
    SymmetryBlockedTensor<T> Ca_occ(gettmp<SymmetryBlockedTensor<T> >("Ca"),
                                    vec(zero,zero), vec(norb,occ_alpha));
    SymmetryBlockedTensor<T> Cb_occ(gettmp<SymmetryBlockedTensor<T> >("Cb"),
                                    vec(zero,zero), vec(norb,occ_beta));

    /*
     * D[ab] = C[ai]*C[bi]
     */
    dDa["ab"]  = Da["ab"];
    dDb["ab"]  = Db["ab"];
     Da = 0;
     Db = 0;
     Da["ab"] += Ca_occ["ai"]*Ca_occ["bi"];
     Db["ab"] += Cb_occ["ai"]*Cb_occ["bi"];
    dDa["ab"] -= Da["ab"];
    dDb["ab"] -= Db["ab"];

    for (int i = 0;i < molecule.getGroup().getNumIrreps();i++)
    {
        //cout << "D " << (i+1) << endl;
        //vector<T> vals;
        //Da(vec(i,i)).getAllData(vals);
        //scal(vals.size(), 2, vals.data(), 1);
        //printmatrix(norb[i], norb[i], vals.data(), 6, 3, 108);
    }

    if (damping > 0.0)
    {
        Da["ab"] += damping*dDa["ab"];
        Db["ab"] += damping*dDb["ab"];
    }
}

template <typename T>
void UHF<T>::DIISExtrap()
{
    SymmetryBlockedTensor<T>& S = get<SymmetryBlockedTensor<T> >("S");
    SymmetryBlockedTensor<T>& Smhalf = gettmp<SymmetryBlockedTensor<T> >("S^-1/2");
    SymmetryBlockedTensor<T>& dF = gettmp<SymmetryBlockedTensor<T> >("dF");
    SymmetryBlockedTensor<T>& Fa = get<SymmetryBlockedTensor<T> >("Fa");
    SymmetryBlockedTensor<T>& Fb = get<SymmetryBlockedTensor<T> >("Fb");
    SymmetryBlockedTensor<T>& Da = get<SymmetryBlockedTensor<T> >("Da");
    SymmetryBlockedTensor<T>& Db = get<SymmetryBlockedTensor<T> >("Db");

    /*
     * Generate the residual:
     *
     *   R = FDS - SDF
     *
     * Then, convert to the orthonormal basis:
     *
     *   ~    -1/2    -1/2
     *   R = S     R S.
     *
     * In this basis we have
     *
     *   ~    -1/2    -1/2  ~    1/2
     *   F = S     F S    , C = S    C, and
     *
     *   ~   ~ ~T    1/2    T  1/2    1/2    1/2
     *   D = C C  = S    C C  S    = S    D S.
     *
     * And so,
     *
     *   ~    -1/2    -1/2  1/2    1/2    1/2    1/2  -1/2    -1/2
     *   R = S     F S     S    D S    - S    D S    S     F S
     *
     *        ~ ~
     *     = [F,D] = 0 at convergence.
     */
    {
        SymmetryBlockedTensor<T> tmp1(Fa);
        SymmetryBlockedTensor<T> tmp2(Fa);

        tmp1["ab"]  =     Fa["ac"]*    Da["cb"];
        tmp2["ab"]  =   tmp1["ac"]*     S["cb"];
        tmp1["ab"]  =      S["ac"]*    Da["cb"];
        tmp2["ab"] -=   tmp1["ac"]*    Fa["cb"];
        tmp1["ab"]  = Smhalf["ac"]*  tmp2["cb"];
          dF["ab"]  =   tmp1["ac"]*Smhalf["cb"];

        tmp1["ab"]  =     Fb["ac"]*    Db["cb"];
        tmp2["ab"]  =   tmp1["ac"]*     S["cb"];
        tmp1["ab"]  =      S["ac"]*    Db["cb"];
        tmp2["ab"] -=   tmp1["ac"]*    Fb["cb"];
        tmp1["ab"]  = Smhalf["ac"]*  tmp2["cb"];
          dF["ab"] +=   tmp1["ac"]*Smhalf["cb"];
    }

    vector< SymmetryBlockedTensor<T>* > Fab(2);
    Fab[0] = &Fa;
    Fab[1] = &Fb;
    diis.extrapolate(Fab, vector< SymmetryBlockedTensor<T>* >(1, &dF));
}

INSTANTIATE_SPECIALIZATIONS(UHF);
