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

template <typename T>
UHF<T>::UHF(const std::string& type, const std::string& name, const Config& config)
: Iterative(config), Task(type, name),
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

    int norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumAlphaElectrons();

    energy = molecule.getNuclearRepulsion();

    vector<int> shapeNN = vec(NS,NS);
    vector<int> sizenn = vec(norb,norb);
    vector<int> sizenO = vec(norb,nalpha);
    vector<int> sizeno = vec(norb,nbeta);
    vector<int> sizenV = vec(norb,norb-nalpha);
    vector<int> sizenv = vec(norb,norb-nbeta);

    put("Fa", new DistTensor<T>(arena, 2, sizenn, shapeNN, false));
    put("Fb", new DistTensor<T>(arena, 2, sizenn, shapeNN, false));
    put("Da", new DistTensor<T>(arena, 2, sizenn, shapeNN, true));
    put("Db", new DistTensor<T>(arena, 2, sizenn, shapeNN, true));

    puttmp("dF", new DistTensor<T>(arena, 2, sizenn, shapeNN, false));
    puttmp("Ca_occ", new DistTensor<T>(arena, 2, sizenO, shapeNN, false));
    puttmp("Cb_occ", new DistTensor<T>(arena, 2, sizeno, shapeNN, false));
    puttmp("Ca_vrt", new DistTensor<T>(arena, 2, sizenV, shapeNN, false));
    puttmp("Cb_vrt", new DistTensor<T>(arena, 2, sizenv, shapeNN, false));
    puttmp("dDa", new DistTensor<T>(arena, 2, sizenn, shapeNN, false));
    puttmp("dDb", new DistTensor<T>(arena, 2, sizenn, shapeNN, false));
    puttmp("S^-1/2", new DistTensor<T>(arena, 2, sizenn, shapeNN, false));

    calcSMinusHalf();

    tic();
    for (int i = 0;iterate();i++)
    {
        double dt = todouble(toc());
        Logger::log(arena) << "Iteration " << (i+1) << " took " << dt << " s" << endl;
        Logger::log(arena) << "Iteration " << (i+1) <<
                              " energy = " << setprecision(15) << energy <<
                              ", convergence = " << setprecision(8) << conv << endl;
        tic();
    }
    toc();

    if (!isConverged()) throw runtime_error("SCF did not converge");

    if (isUsed("S2") || isUsed("multiplicity"))
    {
        calcS2();
    }

    put("energy", new Scalar(arena, energy));
    put("convergence", new Scalar(arena, conv));

    put("occ", new MOSpace<T>(gettmp<DistTensor<T> >("Ca_occ"),
                              gettmp<DistTensor<T> >("Cb_occ")));
    put("vrt", new MOSpace<T>(gettmp<DistTensor<T> >("Ca_vrt"),
                              gettmp<DistTensor<T> >("Cb_vrt")));
}

template <typename T>
void UHF<T>::_iterate()
{
    buildFock();
    DIISExtrap();
    calcEnergy();
    diagonalizeFock();
    calcDensity();

    const Molecule& molecule = get<Molecule>("molecule");
    int norb = molecule.getNumOrbitals();

    DistTensor<T>& dDa = gettmp<DistTensor<T> >("dDa");
    DistTensor<T>& dDb = gettmp<DistTensor<T> >("dDb");

    switch (convtype)
    {
        case MAX_ABS:
            conv = max(dDa.norm(00), dDb.norm(00));
            break;
        case RMSD:
            conv = sqrt((dDa.norm(2)+dDb.norm(2))/(2*norb*norb));
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

    int norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumAlphaElectrons();

    DistTensor<T>& S = get<DistTensor<T> >("S");

    DistTensor<T>& Ca_occ = gettmp<DistTensor<T> >("Ca_occ");
    DistTensor<T>& Cb_occ = gettmp<DistTensor<T> >("Cb_occ");

    DistTensor<T> Delta(S.arena, 2, vec(nalpha,nbeta), vec(NS,NS), false);
    DistTensor<T> tmp(S.arena, 2, vec(nalpha,norb), vec(NS,NS), false);

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

    int norb = molecule.getNumOrbitals();

    DistTensor<T>& S = get<DistTensor<T> >("S");
    DistTensor<T>& Smhalf = gettmp<DistTensor<T> >("S^-1/2");

    vector<T> Ea(norb), Eb(norb);

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

    vector<T> s;
    vector<T> smhalf(norb*norb);

    S.getAllData(s);
    assert(s.size() == norb*norb);

    int info = heev('V', 'U', norb, s.data(), norb, Ea.data());
    assert(info == 0);

    fill(smhalf.begin(), smhalf.end(), 0.0);
    for (int i = 0;i < norb;i++)
    {
        ger(norb, norb, 1/sqrt(Ea[i]), &s[i*norb], 1, &s[i*norb], 1, smhalf.data(), norb);
    }

    if (S.arena.rank == 0)
    {
        vector< tkv_pair<T> > pairs;

        for (int i = 0;i < norb;i++)
        {
            for (int j = 0;j < norb;j++)
            {
                pairs.push_back(tkv_pair<T>(i+j*norb, smhalf[i+j*norb]));
            }
        }

        Smhalf.writeRemoteData(pairs);
    }
    else
    {
        Smhalf.writeRemoteData();
    }

    #endif
}

template <typename T>
void UHF<T>::diagonalizeFock()
{
    const Molecule& molecule = get<Molecule>("molecule");

    int norb = molecule.getNumOrbitals();
    int nalpha = molecule.getNumAlphaElectrons();
    int nbeta = molecule.getNumAlphaElectrons();

    DistTensor<T>& S = get<DistTensor<T> >("S");
    DistTensor<T>& Fa = get<DistTensor<T> >("Fa");
    DistTensor<T>& Fb = get<DistTensor<T> >("Fb");
    DistTensor<T>& Ca_occ = gettmp<DistTensor<T> >("Ca_occ");
    DistTensor<T>& Cb_occ = gettmp<DistTensor<T> >("Cb_occ");
    DistTensor<T>& Ca_vrt = gettmp<DistTensor<T> >("Ca_vrt");
    DistTensor<T>& Cb_vrt = gettmp<DistTensor<T> >("Cb_vrt");

    vector<T> Ea(norb), Eb(norb);
    vector< tkv_pair<T> > pairs;

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

    int info;
    vector<T> fock, s;

    S.getAllData(s);
    assert(s.size() == norb*norb);
    Fa.getAllData(fock);
    assert(fock.size() == norb*norb);
    info = hegv(AXBX, 'V', 'U', norb, fock.data(), norb, s.data(), norb, Ea.data());
    assert(info == 0);

    if (S.arena.rank == 0)
    {
        for (int i = 0;i < nalpha;i++)
        {
            for (int p = 0;p < norb;p++)
            {
                pairs.push_back(tkv_pair<T>(p+i*norb, fock[p+i*norb]));
            }
        }

        Ca_occ.writeRemoteData(pairs);
        pairs.clear();

        for (int i = 0;i < norb-nalpha;i++)
        {
            for (int p = 0;p < norb;p++)
            {
                pairs.push_back(tkv_pair<T>(p+i*norb, fock[p+(i+nalpha)*norb]));
            }
        }

        Ca_vrt.writeRemoteData(pairs);
        pairs.clear();
    }
    else
    {
        Ca_occ.writeRemoteData();
        Ca_vrt.writeRemoteData();
    }

    S.getAllData(s);
    assert(s.size() == norb*norb);
    Fb.getAllData(fock);
    assert(fock.size() == norb*norb);
    info = hegv(AXBX, 'V', 'U', norb, fock.data(), norb, s.data(), norb, Eb.data());
    assert(info == 0);

    if (S.arena.rank == 0)
    {
        for (int i = 0;i < nbeta;i++)
        {
            for (int p = 0;p < norb;p++)
            {
                pairs.push_back(tkv_pair<T>(p+i*norb, fock[p+i*norb]));
            }
        }

        Cb_occ.writeRemoteData(pairs);
        pairs.clear();

        for (int i = 0;i < norb-nbeta;i++)
        {
            for (int p = 0;p < norb;p++)
            {
                pairs.push_back(tkv_pair<T>(p+i*norb, fock[p+(i+nbeta)*norb]));
            }
        }

        Cb_vrt.writeRemoteData(pairs);
    }
    else
    {
        Cb_occ.writeRemoteData();
        Cb_vrt.writeRemoteData();
    }

    #endif

    fixPhase(Ca_occ);
    fixPhase(Cb_occ);
    fixPhase(Ca_vrt);
    fixPhase(Cb_vrt);
}

template <typename T>
void UHF<T>::fixPhase(DistTensor<T>& C)
{
    int norb = C.getLengths()[0];
    int nr = C.getLengths()[1];

    int nproc = C.arena.nproc;
    int rank = C.arena.rank;

    vector< tkv_pair<T> > pairs(norb);

    for (int b = 0;b*nproc < nr;b++)
    {
        int r = b*nproc+rank;

        if (r < nr)
        {
            for (int i = 0;i < norb;i++) pairs[i].k = i+r*norb;

            C.getRemoteData(pairs);

            sort(pairs.begin(), pairs.end());
            int sign = 0;
            for (int i = 0;i < norb;i++)
            {
                if (sign == 0 && abs(pairs[i].d) > 1e-10)
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

            C.writeRemoteData(pairs);
        }
        else
        {
            C.getRemoteData();
            C.writeRemoteData();
        }
    }
}

template <typename T>
void UHF<T>::calcEnergy()
{
    const Molecule& molecule = get<Molecule>("molecule");

    DistTensor<T>& H = get<DistTensor<T> >("H");
    DistTensor<T>& Fa = get<DistTensor<T> >("Fa");
    DistTensor<T>& Fb = get<DistTensor<T> >("Fb");
    DistTensor<T>& Da = get<DistTensor<T> >("Da");
    DistTensor<T>& Db = get<DistTensor<T> >("Db");

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
    DistTensor<T>& dDa = gettmp<DistTensor<T> >("dDa");
    DistTensor<T>& dDb = gettmp<DistTensor<T> >("dDb");
    DistTensor<T>& Da = get<DistTensor<T> >("Da");
    DistTensor<T>& Db = get<DistTensor<T> >("Db");
    DistTensor<T>& Ca_occ = gettmp<DistTensor<T> >("Ca_occ");
    DistTensor<T>& Cb_occ = gettmp<DistTensor<T> >("Cb_occ");

    /*
     * D[ab] = C[ai]*C[bi]
     */
    dDa["ab"]  = Da["ab"];
    dDb["ab"]  = Db["ab"];
     Da["ab"]  = Ca_occ["ai"]*Ca_occ["bi"];
     Db["ab"]  = Cb_occ["ai"]*Cb_occ["bi"];
    dDa["ab"] -= Da["ab"];
    dDb["ab"] -= Db["ab"];

    if (damping > 0.0)
    {
        Da["ab"] += damping*dDa["ab"];
        Db["ab"] += damping*dDb["ab"];
    }
}

template <typename T>
void UHF<T>::DIISExtrap()
{
    DistTensor<T>& S = get<DistTensor<T> >("S");
    DistTensor<T>& Smhalf = gettmp<DistTensor<T> >("S^-1/2");
    DistTensor<T>& dF = gettmp<DistTensor<T> >("dF");
    DistTensor<T>& Fa = get<DistTensor<T> >("Fa");
    DistTensor<T>& Fb = get<DistTensor<T> >("Fb");
    DistTensor<T>& Da = get<DistTensor<T> >("Da");
    DistTensor<T>& Db = get<DistTensor<T> >("Db");

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
        DistTensor<T> tmp1(Fa);
        DistTensor<T> tmp2(Fa);

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

    vector< DistTensor<T>* > Fab(2);
    Fab[0] = &Fa;
    Fab[1] = &Fb;
    diis.extrapolate(Fab, vector< DistTensor<T>* >(1, &dF));
}

INSTANTIATE_SPECIALIZATIONS(UHF);
