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

#include "spinorbital_tensor.hpp"

using namespace std;
using namespace aquarius::tensor;
using namespace aquarius::autocc;

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const SpinorbitalTensor<T>& t, const T val)
: IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T >(0, 0)
{
    nA = 0;
    nM = 0;
    nE = 0;
    nI = 0;
    spin = 0;

    addSpinCase(new DistTensor<T>(t(0), val), ",", "");
}

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const SpinorbitalTensor<T>& other)
: IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T>(other.ndim_, 0)
{
    logical = other.logical;
    nA = other.nA;
    nM = other.nM;
    nE = other.nE;
    nI = other.nI;
    spin = other.spin;

    for (typename vector<SpinCase>::const_iterator sc = other.cases.begin();sc != other.cases.end();++sc)
    {
        SpinCase newsc(*(new DistTensor<T>(*sc->tensor)));
        newsc.logical = sc->logical;
        newsc.log_to_phys = sc->log_to_phys;
        newsc.nA = sc->nA;
        newsc.nM = sc->nM;
        newsc.nE = sc->nE;
        newsc.nI = sc->nI;
        newsc.permFactor = sc->permFactor;
        cases.push_back(newsc);
        tensors_.push_back(typename IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T>::TensorRef(newsc.tensor, true));
    }
}

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const Manifold& left, const Manifold& right, const int spin)
: IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T>(left.np+left.nh+right.np+right.nh, 0), spin(spin)
{
    vector<Line> out_;
    vector<Line> in_;

    for (int i = 0;i <  left.np;i++) out_ += Line(i         , PARTICLE+EXTERNAL+ALPHA);
    for (int i = 0;i < right.nh;i++) out_ += Line(i         , HOLE    +EXTERNAL+ALPHA);
    for (int i = 0;i < right.np;i++)  in_ += Line(i+left.np , PARTICLE+EXTERNAL+ALPHA);
    for (int i = 0;i <  left.nh;i++)  in_ += Line(i+right.nh, HOLE    +EXTERNAL+ALPHA);

    this->logical = out_+in_;

    vector<Line> out(out_);
    vector<Line> in(in_);

    sort(out.begin(), out.end());
    sort(in.begin(), in.end());

    //permFactor = relativeSign(out_, out) *
    //             relativeSign(in_, in);

    nA = count_if(out.begin(), out.end(), isParticle());
    nM = count_if(out.begin(), out.end(), isHole());
    nE = count_if(in.begin(), in.end(), isParticle());
    nI = count_if(in.begin(), in.end(), isHole());
}

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const string& logical, const int spin)
: IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T>(logical.size()-1, 0), spin(spin)
{
    int comma = logical.find(',');
    if (comma == string::npos) throw logic_error("index string is malformed: " + logical);

    if (logical != tolower(logical)) throw logic_error("spinorbital indices must be lowercase");

    this->logical = Line::parse(logical.substr(0,comma) + logical.substr(comma+1));

    vector<Line> out_ = Line::parse(logical.substr(0,comma));
    vector<Line> in_  = Line::parse(logical.substr(comma+1));

    vector<Line> out(out_);
    vector<Line> in(in_);

    sort(out.begin(), out.end());
    sort(in.begin(), in.end());

    //permFactor = relativeSign(out_, out) *
    //             relativeSign(in_, in);

    nA = count_if(out.begin(), out.end(), isParticle());
    nM = count_if(out.begin(), out.end(), isHole());
    nE = count_if(in.begin(), in.end(), isParticle());
    nI = count_if(in.begin(), in.end(), isHole());
}

template<class T>
void SpinorbitalTensor<T>::set_name(char const * name_){
    int i;
    for (i=0; i<cases.size(); i++){
        cases[i].tensor->set_name(name_);
    }
}


template<class T>
void SpinorbitalTensor<T>::addSpinCase(DistTensor<T>* tensor, string logical, string physical, double factor, bool isAlloced)
{
    addSpinCase(*tensor, logical, physical, factor, isAlloced);
}

template<class T>
void SpinorbitalTensor<T>::addSpinCase(DistTensor<T>& tensor, string logical, string physical, double factor, bool isAlloced)
{
    addSpinCase(tensor, Line::parse(logical.substr(0,nA+nM)+logical.substr(nA+nM+1,nE+nI)),
                Line::parse(physical), factor, isAlloced);
}

template<class T>
void SpinorbitalTensor<T>::addSpinCase(DistTensor<T>* tensor, const Manifold& alpha_left,
                 const Manifold& alpha_right, double factor, bool isAlloced)
{
    addSpinCase(*tensor, alpha_left, alpha_right, factor, isAlloced);
}

template<class T>
void SpinorbitalTensor<T>::addSpinCase(DistTensor<T>& tensor, const Manifold& alpha_left,
                 const Manifold& alpha_right, double factor, bool isAlloced)
{
    vector<Line> logical;

    for (int i = 0;i <     alpha_left.np;i++) logical += Line(i                  , PARTICLE+EXTERNAL+ALPHA);
    for (int i = 0;i <  nA-alpha_left.np;i++) logical += Line(i+alpha_left.np    , PARTICLE+EXTERNAL+BETA);
    for (int i = 0;i <    alpha_right.nh;i++) logical += Line(i                  , HOLE    +EXTERNAL+ALPHA);
    for (int i = 0;i < nM-alpha_right.nh;i++) logical += Line(i+alpha_right.nh   , HOLE    +EXTERNAL+BETA);
    for (int i = 0;i <    alpha_right.np;i++) logical += Line(i+nA               , PARTICLE+EXTERNAL+ALPHA);
    for (int i = 0;i < nE-alpha_right.np;i++) logical += Line(i+nA+alpha_right.np, PARTICLE+EXTERNAL+BETA);
    for (int i = 0;i <     alpha_left.nh;i++) logical += Line(i+nM               , HOLE    +EXTERNAL+ALPHA);
    for (int i = 0;i <  nI-alpha_left.nh;i++) logical += Line(i+nM+alpha_left.nh , HOLE    +EXTERNAL+BETA);

    addSpinCase(tensor, logical, logical, factor, isAlloced);
}

template<class T>
void SpinorbitalTensor<T>::addSpinCase(DistTensor<T>& tensor, vector<Line> logical, vector<Line> physical, double factor, bool isAlloced)
{
    SpinCase sc(tensor);
    tensors_.push_back(typename IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T>::TensorRef(&tensor, isAlloced));

    sc.logical = logical;

    if (this->logical.size() != sc.logical.size()) throw logic_error("wrong number of indices");

    for (vector<Line>::const_iterator l1 = sc.logical.begin(), l2 = this->logical.begin();l1 != sc.logical.end();++l1, ++l2)
    {
        if (l1->toAlpha() != l2->toAlpha()) throw logic_error("UHF indices do not match spinorbital");
    }

    vector<Line> out_ = slice(logical, 0, nA+nM);
    vector<Line> in_  = slice(logical, nA+nM, nA+nM+nE+nI);
    vector<Line> phys = physical;

    vector<Line> out(out_);
    vector<Line> in(in_);

    sort(out.begin(), out.end());
    sort(in.begin(), in.end());

    //cout << logical << ' ' << out << ' ' << in << endl;

    if (count_if(out.begin(), out.end(), isParticle()) != nA ||
        count_if(out.begin(), out.end(),     isHole()) != nM ||
        count_if( in.begin(),  in.end(), isParticle()) != nE ||
        count_if( in.begin(),  in.end(),     isHole()) != nI)
    {
        throw logic_error("tensor has the wrong shape");
    }

    sc.permFactor = relativeSign(out_, out) *
                    relativeSign(in_, in) * factor;

    //cout << (out_ + in_) << " = " << sc.permFactor << '*' <<
    //             (out  + in ) << endl << endl;

    if (out.size()+in.size() != phys.size())
        throw logic_error("logical and physical dimensions do not match");

    sc.nA = 0;
    sc.nM = 0;
    sc.nE = 0;
    sc.nI = 0;

    for (int i = 0;i < out.size();i++)
    {
        if (out[i].isAlpha())
        {
            if (out[i].isParticle()) sc.nA++;
            else                     sc.nM++;
        }

        sc.log_to_phys.push_back(find(phys.begin(), phys.end(), out[i])-phys.begin());
        if (sc.log_to_phys.back() == phys.size())
            throw logic_error("logical index not found in physical map: " + str(logical[i]));
    }

    for (int i = 0;i < in.size();i++)
    {
        if (in[i].isAlpha())
        {
            if (in[i].isParticle()) sc.nE++;
            else                    sc.nI++;
        }

        sc.log_to_phys.push_back(find(phys.begin(), phys.end(), in[i])-phys.begin());
        if (sc.log_to_phys.back() == phys.size())
            throw logic_error("logical index not found in physical map: " + str(logical[i+out.size()]));
    }

    if (abs(sc.nA+sc.nM-sc.nE-sc.nI) > abs((int)(out.size()-in.size())))
        throw logic_error("spin case is not valid");

    cases.push_back(sc);

    if (2*(sc.nA+sc.nM-sc.nE-sc.nI)-(nA+nM-nE-nI) != spin)
        throw logic_error("spin is not compatible");

    //cout << "Adding spin case " << sc.nA << ' ' << sc.nM << ' ' << sc.nE << ' ' << sc.nI << endl;
}

template<class T>
void SpinorbitalTensor<T>::matchTypes(const int nin_A, const int nout_A, const vector<Line>& log_A, const int* idx_A,
                                      const int nin_B, const int nout_B, const vector<Line>& log_B, const int* idx_B,
                                      vector<Line>& out_A, vector<Line>& in_A,
                                      vector<Line>& pout_A, vector<Line>& hout_A,
                                      vector<Line>& pin_A, vector<Line>& hin_A,
                                      vector<Line>& sum_A)
{
    for (int i = 0;i < nout_A;i++)
    {
        int j;
        for (j = 0;j < nout_B+nin_B;j++)
        {
            if (idx_A[i] == idx_B[j])
            {
                bool aisp = log_A[i].isParticle();
                bool bisp = log_B[j].isParticle();
                if (aisp != bisp)
                    throw logic_error("types do not match");
                out_A[i] = Line(idx_A[i], log_B[j].getType());
                if (log_B[j].isParticle())
                {
                    pin_A += Line(idx_A[i], log_B[j].getType());
                }
                else
                {
                    hin_A += Line(idx_A[i], log_B[j].getType());
                }
                break;
            }
        }
        if (j == nout_B+nin_B)
        {
            out_A[i] = Line(idx_A[i], log_A[i].getType());
            sum_A += Line(idx_A[i], log_A[i].getType());
        }
    }

    for (int i = nout_A;i < nout_A+nin_A;i++)
    {
        int j;
        for (j = 0;j < nout_B+nin_B;j++)
        {
            if (idx_A[i] == idx_B[j])
            {
                bool aisp = log_A[i].isParticle();
                bool bisp = log_B[j].isParticle();
                if (aisp != bisp)
                    throw logic_error("types do not match");
                in_A[i-nout_A] = Line(idx_A[i], log_B[j].getType());
                if (log_B[j].isParticle())
                {
                    pout_A += Line(idx_A[i], log_B[j].getType());
                }
                else
                {
                    hout_A += Line(idx_A[i], log_B[j].getType());
                }
                break;
            }
        }
        if (j == nout_B+nin_B)
        {
            in_A[i-nout_A] = Line(idx_A[i], log_A[i].getType());
            sum_A += Line(idx_A[i], log_A[i].getType());
        }
    }
}

template<class T>
DistTensor<T>& SpinorbitalTensor<T>::operator()(int nA, int nM, int nE, int nI)
{
    return const_cast<DistTensor<T>&>(const_cast<const SpinorbitalTensor<T>&>(*this)(nA, nM, nE, nI));
}

template<class T>
const DistTensor<T>& SpinorbitalTensor<T>::operator()(int nA, int nM, int nE, int nI) const
{
    for (typename vector<SpinCase>::const_iterator sc = cases.begin();sc != cases.end();++sc)
    {
        if (sc->nA == nA && sc->nM == nM &&
            sc->nE == nE && sc->nI == nI) return *sc->tensor;
    }

    throw logic_error("spin case not found");
}

template<class T>
void SpinorbitalTensor<T>::mult(const T alpha, bool conja, const SpinorbitalTensor<T>& A_, const int* idx_A,
                                               bool conjb, const SpinorbitalTensor<T>& B_, const int* idx_B,
                                const T beta_,                                             const int* idx_C)
{
    const SpinorbitalTensor<T>& A = A_.getDerived();
    const SpinorbitalTensor<T>& B = B_.getDerived();

    int *idx_A_ = new int[A.ndim_];
    int *idx_B_ = new int[B.ndim_];
    int *idx_C_ = new int[ndim_];

    /*
    cout << "contracting: " << alpha << " * " << A.logical << "[";
    for (int i = 0;i < A.ndim_;i++) cout << idx_A[i] << ' ';
    cout << "] " << B.logical << "[";
    for (int i = 0;i < B.ndim_;i++) cout << idx_B[i] << ' ';
    cout << "] -> " << beta_ << " " << logical << "[";
    for (int i = 0;i <   ndim_;i++) cout << idx_C[i] << ' ';
    cout << "]\n" << endl;
    */

    vector<T> beta(cases.size(), beta_);

    vector<Line> sAin(A.nA+A.nM);
    vector<Line> sAout(A.nE+A.nI);

    vector<Line> sBin(B.nA+B.nM);
    vector<Line> sBout(B.nE+B.nI);

    vector<Line> sCin(nA+nM);
    vector<Line> sCout(nE+nI);

    for (typename vector<SpinCase>::const_iterator scC = cases.begin();scC != cases.end();++scC)
    {
        vector<Line> sum;
        vector<Line> sum2;
        vector<Line> apo, bpo;
        vector<Line> aho, bho;
        vector<Line> api, bpi;
        vector<Line> ahi, bhi;

        matchTypes(A.nA+A.nM, A.nE+A.nI,    A.logical, idx_A,
                     nA+  nM,   nE+  nI, scC->logical, idx_C,
                   sAout, sAin, apo, aho, api, ahi, sum);

        matchTypes(B.nA+B.nM, B.nE+B.nI,    B.logical, idx_B,
                     nA+  nM,   nE+  nI, scC->logical, idx_C,
                   sBout, sBin, bpo, bho, bpi, bhi, sum);

        uniq(sum);
        for (int i = 0;i < sum.size();i++)
        {
            sum2 += sum[i].toAlpha();
            sum2 += sum[i].toBeta();
        }

        {
            int i;
            for (i = 0;i < A.nA+A.nM;i++)
            {
                int j;
                for (j = 0;j < B.nA+B.nM;j++)
                {
                    if (idx_A[i] == idx_B[j])
                    {
                        bool aisp = A.logical[i].isParticle();
                        bool bisp = B.logical[j].isParticle();
                        if (aisp != bisp)
                                throw logic_error("types do not match");
                        sAout[i] = sBout[j];
                    }
                }
                for (;j < B.ndim_;j++)
                {
                    if (idx_A[i] == idx_B[j])
                    {
                        bool aisp = A.logical[i].isParticle();
                        bool bisp = B.logical[j].isParticle();
                        if (aisp != bisp)
                                throw logic_error("types do not match");
                        sAout[i] = sBin[j-B.nA-B.nM];
                    }
                }
            }
            for (;i < A.ndim_;i++)
            {
                int j;
                for (j = 0;j < B.nA+B.nM;j++)
                {
                    if (idx_A[i] == idx_B[j])
                    {
                        bool aisp = A.logical[i].isParticle();
                        bool bisp = B.logical[j].isParticle();
                        if (aisp != bisp)
                                throw logic_error("types do not match");
                        sAin[i-A.nA-A.nM] = sBout[j];
                    }
                }
                for (;j < B.ndim_;j++)
                {
                    if (idx_A[i] == idx_B[j])
                    {
                        bool aisp = A.logical[i].isParticle();
                        bool bisp = B.logical[j].isParticle();
                        if (aisp != bisp)
                                throw logic_error("types do not match");
                        sAin[i-A.nA-A.nM] = sBin[j-B.nA-B.nM];
                    }
                }
            }
        }

        for (int i = 0;i < nA+nM;i++)
        {
            sCout[i] = Line(idx_C[i], scC->logical[i].getType());
        }

        for (int i = 0;i < nE+nI;i++)
        {
            sCin[i] = Line(idx_C[nA+nM+i], scC->logical[nA+nM+i].getType());
        }

        Fragment uhfsA("T", sAout, sAin);
        Fragment uhfsB("U", sBout, sBin);
        Fragment uhfsC("V", sCout, sCin);

        vector<Line> isect;

        uniq(apo);
        uniq(bpo);
        isect = intersection(apo, bpo);
        exclude(apo, isect);
        exclude(bpo, isect);
        uniq(aho);
        uniq(bho);
        isect = intersection(aho, bho);
        exclude(aho, isect);
        exclude(bho, isect);
        uniq(api);
        uniq(bpi);
        isect = intersection(api, bpi);
        exclude(api, isect);
        exclude(bpi, isect);
        uniq(ahi);
        uniq(bhi);
        isect = intersection(ahi, bhi);
        exclude(ahi, isect);
        exclude(bhi, isect);

        //cout << apo << "," << bpo << "|" << aho << "," << bho << "|" <<
        //             api << "," << bpi << "|" << ahi << "," << bhi << endl;

        vector< vector<Line> > assym(2);

        Diagram d(Diagram::UHF);
        d += Term(Diagram::UHF)*uhfsA*uhfsB;
        assym[0] = apo;
        assym[1] = bpo;
        if (!apo.empty() && !bpo.empty()) d.antisymmetrize(assym);
        assym[0] = aho;
        assym[1] = bho;
        if (!aho.empty() && !bho.empty()) d.antisymmetrize(assym);
        assym[0] = api;
        assym[1] = bpi;
        if (!api.empty() && !bpi.empty()) d.antisymmetrize(assym);
        assym[0] = ahi;
        assym[1] = bhi;
        if (!ahi.empty() && !bhi.empty()) d.antisymmetrize(assym);
        d.sum(sum);
        d.fixorder(sum2);

        /*
         * Remove terms which are antisymmetrizations of same-spin groups
         */
        vector<Term> terms = d.getTerms();
        for (vector<Term>::iterator t1 = terms.begin();t1 != terms.end();++t1)
        {
            for (vector<Term>::iterator t2 = t1+1;t2 != terms.end();++t2)
            {
                if (Term(*t1).fixorder(true) == Term(*t2).fixorder(true))
                {
                    d -= *t1;
                    break;
                }
            }
        }

        terms = d.getTerms();
        for (vector<Term>::iterator t = terms.begin();t != terms.end();++t)
        {
            *t *= uhfsC;

            //cout << *t << endl;

            double diagFactor = t->getFactor();

            vector<Fragment>::iterator fA, fB, fC;
            for (vector<Fragment>::iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
            {
                if (f->getOp() == "T") fA = f;
                if (f->getOp() == "U") fB = f;
                if (f->getOp() == "V") fC = f;
            }

            vector<Line> out_A = fA->getIndicesOut();
            vector<Line>  in_A = fA->getIndicesIn();
            vector<Line> out_B = fB->getIndicesOut();
            vector<Line>  in_B = fB->getIndicesIn();
            vector<Line> out_C = fC->getIndicesOut();
            vector<Line>  in_C = fC->getIndicesIn();

            int nA_A = count_if(out_A.begin(), out_A.end(), isType<PARTICLE+ALPHA>());
            int nM_A = count_if(out_A.begin(), out_A.end(), isType<    HOLE+ALPHA>());
            int nE_A = count_if( in_A.begin(),  in_A.end(), isType<PARTICLE+ALPHA>());
            int nI_A = count_if( in_A.begin(),  in_A.end(), isType<    HOLE+ALPHA>());

            int nA_B = count_if(out_B.begin(), out_B.end(), isType<PARTICLE+ALPHA>());
            int nM_B = count_if(out_B.begin(), out_B.end(), isType<    HOLE+ALPHA>());
            int nE_B = count_if( in_B.begin(),  in_B.end(), isType<PARTICLE+ALPHA>());
            int nI_B = count_if( in_B.begin(),  in_B.end(), isType<    HOLE+ALPHA>());

            typename vector<SpinCase>::const_iterator scA = A.cases.end();
            for (typename vector<SpinCase>::const_iterator sc = A.cases.begin();sc != A.cases.end();++sc)
            {
                if (nA_A == sc->nA && nM_A == sc->nM &&
                    nE_A == sc->nE && nI_A == sc->nI) scA = sc;
            }
            if (scA == A.cases.end()) throw logic_error("no matching spin case for tensor A");

            typename vector<SpinCase>::const_iterator scB = B.cases.end();
            for (typename vector<SpinCase>::const_iterator sc = B.cases.begin();sc != B.cases.end();++sc)
            {
                if (nA_B == sc->nA && nM_B == sc->nM &&
                    nE_B == sc->nE && nI_B == sc->nI) scB = sc;
            }
            if (scB == B.cases.end()) throw logic_error("no matching spin case for tensor B");

            vector<Line> lA(A.ndim_);
            for (int i = 0;i < A.nA+A.nM;i++)
            {
                idx_A_[scA->log_to_phys[i]] = out_A[i].toInt();
                lA[scA->log_to_phys[i]] = out_A[i];
            }

            for (int i = A.nA+A.nM;i < A.ndim_;i++)
            {
                idx_A_[scA->log_to_phys[i]] = in_A[i-A.nA-A.nM].toInt();
                lA[scA->log_to_phys[i]] = in_A[i-A.nA-A.nM];
            }

            vector<Line> lB(B.ndim_);
            for (int i = 0;i < B.nA+B.nM;i++)
            {
                idx_B_[scB->log_to_phys[i]] = out_B[i].toInt();
                lB[scB->log_to_phys[i]] = out_B[i];
            }

            for (int i = B.nA+B.nM;i < B.ndim_;i++)
            {
                idx_B_[scB->log_to_phys[i]] = in_B[i-B.nA-B.nM].toInt();
                lB[scB->log_to_phys[i]] = in_B[i-B.nA-B.nM];
            }

            vector<Line> lC(ndim_);
            for (int i = 0;i < nA+nM;i++)
            {
                idx_C_[scC->log_to_phys[i]] = out_C[i].toInt();
                lC[scC->log_to_phys[i]] = out_C[i];
            }

            for (int i = nA+nM;i < ndim_;i++)
            {
                idx_C_[scC->log_to_phys[i]] = in_C[i-nA-nM].toInt();
                lC[scC->log_to_phys[i]] = in_C[i-nA-nM];
            }

            /*
            //if (A.logical == "aijk")
            //{
                cout << scA->log_to_phys << ' ' << scB->log_to_phys <<
                        ' ' << scC->log_to_phys << endl;
            cout << "A(" << lA << ") B(" << lB << ") C(" << lC << ") " <<
                    alpha*scA->permFactor*scB->permFactor*scC->permFactor*diagFactor <<
                    ' ' << beta[(int)(scC-cases.begin())] << endl << endl;
            //}
            */

            scC->tensor->mult(alpha*(T)(scA->permFactor*scB->permFactor*scC->permFactor*diagFactor),
            //scC->tensor.mult(alpha*diagFactor,
                              conja, *scA->tensor, idx_A_, conjb, *scB->tensor, idx_B_,
                              beta[(int)(scC-cases.begin())], idx_C_);

            beta[(int)(scC-cases.begin())] = 1.0;
        }
    }

    delete[] idx_A_;
    delete[] idx_B_;
    delete[] idx_C_;
}

template<class T>
void SpinorbitalTensor<T>::sum(const T alpha, bool conja, const SpinorbitalTensor<T>& A_, const int* idx_A,
                               const T beta_,                                             const int* idx_B)
{
    const SpinorbitalTensor<T>& A = A_.getDerived();

    int *idx_A_ = new int[A.ndim_];
    int *idx_B_ = new int[ndim_];

    //cout << "summing [";
    //for (int i = 0;i < A.ndim_;i++) cout << idx_A[i] << ' ';
    //cout << "] [";
    //for (int i = 0;i <   ndim_;i++) cout << idx_B[i] << ' ';
    //cout << "]" << endl;

    vector<T> beta(cases.size(), beta_);

    vector<Line> sAin(A.nA+A.nM);
    vector<Line> sAout(A.nE+A.nI);

    vector<Line> sBin(nA+nM);
    vector<Line> sBout(nE+nI);

    for (typename vector<SpinCase>::const_iterator scB = cases.begin();scB != cases.end();++scB)
    {
        vector<Line> sum;
        vector<Line> sum2;
        vector<Line> apo, aho, api, ahi;

        matchTypes(A.nA+A.nM, A.nE+A.nI,    A.logical, idx_A,
                     nA+  nM,   nE+  nI, scB->logical, idx_B,
                   sAout, sAin, apo, aho, api, ahi, sum);

        for (int i = 0;i < sum.size();i++)
        {
            sum2 += sum[i].toAlpha();
            sum2 += sum[i].toBeta();
        }

        for (int i = 0;i < nA+nM;i++)
        {
            sBout[i] = Line(idx_B[i], scB->logical[i].getType());
        }

        for (int i = 0;i < nE+nI;i++)
        {
            sBin[i] = Line(idx_B[nA+nM+i], scB->logical[nA+nM+i].getType());
        }

        Fragment uhfsA("T", sAout, sAin);
        Fragment uhfsB("U", sBout, sBin);

        Diagram d(Diagram::UHF);
        d += Term(Diagram::UHF)*uhfsA;
        d.sum(sum);
        d.fixorder(sum2);

        vector<Term> terms = d.getTerms();
        for (vector<Term>::iterator t = terms.begin();t != terms.end();++t)
        {
            *t *= uhfsB;

            double diagFactor = t->getFactor();

            //cout << *t << endl;

            vector<Fragment>::iterator fA, fB;
            for (vector<Fragment>::iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
            {
                if (f->getOp() == "T") fA = f;
                if (f->getOp() == "U") fB = f;
            }

            vector<Line> out_A = fA->getIndicesOut();
            vector<Line>  in_A = fA->getIndicesIn();
            vector<Line> out_B = fB->getIndicesOut();
            vector<Line>  in_B = fB->getIndicesIn();

            int nA_A = count_if(out_A.begin(), out_A.end(), isType<PARTICLE+ALPHA>());
            int nM_A = count_if(out_A.begin(), out_A.end(), isType<    HOLE+ALPHA>());
            int nE_A = count_if( in_A.begin(),  in_A.end(), isType<PARTICLE+ALPHA>());
            int nI_A = count_if( in_A.begin(),  in_A.end(), isType<    HOLE+ALPHA>());

            typename vector<SpinCase>::const_iterator scA = A.cases.end();
            for (typename vector<SpinCase>::const_iterator sc = A.cases.begin();sc != A.cases.end();++sc)
            {
                //printf("%d %d %d %d - %d %d %d %d\n", nA_A, nM_A, nE_A, nI_A,
                //       sc->nA, sc->nM, sc->nE, sc->nI);
                if (nA_A == sc->nA && nM_A == sc->nM &&
                    nE_A == sc->nE && nI_A == sc->nI) scA = sc;
            }
            if (scA == A.cases.end()) throw logic_error("no matching spin case for tensor A");

            vector<Line> lA(A.ndim_);
            for (int i = 0;i < A.nA+A.nM;i++)
            {
                idx_A_[scA->log_to_phys[i]] = out_A[i].toInt();
                lA[scA->log_to_phys[i]] = out_A[i];
            }

            for (int i = A.nA+A.nM;i < A.ndim_;i++)
            {
                idx_A_[scA->log_to_phys[i]] = in_A[i-A.nA-A.nM].toInt();
                lA[scA->log_to_phys[i]] = in_A[i-A.nA-A.nM];
            }

            vector<Line> lB(ndim_);
            for (int i = 0;i < nA+nM;i++)
            {
                idx_B_[scB->log_to_phys[i]] = out_B[i].toInt();
                lB[scB->log_to_phys[i]] = out_B[i];
            }

            for (int i = nA+nM;i < ndim_;i++)
            {
                idx_B_[scB->log_to_phys[i]] = in_B[i-nA-nM].toInt();
                lB[scB->log_to_phys[i]] = in_B[i-nA-nM];
            }

            //cout << scA->log_to_phys << ' ' << scB->log_to_phys << endl;

            //if (logical == "aijk")
            //{
            //cout << Fragment("A", lA, vector<Line>()) << ' ' <<
            //             Fragment("B", lB, vector<Line>()) << ' ' <<
            //             alpha << '*' << scA->permFactor << '*' <<
            //             scB->permFactor << '*' << diagFactor <<
            //             //alpha << '*' << diagFactor <<
            //             ' ' << beta[(int)(scB-cases.begin())] << endl;
            //}

            scB->tensor->sum(alpha*(T)(scA->permFactor*scB->permFactor*diagFactor),
            //scB->tensor.sum(alpha*diagFactor,
                              conja, *scA->tensor, idx_A_,
                              beta[(int)(scB-cases.begin())], idx_B_);

            beta[(int)(scB-cases.begin())] = 1.0;
        }
    }

    delete[] idx_A_;
    delete[] idx_B_;
}

template<class T>
void SpinorbitalTensor<T>::scale(const T alpha, const int* idx_A)
{
    int *idx_A_ = new int[ndim_];

    vector<Line> sAin(nA+nM);
    vector<Line> sAout(nE+nI);

    for (typename vector<SpinCase>::const_iterator scA = cases.begin();scA != cases.end();++scA)
    {
        for (int i = 0;i < nA+nM;i++)
        {
            sAout[i] = Line(idx_A[i], scA->logical[i].getType());
        }

        for (int i = 0;i < nE+nI;i++)
        {
            sAin[i] = Line(idx_A[nA+nM+i], scA->logical[nA+nM+i].getType());
        }

        Term t(Diagram::UHF);
        t *= Fragment("T", sAout, sAin);

        Fragment fA = t.getFragments()[0];

        vector<Line> out_A = fA.getIndicesOut();
        vector<Line>  in_A = fA.getIndicesIn();

        vector<Line> lA(ndim_);
        for (int i = 0;i < nA+nM;i++)
        {
            idx_A_[scA->log_to_phys[i]] = out_A[i].toInt();
            lA[scA->log_to_phys[i]] = out_A[i];
        }

        for (int i = nA+nM;i < ndim_;i++)
        {
            idx_A_[scA->log_to_phys[i]] = in_A[i-nA-nM].toInt();
            lA[scA->log_to_phys[i]] = in_A[i-nA-nM];
        }

        //cout << Fragment("A", lA, vector<Line>()) << ' ' <<
        //        alpha*scA->permFactor << endl;
        //          alpha << endl;

        //scA->tensor->scale(alpha*scA->permFactor, idx_A_);
        scA->tensor->scale(alpha, idx_A_);
    }

    delete[] idx_A_;
}

template<class T>
void SpinorbitalTensor<T>::div(const T alpha, bool conja, const SpinorbitalTensor<T>& A,
                                              bool conjb, const SpinorbitalTensor<T>& B, const T beta)
{
    #ifdef VALIDATE_INPUTS
    if (ndim_ != A.getDerived().ndim_ ||
        ndim_ != B.getDerived().ndim_) throw InvalidNdimError();
    for (int i = 0;i < ndim_;i++)
    {
        if (logical[i].isParticle() != A.getDerived().logical[i].isParticle() ||
            logical[i].isParticle() != B.getDerived().logical[i].isParticle())
            throw logic_error("types do not match");
    }
    #endif //VALIDATE_INPUTS

    for (int i = 0;i < cases.size();i++)
    {
        #ifdef VALIDATE_INPUTS
        for (int j = 0;j < ndim_;j++)
        {
            if (cases[i].log_to_phys[j] != A.getDerived().cases[i].log_to_phys[j] ||
                cases[i].log_to_phys[j] != B.getDerived().cases[i].log_to_phys[j])
                throw logic_error("types do not match");
        }
        #endif //VALIDATE_INPUTS

        beta*(*cases[i].tensor) += alpha*(*A.getDerived().cases[i].tensor)/
                                         (*B.getDerived().cases[i].tensor);
    }
}

template<class T>
void SpinorbitalTensor<T>::invert(const T alpha, bool conja, const SpinorbitalTensor<T>& A, const T beta)
{
    #ifdef VALIDATE_INPUTS
    if (ndim_ != A.getDerived().ndim_ ||
        ndim_ != B.getDerived().ndim_) throw InvalidNdimError();
    for (int i = 0;i < ndim_;i++)
    {
        if (logical[i].isParticle() != A.getDerived().logical[i].isParticle())
            throw logic_error("types do not match");
    }
    #endif //VALIDATE_INPUTS

    for (int i = 0;i < tensors_.size();i++)
    {
        #ifdef VALIDATE_INPUTS
        for (int j = 0;j < ndim_;j++)
        {
            if (cases[i].log_to_phys[j] != A.getDerived().logical[i].log_to_phys[j])
                throw logic_error("types do not match");
        }
        #endif //VALIDATE_INPUTS

        beta*(*cases[i].tensor) += alpha/(*A.getDerived().cases[i].tensor);
    }
}

template<class T>
T SpinorbitalTensor<T>::dot(bool conja, const SpinorbitalTensor<T>& A, const int* idx_A,
                            bool conjb,                                const int* idx_B) const
{
    DistTensor<T> one(A(0), (T)1);
    DistTensor<T> dt(A(0), (T)0);
    SpinorbitalTensor<T> sodt(",");
    sodt.addSpinCase(dt, ",", "");
    sodt.mult(1, conja,     A, idx_A,
                 conjb, *this, idx_B,
              0,                NULL);
    return scalar(dt*one);
}

INSTANTIATE_SPECIALIZATIONS(SpinorbitalTensor);
