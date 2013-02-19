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

#ifndef _AQUARIUS_AUTOCC_SPINORBITAL_HPP_
#define _AQUARIUS_AUTOCC_SPINORBITAL_HPP_

#include <cassert>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cctype>
#include <iostream>
#include <cstdio>

#include "tensor/tensor.hpp"
#include "tensor/dist_tensor.hpp"
#include "stl_ext/stl_ext.hpp"
#include "autocc/autocc.hpp"

namespace aquarius
{
namespace tensor
{

template<class Base>
class SpinorbitalTensor : public IndexableTensor< SpinorbitalTensor<Base> >
{
    INHERIT_FROM_INDEXABLE_TENSOR(SpinorbitalTensor<Base>)

    protected:
        struct SpinCase
        {
            Base* tensor;
            std::string logical;
            int nA, nE, nM, nI;
            std::vector<int> log_to_phys;
            bool isAlloced;
            double permFactor;

            SpinCase(Base& tensor, bool isAlloced = false) : tensor(&tensor), isAlloced(isAlloced) {}
            SpinCase(Base* tensor, bool isAlloced = true) : tensor(tensor), isAlloced(isAlloced) {}
        };

        std::string logical;
        int nA, nE, nM, nI;
        //double permFactor;
        std::vector<SpinCase> cases;

    public:
        SpinorbitalTensor(const SpinorbitalTensor<Base>& other)
        : IndexableTensor< SpinorbitalTensor<Base> >(other.ndim_)
        {
            logical = other.logical;
            nA = other.nA;
            nM = other.nM;
            nE = other.nE;
            nI = other.nI;

            for (typename std::vector<SpinCase>::const_iterator sc = other.cases.begin();sc != other.cases.end();++sc)
            {
                SpinCase newsc(new Base(*(sc->tensor)));
                newsc.logical = sc->logical;
                newsc.log_to_phys = sc->log_to_phys;
                newsc.nA = sc->nA;
                newsc.nM = sc->nM;
                newsc.nE = sc->nE;
                newsc.nI = sc->nI;
                newsc.permFactor = sc->permFactor;
                cases.push_back(newsc);
            }
        }

        SpinorbitalTensor(std::string logical)
        : IndexableTensor< SpinorbitalTensor<Base> >(logical.size()-1)
        {
            int comma = logical.find(',');
            if (comma == std::string::npos) throw std::logic_error("index std::string is malformed: " + logical);

            if (logical != tolower(logical)) throw std::logic_error("spinorbital indices must be lowercase");

            this->logical = logical.substr(0,comma) + logical.substr(comma+1);

            std::vector<autocc::Line> out_ = autocc::Line::parse(logical.substr(0,comma));
            std::vector<autocc::Line> in_  = autocc::Line::parse(logical.substr(comma+1));

            std::vector<autocc::Line> out(out_);
            std::vector<autocc::Line> in(in_);

            sort(out.begin(), out.end());
            sort(in.begin(), in.end());

            //permFactor = relativeSign(out_, out) *
            //             relativeSign(in_, in);

            nA = count_if(out.begin(), out.end(), autocc::isParticle());
            nM = count_if(out.begin(), out.end(), autocc::isHole());
            nE = count_if(in.begin(), in.end(), autocc::isParticle());
            nI = count_if(in.begin(), in.end(), autocc::isHole());
        }

        virtual ~SpinorbitalTensor()
        {
            for (typename std::vector<SpinCase>::iterator i = cases.begin();i != cases.end();++i)
            {
                if (i->isAlloced) delete i->tensor;
            }
        }

        void addSpinCase(Base* tensor, std::string logical, std::string physical, double factor = 1.0, bool isAlloced = true)
        {
            addSpinCase(*tensor, logical, physical, factor, isAlloced);
        }

        void addSpinCase(Base& tensor, std::string logical, std::string physical, double factor = 1.0, bool isAlloced = false)
        {
            SpinCase sc(tensor, isAlloced);

            int comma = logical.find(',');
            if (comma == std::string::npos) throw std::logic_error("index std::string is malformed: " + logical);

            sc.logical = logical.substr(0,comma) + logical.substr(comma+1);

            if (this->logical.size() != sc.logical.size()) throw std::logic_error("wrong number of indices");

            for (std::string::const_iterator l1 = sc.logical.begin(), l2 = this->logical.begin();l1 != sc.logical.end();++l1, ++l2)
            {
                if (tolower(*l1) != tolower(*l2)) throw std::logic_error("UHF indices do not match spinorbital");
            }

            std::vector<autocc::Line> out_ = autocc::Line::parse(logical.substr(0,comma));
            std::vector<autocc::Line> in_  = autocc::Line::parse(logical.substr(comma+1));
            std::vector<autocc::Line> phys = autocc::Line::parse(physical);

            std::vector<autocc::Line> out(out_);
            std::vector<autocc::Line> in(in_);

            std::sort(out.begin(), out.end());
            std::sort(in.begin(), in.end());

            //std::cout << logical << ' ' << out << ' ' << in << std::endl;

            if (std::count_if(out.begin(), out.end(), autocc::isParticle()) != nA ||
                std::count_if(out.begin(), out.end(),     autocc::isHole()) != nM ||
                std::count_if( in.begin(),  in.end(), autocc::isParticle()) != nE ||
                std::count_if( in.begin(),  in.end(),     autocc::isHole()) != nI)
            {
                throw std::logic_error("tensor has the wrong shape");
            }

            sc.permFactor = relativeSign(out_, out) *
                            relativeSign(in_, in) * factor;

            //std::cout << (out_ + in_) << " = " << sc.permFactor << '*' <<
            //             (out  + in ) << std::endl << std::endl;

            if (out.size()+in.size() != phys.size())
                throw std::logic_error("logical and physical dimensions do not match");

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

                sc.log_to_phys.push_back(std::find(phys.begin(), phys.end(), out[i])-phys.begin());
                if (sc.log_to_phys.back() == phys.size())
                    throw std::logic_error("logical index not found in physical map: " + logical[i]);
            }

            for (int i = 0;i < in.size();i++)
            {
                if (in[i].isAlpha())
                {
                    if (in[i].isParticle()) sc.nE++;
                    else                    sc.nI++;
                }

                sc.log_to_phys.push_back(std::find(phys.begin(), phys.end(), in[i])-phys.begin());
                if (sc.log_to_phys.back() == phys.size())
                    throw std::logic_error("logical index not found in physical map: " + logical[i+out.size()]);
            }

            if (abs(sc.nA+sc.nM-sc.nE-sc.nI) > abs((int)(out.size()-in.size())))
                throw std::logic_error("spin case is not valid");

            cases.push_back(sc);

            //std::cout << "Adding spin case " << sc.nA << ' ' << sc.nM << ' ' << sc.nE << ' ' << sc.nI << std::endl;
        }

        int getNumSpinCases()
        {
            return cases.size();
        }

        Base& getSpinCase(int sc)
        {
            return *(cases[sc].tensor);
        }

        const Base& getSpinCase(int sc) const
        {
            return *(cases[sc].tensor);
        }

        virtual SpinorbitalTensor<Base>& operator=(const double val)
        {
            for (typename std::vector<SpinCase>::iterator i = cases.begin();i != cases.end();++i)
            {
                i->tensor->operator=(val);
            }

            return *this;
        }

    protected:
        static void matchTypes(const int nin_A, const int nout_A, const std::string& log_A, const int* idx_A,
                               const int nin_B, const int nout_B, const std::string& log_B, const int* idx_B,
                               std::vector<autocc::Line>& out_A, std::vector<autocc::Line>& in_A,
                               std::vector<autocc::Line>& pout_A, std::vector<autocc::Line>& hout_A,
                               std::vector<autocc::Line>& pin_A, std::vector<autocc::Line>& hin_A,
                               std::vector<autocc::Line>& sum_A)
        {
            for (int i = 0;i < nout_A;i++)
            {
                int j;
                for (j = 0;j < nout_B+nin_B;j++)
                {
                    if (idx_A[i] == idx_B[j])
                    {
                        bool aisp = autocc::Line(log_A[i]).isParticle();
                        bool bisp = autocc::Line(log_B[j]).isParticle();
                        if (aisp != bisp)
                            throw std::logic_error("types do not match");
                        out_A[i] = autocc::Line(idx_A[i], autocc::Line(log_B[j]).getType());
                        if (autocc::Line(log_B[j]).isParticle())
                        {
                            pin_A += autocc::Line(idx_A[i], autocc::Line(log_B[j]).getType());
                        }
                        else
                        {
                            hin_A += autocc::Line(idx_A[i], autocc::Line(log_B[j]).getType());
                        }
                        break;
                    }
                }
                if (j == nout_B+nin_B)
                {
                    out_A[i] = autocc::Line(idx_A[i], autocc::Line(log_A[i]).getType());
                    sum_A += autocc::Line(idx_A[i], autocc::Line(log_A[i]).getType());
                }
            }

            for (int i = nout_A;i < nout_A+nin_A;i++)
            {
                int j;
                for (j = 0;j < nout_B+nin_B;j++)
                {
                    if (idx_A[i] == idx_B[j])
                    {
                        bool aisp = autocc::Line(log_A[i]).isParticle();
                        bool bisp = autocc::Line(log_B[j]).isParticle();
                        if (aisp != bisp)
                            throw std::logic_error("types do not match");
                        in_A[i-nout_A] = autocc::Line(idx_A[i], autocc::Line(log_B[j]).getType());
                        if (autocc::Line(log_B[j]).isParticle())
                        {
                            pout_A += autocc::Line(idx_A[i], autocc::Line(log_B[j]).getType());
                        }
                        else
                        {
                            hout_A += autocc::Line(idx_A[i], autocc::Line(log_B[j]).getType());
                        }
                        break;
                    }
                }
                if (j == nout_B+nin_B)
                {
                    in_A[i-nout_A] = autocc::Line(idx_A[i], autocc::Line(log_A[i]).getType());
                    sum_A += autocc::Line(idx_A[i], autocc::Line(log_A[i]).getType());
                }
            }
        }

        virtual void mult(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                                              const SpinorbitalTensor<Base>& B, const int* idx_B,
                          const double beta_,                                   const int* idx_C)
        {
            int *idx_A_ = new int[A.ndim_];
            int *idx_B_ = new int[B.ndim_];
            int *idx_C_ = new int[ndim_];

            //std::cout << "summing " << A.logical << "[";
            //for (int i = 0;i < A.ndim_;i++) std::cout << idx_A[i] << ' ';
            //std::cout << "] " << B.logical << "[";
            //for (int i = 0;i < B.ndim_;i++) std::cout << idx_B[i] << ' ';
            //std::cout << "] " << logical << "[";
            //for (int i = 0;i <   ndim_;i++) std::cout << idx_C[i] << ' ';
            //std::cout << "]" << std::endl;

            std::vector<double> beta(cases.size(), beta_);

            std::vector<autocc::Line> sAin(A.nA+A.nM);
            std::vector<autocc::Line> sAout(A.nE+A.nI);

            std::vector<autocc::Line> sBin(B.nA+B.nM);
            std::vector<autocc::Line> sBout(B.nE+B.nI);

            std::vector<autocc::Line> sCin(nA+nM);
            std::vector<autocc::Line> sCout(nE+nI);

            for (typename std::vector<SpinCase>::const_iterator scC = cases.begin();scC != cases.end();++scC)
            {
                std::vector<autocc::Line> sum;
                std::vector<autocc::Line> sum2;
                std::vector<autocc::Line> apo, bpo;
                std::vector<autocc::Line> aho, bho;
                std::vector<autocc::Line> api, bpi;
                std::vector<autocc::Line> ahi, bhi;

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
                                bool aisp = autocc::Line(A.logical[i]).isParticle();
                                bool bisp = autocc::Line(B.logical[j]).isParticle();
                                if (aisp != bisp)
                                        throw std::logic_error("types do not match");
                                sAout[i] = sBout[j];
                            }
                        }
                        for (;j < B.ndim_;j++)
                        {
                            if (idx_A[i] == idx_B[j])
                            {
                                bool aisp = autocc::Line(A.logical[i]).isParticle();
                                bool bisp = autocc::Line(B.logical[j]).isParticle();
                                if (aisp != bisp)
                                        throw std::logic_error("types do not match");
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
                                bool aisp = autocc::Line(A.logical[i]).isParticle();
                                bool bisp = autocc::Line(B.logical[j]).isParticle();
                                if (aisp != bisp)
                                        throw std::logic_error("types do not match");
                                sAin[i-A.nA-A.nM] = sBout[j];
                            }
                        }
                        for (;j < B.ndim_;j++)
                        {
                            if (idx_A[i] == idx_B[j])
                            {
                                bool aisp = autocc::Line(A.logical[i]).isParticle();
                                bool bisp = autocc::Line(B.logical[j]).isParticle();
                                if (aisp != bisp)
                                        throw std::logic_error("types do not match");
                                sAin[i-A.nA-A.nM] = sBin[j-B.nA-B.nM];
                            }
                        }
                    }
                }

                for (int i = 0;i < nA+nM;i++)
                {
                    sCout[i] = autocc::Line(idx_C[i], autocc::Line(scC->logical[i]).getType());
                }

                for (int i = 0;i < nE+nI;i++)
                {
                    sCin[i] = autocc::Line(idx_C[nA+nM+i], autocc::Line(scC->logical[nA+nM+i]).getType());
                }

                autocc::Fragment uhfsA("T", sAout, sAin);
                autocc::Fragment uhfsB("U", sBout, sBin);
                autocc::Fragment uhfsC("V", sCout, sCin);

                std::vector<autocc::Line> isect;

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

                //std::cout << apo << "," << bpo << "|" << aho << "," << bho << "|" <<
                //             api << "," << bpi << "|" << ahi << "," << bhi << std::endl;

                std::vector< std::vector<autocc::Line> > assym(2);

                autocc::Diagram d(autocc::Diagram::UHF);
                d += autocc::Term(autocc::Diagram::UHF)*uhfsA*uhfsB;
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
                 * Remove autocc::Terms which are antisymmetrizations of same-spin groups
                 */
                std::vector<autocc::Term> terms = d.getTerms();
                for (std::vector<autocc::Term>::iterator t1 = terms.begin();t1 != terms.end();++t1)
                {
                    for (std::vector<autocc::Term>::iterator t2 = t1+1;t2 != terms.end();++t2)
                    {
                        if (autocc::Term(*t1).fixorder(true) == autocc::Term(*t2).fixorder(true))
                        {
                            d -= *t1;
                            break;
                        }
                    }
                }

                terms = d.getTerms();
                for (std::vector<autocc::Term>::iterator t = terms.begin();t != terms.end();++t)
                {
                    *t *= uhfsC;

                    //std::cout << *t << std::endl;

                    double diagFactor = t->getFactor();

                    std::vector<autocc::Fragment>::iterator fA, fB, fC;
                    for (std::vector<autocc::Fragment>::iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
                    {
                        if (f->getOp() == "T") fA = f;
                        if (f->getOp() == "U") fB = f;
                        if (f->getOp() == "V") fC = f;
                    }

                    std::vector<autocc::Line> out_A = fA->getIndicesOut();
                    std::vector<autocc::Line>  in_A = fA->getIndicesIn();
                    std::vector<autocc::Line> out_B = fB->getIndicesOut();
                    std::vector<autocc::Line>  in_B = fB->getIndicesIn();
                    std::vector<autocc::Line> out_C = fC->getIndicesOut();
                    std::vector<autocc::Line>  in_C = fC->getIndicesIn();

                    int nA_A = count_if(out_A.begin(), out_A.end(), autocc::isType<autocc::PARTICLE+autocc::ALPHA>());
                    int nM_A = count_if(out_A.begin(), out_A.end(), autocc::isType<    autocc::HOLE+autocc::ALPHA>());
                    int nE_A = count_if( in_A.begin(),  in_A.end(), autocc::isType<autocc::PARTICLE+autocc::ALPHA>());
                    int nI_A = count_if( in_A.begin(),  in_A.end(), autocc::isType<    autocc::HOLE+autocc::ALPHA>());

                    int nA_B = count_if(out_B.begin(), out_B.end(), autocc::isType<autocc::PARTICLE+autocc::ALPHA>());
                    int nM_B = count_if(out_B.begin(), out_B.end(), autocc::isType<    autocc::HOLE+autocc::ALPHA>());
                    int nE_B = count_if( in_B.begin(),  in_B.end(), autocc::isType<autocc::PARTICLE+autocc::ALPHA>());
                    int nI_B = count_if( in_B.begin(),  in_B.end(), autocc::isType<    autocc::HOLE+autocc::ALPHA>());

                    typename std::vector<SpinCase>::const_iterator scA = A.cases.end();
                    for (typename std::vector<SpinCase>::const_iterator sc = A.cases.begin();sc != A.cases.end();++sc)
                    {
                        if (nA_A == sc->nA && nM_A == sc->nM &&
                            nE_A == sc->nE && nI_A == sc->nI) scA = sc;
                    }
                    if (scA == A.cases.end()) throw std::logic_error("no matching spin case for tensor A");

                    typename std::vector<SpinCase>::const_iterator scB = B.cases.end();
                    for (typename std::vector<SpinCase>::const_iterator sc = B.cases.begin();sc != B.cases.end();++sc)
                    {
                        if (nA_B == sc->nA && nM_B == sc->nM &&
                            nE_B == sc->nE && nI_B == sc->nI) scB = sc;
                    }
                    if (scB == B.cases.end()) throw std::logic_error("no matching spin case for tensor B");

                    std::vector<autocc::Line> lA(A.ndim_);
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

                    std::vector<autocc::Line> lB(B.ndim_);
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

                    std::vector<autocc::Line> lC(ndim_);
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

                    //if (A.logical == "aijk")
                    //{
                    //    std::cout << scA->log_to_phys << ' ' << scB->log_to_phys <<
                    //            ' ' << scC->log_to_phys << std::endl;
                    //std::cout << autocc::Fragment("A", lA, std::vector<autocc::Line>()) << ' ' <<
                    //        autocc::Fragment("B", lB, std::vector<autocc::Line>()) << ' ' <<
                    //        autocc::Fragment("C", lC, std::vector<autocc::Line>()) << ' ' <<
                    //        alpha*diagFactor <<
                    //        ' ' << beta[(int)(scC-cases.begin())] << std::endl;
                    //}

                    scC->tensor->mult(alpha*scA->permFactor*scB->permFactor*scC->permFactor*diagFactor,
                    //scC->tensor->mult(alpha*diagFactor,
                                      *(scA->tensor), idx_A_, *(scB->tensor), idx_B_,
                                      beta[(int)(scC-cases.begin())], idx_C_);

                    beta[(int)(scC-cases.begin())] = 1.0;
                }
            }

            delete[] idx_A_;
            delete[] idx_B_;
            delete[] idx_C_;
        }

        virtual void sum(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                         const double beta_,                                   const int* idx_B)
        {
            int *idx_A_ = new int[A.ndim_];
            int *idx_B_ = new int[ndim_];

            //std::cout << "summing [";
            //for (int i = 0;i < A.ndim_;i++) std::cout << idx_A[i] << ' ';
            //std::cout << "] [";
            //for (int i = 0;i <   ndim_;i++) std::cout << idx_B[i] << ' ';
            //std::cout << "]" << std::endl;

            std::vector<double> beta(cases.size(), beta_);

            std::vector<autocc::Line> sAin(A.nA+A.nM);
            std::vector<autocc::Line> sAout(A.nE+A.nI);

            std::vector<autocc::Line> sBin(nA+nM);
            std::vector<autocc::Line> sBout(nE+nI);

            for (typename std::vector<SpinCase>::const_iterator scB = cases.begin();scB != cases.end();++scB)
            {
                std::vector<autocc::Line> sum;
                std::vector<autocc::Line> sum2;
                std::vector<autocc::Line> apo, aho, api, ahi;

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
                    sBout[i] = autocc::Line(idx_B[i], autocc::Line(scB->logical[i]).getType());
                }

                for (int i = 0;i < nE+nI;i++)
                {
                    sBin[i] = autocc::Line(idx_B[nA+nM+i], autocc::Line(scB->logical[nA+nM+i]).getType());
                }

                autocc::Fragment uhfsA("T", sAout, sAin);
                autocc::Fragment uhfsB("U", sBout, sBin);

                autocc::Diagram d(autocc::Diagram::UHF);
                d += autocc::Term(autocc::Diagram::UHF)*uhfsA;
                d.sum(sum);
                d.fixorder(sum2);

                std::vector<autocc::Term> terms = d.getTerms();
                for (std::vector<autocc::Term>::iterator t = terms.begin();t != terms.end();++t)
                {
                    *t *= uhfsB;

                    double diagFactor = t->getFactor();

                    //std::cout << *t << std::endl;

                    std::vector<autocc::Fragment>::iterator fA, fB;
                    for (std::vector<autocc::Fragment>::iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
                    {
                        if (f->getOp() == "T") fA = f;
                        if (f->getOp() == "U") fB = f;
                    }

                    std::vector<autocc::Line> out_A = fA->getIndicesOut();
                    std::vector<autocc::Line>  in_A = fA->getIndicesIn();
                    std::vector<autocc::Line> out_B = fB->getIndicesOut();
                    std::vector<autocc::Line>  in_B = fB->getIndicesIn();

                    int nA_A = count_if(out_A.begin(), out_A.end(), autocc::isType<autocc::PARTICLE+autocc::ALPHA>());
                    int nM_A = count_if(out_A.begin(), out_A.end(), autocc::isType<    autocc::HOLE+autocc::ALPHA>());
                    int nE_A = count_if( in_A.begin(),  in_A.end(), autocc::isType<autocc::PARTICLE+autocc::ALPHA>());
                    int nI_A = count_if( in_A.begin(),  in_A.end(), autocc::isType<    autocc::HOLE+autocc::ALPHA>());

                    typename std::vector<SpinCase>::const_iterator scA = A.cases.end();
                    for (typename std::vector<SpinCase>::const_iterator sc = A.cases.begin();sc != A.cases.end();++sc)
                    {
                        //printf("%d %d %d %d - %d %d %d %d\n", nA_A, nM_A, nE_A, nI_A,
                        //       sc->nA, sc->nM, sc->nE, sc->nI);
                        if (nA_A == sc->nA && nM_A == sc->nM &&
                            nE_A == sc->nE && nI_A == sc->nI) scA = sc;
                    }
                    if (scA == A.cases.end()) throw std::logic_error("no matching spin case for tensor A");

                    std::vector<autocc::Line> lA(A.ndim_);
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

                    std::vector<autocc::Line> lB(ndim_);
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

                    //std::cout << scA->log_to_phys << ' ' << scB->log_to_phys << std::endl;

                    //if (logical == "aijk")
                    //{
                    //std::cout << autocc::Fragment("A", lA, std::vector<autocc::Line>()) << ' ' <<
                    //             autocc::Fragment("B", lB, std::vector<autocc::Line>()) << ' ' <<
                    //             alpha << '*' << scA->permFactor << '*' <<
                    //             scB->permFactor << '*' << diagFactor <<
                    //             //alpha << '*' << diagFactor <<
                    //             ' ' << beta[(int)(scB-cases.begin())] << std::endl;
                    //}

                    scB->tensor->sum(alpha*scA->permFactor*scB->permFactor*diagFactor,
                    //scB->tensor->sum(alpha*diagFactor,
                                      *(scA->tensor), idx_A_,
                                      beta[(int)(scB-cases.begin())], idx_B_);

                    beta[(int)(scB-cases.begin())] = 1.0;
                }
            }

            delete[] idx_A_;
            delete[] idx_B_;
        }

        virtual void scale(const double alpha, const int* idx_A)
        {
            int *idx_A_ = new int[ndim_];

            std::vector<autocc::Line> sAin(nA+nM);
            std::vector<autocc::Line> sAout(nE+nI);

            for (typename std::vector<SpinCase>::const_iterator scA = cases.begin();scA != cases.end();++scA)
            {
                for (int i = 0;i < nA+nM;i++)
                {
                    sAout[i] = autocc::Line(idx_A[i], autocc::Line(scA->logical[i]).getType());
                }

                for (int i = 0;i < nE+nI;i++)
                {
                    sAin[i] = autocc::Line(idx_A[nA+nM+i], autocc::Line(scA->logical[nA+nM+i]).getType());
                }

                autocc::Term t(autocc::Diagram::UHF);
                t *= autocc::Fragment("T", sAout, sAin);

                autocc::Fragment fA = t.getFragments()[0];

                std::vector<autocc::Line> out_A = fA.getIndicesOut();
                std::vector<autocc::Line>  in_A = fA.getIndicesIn();

                std::vector<autocc::Line> lA(ndim_);
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

                //std::cout << autocc::Fragment("A", lA, std::vector<autocc::Line>()) << ' ' <<
                //        alpha*scA->permFactor << std::endl;
                //          alpha << std::endl;

                //scA->tensor->scale(alpha*scA->permFactor, idx_A_);
                scA->tensor->scale(alpha, idx_A_);
            }

            delete[] idx_A_;
        }
};

template <typename T>
double scalar(const IndexedTensor< SpinorbitalTensor< DistTensor<T> > >& other)
{
    DistTensor<T> dt(other.tensor_.getSpinCase(0).getCTF());
    SpinorbitalTensor< DistTensor<T> > sodt(",");
    sodt.addSpinCase(dt, ",", "");
    size_t n;
    double ret, *val;
    sodt[""] = other;
    dt.getAllData(n, val);
    assert(n==1);
    ret = val[0];
    free(val);
    return ret;
}

template <typename  T>
double scalar(const IndexedTensorMult< SpinorbitalTensor< DistTensor<T> > >& other)
{
    DistTensor<T> dt(other.A_.tensor_.getSpinCase(0).ctf);
    SpinorbitalTensor< DistTensor<T> > sodt(",");
    sodt.addSpinCase(dt, ",", "");
    int64_t n;
    double ret, *val;
    sodt[""] = other;
    dt.getAllData(n, val);
    assert(n==1);
    ret = val[0];
    free(val);
    return ret;
}

}
}

#endif
