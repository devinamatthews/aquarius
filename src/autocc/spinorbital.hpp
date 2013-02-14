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
 * ARE DISCLAIMED. IN NO EVENT SHALL EDGAR SOLOMONIK BE LIABLE FOR ANY
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

#include "tensor.hpp"
#include "dist_tensor.hpp"

#include "stl_ext/stl_ext.hpp"

#include "line.hpp"
#include "diagram.hpp"
#include "fragment.hpp"
#include "term.hpp"

namespace aquarius
{
namespace autocc
{

template<class Base>
class SpinorbitalTensor : public libtensor::Tensor< SpinorbitalTensor<Base> >
{
    friend class libtensor::Tensor< SpinorbitalTensor<Base> >;

    protected:
        using libtensor::Tensor< SpinorbitalTensor<Base> >::ndim_;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::len_;
    public:
        using libtensor::Tensor< SpinorbitalTensor<Base> >::mult;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::contract;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::weight;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::outerProduct;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::transpose;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::diagonal;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::replicate;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::trace;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::sum;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::scale;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::operator=;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::operator+=;
        using libtensor::Tensor< SpinorbitalTensor<Base> >::operator-=;

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
        : libtensor::Tensor< SpinorbitalTensor<Base> >(other.ndim_, other.len_)
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
        : libtensor::Tensor< SpinorbitalTensor<Base> >(logical.size()-1, std::vector<int>(logical.size()-1, 0).data())
        {
            int comma = logical.find(',');
            if (comma == std::string::npos) throw std::logic_error("index std::string is malformed: " + logical);

            if (logical != tolower(logical)) throw std::logic_error("spinorbital indices must be lowercase");

            this->logical = logical.substr(0,comma) + logical.substr(comma+1);

            std::vector<Line> out_ = Line::parse(logical.substr(0,comma));
            std::vector<Line> in_  = Line::parse(logical.substr(comma+1));

            std::vector<Line> out(out_);
            std::vector<Line> in(in_);

            sort(out.begin(), out.end());
            sort(in.begin(), in.end());

            //permFactor = relativeSign(out_, out) *
            //             relativeSign(in_, in);

            nA = count_if(out.begin(), out.end(), isParticle());
            nM = count_if(out.begin(), out.end(), isHole());
            nE = count_if(in.begin(), in.end(), isParticle());
            nI = count_if(in.begin(), in.end(), isHole());
        }

        ~SpinorbitalTensor()
        {
            for (typename std::vector<SpinCase>::iterator i = cases.begin();i != cases.end();++i)
            {
                if (i->isAlloced) delete i->tensor;
            }
        }

        libtensor::IndexedTensor< SpinorbitalTensor<Base> > operator*(const double factor)
        {
            libtensor::IndexedTensor< SpinorbitalTensor<Base> > it(*this, logical.c_str());
            return factor*it;
        }

        const libtensor::IndexedTensor< SpinorbitalTensor<Base> > operator*(const double factor) const
        {
            libtensor::IndexedTensor< SpinorbitalTensor<Base> > it(*this, logical.c_str());
            return factor*it;
        }

        SpinorbitalTensor<Base>& operator=(const SpinorbitalTensor<Base>& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            for (int i = 0;i < ndim_;i++)
            {
                if (Line(logical[i]).getType() != Line(other.logical[i]).getType())
                    throw IndexMismatchError();
            }
            #endif //VALIDATE_INPUTS

            (*this)[logical.c_str()] = other[logical.c_str()];

            return *this;
        }

        SpinorbitalTensor<Base>& operator+=(const SpinorbitalTensor<Base>& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            for (int i = 0;i < ndim_;i++)
            {
                if (Line(logical[i]).getType() != Line(other.logical[i]).getType())
                    throw IndexMismatchError();
            }
            #endif //VALIDATE_INPUTS

            (*this)[logical.c_str()] += other[logical.c_str()];

            return *this;
        }

        SpinorbitalTensor<Base>& operator-=(const SpinorbitalTensor<Base>& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            for (int i = 0;i < ndim_;i++)
            {
                if (Line(logical[i]).getType() != Line(other.logical[i]).getType())
                    throw IndexMismatchError();
            }
            #endif //VALIDATE_INPUTS

            (*this)[logical.c_str()] -= other[logical.c_str()];

            return *this;
        }

        SpinorbitalTensor<Base>& operator=(const libtensor::IndexedTensor< SpinorbitalTensor<Base> >& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            for (int i = 0;i < ndim_;i++)
            {
                if (Line(logical[i]).getType() != Line(other.logical[i]).getType())
                    throw IndexMismatchError();
            }
            if (strcmp(logical.c_str(), other.idx_) != 0) throw IndexMismatchError();
            #endif //VALIDATE_INPUTS

            (*this)[logical.c_str()] = other;

            return *this;
        }

        SpinorbitalTensor<Base>& operator+=(const libtensor::IndexedTensor< SpinorbitalTensor<Base> >& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            for (int i = 0;i < ndim_;i++)
            {
                if (Line(logical[i]).getType() != Line(other.logical[i]).getType())
                    throw IndexMismatchError();
            }
            if (strcmp(logical.c_str(), other.idx_) != 0) throw IndexMismatchError();
            #endif //VALIDATE_INPUTS

            (*this)[logical.c_str()] += other;

            return *this;
        }

        SpinorbitalTensor<Base>& operator-=(const libtensor::IndexedTensor< SpinorbitalTensor<Base> >& other)
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            for (int i = 0;i < ndim_;i++)
            {
                if (Line(logical[i]).getType() != Line(other.logical[i]).getType())
                    throw IndexMismatchError();
            }
            if (strcmp(logical.c_str(), other.idx_) != 0) throw IndexMismatchError();
            #endif //VALIDATE_INPUTS

            (*this)[logical.c_str()] -= other;

            return *this;
        }

        libtensor::IndexedTensorMult< SpinorbitalTensor<Base> > operator*(const SpinorbitalTensor<Base>& other) const
        {
            #ifdef VALIDATE_INPUTS
            if (ndim_ != other.ndim_) throw InvalidNdimError();
            for (int i = 0;i < ndim_;i++)
            {
                if (Line(logical[i]).getType() != Line(other.logical[i]).getType())
                    throw IndexMismatchError();
            }
            #endif //VALIDATE_INPUTS

            return (*this)[logical.c_str()]*other[logical.c_str()];
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

            std::vector<Line> out_ = Line::parse(logical.substr(0,comma));
            std::vector<Line> in_  = Line::parse(logical.substr(comma+1));
            std::vector<Line> phys = Line::parse(physical);

            std::vector<Line> out(out_);
            std::vector<Line> in(in_);

            std::sort(out.begin(), out.end());
            std::sort(in.begin(), in.end());

            //std::cout << logical << ' ' << out << ' ' << in << std::endl;

            if (std::count_if(out.begin(), out.end(), isParticle()) != nA ||
                std::count_if(out.begin(), out.end(),     isHole()) != nM ||
                std::count_if( in.begin(),  in.end(), isParticle()) != nE ||
                std::count_if( in.begin(),  in.end(),     isHole()) != nI)
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

        virtual void mult(const double alpha, const SpinorbitalTensor<Base>& A, const char* idx_A,
                                              const SpinorbitalTensor<Base>& B, const char* idx_B,
                          const double beta_,                                   const char* idx_C)
        {
            int *idx_A_ = new int[A.ndim_];
            int *idx_B_ = new int[B.ndim_];
            int *idx_C_ = new int[ndim_];

            std::vector<double> beta(cases.size(), beta_);

            std::string sA = "T(" + std::string(idx_A).substr(0, A.nA+A.nM) + "," +
                                    std::string(idx_A).substr(A.nA+A.nM) + ")";
            std::string sB = "U(" + std::string(idx_B).substr(0, B.nA+B.nM) + "," +
                                    std::string(idx_B).substr(B.nA+B.nM) + ")";
            std::string sC = "V(" + std::string(idx_C).substr(0, nA+nM) + "," +
                                    std::string(idx_C).substr(nA+nM) + ")";

            std::string ape, ahe, bpe, bhe, ssum;

            for (int i = 0;i < strlen(idx_A);i++)
            {
                if (idx_A[i] != tolower(idx_A[i])) throw std::logic_error("spinorbital indices must be lowercase");
                if (idx_A[i] >= 'a' && idx_A[i] <= 'd' &&
                    sB.find(idx_A[i]) == std::string::npos) ape += idx_A[i];
                if (idx_A[i] >= 'i' && idx_A[i] <= 'l' &&
                    sB.find(idx_A[i]) == std::string::npos) ahe += idx_A[i];
                if (idx_A[i] >= 'e' && idx_A[i] <= 'h' &&
                    sC.find(idx_A[i]) == std::string::npos) ssum += idx_A[i];
                if (idx_A[i] >= 'm' && idx_A[i] <= 'p' &&
                    sC.find(idx_A[i]) == std::string::npos) ssum += idx_A[i];
            }

            for (int i = 0;i < strlen(idx_B);i++)
            {
                if (idx_B[i] != tolower(idx_B[i])) throw std::logic_error("spinorbital indices must be lowercase");
                if (idx_B[i] >= 'a' && idx_B[i] <= 'd' &&
                    sA.find(idx_B[i]) == std::string::npos) bpe += idx_B[i];
                if (idx_B[i] >= 'i' && idx_B[i] <= 'l' &&
                    sA.find(idx_B[i]) == std::string::npos) bhe += idx_B[i];
                if (idx_B[i] >= 'e' && idx_B[i] <= 'h' &&
                    sA.find(idx_B[i]) == std::string::npos &&
                    sC.find(idx_B[i]) == std::string::npos) ssum += idx_B[i];
                if (idx_B[i] >= 'm' && idx_B[i] <= 'p' &&
                    sA.find(idx_B[i]) == std::string::npos &&
                    sC.find(idx_B[i]) == std::string::npos) ssum += idx_B[i];
            }

            std::vector<Line> sum = Line::parse(ssum);
            std::vector<Line> sum2 = sum + Line::parse(toupper(ssum));

            for (int i = 0;i < strlen(idx_C);i++)
            {
                if (idx_C[i] != tolower(idx_C[i])) throw std::logic_error("spinorbital indices must be lowercase");
            }

            std::string term;
            if (!ape.empty() && !bpe.empty()) term += "P(" + ape + "|" + bpe + ") ";
            if (!ahe.empty() && !bhe.empty()) term += "P(" + ahe + "|" + bhe + ") ";
            term += sA + ' ' + sB;

            for (typename std::vector<SpinCase>::const_iterator scC = cases.begin();scC != cases.end();++scC)
            {
                std::string uhfterm(term), uhfsC(sC);
                std::string from(logical), to(scC->logical);

                translate(from, logical+toupper(logical), std::string(idx_C)+toupper(std::string(idx_C)));
                translate(  to, logical+toupper(logical), std::string(idx_C)+toupper(std::string(idx_C)));

                translate(uhfterm, from, to);
                translate(  uhfsC, from, to);

                Diagram d(Diagram::UHF, std::vector<std::string>(1, uhfterm));
                d.sum(sum);
                d.fixorder(sum2);

                std::vector<Term> terms = d.getTerms();
                for (std::vector<Term>::iterator t1 = terms.begin();t1 != terms.end();++t1)
                {
                    for (std::vector<Term>::iterator t2 = t1+1;t2 != terms.end();++t2)
                    {
                        if (Term(*t1).fixorder(true) == Term(*t2).fixorder(true))
                        {
                            d -= *t1;
                            break;
                        }
                    }
                }

                //if (A.logical == "aijk")
                //{
                //    std::cout << std::endl<< uhfterm << ' ' << uhfsC << ' ' << ssum << std::endl;
                //    std::cout << d;
                //}

                terms = d.getTerms();
                for (std::vector<Term>::iterator t = terms.begin();t != terms.end();++t)
                {
                    *t *= Fragment(uhfsC);

                    double diagFactor = t->getFactor();

                    std::vector<Fragment>::iterator fA, fB, fC;
                    for (std::vector<Fragment>::iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
                    {
                        if (f->getOp() == "T") fA = f;
                        if (f->getOp() == "U") fB = f;
                        if (f->getOp() == "V") fC = f;
                    }

                    std::vector<Line> out_A = fA->getIndicesOut();
                    std::vector<Line>  in_A = fA->getIndicesIn();
                    std::vector<Line> out_B = fB->getIndicesOut();
                    std::vector<Line>  in_B = fB->getIndicesIn();
                    std::vector<Line> out_C = fC->getIndicesOut();
                    std::vector<Line>  in_C = fC->getIndicesIn();

                    int nA_A = count_if(out_A.begin(), out_A.end(), isType<PARTICLE+ALPHA>());
                    int nM_A = count_if(out_A.begin(), out_A.end(), isType<    HOLE+ALPHA>());
                    int nE_A = count_if( in_A.begin(),  in_A.end(), isType<PARTICLE+ALPHA>());
                    int nI_A = count_if( in_A.begin(),  in_A.end(), isType<    HOLE+ALPHA>());

                    int nA_B = count_if(out_B.begin(), out_B.end(), isType<PARTICLE+ALPHA>());
                    int nM_B = count_if(out_B.begin(), out_B.end(), isType<    HOLE+ALPHA>());
                    int nE_B = count_if( in_B.begin(),  in_B.end(), isType<PARTICLE+ALPHA>());
                    int nI_B = count_if( in_B.begin(),  in_B.end(), isType<    HOLE+ALPHA>());

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

                    std::vector<Line> lA(A.ndim_);
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

                    std::vector<Line> lB(B.ndim_);
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

                    std::vector<Line> lC(ndim_);
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
                    //std::cout << Fragment("A", lA, std::vector<Line>()) << ' ' <<
                    //        Fragment("B", lB, std::vector<Line>()) << ' ' <<
                    //        Fragment("C", lC, std::vector<Line>()) << ' ' <<
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

        virtual void contract(const double alpha, const SpinorbitalTensor<Base>& A, const char* idx_A,
                                                  const SpinorbitalTensor<Base>& B, const char* idx_B,
                              const double beta,                                    const char* idx_C)
        {
            assert(0);
        }

        virtual void weight(const double alpha, const SpinorbitalTensor<Base>& A, const char* idx_A,
                                                const SpinorbitalTensor<Base>& B, const char* idx_B,
                            const double beta,                                    const char* idx_C)
        {
            assert(0);
        }

        virtual void outerProduct(const double alpha, const SpinorbitalTensor<Base>& A, const char* idx_A,
                                                      const SpinorbitalTensor<Base>& B, const char* idx_B,
                                  const double beta,                                    const char* idx_C)
        {
            assert(0);
        }

        virtual void transpose(const double alpha, const SpinorbitalTensor<Base>& A, const char* idx_A,
                               const double beta,                                    const char* idx_B)
        {
            assert(0);
        }

        virtual void diagonal(const double alpha, const SpinorbitalTensor<Base>& A, const char* idx_A,
                              const double beta,                                    const char* idx_B)
        {
            assert(0);
        }

        virtual void replicate(const double alpha, const SpinorbitalTensor<Base>& A, const char* idx_A,
                               const double beta,                                    const char* idx_B)
        {
            assert(0);
        }

        virtual void sum(const double alpha, const SpinorbitalTensor<Base>& A, const char* idx_A,
                         const double beta_,                                   const char* idx_B)
        {
            int *idx_A_ = new int[A.ndim_];
            int *idx_B_ = new int[ndim_];

            std::vector<double> beta(cases.size(), beta_);

            std::string sA = "T(" + std::string(idx_A).substr(0, A.nA+A.nM) + "," +
                                    std::string(idx_A).substr(A.nA+A.nM) + ")";
            std::string sB = "U(" + std::string(idx_B).substr(0, nA+nM) + "," +
                                    std::string(idx_B).substr(nA+nM) + ")";

            std::string ssum;

            for (int i = 0;i < strlen(idx_A);i++)
            {
                if (idx_A[i] != tolower(idx_A[i])) throw std::logic_error("spinorbital indices must be lowercase");
                if (idx_A[i] >= 'e' && idx_A[i] <= 'h' &&
                    sB.find(idx_A[i]) == std::string::npos) ssum += idx_A[i];
                if (idx_A[i] >= 'm' && idx_A[i] <= 'p' &&
                    sB.find(idx_A[i]) == std::string::npos) ssum += idx_A[i];
            }

            std::vector<Line> sum = Line::parse(ssum);
            std::vector<Line> sum2 = sum + Line::parse(toupper(ssum));

            for (int i = 0;i < strlen(idx_B);i++)
            {
                if (idx_B[i] != tolower(idx_B[i])) throw std::logic_error("spinorbital indices must be lowercase");
            }

            std::string term = sA;

            for (typename std::vector<SpinCase>::const_iterator scB = cases.begin();scB != cases.end();++scB)
            {
                std::string uhfterm(term), uhfsB(sB);
                std::string from(logical), to(scB->logical);

                translate(from, logical+toupper(logical), std::string(idx_B)+toupper(std::string(idx_B)));
                translate(  to, logical+toupper(logical), std::string(idx_B)+toupper(std::string(idx_B)));

                translate(uhfterm, from, to);
                translate(  uhfsB, from, to);

                Diagram d(Diagram::UHF);
                d += Term(Diagram::UHF, uhfterm);
                d.sum(sum);
                d.fixorder(sum2);

                //if (logical == "aijk")
                //{
                //std::cout << std::endl<< uhfterm << ' ' << uhfsB << std::endl;
                //std::cout << d;
                //}

                std::vector<Term> terms = d.getTerms();
                for (std::vector<Term>::iterator t = terms.begin();t != terms.end();++t)
                {
                    *t *= Fragment(uhfsB);

                    double diagFactor = t->getFactor();

                    //std::cout << *t << std::endl;

                    std::vector<Fragment>::iterator fA, fB;
                    for (std::vector<Fragment>::iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
                    {
                        if (f->getOp() == "T") fA = f;
                        if (f->getOp() == "U") fB = f;
                    }

                    std::vector<Line> out_A = fA->getIndicesOut();
                    std::vector<Line>  in_A = fA->getIndicesIn();
                    std::vector<Line> out_B = fB->getIndicesOut();
                    std::vector<Line>  in_B = fB->getIndicesIn();

                    int nA_A = count_if(out_A.begin(), out_A.end(), isType<PARTICLE+ALPHA>());
                    int nM_A = count_if(out_A.begin(), out_A.end(), isType<    HOLE+ALPHA>());
                    int nE_A = count_if( in_A.begin(),  in_A.end(), isType<PARTICLE+ALPHA>());
                    int nI_A = count_if( in_A.begin(),  in_A.end(), isType<    HOLE+ALPHA>());

                    typename std::vector<SpinCase>::const_iterator scA = A.cases.end();
                    for (typename std::vector<SpinCase>::const_iterator sc = A.cases.begin();sc != A.cases.end();++sc)
                    {
                        //printf("%d %d %d %d - %d %d %d %d\n", nA_A, nM_A, nE_A, nI_A,
                        //       sc->nA, sc->nM, sc->nE, sc->nI);
                        if (nA_A == sc->nA && nM_A == sc->nM &&
                            nE_A == sc->nE && nI_A == sc->nI) scA = sc;
                    }
                    if (scA == A.cases.end()) throw std::logic_error("no matching spin case for tensor A");

                    std::vector<Line> lA(A.ndim_);
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

                    std::vector<Line> lB(ndim_);
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
                    //std::cout << Fragment("A", lA, std::vector<Line>()) << ' ' <<
                    //             Fragment("B", lB, std::vector<Line>()) << ' ' <<
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

        virtual void trace(const double alpha, const SpinorbitalTensor<Base>& A, const char* idx_A,
                           const double beta,                                    const char* idx_B)
        {
            assert(0);
        }

        virtual void scale(const double alpha, const char* idx_A)
        {
            int *idx_A_ = new int[ndim_];

            std::string sA = "T(" + std::string(idx_A).substr(0, nA+nM) + "," +
                                    std::string(idx_A).substr(nA+nM) + ")";

            for (int i = 0;i < strlen(idx_A);i++)
            {
                if (idx_A[i] != tolower(idx_A[i])) throw std::logic_error("spinorbital indices must be lowercase");
            }

            for (typename std::vector<SpinCase>::const_iterator scA = cases.begin();scA != cases.end();++scA)
            {
                std::string uhfsA(sA);
                std::string from(logical), to(scA->logical);

                translate(from, logical+toupper(logical), std::string(idx_A)+toupper(std::string(idx_A)));
                translate(  to, logical+toupper(logical), std::string(idx_A)+toupper(std::string(idx_A)));

                translate(  uhfsA, from, to);

                Term t(Diagram::UHF);
                t *= Fragment(uhfsA);

                Fragment fA = t.getFragments()[0];

                std::vector<Line> out_A = fA.getIndicesOut();
                std::vector<Line>  in_A = fA.getIndicesIn();

                std::vector<Line> lA(ndim_);
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

                //std::cout << Fragment("A", lA, std::vector<Line>()) << ' ' <<
                //        alpha*scA->permFactor << std::endl;
                //          alpha << std::endl;

                //scA->tensor->scale(alpha*scA->permFactor, idx_A_);
                scA->tensor->scale(alpha, idx_A_);
            }

            delete[] idx_A_;
        }

        virtual void mult(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                                              const SpinorbitalTensor<Base>& B, const int* idx_B,
                          const double beta,                                    const int* idx_C)
        {
            assert(0);
        }

        virtual void contract(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                                                  const SpinorbitalTensor<Base>& B, const int* idx_B,
                              const double beta,                                    const int* idx_C)
        {
            assert(0);
        }

        virtual void weight(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                                                const SpinorbitalTensor<Base>& B, const int* idx_B,
                            const double beta,                                    const int* idx_C)
        {
            assert(0);
        }

        virtual void outerProduct(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                                                      const SpinorbitalTensor<Base>& B, const int* idx_B,
                                  const double beta,                                    const int* idx_C)
        {
            assert(0);
        }

        virtual void transpose(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                               const double beta,                                    const int* idx_B)
        {
            assert(0);
        }

        virtual void diagonal(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                              const double beta,                                    const int* idx_B)
        {
            assert(0);
        }

        virtual void replicate(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                               const double beta,                                    const int* idx_B)
        {
            assert(0);
        }

        virtual void sum(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                         const double beta,                                    const int* idx_B)
        {
            assert(0);
        }

        virtual void trace(const double alpha, const SpinorbitalTensor<Base>& A, const int* idx_A,
                           const double beta,                                    const int* idx_B)
        {
            assert(0);
        }

        virtual void scale(const double alpha, const int* idx_A)
        {
            assert(0);
        }
};

}
}

namespace libtensor
{
    template<>
    double scalar(const IndexedTensor< aquarius::autocc::SpinorbitalTensor<DistTensor> >& other);

    template<>
    double scalar(const IndexedTensorMult< aquarius::autocc::SpinorbitalTensor<DistTensor> >& other);
}

#endif
