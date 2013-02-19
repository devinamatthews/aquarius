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

#ifndef _AQUARIUS_SLIDE_HPP_
#define _AQUARIUS_SLIDE_HPP_

#include <cstddef>
#include <string>
#include <vector>
#include <stdexcept>

#include "slide.h"

namespace aquarius
{
namespace slide
{

class SLIDE;
class Context;
class Center;
class Element;
class Shell;

class SLIDEError;

void swap(Context& first, Context& second);
void swap(Center& first, Center& second);
void swap(Element& first, Element& second);
void swap(Shell& first, Shell& second);

class SLIDE
{
    private:
        static SLIDE* slide;

        SLIDE() {}

        ~SLIDE() {}

    public:
        static SLIDE& init();

        static void finish();

        static void setCartesianOrdering(const int L, const int* ordering);

        static const int* getCartesianOrdering(const int L);

        static void setSphericalOrdering(const int L, const int* ordering);

        static const int* getSphericalOrdering(const int L);

        static void setGroup(const int group);

        static int getGroup();

        static int getGroupOrder();

        static void setOrdering(const int ordering);

        static int getOrdering();

        static void setGeometryTolerance(const double tolerance);

        static double getGeometryTolerance();

        static const char* getGroupLabel(const int group);

        static const char* getLLabel(const int L);

        static const char* getOpLabel(const int op);

        static const char* getIrrepLabel(const int irrep);

        static const char* getFuncLabel(const int L, const bool spherical, const int func);
};

class Context
{
    friend void swap(Context& first, Context& second)
    {
        using std::swap;
        swap(first.context, second.context);
    }

    protected:
        context_t* context;
        const Shell *a, *b, *c, *d;

    public:
        Context();

        Context(const Context& other);

        ~Context();

        Context& operator=(Context other);

        size_t process2eInts(const size_t nprocess, double* integrals, idx4_t* indices, const double cutoff);

        size_t process1eInts(const size_t nprocess, double* integrals, idx2_t* indices, const double cutoff);

        void setNumProcessed(size_t num_processed);

        size_t getNumProcessed() const;

        size_t getNumIntegrals() const;

        double* getIntegrals();

        const double* getIntegrals() const;

        void calcERI(const double alpha, const double beta,
                     const Shell& a, const Shell& b, const Shell& c, const Shell& d);

        void calcOVI(const double alpha, const double beta,
                     const Shell& a, const Shell& b);

        void calcKEI(const double alpha, const double beta,
                     const Shell& a, const Shell& b);

        void calcNAI(const double alpha, const double beta,
                     const Shell& a, const Shell& b, const Center centers[], const int ncenters);

        const Shell& getA() const;
        const Shell& getB() const;
        const Shell& getC() const;
        const Shell& getD() const;


};

class Element
{
    friend class Center;
    friend void swap(Element& first, Element& second)
    {
        using std::swap;
        swap(first.element, second.element);
    }

    protected:
        element_t* element;

    private:
        Element(element_t* element) : element(element) {}

    public:
        Element(const Element& other);

        ~Element();

        Element& operator=(Element other);

        static Element getElement(const std::string& symbol);

        void setCharge(const double charge);

        double getCharge() const;

        int getAtomicNumber() const;

        int getNucleonNumber() const;

        double getMass() const;

        int getSpin() const;

        const char* getName() const;

        const char* getSymbol() const;

        bool operator==(const Element& other) const;
};

class Center
{
    friend class Context;
    friend class Shell;
    friend void swap(Center& first, Center& second)
    {
        using std::swap;
        swap(first.center, second.center);
        swap(first.element, second.element);
    }

    protected:
        center_t* center;
        Element element;

    public:
        Center(const Center& other);

        Center(const double pos[3], const Element& element);

        ~Center();

        Center& operator=(Center other);

        static double nuclearRepulsion(const std::vector<Center>& centers);

        const int* getStabilizer() const;

        int getDegeneracy() const;

        const Element& getElement() const;

        const double* getCenter(const int degen) const;

        const double* getCenterAfterOp(const int op) const;
};

class Shell
{
    friend class Context;
    friend void swap(Shell& first, Shell& second)
    {
        using std::swap;
        swap(first.center, second.center);
        swap(first.shell, second.shell);
    }

    protected:
        shell_t* shell;
        Center center;

    public:
        Shell(const Shell& other);

        Shell(const Shell& other, const bool spherical, const bool contaminants);

        Shell(const Center& pos, const int L, const int nprim, const int ncontr, const bool spherical, const bool contaminants,
              const double* exponents, const double* coefficients, const int idx[8] = (int[]){0,0,0,0,0,0,0,0});

        ~Shell();

        Shell& operator=(Shell other);

        void funcCenters(const int func, const int irrep, int* proj) const;

        void aoToSo(const int primitive_ordering, double* aoso, const int ld) const;

        int getL() const;

        int getNPrim() const;

        int getNContr() const;

        int getNFunc() const;

        int getNFuncInIrrep(const int irrep) const;

        const int* getNFuncInEachIrrep() const;

        const int* getIdx() const;

        int* getIdx();

        bool getSpherical() const;

        bool getContaminants() const;

        int getIrrepOfFunc(const int func, const int degen) const;

        const int* getIrrepsOfFunc(const int func) const;

        const double* getExponents() const;

        const double* getCoefficients() const;

        const Center& getCenter() const;
};

class SLIDEError : public std::runtime_error
{
    public:
        SLIDEError(const std::string& what_string) : runtime_error(what_string) {}
};

}
}

#endif
