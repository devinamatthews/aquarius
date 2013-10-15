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

#ifndef _AQUARIUS_INTEGRALS_ELEMENT_HPP_
#define _AQUARIUS_INTEGRALS_ELEMENT_HPP_

#include <cstddef>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <strings.h>

#include "memory/memory.h"
#include "util/math_ext.h"
#include "util/blas.h"
#include "util/stl_ext.hpp"
#include "symmetry/symmetry.hpp"

namespace aquarius
{
namespace input
{
class Molecule;
}

namespace integrals
{

class Element
{
    protected:
        double charge;
        int Z;
        int A;
        int S;
        std::string symbol;
        std::string name;
        double mass;

        Element(double charge, int Z, int A, int S, const char* symbol,
                const char* name, double mass)
        : charge(charge), Z(Z), A(A), S(S), symbol(symbol), name(name), mass(mass) {}

    public:
        static Element getElement(const std::string& symbol, int A);

        static Element getElement(const std::string& symbol);

        void setCharge(double charge) { this->charge = charge; }

        double getCharge() const { return charge; }

        int getAtomicNumber() const { return Z; }

        int getNucleonNumber() const { return A; }

        void setMass(double mass) { this->mass = mass; }

        double getMass() const { return mass; }

        int getSpin() const { return S; }

        const std::string& getName() const { return name; }

        const std::string& getSymbol() const { return symbol; }

        bool operator==(const Element& other) const
        {
            return Z == other.Z && A == other.A;
        }
};

}
}

#endif
