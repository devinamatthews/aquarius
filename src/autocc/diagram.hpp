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

#ifndef _AQUARIUS_AUTOCC_DIAGRAM_HPP_
#define _AQUARIUS_AUTOCC_DIAGRAM_HPP_

#include <ostream>
#include <vector>
#include <string>

namespace aquarius
{
namespace autocc
{

class Diagram;
class Term;
class Fragment;
class Fraction;
class Line;
class Manifold;

std::ostream& operator<<(std::ostream& out, const Diagram& d);

class Diagram
{
    friend class Term;
    friend class Fragment;

    public:
        enum Type
        {
            SPINORBITAL,
            UHF,
            RHF,
            SKELETON
        };

    protected:
        Type type;
        std::vector<Term> terms;

    public:
        Diagram(const Type type) : type(type) {}

        Diagram(const Type type, const std::vector<Term>& terms);

        Diagram(const Type type, const std::vector<std::string>& strings);

        Diagram& convert(const Type type);

        Diagram& fixorder(bool all = false);

        Diagram& fixorder(const std::vector<Line>& which);

        Diagram& fixexternal();

        Diagram& negate();

        Diagram operator-() const;

        Diagram& operator*=(const Diagram& other);

        Diagram& operator*=(const Term& other);

        Diagram& operator*=(const Fraction& other);

        Diagram& operator*=(const int other);

        Diagram operator*(const Diagram& other) const;

        Diagram operator*(const Term& other) const;

        Diagram operator*(const Fraction& other) const;

        Diagram operator*(const int other) const;

        Diagram& operator+=(const Diagram& other);

        Diagram& operator+=(const Term& other);

        Diagram& operator-=(const Diagram& other);

        Diagram& operator-=(const Term& other);

        Diagram operator+(const Diagram& other) const;

        Diagram operator-(const Diagram& other) const;

        bool operator==(const Diagram& other) const;

        const Term& operator[](const int i) const;

        int numTerms() const;

        const std::vector<Term>& getTerms() const;

        Diagram& sumall();

        Diagram& sum(const std::vector<Line>& which);

    protected:
        Diagram& suminternal();

        Diagram& translate(const std::vector<Line>& from, const std::vector<Line>& to);

        Diagram& expanduhf();

        Diagram& expandrhf();

        Diagram& symmetrize();

        Diagram& antisymmetrize(const std::vector< std::vector<Line> >& assym);
};

}
}

#endif
