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

#ifndef _AQUARIUS_AUTOCC_TERM_HPP_
#define _AQUARIUS_AUTOCC_TERM_HPP_

#include <ostream>
#include <vector>

#include "fraction.hpp"
#include "diagram.hpp"

namespace aquarius
{
namespace autocc
{

class Term;
class Fragment;
class Line;
class Manifold;

std::ostream& operator<<(std::ostream& out, const aquarius::autocc::Term& t);

class Term
{
    friend class Diagram;
    friend class Fragment;

    friend std::ostream& operator<<(std::ostream& out, const Term& t);

    private:
        Diagram::Type type;
        Fraction factor;
        std::vector<Fragment> fragments;

    public:
        Term(const Term& other)
        : type(other.type), factor(other.factor), fragments(other.fragments) {}

        Term(const Diagram::Type type, const std::string s);

        Term(const Diagram::Type type = Diagram::SPINORBITAL, const Fraction factor = Fraction()) : type(type), factor(factor) {}

        Term(const Diagram::Type type, const Fraction factor, const std::vector<Fragment>& fragments);

        std::vector<Line> indices() const;

        Term& fixorder(bool all = false);

        Term& fixorder(const std::vector<Line>& which);

        Term& fixexternal();

        Term& negate();

        Term operator-() const;

        Term& operator*=(const Term& other);

        Term& operator*=(const Fragment& other);

        Term& operator*=(const Fraction& other);

        Term& operator*=(const int other);

        Diagram operator*(const Diagram& other) const;

        Term operator*(const Term& other) const;

        Term operator*(const Fragment& other) const;

        Term operator*(const Fraction& other) const;

        Term operator*(const int other) const;

        bool operator<(const Term& other) const;

        bool operator==(const Term& other) const;

        const Fraction& getFactor() const;

        std::pair<Manifold,Manifold> getShape() const;

        void getShape(Manifold& left, Manifold& right) const;

        Term& translate(const std::vector<Line>& from, const std::vector<Line>& to);

        const std::vector<Fragment>& getFragments() const;

        std::vector<Fragment>& getFragments();

    protected:

        Diagram sumall() const;

        Diagram suminternal() const;

        Diagram sum(const std::vector<Line>& which_) const;

        bool checkspin() const;

        Diagram expanduhf() const;

        Diagram expandrhf() const;

        static Diagram doexpand(const Term& term, std::vector<Diagram>::iterator first, std::vector<Diagram>::iterator last);

        Diagram symmetrize() const;

        Diagram antisymmetrize(const std::vector< std::vector<Line> >& assym) const;

        void canonicalize();
};

}
}

#endif
