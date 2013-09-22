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

#ifndef _AQUARIUS_AUTOCC_FRAGMENT_HPP_
#define _AQUARIUS_AUTOCC_FRAGMENT_HPP_

#include <ostream>
#include <string>
#include <vector>
#include <cassert>

#include "diagram.hpp"

namespace aquarius
{
namespace autocc
{

class Fragment;
class Line;
class Term;
class Manifold;

std::ostream& operator<<(std::ostream& out, const Fragment& f);

template <typename T>
int relativeSign(const T& s1, const T& s2)
{
    return relativeSign(s1.begin(), s1.end(), s2.begin(), s2.end());
}

template <typename iterator>
int relativeSign(iterator b1, iterator e1, iterator b2, iterator e2)
{
    assert((int)(e2-b2) == (int)(e1-b1));

    int sign = 1;
    int n = (int)(e1-b1);
    std::vector<bool> seen(n, false);

    for (int i = 0;i < n;i++)
    {
        if (seen[i]) continue;
        int j = i;
        while (true)
        {
            int k;
            for (k = 0;k < n && (!(*(b1+k) == *(b2+j)) || seen[k]);k++);
            assert(k < n);
            j = k;
            seen[j] = true;
            if (j == i) break;
            sign = -sign;
        }
    }

    return sign;
}

class Fragment
{
    friend class Term;
    friend class Diagram;

    friend std::ostream& operator<<(std::ostream& out, const Fragment& f);

    protected:
        std::string op;
        std::vector<Line> out;
        std::vector<Line> in;

    public:
        Fragment(const std::string& s);

        Fragment(const std::string& op, const std::vector<Line>& out, const std::vector<Line>& in)
        : op(op), out(out), in(in) {}

        bool operator<(const Fragment& other) const;

        bool operator==(const Fragment& other) const;

        std::pair<Manifold,Manifold> getShape() const;

        void getShape(Manifold& left, Manifold& right) const;

        std::vector<Line> indices() const;

        const std::string& getOp() const;

        const std::vector<Line>& getIndicesOut() const;

        const std::vector<Line>& getIndicesIn() const;

    protected:
        Fragment& translate(const std::vector<Line>& from, const std::vector<Line>& to);

        int canonicalize(const Diagram::Type type);
};

}
}

#endif
