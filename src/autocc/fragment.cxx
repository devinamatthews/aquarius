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

#include <stdexcept>
#include <algorithm>
#include <utility>
#include <cassert>
#include <iostream>

#include "autocc.hpp"

using namespace std;

namespace aquarius
{
namespace autocc
{

int relativeSign(const vector<Line>& s1, const vector<Line>& s2)
{
    int i, j, k;
    int sign = 1;
    bool *seen = new bool[s1.size()];

    for (i = 0;i < s1.size();i++) seen[i] = false;

    for (i = 0;i < s1.size();i++)
    {
        if (seen[i]) continue;
        j = i;
        while (true)
        {
            for (k = 0;k < s1.size() && (!(s1[k] == s2[j]) || seen[k]);k++);
            assert(k < s1.size());
            j = k;
            seen[j] = true;
            if (j == i) break;
            sign = -sign;
        }
    }

    delete[] seen;

    return sign;
}

ostream& operator<<(ostream& out, const Fragment& f)
{
    out << f.op << "(" << f.out << "," << f.in << ")";
    return out;
}

Fragment::Fragment(const std::string& s)
{
    int lp = s.find('(');
    int comma = s.find(',');
    int rp = s.find(')');

    if (lp == string::npos || rp == string::npos)
    {
        throw logic_error("malformed fragment: " + s);
    }

    if (rp-lp == 3)
    {
        op = s.substr(0,lp);
        out = Line::parse(s.substr(lp+1, 1));
        in = Line::parse(s.substr(lp+2, 1));
    }
    else
    {
        if (comma == string::npos)
        {
            throw logic_error("malformed fragment: " + s);
        }

        op = s.substr(0,lp);
        out = Line::parse(s.substr(lp+1, comma-lp-1));
        in = Line::parse(s.substr(comma+1, rp-comma-1));
    }
}

bool Fragment::operator<(const Fragment& other) const
{
    int i = op.compare(other.op);

    if (i < 0)
    {
        return true;
    }
    else if (i > 0)
    {
        return false;
    }

    if (getShape() < other.getShape())
    {
        return true;
    }
    else if (other.getShape() < getShape())
    {
        return false;
    }

    int npext1 = filter_copy(in+out, isType<PARTICLE+EXTERNAL>()).size();
    int npext2 = filter_copy(other.in+other.out, isType<PARTICLE+EXTERNAL>()).size();

    if (npext1 < npext2)
    {
        return true;
    }
    else if (npext2 < npext1)
    {
        return false;
    }

    int nhext1 = filter_copy(in+out, isType<HOLE+EXTERNAL>()).size();
    int nhext2 = filter_copy(other.in+other.out, isType<HOLE+EXTERNAL>()).size();

    if (nhext1 < nhext2)
    {
        return true;
    }
    else if (nhext2 < nhext1)
    {
        return false;
    }

    if (lexicographical_compare(out.begin(), out.end(), other.out.begin(), other.out.end()))
    {
        return true;
    }
    else if (lexicographical_compare(other.out.begin(), other.out.end(), out.begin(), out.end()))
    {
        return false;
    }

    if (lexicographical_compare(in.begin(), in.end(), other.in.begin(), other.in.end()))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Fragment::operator==(const Fragment& other) const
{
    return op.compare(other.op) == 0 &&
           equal(out.begin(), out.end(), other.out.begin()) &&
           equal(in.begin(), in.end(), other.in.begin());
}

const string& Fragment::getOp() const
{
    return op;
}

const vector<Line>& Fragment::getIndicesOut() const
{
    return out;
}

const vector<Line>& Fragment::getIndicesIn() const
{
    return in;
}

Fragment& Fragment::translate(const vector<Line>& from, const vector<Line>& to)
{
    std::translate<Line>(out, from, to);
    std::translate<Line>(in, from, to);

    return *this;
}

int Fragment::canonicalize(const Diagram::Type type)
{
    int sign = 1;

    if (type == Diagram::UHF || type == Diagram::SPINORBITAL)
    {
        vector<Line> newout(out);
        vector<Line> newin(in);
        sort(newout.begin(), newout.end());
        sort(newin.begin(), newin.end());
        //cout << out << "," << in << " -> " << newout << "," << newin << " " << relativeSign(out, newout) << " " << relativeSign(in, newin) << endl;
        sign = relativeSign(out, newout)*relativeSign(in, newin);
        out = newout;
        in = newin;
    }
    else
    {
        vector< pair<Line,Line> > particles;
        for (int i = 0;i < out.size();i++) particles.push_back(make_pair(out[i], in[i]));
        sort(particles.begin(), particles.end(), compareFirst<Line,Line>);
        for (int i = 0;i < out.size();i++)
        {
            out[i] = particles[i].first;
            in[i] = particles[i].second;
        }
    }

    return sign;
}

vector<Line> Fragment::indices() const
{
    vector<Line> indices = out + in;
    sort(indices.begin(), indices.end());
    vector<Line>::iterator end = unique(indices.begin(), indices.end());
    indices.resize((int)(end - indices.begin()));

    return indices;
}

pair<Manifold,Manifold> Fragment::getShape() const
{
    return make_pair(Manifold(count_if(out.begin(), out.end(), isType<PARTICLE+EXTERNAL>()),
                              count_if( in.begin(),  in.end(), isType<HOLE+EXTERNAL>())),
                     Manifold(count_if( in.begin(),  in.end(), isType<PARTICLE+EXTERNAL>()),
                              count_if(out.begin(), out.end(), isType<HOLE+EXTERNAL>())));
}

void Fragment::getShape(Manifold& left, Manifold& right) const
{
    left.np  = count_if(out.begin(), out.end(), isType<PARTICLE+EXTERNAL>());
    left.nh  = count_if( in.begin(),  in.end(), isType<HOLE+EXTERNAL>());
    right.np = count_if( in.begin(),  in.end(), isType<PARTICLE+EXTERNAL>());
    right.nh = count_if(out.begin(), out.end(), isType<HOLE+EXTERNAL>());
}

}
}
