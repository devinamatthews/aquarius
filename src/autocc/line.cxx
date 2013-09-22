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

#include <ostream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <stdexcept>
#include <cstdio>
#include <climits>

#include "autocc.hpp"

using namespace std;

namespace aquarius
{
namespace autocc
{

ostream& operator<<(ostream& out, const Line& l)
{
    char c = (char)((int)'A'+l.getType()+((int)'a'-(int)'A')*l.getSpin());
    out << c << (l.isOccupied() ? 'o' : 'v') << l.getIndex();
    return out;
}

ostream& operator<<(ostream& out, const vector<Line>& v)
{
    for (vector<Line>::const_iterator l = v.begin();l != v.end();)
    {
        out << *l;
        if (++l != v.end()) out << " ";
    }

    return out;
}

Manifold& Manifold::operator+=(const Manifold& other)
{
    for (map<int,int>::const_iterator i = other.np.begin();i != other.np.end();++i)
    {
        np[i->first] += i->second;
    }
    for (map<int,int>::const_iterator i = other.nh.begin();i != other.nh.end();++i)
    {
        nh[i->first] += i->second;
    }
    return *this;
}

Manifold& Manifold::operator-=(const Manifold& other)
{
    for (map<int,int>::const_iterator i = other.np.begin();i != other.np.end();++i)
    {
        np[i->first] -= i->second;
    }
    for (map<int,int>::const_iterator i = other.nh.begin();i != other.nh.end();++i)
    {
        nh[i->first] -= i->second;
    }
    return *this;
}

Manifold Manifold::operator+(const Manifold& other) const
{
    return Manifold(*this) += other;
}

Manifold Manifold::operator-(const Manifold& other) const
{
    return Manifold(*this) -= other;
}


bool Manifold::operator==(const Manifold& other) const
{
    for (map<int,int>::const_iterator i1  = np.begin(), i2  = other.np.begin();
                                      i1 != np.end() || i2 != other.np.end();)
    {
        int v1, v2;

        if (i2 == other.np.end() || i1->first < i2->first)
        {
            v1 = i1->second;
            v2 = 0;
            ++i1;
        }
        else if (i1 == np.end() || i2->first < i1->first)
        {
            v1 = 0;
            v2 = i2->second;
            ++i2;
        }
        else
        {
            v1 = i1->second;
            v2 = i2->second;
            ++i1;
            ++i2;
        }

        if (v1 != v2) return false;
    }

    for (map<int,int>::const_iterator i1  = nh.begin(), i2  = other.nh.begin();
                                      i1 != nh.end() || i2 != other.nh.end();)
    {
        int v1, v2;

        if (i2 == other.nh.end() || i1->first < i2->first)
        {
            v1 = i1->second;
            v2 = 0;
            ++i1;
        }
        else if (i1 == nh.end() || i2->first < i1->first)
        {
            v1 = 0;
            v2 = i2->second;
            ++i2;
        }
        else
        {
            v1 = i1->second;
            v2 = i2->second;
            ++i1;
            ++i2;
        }

        if (v1 != v2) return false;
    }

    return true;
}

bool Manifold::operator<(const Manifold& other) const
{
    for (map<int,int>::const_iterator i1  = np.begin(), i2  = other.np.begin();
                                      i1 != np.end() || i2 != other.np.end();)
    {
        int v1, v2;

        if (i2 == other.np.end() || i1->first < i2->first)
        {
            v1 = i1->second;
            v2 = 0;
            ++i1;
        }
        else if (i1 == np.end() || i2->first < i1->first)
        {
            v1 = 0;
            v2 = i2->second;
            ++i2;
        }
        else
        {
            v1 = i1->second;
            v2 = i2->second;
            ++i1;
            ++i2;
        }

        if (v1 < v2)
        {
            return true;
        }
        else if (v2 < v1)
        {
            return false;
        }
    }

    for (map<int,int>::const_iterator i1  = nh.begin(), i2  = other.nh.begin();
                                      i1 != nh.end() || i2 != other.nh.end();)
    {
        int v1, v2;

        if (i2 == other.nh.end() || i1->first < i2->first)
        {
            v1 = i1->second;
            v2 = 0;
            ++i1;
        }
        else if (i1 == nh.end() || i2->first < i1->first)
        {
            v1 = 0;
            v2 = i2->second;
            ++i2;
        }
        else
        {
            v1 = i1->second;
            v2 = i2->second;
            ++i1;
            ++i2;
        }

        if (v1 < v2)
        {
            return true;
        }
        else if (v2 < v1)
        {
            return false;
        }
    }

    return false;
}

Line::Line(const string& label)
{
    istringstream iss(label);
    char c, o;
    unsigned int i;

    iss >> c >> o >> i;

    assert(isalpha(c));
    assert(o == 'o' || o == 'v');
    assert(i <= INDEX_MAX);

    index = (islower(c) ? _BETA : 0)|(o == 'o' ? _OCCUPIED : 0)|((toupper(c)-(int)'A')<<TYPE_SHIFT)|i;
}

/*
Line::Line(char label)
{
    assert(isalpha(c));
    index = (islower(c) ? _BETA : 0)|((toupper(c)-(int)'A')<<TYPE_SHIFT);
}
*/

Line::Line(int index, int type, int occupancy, int spin)
{
    assert(index >= 0 && index < INDEX_MAX);
    assert(type >= 0 && type < TYPE_MAX);
    assert(occupancy == OCCUPIED || occupancy == VIRTUAL);
    assert(spin == ALPHA || spin == BETA);

    this->index = (spin<<SPIN_SHIFT)|(occupancy<<OCC_SHIFT)|(type<<TYPE_SHIFT)|index;
}

Line::Line(int index, char type_, int occupancy, int spin)
{
    int type = toupper(type_)-(int)'A';

    assert(index >= 0 && index < INDEX_MAX);
    assert(type >= 0 && type < TYPE_MAX);
    assert(occupancy == OCCUPIED || occupancy == VIRTUAL);
    assert(spin == ALPHA || spin == BETA);

    this->index = (spin<<SPIN_SHIFT)|(occupancy<<OCC_SHIFT)|(type<<TYPE_SHIFT)|index;
}

vector<Line> Line::parse(const string& s)
{
    vector<Line> v;
    string l;
    istringstream iss(s);
    while (iss >> l)
    {
        v.push_back(Line(l));
    }
    return v;
}

}
}
