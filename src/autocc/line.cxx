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

static char typechars[] = {'a', 'e', 'i', 'm', 'A', 'E', 'I', 'M'};

ostream& operator<<(ostream& out, const Line& l)
{
    out << typechars[((l.index)&Line::TYPE_MASK) >> 13] << ((l.index)&Line::INDEX_MASK);
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

std::ostream& operator<<(std::ostream& out, const Manifold& m)
{
    return out << m.np << 'p' << m.nh << 'h';
}

const Manifold Manifold::MIN_VALUE = Manifold(0, 0);
const Manifold Manifold::MAX_VALUE = Manifold(INT_MAX, INT_MAX);

Manifold::Manifold(const char *_s)
: np(0), nh(0)
{
    parse(string(_s));
}

Manifold::Manifold(const string& s)
: np(0), nh(0)
{
    parse(s);
}

void Manifold::parse(const string& s)
{
    if (s.length() == 0) return;

    if (s.find('p') != string::npos)
    {
        if (s.find('h') != string::npos)
        {
            if (sscanf(s.c_str(), "%dp%dh", &np, &nh) != 2)
            {
                throw logic_error("malformed manifold: " + s);
            }
        }
        else
        {

            if (sscanf(s.c_str(), "%dp", &np) != 1)
            {
                throw logic_error("malformed manifold: " + s);
            }
        }
    }
    else
    {
        if (s.find('h') != string::npos)
        {
            if (sscanf(s.c_str(), "%dh", &nh) != 1)
            {
                throw logic_error("malformed manifold: " + s);
            }
        }
        else
        {
            throw logic_error("malformed manifold: " + s);
        }
    }
}

Manifold& Manifold::operator+=(const Manifold& other)
{
    np += other.np; nh += other.nh; return *this;
}

Manifold& Manifold::operator-=(const Manifold& other)
{
    np -= other.np; nh -= other.nh; return *this;
}

Manifold Manifold::operator+(const Manifold& other) const
{
    return Manifold(*this) += other;
}

Manifold Manifold::operator-(const Manifold& other) const
{
    return Manifold(*this) -= other;
}

bool Manifold::operator<(const Manifold& other) const
{
   if (np < other.np)
   {
       return true;
   }
   else if (other.np < np)
   {
       return false;
   }
   else
   {
       return nh < other.nh;
   }
}

bool Manifold::operator==(const Manifold& other) const
{
    return np == other.np && nh == other.nh;
}

Line::Line(const string& label)
{
    istringstream iss(label);
    char c;
    unsigned int i;

    iss >> c >> i;

    assert(i <= INDEX_MAX);

    switch (c)
    {
        case 'a':
            index = i;
            break;
        case 'i':
            index = i|_HOLE;
            break;
        case 'e':
            index = i|_INTERNAL;
            break;
        case 'm':
            index = i|_INTERNAL|_HOLE;
            break;
        case 'A':
            index = i|_ALPHA;
            break;
        case 'I':
            index = i|_ALPHA|_HOLE;
            break;
        case 'E':
            index = i|_ALPHA|_INTERNAL;
            break;
        case 'M':
            index = i|_ALPHA|_INTERNAL|_HOLE;
            break;
        default:
            throw logic_error("invalid line type: " + c);
            break;
    }
}

Line::Line(char label)
{
    index = 0;

    if (label >= 'A' && label <= 'P')
    {
        index |= _ALPHA;

        if (label >= 'I')
        {
            index |= _HOLE;

            if (label >= 'M')
            {
                index |= _INTERNAL;
                index |= (unsigned int)(label-'M');
            }
            else
            {
                index |= (unsigned int)(label-'I');
            }
        }
        else
        {
            if (label >= 'E')
            {
                index |= _INTERNAL;
                index |= (unsigned int)(label-'E');
            }
            else
            {
                index |= (unsigned int)(label-'A');
            }
        }
    }
    else if (label >= 'a' && label <= 'p')
    {
        if (label >= 'i')
        {
            index |= _HOLE;

            if (label >= 'm')
            {
                index |= _INTERNAL;
                index |= (unsigned int)(label-'m');
            }
            else
            {
                index |= (unsigned int)(label-'i');
            }
        }
        else
        {
            if (label >= 'e')
            {
                index |= _INTERNAL;
                index |= (unsigned int)(label-'e');
            }
            else
            {
                index |= (unsigned int)(label-'a');
            }
        }
    }
    else
    {
        throw logic_error("invalid label: " + label);
    }
}

Line::Line(unsigned int index, unsigned int type)
{
    this->index = index;

    if (type&INTERNAL)
    {
        if (type&EXTERNAL) throw logic_error("index cannot be both INTERNAL and EXTERNAL");
        this->index |= _INTERNAL;
    }
    else if (!type&EXTERNAL)
    {
        throw logic_error("index must be either INTERNAL or EXTERNAL");
    }

    if (type&HOLE)
    {
        if (type&PARTICLE) throw logic_error("index cannot be both HOLE and PARTICLE");
        this->index |= _HOLE;
    }
    else if (!type&PARTICLE)
    {
        throw logic_error("index must be either HOLE or PARTICLE");
    }

    if (type&ALPHA)
    {
        if (type&BETA) throw logic_error("index cannot be both ALPHA and BETA");
        this->index |= _ALPHA;
    }
    else if (!type&BETA)
    {
        throw logic_error("index must be either ALPHA or BETA");
    }
}

vector<Line> Line::parse(string s)
{
    vector<Line> v;
    if (s.find_first_of("0123456789") != string::npos)
    {
        string l;
        istringstream iss(s);
        while (iss >> l)
        {
            v.push_back(Line(l));
        }
    }
    else
    {
        for (int i = 0;i < s.size();i++) v.push_back(Line(s[i]));
    }
    return v;
}

//Line::operator unsigned int() const
//{
//    return index&INDEX_MASK;
//}

Line& Line::operator=(const Line& line)
{
    this->index = line.index;
    return *this;
}

//Line& Line::operator=(unsigned int index)
//{
//    this->index = (this->index&TYPE_MASK) | index;
//}

bool Line::operator<(const Line& other) const
{
    return index < other.index;
}

bool Line::operator==(const Line& other) const
{
    return index == other.index;
}

int Line::getType() const
{
    int type = 0;

    if (index&_INTERNAL)
    {
        type |= INTERNAL;
    }
    else
    {
        type |= EXTERNAL;
    }

    if (index&_HOLE)
    {
        type |= HOLE;
    }
    else
    {
        type |= PARTICLE;
    }

    if (index&_ALPHA)
    {
        type |= ALPHA;
    }
    else
    {
        type |= BETA;
    }

    return type;
}

int Line::getIndex() const
{
    return index&INDEX_MASK;
}

int Line::toInt() const
{
    return index;
}

bool Line::isInternal() const
{
    return index&_INTERNAL;
}

bool Line::isExternal() const
{
    return !(index&_INTERNAL);
}

bool Line::isHole() const
{
    return index&_HOLE;
}

bool Line::isParticle() const
{
    return !(index&_HOLE);
}

bool Line::isAlpha() const
{
    return index&_ALPHA;
}

bool Line::isBeta() const
{
    return !(index&_ALPHA);
}

Line Line::toAlpha() const
{
    return Line((uint16_t)(index|_ALPHA));
}

Line Line::toBeta() const
{
    return Line((uint16_t)(index&((~_ALPHA)&0xFFFF)));
}

}
}
