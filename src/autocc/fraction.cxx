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

#include <sstream>
#include <stdexcept>
#include <iostream>

#include "autocc.hpp"

using namespace std;

namespace aquarius
{
namespace autocc
{

ostream& operator<<(ostream& out, const Fraction& f)
{
    if (f.den == 1)
    {
        return out << f.num;
    }
    else
    {
        return out << f.num << "/" << f.den;
    }
}

bool operator<(const Fraction& f, const double d)
{
    return f.num < d*f.den;
}

bool operator<(const double d, const Fraction& f)
{

    return d*f.den < f.num;
}

Fraction::Fraction(int num, int den)
: num(num), den(den)
{
    reduce();
}

Fraction::Fraction(const string& s)
{
    int i = s.find('/');

    if (i == string::npos)
    {
        istringstream iss(s);
        iss >> num;
        den = 1;
    }
    else
    {
        istringstream iss(s.substr(i+1));
        iss >> den;
        istringstream iss2(s.substr(0,i));
        iss2 >> num;
    }

    reduce();
}

void Fraction::reduce()
{
    if (den == 0)
    {
        throw logic_error("division by zero");
    }
    else if (den < 0)
    {
        num = -num;
        den = -den;
    }

    int a = (num < 0 ? -num : num);
    int b = den;
    int c;

    while (a != 0)
    {
       c = a; a = b%a;  b = c;
    }

    num /= b;
    den /= b;
}

Fraction::operator double() const
{
    return double(num)/double(den);
}

Fraction Fraction::operator-() const
{
    return Fraction(-num, den);
}

Fraction& Fraction::negate()
{
    num = -num;
    return *this;
}

Fraction Fraction::inverse() const
{
    return Fraction(den, num);
}

Fraction& Fraction::invert()
{
    int tmp = num;
    num = den;
    den = tmp;
    return *this;
}

Fraction& Fraction::operator+=(const Fraction& other)
{
    num = num*other.den + other.num*den;
    den *= other.den;
    reduce();
    return *this;
}

Fraction& Fraction::operator-=(const Fraction& other)
{
    num = num*other.den - other.num*den;
    den *= other.den;
    reduce();
    return *this;
}

Fraction& Fraction::operator*=(const Fraction& other)
{
    num *= other.num;
    den *= other.den;
    reduce();
    return *this;
}

Fraction& Fraction::operator*=(const int other)
{
    num *= other;
    reduce();
    return *this;
}

Fraction& Fraction::operator/=(const Fraction& other)
{
    num *= other.den;
    den *= other.num;
    reduce();
    return *this;
}

Fraction& Fraction::operator/=(const int other)
{
    den *= other;
    reduce();
    return *this;
}

Fraction Fraction::operator+(const Fraction& other) const
{
    return Fraction(num*other.den + other.num*den, den*other.den);
}

Fraction Fraction::operator-(const Fraction& other) const
{
    return Fraction(num*other.den - other.num*den, den*other.den);
}

Fraction Fraction::operator*(const Fraction& other) const
{
    return Fraction(num*other.num, den*other.den);
}

Fraction Fraction::operator/(const Fraction& other) const
{
    return Fraction(num*other.den, den*other.num);
}

bool Fraction::operator<(const Fraction& other) const
{
    return num*other.den < other.num*den;
}

bool Fraction::operator==(const int other) const
{
    return num == other && den == 1;
}

bool Fraction::operator==(const Fraction& other) const
{
    return num == other.num && den == other.den;
}

}
}
