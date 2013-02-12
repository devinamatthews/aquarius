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

#ifndef _AQUARIUS_AUTOCC_FRACTION_HPP_
#define _AQUARIUS_AUTOCC_FRACTION_HPP_

#include <ostream>
#include <string>

namespace aquarius
{
namespace autocc
{

class Fraction;

std::ostream& operator<<(std::ostream& out, const Fraction& f);

bool operator<(const aquarius::autocc::Fraction& f, const double d);
bool operator<(const double d, const aquarius::autocc::Fraction& f);

class Fraction
{
    friend std::ostream& operator<<(std::ostream& out, const Fraction& f);
    friend bool operator<(const Fraction& f, const double d);
    friend bool operator<(const double d, const Fraction& f);

    protected:
        int num, den;

    public:
        Fraction(int num = 1, int den = 1);

        Fraction(const std::string& s);

        operator double() const;

        Fraction operator-() const;

        Fraction& negate();

        Fraction inverse() const;

        Fraction& invert();

        Fraction& operator+=(const Fraction& other);

        Fraction& operator-=(const Fraction& other);

        Fraction& operator*=(const Fraction& other);

        Fraction& operator*=(const int other);

        Fraction& operator/=(const Fraction& other);

        Fraction& operator/=(const int other);

        Fraction operator+(const Fraction& other) const;

        Fraction operator-(const Fraction& other) const;

        Fraction operator*(const Fraction& other) const;

        Fraction operator/(const Fraction& other) const;

        bool operator==(const int other) const;

        bool operator==(const Fraction& other) const;

        bool operator<(const Fraction& other) const;

    protected:
        void reduce();
};

}
}

#endif
