#ifndef _AQUARIUS_AUTOCC_FRACTION_HPP_
#define _AQUARIUS_AUTOCC_FRACTION_HPP_

#include "util/global.hpp"

namespace aquarius
{
namespace autocc
{

class Fraction;

ostream& operator<<(ostream& out, const Fraction& f);

bool operator<(const aquarius::autocc::Fraction& f, const double d);
bool operator<(const double d, const aquarius::autocc::Fraction& f);

class Fraction
{
    friend ostream& operator<<(ostream& out, const Fraction& f);
    friend bool operator<(const Fraction& f, const double d);
    friend bool operator<(const double d, const Fraction& f);

    protected:
        int num, den;

    public:
        Fraction(int num = 1, int den = 1);

        Fraction(const string& s);

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
