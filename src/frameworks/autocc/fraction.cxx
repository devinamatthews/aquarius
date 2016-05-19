#include "autocc.hpp"

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
