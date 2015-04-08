#ifndef _AQUARIUS_AUTOCC_DIAGRAM_HPP_
#define _AQUARIUS_AUTOCC_DIAGRAM_HPP_

#include "util/global.hpp"

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

ostream& operator<<(ostream& out, const Diagram& d);

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
        vector<Term> terms;

    public:
        Diagram(const Type type) : type(type) {}

        Diagram(const Type type, const vector<Term>& terms);

        Diagram(const Type type, const vector<string>& strings);

        Diagram& convert(const Type type);

        Diagram& fixorder(bool all = false);

        Diagram& fixorder(const vector<Line>& which);

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

        const vector<Term>& getTerms() const;

        Diagram& sumall();

        Diagram& sum(const vector<Line>& which);

        Diagram& antisymmetrize(const vector<vector<Line>>& assym);

    protected:
        Diagram& suminternal();

        Diagram& translate(const vector<Line>& from, const vector<Line>& to);

        Diagram& expanduhf();

        Diagram& expandrhf();

        Diagram& symmetrize();
};

}
}

#endif
