#ifndef _AQUARIUS_AUTOCC_TERM_HPP_
#define _AQUARIUS_AUTOCC_TERM_HPP_

#include "util/global.hpp"

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

ostream& operator<<(ostream& out, const aquarius::autocc::Term& t);

class Term
{
    friend class Diagram;
    friend class Fragment;

    friend ostream& operator<<(ostream& out, const Term& t);

    private:
        Diagram::Type type;
        Fraction factor;
        vector<Fragment> fragments;

    public:
        Term(const Term& other)
        : type(other.type), factor(other.factor), fragments(other.fragments) {}

        Term(const Diagram::Type type, const string s);

        Term(const Diagram::Type type = Diagram::SPINORBITAL, const Fraction factor = Fraction()) : type(type), factor(factor) {}

        Term(const Diagram::Type type, const Fraction factor, const vector<Fragment>& fragments);

        vector<Line> indices() const;

        vector<Line> internal() const;

        vector<Line> external() const;

        Term& fixorder(bool all = false);

        Term& fixorder(const vector<Line>& which);

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

        pair<Manifold,Manifold> getShape() const;

        void getShape(Manifold& left, Manifold& right) const;

        Term& translate(const vector<Line>& from, const vector<Line>& to);

        const vector<Fragment>& getFragments() const;

        vector<Fragment>& getFragments();

    protected:

        Diagram sumall() const;

        Diagram suminternal() const;

        Diagram sum(const vector<Line>& which_) const;

        bool checkspin() const;

        Diagram expanduhf() const;

        Diagram expandrhf() const;

        static Diagram doexpand(const Term& term, vector<Diagram>::iterator first, vector<Diagram>::iterator last);

        Diagram symmetrize() const;

        Diagram antisymmetrize(const vector< vector<Line> >& assym) const;

        void canonicalize();
};

}
}

#endif
