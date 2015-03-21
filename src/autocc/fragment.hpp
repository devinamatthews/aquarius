#ifndef _AQUARIUS_AUTOCC_FRAGMENT_HPP_
#define _AQUARIUS_AUTOCC_FRAGMENT_HPP_

#include "util/global.hpp"

#include "diagram.hpp"

namespace aquarius
{
namespace autocc
{

class Fragment;
class Line;
class Term;
class Manifold;

ostream& operator<<(ostream& out, const Fragment& f);

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
    vector<bool> seen(n, false);

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

    friend ostream& operator<<(ostream& out, const Fragment& f);

    protected:
        string op;
        vector<Line> out;
        vector<Line> in;

    public:
        Fragment(const string& s);

        Fragment(const string& op, const vector<Line>& out, const vector<Line>& in)
        : op(op), out(out), in(in) {}

        bool operator<(const Fragment& other) const;

        bool operator==(const Fragment& other) const;

        pair<Manifold,Manifold> getShape() const;

        void getShape(Manifold& left, Manifold& right) const;

        vector<Line> indices() const;

        const string& getOp() const;

        const vector<Line>& getIndicesOut() const;

        const vector<Line>& getIndicesIn() const;

    protected:
        Fragment& translate(const vector<Line>& from, const vector<Line>& to);

        int canonicalize(const Diagram::Type type);
};

}
}

#endif
