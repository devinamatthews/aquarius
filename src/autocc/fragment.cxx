#include "autocc.hpp"

namespace aquarius
{
namespace autocc
{

ostream& operator<<(ostream& out, const Fragment& f)
{
    out << f.op << "(" << f.out << "," << f.in << ")";
    return out;
}

Fragment::Fragment(const string& s)
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

    /*
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
    */

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
    aquarius::translate(out, from, to);
    aquarius::translate(in, from, to);

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
        sort(particles.begin(), particles.end(),
             [](const pair<Line,Line>& p1, const pair<Line,Line>& p2) { return p1.first < p2.first; });
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
    return uniq_copy(out+in);
}

pair<Manifold,Manifold> Fragment::getShape() const
{
    pair<Manifold,Manifold> p;
    getShape(p.first, p.second);
    return p;
}

void Fragment::getShape(Manifold& left, Manifold& right) const
{
    for (vector<Line>::const_iterator i = out.begin();i != out.end();++i)
    {
        if (i->isVirtual())
        {
            left.np[i->getType()]++;
        }
        else
        {
            right.nh[i->getType()]++;
        }
    }

    for (vector<Line>::const_iterator i = in.begin();i != in.end();++i)
    {
        if (i->isVirtual())
        {
            right.np[i->getType()]++;
        }
        else
        {
            left.nh[i->getType()]++;
        }
    }
}

}
}
