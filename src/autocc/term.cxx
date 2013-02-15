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

#include <algorithm>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "autocc.hpp"

using namespace std;

namespace aquarius
{
namespace autocc
{

ostream& operator<<(ostream& out, const Term& t)
{
    out << t.factor;
    for (vector<Fragment>::const_iterator it = t.fragments.begin();it != t.fragments.end();++it) out << " " << *it;
    return out;
}


Term::Term(const Diagram::Type type, const std::string s)
: type(type)
{
    istringstream iss(s);
    string token;

    while (iss >> token)
    {
        if (token == "-")
        {
            if (factor == 1)
            {
                factor = -1;
            }
            else
            {
                throw logic_error("misplaced minus sign: " + s);
            }
        }
        else if (isdigit(token[0]) ||
                 token.find('/') != string::npos ||
                 token[0] == '-')
        {
            if (factor == -1)
            {
                factor = -Fraction(token);
            }
            else if (factor == 1)
            {
                factor = Fraction(token);
            }
            else
            {
                throw logic_error("multiple factors given: " + s);
            }
        }
        else
        {
            fragments.push_back(Fragment(token));
        }
    }

    canonicalize();
}

Term::Term(const Diagram::Type type, const Fraction factor, const vector<Fragment>& fragments)
: type(type), factor(factor), fragments(fragments)
{
    canonicalize();
}

Term& Term::fixorder(const vector<Line>& which)
{
    vector<Line> oldah, oldbh, oldap, oldbp;
    vector<Line> newah, newbh, newap, newbp;

    int ahidx = 0, bhidx = 0, apidx = 0, bpidx = 0;
    vector<Line> inds = indices();
    for (vector<Line>::iterator l = inds.begin();l != inds.end();++l)
    {
        if (l->isInternal() && find(which.begin(), which.end(), *l) == which.end())
        {
            if (l->isParticle())
            {
                if (l->isAlpha())
                {
                    apidx++;
                }
                else
                {
                    bpidx++;
                }
            }
            else
            {
                if (l->isAlpha())
                {
                    ahidx++;
                }
                else
                {
                    bhidx++;
                }
            }
        }
    }

    for (vector<Fragment>::iterator f = getFragments().begin();f != getFragments().end();++f)
    {
        vector<Line> internal;

        for (vector<Line>::iterator i = f->out.begin();i != f->out.end();++i)
        {
            if (find(which.begin(), which.end(), *i) != which.end())
            {
                internal.push_back(*i);
            }
        }

        for (vector<Line>::iterator i = f->in.begin();i != f->in.end();++i)
        {
            if (find(which.begin(), which.end(), *i) != which.end())
            {
                internal.push_back(*i);
            }
        }

        for (vector<Line>::iterator l = internal.begin();l != internal.end();++l)
        {
            if (l->isParticle())
            {
                if (l->isAlpha())
                {
                    if (find(oldap.begin(), oldap.end(), *l) != oldap.end()) continue;
                    oldap.push_back(*l);
                    newap.push_back(Line(apidx++, INTERNAL+PARTICLE+ALPHA));
                }
                else
                {
                    if (find(oldbp.begin(), oldbp.end(), *l) != oldbp.end()) continue;
                    oldbp.push_back(*l);
                    newbp.push_back(Line(bpidx++, INTERNAL+PARTICLE+BETA));
                }
            }
            else
            {
                if (l->isAlpha())
                {
                    if (find(oldah.begin(), oldah.end(), *l) != oldah.end()) continue;
                    oldah.push_back(*l);
                    newah.push_back(Line(ahidx++, INTERNAL+HOLE+ALPHA));
                }
                else
                {
                    if (find(oldbh.begin(), oldbh.end(), *l) != oldbh.end()) continue;
                    oldbh.push_back(*l);
                    newbh.push_back(Line(bhidx++, INTERNAL+HOLE+BETA));
                }
            }
        }
    }

    translate(oldah+oldbh+oldap+oldbp, newah+newbh+newap+newbp);

    return *this;
}

Term& Term::fixorder(bool all)
{
    vector<Line> oldah, oldbh, oldap, oldbp;
    vector<Line> newah, newbh, newap, newbp;

    int ahidx = 0, bhidx = 0, apidx = 0, bpidx = 0;
    for (vector<Fragment>::iterator f = getFragments().begin();f != getFragments().end();++f)
    {
        vector<Line> internal;

        if (all)
        {
            internal += f->out + f->in;
        }
        else
        {
            internal += filter_copy(f->out + f->in, isInternal());
        }

        for (vector<Line>::iterator l = internal.begin();l != internal.end();++l)
        {
            if (l->isParticle())
            {
                if (l->isAlpha())
                {
                    if (find(oldap.begin(), oldap.end(), *l) != oldap.end()) continue;
                    oldap.push_back(*l);
                    newap.push_back(Line(apidx++, INTERNAL+PARTICLE+ALPHA));
                }
                else
                {
                    if (find(oldbp.begin(), oldbp.end(), *l) != oldbp.end()) continue;
                    oldbp.push_back(*l);
                    newbp.push_back(Line(bpidx++, INTERNAL+PARTICLE+BETA));
                }
            }
            else
            {
                if (l->isAlpha())
                {
                    if (find(oldah.begin(), oldah.end(), *l) != oldah.end()) continue;
                    oldah.push_back(*l);
                    newah.push_back(Line(ahidx++, INTERNAL+HOLE+ALPHA));
                }
                else
                {
                    if (find(oldbh.begin(), oldbh.end(), *l) != oldbh.end()) continue;
                    oldbh.push_back(*l);
                    newbh.push_back(Line(bhidx++, INTERNAL+HOLE+BETA));
                }
            }
        }
    }

    translate(oldah+oldbh+oldap+oldbp, newah+newbh+newap+newbp);

    /*
    vector<Line> old(newah+newbh+newap+newbp);
    vector<Line> mininds(old);
    string minstr = str(*this);
    do
    {
        do
        {
            do
            {
                do
                {
                    Term tmp(*this);
                    tmp.translate(old, newah+newbh+newap+newbp);
                    string tmpstr = str(tmp);
                    if (tmpstr < minstr)
                    {
                        minstr = tmpstr;
                        mininds = newah+newbh+newap+newbp;
                    }
                }
                while (next_permutation(newbh.begin(), newbh.end()));
            }
            while (next_permutation(newah.begin(), newah.end()));
        }
        while (next_permutation(newbp.begin(), newbp.end()));
    }
    while (next_permutation(newap.begin(), newap.end()));

    translate(old, mininds);
    */

    return *this;
}

Term& Term::fixexternal()
{
    vector<Line> oldah, oldbh, oldap, oldbp;
    vector<Line> newah, newbh, newap, newbp;

    int ahidx = 0, bhidx = 0, apidx = 0, bpidx = 0;
    for (vector<Fragment>::iterator f = getFragments().begin();f != getFragments().end();++f)
    {
        vector<Line> ext = filter_copy(f->out + f->in, isExternal());
        for (vector<Line>::iterator l = ext.begin();l != ext.end();++l)
        {
            if (l->isParticle())
            {
                if (l->isAlpha())
                {
                    if (find(oldap.begin(), oldap.end(), *l) != oldap.end()) continue;
                    oldap.push_back(*l);
                    newap.push_back(Line(apidx++, EXTERNAL+PARTICLE+ALPHA));
                }
                else
                {
                    if (find(oldbp.begin(), oldbp.end(), *l) != oldbp.end()) continue;
                    oldbp.push_back(*l);
                    newbp.push_back(Line(bpidx++, EXTERNAL+PARTICLE+BETA));
                }
            }
            else
            {
                if (l->isAlpha())
                {
                    if (find(oldah.begin(), oldah.end(), *l) != oldah.end()) continue;
                    oldah.push_back(*l);
                    newah.push_back(Line(ahidx++, EXTERNAL+HOLE+ALPHA));
                }
                else
                {
                    if (find(oldbh.begin(), oldbh.end(), *l) != oldbh.end()) continue;
                    oldbh.push_back(*l);
                    newbh.push_back(Line(bhidx++, EXTERNAL+HOLE+BETA));
                }
            }
        }
    }

    translate(oldah+oldbh+oldap+oldbp, newah+newbh+newap+newbp);

    *this *= relativeSign(oldah, newah);
    *this *= relativeSign(oldbh, newbh);
    *this *= relativeSign(oldap, newap);
    *this *= relativeSign(oldbp, newbp);

    return *this;
}

Term& Term::operator*=(const Term& other)
{
    if (type != other.type) throw logic_error("only terms of the same type may be multiplied");
    factor *= other.factor;
    fragments.insert(fragments.end(), other.fragments.begin(), other.fragments.end());
    //sort(fragments.begin(), fragments.end());
    canonicalize();
    return *this;
}

Term& Term::operator*=(const Fragment& other)
{
    fragments.push_back(other);
    //sort(fragments.begin(), fragments.end());
    canonicalize();
    return *this;
}

Term& Term::operator*=(const Fraction& other)
{
    factor *= other;
    return *this;
}

Term& Term::operator*=(const int other)
{
    factor *= other;
    return *this;
}

Diagram Term::operator*(const Diagram& other) const
{
    return other*(*this);
}

Term Term::operator*(const Term& other) const
{
    return Term(*this) *= other;
}

Term Term::operator*(const Fragment& other) const
{
    return Term(*this) *= other;
}

Term Term::operator*(const Fraction& other) const
{
    return Term(*this) *= other;
}

Term Term::operator*(const int other) const
{
    return Term(*this) *= other;
}

Term& Term::negate()
{
    factor.negate();
    return *this;
}

Term Term::operator-() const
{
    return Term(type, -factor, fragments);
}

bool Term::operator<(const Term& other) const
{
    if (fragments.size() < other.fragments.size())
    {
        return true;
    }
    else if (fragments.size() > other.fragments.size())
    {
        return false;
    }

    for (int i = 0;i < fragments.size();i++)
    {
        if (fragments[i] < other.fragments[i])
        {
            return true;
        }
        else if (other.fragments[i] < fragments[i])
        {
            return false;
        }
    }

    return false;
}

bool Term::operator==(const Term& other) const
{
    if (fragments.size() != other.fragments.size()) return false;

    for (int i = 0;i < fragments.size();i++)
    {
        if (!(fragments[i] == other.fragments[i])) return false;
    }

    return true;
}

const Fraction& Term::getFactor() const
{
    return factor;
}

Diagram Term::sumall() const
{
    Diagram diagram(type);

    vector<Line> indices(this->indices());
    vector<Line> old(indices);

    vector<Line> uhf(old);
    for (int x = 0;x < 1<<old.size();x++)
    {
        for (int i = 0;i < old.size();i++)
        {
            if ((x>>i) & 1)
            {
                uhf[i] = old[i].toAlpha();
            }
            else
            {
                uhf[i] = old[i].toBeta();
            }
        }

        Term newterm(*this);
        newterm.translate(old, uhf);

        if (newterm.checkspin())
        {
            //cout << "before:" << endl << diagram;
            //cout << "trans: " << old << "|" << uhf << endl;
            //cout << "newterm: " << newterm << endl;
            diagram += newterm;
            //cout << "after:" << endl << diagram << endl;
        }
    }

    return diagram;
}

Diagram Term::suminternal() const
{
    Diagram diagram(type);

    vector<Line> indices(this->indices());
    vector<Line> old(indices);
    filter(old, isInternal());

    //cout << "inds: " << indices << endl;
    //cout << "old: " << old << endl;

    vector<Line> uhf(old);
    for (int x = 0;x < 1<<old.size();x++)
    {
        for (int i = 0;i < old.size();i++)
        {
            if ((x>>i) & 1)
            {
                uhf[i] = old[i].toAlpha();
            }
            else
            {
                uhf[i] = old[i].toBeta();
            }
        }

        //cout << "trans: " << old << "|" << uhf << endl;

        Term newterm(*this);
        newterm.translate(old, uhf);

        //cout << "newterm: " << newterm << endl;

        if (newterm.checkspin())
        {
            //print 'old:',self
            //print 'new:',newterm
            diagram += newterm;
        }
    }

    return diagram;
}

Diagram Term::sum(const vector<Line>& which_) const
{
    Diagram diagram(type);

    vector<Line> which(which_);
    vector<Line> indices(this->indices());
    vector<Line> old(indices);

    sort(which.begin(), which.end());
    vector<Line>::iterator end = set_intersection(indices.begin(), indices.end(),
                                                  which.begin(), which.end(), old.begin());
    old.resize((int)(end-old.begin()));

    //cout << "inds: " << indices << endl;
    //cout << "old: " << old << endl;

    vector<Line> uhf(old);
    for (int x = 0;x < 1<<old.size();x++)
    {
        for (int i = 0;i < old.size();i++)
        {
            if ((x>>i) & 1)
            {
                uhf[i] = old[i].toAlpha();
            }
            else
            {
                uhf[i] = old[i].toBeta();
            }
        }

        //cout << "trans: " << old << "|" << uhf << endl;

        Term newterm(*this);
        newterm.translate(old, uhf);

        //cout << "newterm: " << newterm << endl;

        if (newterm.checkspin())
        {
            //print 'old:',self
            //print 'new:',newterm
            diagram += newterm;
        }
    }

    return diagram;
}

bool Term::checkspin() const
{
    for (vector<Fragment>::const_iterator it = fragments.begin();it != fragments.end();++it)
    {
        int nl = 0;
        for (int i = 0;i < it->out.size();i++) if (it->out[i].isAlpha()) nl++;
        int nr = 0;
        for (int i = 0;i < it->in.size();i++) if (it->in[i].isAlpha()) nr++;
        if (abs(nl-nr) > abs((int)(it->out.size()-it->in.size()))) return false;
    }

    return true;
}

Term& Term::translate(const vector<Line>& from, const vector<Line>& to)
{
    for (vector<Fragment>::iterator it = fragments.begin();it != fragments.end();++it)
    {
        it->translate(from, to);
    }

    canonicalize();

    return *this;
}

vector<Line> Term::indices() const
{
    vector<Line> indices;

    for (vector<Fragment>::const_iterator f = fragments.begin();f != fragments.end();++f)
    {
        vector<Line> uniq = f->out + f->in;
        sort(uniq.begin(), uniq.end());
        vector<Line>::iterator end = unique(uniq.begin(), uniq.end());
        uniq.resize((int)(end - uniq.begin()));

        indices += uniq;
    }

    sort(indices.begin(), indices.end());
    vector<Line>::iterator end = unique(indices.begin(), indices.end());
    indices.resize((int)(end - indices.begin()));

    return indices;
}

const vector<Fragment>& Term::getFragments() const
{
    return fragments;
}

vector<Fragment>& Term::getFragments()
{
    return fragments;
}

Diagram Term::expanduhf() const
{
    vector<Diagram> choices;

    for (vector<Fragment>::const_iterator it = fragments.begin();it != fragments.end();++it)
    {
        Fragment f(*it);

        vector< vector<Line> > alpha;
        vector< vector<Line> > beta;
        for (int i = 0;i < f.out.size();i++) f.out[i] = f.out[i].toBeta();
        for (int i = 0;i < f.in.size();i++)
        {
            if (f.in[i].isAlpha())
            {
                f.in[i] = f.in[i].toBeta();
                alpha.push_back(vector<Line>(1, f.in[i]));
            }
            else
            {
                beta.push_back(vector<Line>(1, f.in[i]));
            }
        }

        choices.push_back(Term(type, Fraction(1,1), vector<Fragment>(1, f)).antisymmetrize(alpha).antisymmetrize(beta));
    }

    return doexpand(Term(type, factor), choices.begin(), choices.end());
}

Diagram Term::expandrhf() const
{
    vector<Diagram> choices;

    for (vector<Fragment>::const_iterator it = fragments.begin();it != fragments.end();++it)
    {
        Fragment f(*it);

        vector<Line> f1;
        vector<Line> f2;
        for (int i = 0;i < f.out.size();i++)
        {
            if (f.out[i].isAlpha())
            {
                if (f.out[i].isInternal())
                {
                    f.out[i] = f.out[i].toBeta();
                    f1 += f.out[i];
                }
                else
                {
                    f2 += f.out[i];
                }
            }
            else
            {
                f2 += f.out[i];
            }
        }

        for (int i = 0;i < f.in.size();i++)
        {
            if (f.in[i].isAlpha() && f.in[i].isInternal())
            {
                f.in[i] = f.in[i].toBeta();
            }
        }

        //cout << endl << "f1: " << f1 << endl;
        //cout << "f2: " << f2 << endl;

        Diagram choice(type, vector<Term>(1, Term(type, Fraction(1,1), vector<Fragment>(1, f))));
        for (int i = 0;i < f1.size();i++)
        {
            vector< vector<Line> > assym;
            assym.push_back(slice(f1, i, i+1));
            assym.push_back(slice(f1, i+1) + f2);
            //cout << "assym: " << assym[0] << "|" << assym[1] << endl;
            //cout << "choice: " << endl << choice;
            //cout << "new: " << endl << Diagram(choice).antisymmetrize(assym);
            choice += Diagram(choice).antisymmetrize(assym);
            //cout << "choice: " << endl << choice;
        }
        choices.push_back(choice);
    }

    return doexpand(Term(type, factor), choices.begin(), choices.end());
}

Diagram Term::doexpand(const Term& term, vector<Diagram>::iterator first, vector<Diagram>::iterator last)
{
    if (first == last) return Diagram(term.type, vector<Term>(1, term));

    Diagram diagram(term.type);
    for (vector<Term>::iterator it = first->terms.begin();it != first->terms.end();++it)
    {
        diagram += doexpand(term*(*it), first+1, last);
    }

    return diagram;
}

Diagram Term::symmetrize() const
{
    vector<Line> inds(indices());
    sort(inds.begin(), inds.end());
    vector< vector<Line> > particles;
    vector<Line> occ(inds);
    vector<Line> vrt(inds);

    filter(vrt, isType<PARTICLE+EXTERNAL>());
    filter(occ, isType<HOLE+EXTERNAL>());

    for (int i = 0;i < vrt.size();i++) particles.push_back(slice(vrt, i, i+1) + slice(occ, i, i+1));

    Diagram diagram(type, vector<Term>(1, *this));
    for (int i = 0;i < particles.size();i++)
    {
        Diagram tmp(type);
        for (int j = i+1;j < particles.size();j++)
        {
            tmp += Diagram(diagram).translate(particles[i]+particles[j], particles[j]+particles[i]);
        }
        diagram += tmp;
    }

    vector< vector<Line> > alpha;
    vector< vector<Line> > beta;
    for (int i = 0;i < vrt.size();i++)
    {
        if (vrt[i].isAlpha())
        {
            alpha.push_back(vector<Line>(1,vrt[i]));
        }
        else
        {
            beta.push_back(vector<Line>(1,vrt[i]));
        }
    }

    diagram.antisymmetrize(alpha).antisymmetrize(beta);

    inds = diagram[0].indices();
    vector<Line> alphas(inds);
    filter(alphas, isAlpha());
    vector<Line> betas(alphas);
    for (int i = 0;i < betas.size();i++) betas[i] = betas[i].toBeta();

    return diagram.translate(alphas, betas);
}

Diagram Term::antisymmetrize(const vector< vector<Line> >& _assym) const
{
    vector< vector<Line> > assym(_assym);
    Diagram diagram(type, vector<Term>(1, *this));

    for (int k = assym.size()-1;k > 0;k--)
    {
        Diagram tmp(type);
        for (int n = 1;n <= min(assym[k-1].size(), assym[k].size());n++)
        {
            string p1 = string(n, '0') + string(assym[k-1].size()-n, '1');
            do
            {
                string p2 = string(n, '0') + string(assym[k].size()-n, '1');
                do
                {
                    vector<Line> i1, i2;

                    for (int i = 0;i < p1.size();i++) if (p1[i] == '0') i1 += assym[k-1][i];
                    for (int i = 0;i < p2.size();i++) if (p2[i] == '0') i2 += assym[k][i];

                    //cout << (i1+i2) << " -> " << (i2+i1) << endl;

                    if (n&1)
                    {
                        tmp += (-diagram).translate(i1+i2, i2+i1);
                    }
                    else
                    {
                        tmp += Diagram(diagram).translate(i1+i2, i2+i1);
                    }
                } while (next_permutation(p2.begin(), p2.end()));
            } while (next_permutation(p1.begin(), p1.end()));
        }
        //cout << "diagram: " << endl << diagram;
        //cout << "tmp: " << endl << tmp;
        diagram += tmp;
        assym[k-1] += assym[k];
        assym.pop_back();
    }

    return diagram;
}

void Term::canonicalize()
{
    vector<Fragment>::iterator it;
    //Term old(*this);

    for (it = fragments.begin();it != fragments.end();++it)
    {
        factor *= it->canonicalize(type);
    }

    sort(fragments.begin(), fragments.end());

    //cout << old << " -> " << *this << endl;
}

pair<Manifold,Manifold>& operator+=(pair<Manifold,Manifold>& p1, pair<Manifold,Manifold> p2)
{
    p1.first += p2.first;
    p1.second += p2.second;
    return p1;
}

pair<Manifold,Manifold> Term::getShape() const
{
    pair<Manifold,Manifold> shape;

    for (vector<Fragment>::const_iterator it = fragments.begin();it != fragments.end();++it)
    {
        shape += it->getShape();
    }

    return shape;
}

void Term::getShape(Manifold& left, Manifold& right) const
{
    pair<Manifold,Manifold> shape;

    for (vector<Fragment>::const_iterator it = fragments.begin();it != fragments.end();++it)
    {
        shape += it->getShape();
    }

    left = shape.first;
    right = shape.second;
}

}
}
