#include "autocc.hpp"

namespace aquarius
{
namespace autocc
{

ostream& operator<<(ostream& out, const Term& t)
{
    out << t.factor;
    for (auto& f : t.fragments) out << " " << f;
    return out;
}

Term::Term(const Diagram::Type type, const string s)
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
    int idx = 0;
    vector<Line> inds = indices();
    for (vector<Line>::iterator l = inds.begin();l != inds.end();++l)
    {
        if (find(which.begin(), which.end(), *l) == which.end())
        {
            idx = max(idx, l->getIndex());
        }
    }

    vector<Line> oldlines, newlines;

    for (vector<Fragment>::iterator f = getFragments().begin();f != getFragments().end();++f)
    {
        vector<Line> tofix = intersection(which, filter_copy(f->indices(), and1(isAlpha(),isVirtual())));

        for (vector<Line>::iterator l = tofix.begin();l != tofix.end();++l)
        {
            oldlines += *l;
            newlines += l->asIndex(++idx);
        }
    }

    for (vector<Fragment>::iterator f = getFragments().begin();f != getFragments().end();++f)
    {
        vector<Line> tofix = intersection(which, filter_copy(f->indices(), and1(isAlpha(),isOccupied())));

        for (vector<Line>::iterator l = tofix.begin();l != tofix.end();++l)
        {
            oldlines += *l;
            newlines += l->asIndex(++idx);
        }
    }

    for (vector<Fragment>::iterator f = getFragments().begin();f != getFragments().end();++f)
    {
        vector<Line> tofix = intersection(which, filter_copy(f->indices(), and1(isBeta(),isVirtual())));

        for (vector<Line>::iterator l = tofix.begin();l != tofix.end();++l)
        {
            oldlines += *l;
            newlines += l->asIndex(++idx);
        }
    }

    for (vector<Fragment>::iterator f = getFragments().begin();f != getFragments().end();++f)
    {
        vector<Line> tofix = intersection(which, filter_copy(f->indices(), and1(isBeta(),isOccupied())));

        for (vector<Line>::iterator l = tofix.begin();l != tofix.end();++l)
        {
            oldlines += *l;
            newlines += l->asIndex(++idx);
        }
    }

    translate(oldlines, newlines);

    return *this;
}

Term& Term::fixorder(bool all)
{
    if (all)
    {
        fixorder(indices());
    }
    else
    {
        fixorder(internal());
    }

    return *this;
}

Term& Term::fixexternal()
{
    assert(0);

    /*
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
    */

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
    return sum(indices());
}

Diagram Term::suminternal() const
{
    return sum(internal());
}

Diagram Term::sum(const vector<Line>& which) const
{
    Diagram diagram(type);

    vector<Line> old = intersection(which, indices());

    vector<Line> uhf(old);
    for (int x = 0;x < 1<<old.size();x++)
    {
        for (int i = 0;i < old.size();i++)
        {
            if ((x>>i) & 1)
            {
                uhf[i] = old[i].asAlpha();
            }
            else
            {
                uhf[i] = old[i].asBeta();
            }
        }

        Term newterm(*this);
        newterm.translate(old, uhf);

        if (newterm.checkspin())
        {
            diagram += newterm;
        }
    }

    return diagram;
}

bool Term::checkspin() const
{
    for (vector<Fragment>::const_iterator it = fragments.begin();it != fragments.end();++it)
    {
        int nl = count_if(it->out.begin(), it->out.end(), isAlpha());
        int nr = count_if(it->in.begin(), it->in.end(), isAlpha());
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
        indices += f->indices();
    }

    return uniq_copy(indices);
}

vector<Line> Term::internal() const
{
    vector<Line> internal;

    for (vector<Fragment>::const_iterator f1 = getFragments().begin();f1 != getFragments().end();++f1)
    {
        for (vector<Fragment>::const_iterator f2 = f1+1;f2 != getFragments().end();++f2)
        {
            internal += intersection(f1->indices(), f2->indices());
        }

        internal += exclude_copy(f1->indices(), uniq_copy(f1->indices()));
    }

    return uniq_copy(internal);
}

vector<Line> Term::external() const
{
    return exclude_copy(indices(), internal());
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
        for (int i = 0;i < f.out.size();i++) f.out[i].toBeta();
        for (int i = 0;i < f.in.size();i++)
        {
            if (f.in[i].isAlpha())
            {
                f.in[i].toBeta();
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
    vector<Line> _internal = internal();
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
                if (find(_internal.begin(), _internal.end(), f.out[i]) != _internal.end())
                {
                    f.out[i].toBeta();
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
            if (f.in[i].isAlpha() && find(_internal.begin(), _internal.end(), f.out[i]) != _internal.end())
            {
                f.in[i].toBeta();
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
    vector<Line> _external = external();
    vector< vector<Line> > particles;
    vector< vector<Line> > alpha;
    vector< vector<Line> > beta;

    /*
     * Attempt to trace each external line to another external line.
     * This will either be possible in which case the pair of lines defines
     * a "particle", or not possible in the case that the term
     * contains non-particle-conserving operators.
     */
    for (vector<Line>::iterator l = _external.begin();;)
    {
        if (l == _external.end()) break;

        const Line *l1 = &(*l), *l2;

        do
        {
            l2 = NULL;
            bool found = false;

            for (vector<Fragment>::const_iterator f = fragments.begin();f != fragments.end();++f)
            {
                for (int i = 0;i < f->in.size();i++)
                {
                    if (f->in[i] == *l1)
                    {
                        if (i < f->out.size()) l2 = &f->out[i];
                        found = true;
                        break;
                    }
                }

                if (found) break;

                for (int i = 0;i < f->out.size();i++)
                {
                    if (f->out[i] == *l1)
                    {
                        if (i < f->in.size()) l2 = &f->in[i];
                        found = true;
                        break;
                    }
                }

                if (found) break;
            }

            swap(l1, l2);
        }
        while (l1 != NULL && find(_external.begin(), _external.end(), *l1) == _external.end());

        if (l1 != NULL)
        {
            _external.erase(find(_external.begin(), _external.end(), *l1));
            assert(l->getSpin() == l1->getSpin());
            particles += {*l, *l1};
            (l->isAlpha() ? alpha : beta) += {min(*l, *l1)};
        }

        l = _external.erase(l);
    }

    Diagram diagram(type, {*this});
    for (int i = 0;i < particles.size();i++)
    {
        Diagram tmp(type);
        for (int j = i+1;j < particles.size();j++)
        {
            tmp += Diagram(diagram).translate(particles[i]+particles[j], particles[j]+particles[i]);
        }
        diagram += tmp;
    }

    diagram.antisymmetrize(alpha).antisymmetrize(beta);

    vector<Line> alphas = filter_copy(indices(), isAlpha());
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
    getShape(shape.first, shape.second);
    return shape;
}

void Term::getShape(Manifold& left, Manifold& right) const
{
    vector<Line> ex = external();
    vector<Line> out, in;

    for (vector<Fragment>::const_iterator f = fragments.begin();f != fragments.end();++f)
    {
        operator+=(out, f->out);
        out += f->out;
        in += f->in;
    }

    vector<Line> exout = intersection(ex, uniq_copy(out));
    vector<Line> exin = intersection(ex, uniq_copy(in));

    for (vector<Line>::const_iterator i = exout.begin();i != exout.end();++i)
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

    for (vector<Line>::const_iterator i = exin.begin();i != exin.end();++i)
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
