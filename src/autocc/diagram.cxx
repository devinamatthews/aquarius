#include "autocc.hpp"

namespace aquarius
{
namespace autocc
{

ostream& operator<<(ostream& out, const Diagram& d)
{
    for (vector<Term>::const_iterator it = d.getTerms().begin();it != d.getTerms().end();++it) out << *it << endl;
    return out;
}

Diagram::Diagram(const Type type, const vector<Term>& terms)
: type(type)
{
    for (vector<Term>::const_iterator t = terms.begin();t != terms.end();++t) *this += *t;
}

Diagram::Diagram(const Type type, const vector<string>& strings)
: type(type)
{
    for (vector<string>::const_iterator str = strings.begin();str != strings.end();++str)
    {
        Diagram temp(type);
        vector<string> assyms;
        vector<Fragment> fragments;
        Fraction factor;
        istringstream iss(*str);
        string token;

        while (iss >> token)
        {
            if (token[0] == 'p' || token[0] == 'P')
            {
                assyms.push_back(token);
            }
            else if (token == "-")
            {
                factor = Fraction(-1,1);
            }
            else if (isdigit(token[0]) ||
                     token.find('/') != string::npos ||
                     token[0] == '-')
            {
                if (factor == Fraction(-1,1))
                {
                    factor = -Fraction(token);
                }
                else
                {
                    factor = Fraction(token);
                }
            }
            else
            {
                fragments.push_back(Fragment(token));
            }
        }

        temp += Term(type, factor, fragments);

        vector<Line> inds(temp.terms[0].indices());
        for (vector<string>::iterator s = assyms.begin();s != assyms.end();++s)
        {
            /*
            for (int i = 0;i < s->size();i++)
            {
                if (isalpha((*s)[i]) && find(inds.begin(), inds.end(), Line((*s)[i])) == inds.end())
                {
                    if (isupper((*s)[i]))
                    {
                        replace(s->begin(), s->end(), (*s)[i], (char)tolower((*s)[i]));
                    }
                    else
                    {
                        replace(s->begin(), s->end(), (*s)[i], (char)toupper((*s)[i]));
                    }
                }
            }
            */

            int lp = s->find('(');
            int oldbar = lp;
            if (lp == string::npos) throw logic_error("malformed antisymmetrizer: " + *s);

            vector< vector<Line> > assym;
            for (int bar;(bar = s->find('|', oldbar+1)) != string::npos;oldbar = bar)
            {
                assym.push_back(Line::parse(s->substr(oldbar+1, bar-oldbar-1)));
            }

            int rp = s->find(')', oldbar+1);
            if (rp == string::npos) throw logic_error("malformed antisymmetrizer: " + *s);

            assym.push_back(Line::parse(s->substr(oldbar+1, rp-oldbar-1)));

            if (assym.size() == 1 && assym[0].size() == 2)
            {
                assym.push_back(vector<Line>(1, assym[0][1]));
                assym[0] = vector<Line>(1, assym[0][0]);
            }

            temp.antisymmetrize(assym);
        }

        *this += temp;
    }
}

Diagram& Diagram::convert(const Type type)
{
    switch (type)
    {
        case SPINORBITAL:
            throw logic_error("It it not possible to convert to SPINORBITAL");
            break;
        case UHF:
            if (this->type != SPINORBITAL) throw logic_error("It it only possible to convert to UHF from SPINORBITAL");
            this->type = UHF;
            suminternal();
            break;
        case RHF:
            throw logic_error("It it not possible to convert to RHF");
            break;
        case SKELETON:
            switch (this->type)
            {
                case SPINORBITAL:
                    this->type = UHF;
                    //cout << "a:\n" << *this;
                    suminternal();
                    //cout << "b:\n" << *this;
                    this->type = SKELETON;
                    expanduhf();
                    //cout << "c:\n" << *this;
                    break;
                case UHF:
                    this->type = SKELETON;
                    expanduhf();
                    break;
                case RHF:
                    this->type = SKELETON;
                    //cout << "d:\n" << *this;
                    expandrhf();
                    //cout << "e:\n" << *this;
                    symmetrize();
                    //cout << "f:\n" << *this;
                    break;
                case SKELETON:
                    break;
            }
            break;
    }

    return *this;
}

Diagram& Diagram::fixorder(bool all)
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term) *this += term->fixorder(all);
    return *this;
}

Diagram& Diagram::fixorder(const vector<Line>& which)
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term) *this += term->fixorder(which);
    return *this;
}

Diagram& Diagram::fixexternal()
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term) *this += term->fixexternal();
    return *this;
}

Diagram& Diagram::negate()
{
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term) term->negate();
    return *this;
}

Diagram Diagram::operator-() const
{
    return Diagram(*this).negate();
}

Diagram& Diagram::operator*=(const Diagram& other)
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator t1 = terms.begin();t1 != terms.end();++t1)
    {
        for (vector<Term>::const_iterator t2 = other.terms.begin();t2 != other.terms.end();++t2)
        {
            *this += (*t1)*(*t2);
        }
    }
    return *this;
}

Diagram& Diagram::operator*=(const Term& other)
{
    for (vector<Term>::iterator t1 = terms.begin();t1 != terms.end();++t1)
    {
        *t1 *= other;
    }
    sort(terms.begin(), terms.end());
    return *this;
}

Diagram& Diagram::operator*=(const Fraction& other)
{
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term) *term *= other;
    return *this;
}

Diagram& Diagram::operator*=(const int other)
{
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term) *term *= other;
    return *this;
}

Diagram Diagram::operator*(const Diagram& other) const
{
    return Diagram(*this) *= other;
}

Diagram Diagram::operator*(const Term& other) const
{
    return Diagram(*this) *= other;
}

Diagram Diagram::operator*(const Fraction& other) const
{
    return Diagram(*this) *= other;
}

Diagram Diagram::operator*(const int other) const
{
    return Diagram(*this) *= other;
}

Diagram& Diagram::operator+=(const Diagram& other)
{
    for (vector<Term>::const_iterator term = other.terms.begin();term != other.terms.end();++term) *this += *term;
    return *this;
}

Diagram& Diagram::operator+=(const Term& other)
{
    vector<Term>::iterator i = upper_bound(terms.begin(), terms.end(), other);

    if (i != terms.begin() && *(i-1) == other)
    {
        --i;
        i->factor += other.factor;
        if (i->factor == 0) terms.erase(i);
    }
    else
    {
        terms.insert(i, other);
    }

    return *this;
}

Diagram& Diagram::operator-=(const Diagram& other)
{
    for (vector<Term>::const_iterator term = other.terms.begin();term != other.terms.end();++term) *this -= *term;
    return *this;
}

Diagram& Diagram::operator-=(const Term& other)
{
    *this += -other;
    return *this;
}

Diagram Diagram::operator+(const Diagram& other) const
{
    return Diagram(*this) += other;
}

Diagram Diagram::operator-(const Diagram& other) const
{
    return Diagram(*this) -= other;
}

bool Diagram::operator==(const Diagram& other) const
{
    if (terms.size() != other.terms.size()) return false;

    for (int i = 0;i < terms.size();i++)
    {
        if (!(terms[i].factor == other.terms[i].factor) || !(terms[i] == other.terms[i])) return false;
    }

    return true;
}

const Term& Diagram::operator[](const int i) const
{
    return terms[i];
}

int Diagram::numTerms() const
{
    return terms.size();
}

const vector<Term>& Diagram::getTerms() const
{
    return terms;
}

Diagram& Diagram::sumall()
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term)
    {
        term->type = UHF;
        *this += term->sumall();
    }
    return *this;
}

Diagram& Diagram::sum(const vector<Line>& which)
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term)
    {
        term->type = UHF;
        *this += term->sum(which);
    }
    return *this;
}

Diagram& Diagram::suminternal()
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term)
    {
        term->type = UHF;
        *this += term->suminternal();
    }
    return *this;
}

Diagram& Diagram::translate(const vector<Line>& from, const vector<Line>& to)
{
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term)
    {
        term->translate(from, to);
    }
    sort(terms.begin(), terms.end());
    return *this;
}

Diagram& Diagram::expanduhf()
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term)
    {
        term->type = SKELETON;
        *this += term->expanduhf();
    }
    return *this;
}

Diagram& Diagram::expandrhf()
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term)
    {
        term->type = SKELETON;
        *this += term->expandrhf();
    }
    return *this;
}

Diagram& Diagram::symmetrize()
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term)
    {
        *this += term->symmetrize();
    }
    return *this;
}

Diagram& Diagram::antisymmetrize(const vector< vector<Line> >& assym)
{
    vector<Term> terms;
    terms.swap(this->terms);
    for (vector<Term>::iterator term = terms.begin();term != terms.end();++term)
    {
        *this += term->antisymmetrize(assym);
    }
    return *this;
}

}
}
