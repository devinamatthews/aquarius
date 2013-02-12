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

#include <cassert>
#include <algorithm>
#include <iostream>

#include "autocc.hpp"

using namespace std;

namespace aquarius
{
namespace autocc
{

void Operator::canonicalize(Term& t)
{
    vector<Line> oldinds(t.indices());
    vector<Line> newinds(oldinds);

    int ep = 0;
    int eh = 0;
    int ip = 0;
    int ih = 0;
    for (vector<Line>::iterator l = newinds.begin();l != newinds.end();++l)
    {
        int n = 0;
        for (vector<Fragment>::const_iterator f = t.getFragments().begin();f != t.getFragments().end();++f)
        {
            vector<Line> finds(f->indices());
            if (find(finds.begin(), finds.end(), *l) != finds.end()) n++;
        }

        if (n > 1)
        {
            if (l->isParticle())
            {
                *l = Line(ip++, PARTICLE+INTERNAL+ALPHA);
            }
            else
            {
                *l = Line(ih++, HOLE+INTERNAL+ALPHA);
            }
        }
        else
        {
            if (l->isParticle())
            {
                *l = Line(ep++, PARTICLE+EXTERNAL+ALPHA);
            }
            else
            {
                *l = Line(eh++, HOLE+EXTERNAL+ALPHA);
            }
        }
    }

    t.translate(oldinds, newinds);
}

Operator::Operator()
: ref_count(0)
{
    //cout << "constructing " << this << endl;
}

void Operator::increment()
{
    ref_count++;
    //cout << "incrementing " << this << " to " << ref_count << endl;
}

void Operator::decrement()
{
    assert(ref_count > 0);
    ref_count--;
    //cout << "decrementing " << this << " to " << ref_count << endl;
    if (ref_count == 0) delete this;
}

Operator::~Operator()
{
    //cout << "destructing " << this << endl;
}

OperatorProduct& Operator::operator*(Operator& other)
{
    return *(new OperatorProduct(*this, other));
}

OperatorSum& Operator::operator+(Operator& other)
{
    return *(new OperatorSum(*this, other));
}

BasicOperator::BasicOperator()
{
    ref_count = 1;
}

BasicOperator::BasicOperator(const Diagram& diagram)
: terms(diagram.getTerms())
{
    ref_count = 1;
    canonicalize();
}

BasicOperator::BasicOperator(const string& term)
: terms(1, Term(Diagram::SPINORBITAL, term))
{
    ref_count = 1;
    canonicalize();
}

BasicOperator::BasicOperator(const Term& term)
: terms(1, term)
{
    ref_count = 1;
    canonicalize();
}

BasicOperator::BasicOperator(const vector<Term>& terms)
: terms(terms)
{
    ref_count = 1;
    canonicalize();
}

void BasicOperator::canonicalize()
{
    for (vector<Term>::iterator t = terms.begin();t != terms.end();++t)
    {
        Operator::canonicalize(*t);
    }
}

Diagram BasicOperator::resolve(const Manifold& left, const Manifold& right) const
{
    Diagram matches(Diagram::SPINORBITAL);

    for (vector<Term>::const_iterator it = terms.begin();it != terms.end();++it)
    {
        pair<Manifold,Manifold> shape(it->getShape());
        if (shape.first == left && shape.second == right)
        {
            matches += *it;
        }
    }

    return matches;
}

ManifoldGenerator* BasicOperator::matching(const Manifold& leftMin, const Manifold& leftMax,
                                           const Manifold& rightMin, const Manifold& rightMax) const
{
    return new BasicManifoldGenerator(leftMin, leftMax, rightMin, rightMax, *this);
}

ComplexOperator::ComplexOperator(Operator& other)
: operand(&other)
{
    ref_count = 1;
    operand->increment();
}

ComplexOperator& ComplexOperator::operator=(const ComplexOperator& other)
{
    operand->decrement();
    operand = other.operand;
    operand->increment();
    return *this;
}

ComplexOperator::~ComplexOperator()
{
    operand->decrement();
}

Diagram ComplexOperator::resolve(const Manifold& left, const Manifold& right) const
{
    return operand->resolve(left, right);
}

ManifoldGenerator* ComplexOperator::matching(const Manifold& leftMin, const Manifold& leftMax,
                                             const Manifold& rightMin, const Manifold& rightMax) const
{
    return operand->matching(leftMin, leftMax, rightMin, rightMax);
}

ExponentialOperator::ExponentialOperator(const ExponentialOperator& other)
: expansion(other.expansion)
{
    expansion->increment();
}

ExponentialOperator::ExponentialOperator(Operator& op, int order)
{
    expansion = new BasicOperator("1");

    Operator* prod = NULL;
    for (int i = 1;i <= order;i++)
    {
        if (prod == NULL)
            prod = &op;
        else
            prod = &((*prod)*op);
        expansion = &((*expansion)+(*prod));
    }

    expansion->increment();
}

ExponentialOperator& ExponentialOperator::operator=(const ExponentialOperator& other)
{
    expansion->decrement();
    expansion = other.expansion;
    expansion->increment();
    return *this;
}

ExponentialOperator::~ExponentialOperator()
{
    expansion->decrement();
}

Diagram ExponentialOperator::resolve(const Manifold& left, const Manifold& right) const
{
    return expansion->resolve(left, right);
}

ManifoldGenerator* ExponentialOperator::matching(const Manifold& leftMin, const Manifold& leftMax,
                                                 const Manifold& rightMin, const Manifold& rightMax) const
{
    return expansion->matching(leftMin, leftMax, rightMin, rightMax);
}

OperatorProduct::OperatorProduct(const OperatorProduct& other)
: l(other.l), r(other.r), flags(other.flags)
{
    l->increment();
    r->increment();
}

OperatorProduct::OperatorProduct(Operator& _l, Operator& _r, int flags)
: l(&_l), r(&_r), flags(flags)
{
    l->increment();
    r->increment();
}

OperatorProduct& OperatorProduct::operator=(const OperatorProduct& other)
{
    l->decrement();
    r->decrement();
    l = other.l;
    r = other.r;
    l->increment();
    r->increment();
    return *this;
}

OperatorProduct::~OperatorProduct()
{
    l->decrement();
    r->decrement();
}

bool OperatorProduct::isConnected(const Term& t)
{
    const vector<Fragment>& f(t.getFragments());

    if (f.size() == 1) return true;

    for (vector<Fragment>::const_iterator i1 = f.begin();i1 != f.end();++i1)
    {
        bool found = false;
        vector<Line> inds1(i1->indices());
        for (vector<Line>::iterator l1 = inds1.begin();l1 != inds1.end();++l1)
        {
            if (l1->isInternal())
            {
                for (vector<Fragment>::const_iterator i2 = f.begin();i2 != f.end();++i2)
                {
                    if (i1 == i2) continue;
                    vector<Line> inds2(i2->indices());
                    found |= find(inds2.begin(), inds2.end(), *l1) != inds2.end();
                }
            }
        }
        if (!found) return false;
    }

    return true;
}

bool OperatorProduct::isOpen(const Term& t)
{
    vector<Line> inds(t.indices());
    for (vector<Line>::iterator l = inds.begin();l != inds.end();++l) if (l->isExternal()) return true;
    return false;
}

Diagram OperatorProduct::resolve(const Manifold& left, const Manifold& right) const
{
    Diagram res(Diagram::SPINORBITAL);
    ManifoldGenerator *lGen, *rGen;
    Manifold ll, lr, rl, rr, c;

    /*
     * Find all manifold pairs for l which don't overrun left.
     * The loop over lr may be cut short if it is determined that no higher values can produce
     * a valid contraction.
     */
    lGen = l->matching(Manifold::MIN_VALUE, left, Manifold::MIN_VALUE, Manifold::MAX_VALUE);
    while (lGen->next(ll, lr))
    {
        /*
         * For this value of lr, loop over all possible ways to contract.
         */
        for (c.np = max(0,lr.np-right.np);c.np <= lr.np;c.np++)
        {
            for (c.nh = max(0,lr.nh-right.nh);c.nh <= lr.nh;c.nh++)
            {
                rGen = r->matching(left-ll+c, left-ll+c, right-lr+c, right-lr+c);
                while (rGen->next(rl, rr))
                {
                    if (ll == 0 && lr == c &&
                        rr == 0 && rl == c)
                    {
                        if (!(flags&CLOSED)) continue;
                    }
                    else
                    {
                        if (!(flags&OPEN)) continue;
                    }

                    if (c == 0 && !(ll == 0 && lr == 0) && !(rl == 0 && rr == 0))
                    {
                        if (!(flags&DISCONNECTED)) continue;
                    }
                    else
                    {
                        if (!(flags&CONNECTED)) continue;
                    }

                    vector<Term> lFragments(l->resolve(ll, lr).getTerms());
                    vector<Term> rFragments(r->resolve(rl, rr).getTerms());

                    for (vector<Term>::iterator lTerm = lFragments.begin();lTerm != lFragments.end();++lTerm)
                    {
                        for (vector<Term>::iterator rTerm = rFragments.begin();rTerm != rFragments.end();++rTerm)
                        {
                            doProduct(res, *lTerm, *rTerm, c);
                        }
                    }
                }
                delete rGen;
            }
        }
    }
    delete lGen;

    return res;
}

OperatorProduct& OperatorProduct::operator()(int flags)
{
    this->flags = (flags&0x3 ? flags&0x3 : this->flags&0x3) |
                  (flags&0xC ? flags&0xC : this->flags&0xC);
    return *this;
}

void OperatorProduct::doProduct(Diagram& res, Term l, Term r, const Manifold& c) const
{
    vector<Line> which_p_lr(c.np);
    vector<Line> which_h_lr(c.nh);
    vector<Line> which_p_rl(c.np);
    vector<Line> which_h_rl(c.nh);

    /*
     * Increase indices in r so that it does not conflict with l.
     */
    int maxep = -1;
    int maxeh = -1;
    int maxip = -1;
    int maxih = -1;
    vector<Line> linds = l.indices();
    vector<Line> oldrinds(linds);
    for (vector<Line>::iterator i = linds.begin();i != linds.end();++i)
    {
        if (i->isExternal())
        {
            if (i->isParticle())
                *i = Line(++maxep, EXTERNAL+PARTICLE+ALPHA);
            else
                *i = Line(++maxeh, EXTERNAL+HOLE+ALPHA);
        }
        else
        {
            if (i->isParticle())
                *i = Line(++maxip, INTERNAL+PARTICLE+ALPHA);
            else
                *i = Line(++maxih, INTERNAL+HOLE+ALPHA);
        }
    }

    l.translate(oldrinds, linds);

    vector<Line> rinds = r.indices();
    oldrinds = rinds;
    for (vector<Line>::iterator i = rinds.begin();i != rinds.end();++i)
    {
        if (i->isExternal())
        {
            if (i->isParticle())
                *i = Line(++maxep, EXTERNAL+PARTICLE+ALPHA);
            else
                *i = Line(++maxeh, EXTERNAL+HOLE+ALPHA);
        }
        else
        {
            if (i->isParticle())
                *i = Line(++maxip, INTERNAL+PARTICLE+ALPHA);
            else
                *i = Line(++maxih, INTERNAL+HOLE+ALPHA);
        }
    }

    r.translate(oldrinds, rinds);

    vector<Line> c_p(c.np);
    vector<Line> c_h(c.nh);
    for (vector<Line>::iterator i = c_p.begin();i != c_p.end();++i) *i = Line(++maxip, INTERNAL+PARTICLE+ALPHA);
    for (vector<Line>::iterator i = c_h.begin();i != c_h.end();++i) *i = Line(++maxih, INTERNAL+HOLE+ALPHA);

    ProductGenerator which_lr(c, l, ProductGenerator::RIGHT);
    while (which_lr.next(which_p_lr.data(), which_h_lr.data()))
    {
        ProductGenerator which_rl(c, r, ProductGenerator::LEFT);
        while (which_rl.next(which_p_rl.data(), which_h_rl.data()))
        {
            Term tmp = l*r;
            tmp.translate(which_p_lr, c_p);
            tmp.translate(which_h_lr, c_h);
            tmp.translate(which_p_rl, c_p);
            tmp.translate(which_h_rl, c_h);
            canonicalize(tmp);

            if (isConnected(tmp))
            {
                if (!(flags&CONNECTED)) continue;
            }
            else
            {
                if (!(flags&DISCONNECTED)) continue;
            }

            if (isOpen(tmp))
            {
                if (!(flags&OPEN)) continue;
            }
            else
            {
                if (!(flags&CLOSED)) continue;
            }

            tmp.fixexternal().fixorder();
            tmp *= tmp.getFactor().inverse();

            res += tmp;
        }
    }
}

ManifoldGenerator* OperatorProduct::matching(const Manifold& leftMin, const Manifold& leftMax,
                                             const Manifold& rightMin, const Manifold& rightMax) const
{
    return new ProductManifoldGenerator(leftMin, leftMax, rightMin, rightMax, *this);
}

OperatorProduct::ProductGenerator::ProductGenerator(const Manifold& c, const Term& term, const Side side)
: c(c), term(term), side(side)
{
    vector<Fragment> fragments(term.getFragments());

    for (vector<Fragment>::iterator f = fragments.begin();f != fragments.end();++f)
    {
        vector<Line> finds;

        if (side == LEFT)
            pindices.push_back(filter_copy(f->getIndicesOut(), isParticle()));
        else
            pindices.push_back(filter_copy(f->getIndicesIn(), isParticle()));

        how_many_p.push_back(0);
    }

    for (vector<Fragment>::iterator f = fragments.begin();f != fragments.end();++f)
    {
        vector<Line> finds;

        if (side == LEFT)
            hindices.push_back(filter_copy(f->getIndicesIn(), isHole()));
        else
            hindices.push_back(filter_copy(f->getIndicesOut(), isHole()));

        how_many_h.push_back(0);
    }

    which_bin_p.resize(c.np, 0);
    which_bin_h.resize(c.nh, 0);
}

bool OperatorProduct::ProductGenerator::generator_next(Line which_p[], Line which_h[])
{
    pdone = false;
    while (!pdone)
    {
        how_many_p.assign(how_many_p.size(), 0);
        for (pi = 0;pi < c.np;pi++)
            how_many_p[which_bin_p[pi]]++;

        hdone = false;
        while (!hdone)
        {
            how_many_h.assign(how_many_h.size(), 0);
            for (hi = 0;hi < c.nh;hi++)
                how_many_h[which_bin_h[hi]]++;

            valid = true;
            idx = 0;
            for (pi = 0;pi < how_many_p.size();pi++)
            {
                if (how_many_p[pi] > pindices[pi].size())
                {
                    valid = false;
                    break;
                }
                copy(pindices[pi].begin(), pindices[pi].begin()+how_many_p[pi], which_p+idx);
                idx += how_many_p[pi];
            }

            idx = 0;
            for (hi = 0;hi < how_many_h.size();hi++)
            {
                if (how_many_h[hi] > hindices[hi].size())
                {
                    valid = false;
                    break;
                }
                copy(hindices[hi].begin(), hindices[hi].begin()+how_many_h[hi], which_h+idx);
                idx += how_many_h[hi];
            }

            if (valid)
            {
                generator_yield();
            }

            for (hi = 0;hi < c.nh;hi++)
            {
                which_bin_h[hi]++;
                if (which_bin_h[hi] >= how_many_h.size())
                {
                    which_bin_h[hi] = 0;
                    if (hi == c.nh-1) hdone = true;
                }
                else
                {
                    break;
                }
            }
            if (c.nh == 0) hdone = true;
        }

        for (pi = 0;pi < c.np;pi++)
        {
            which_bin_p[pi]++;
            if (which_bin_p[pi] >= how_many_p.size())
            {
                which_bin_p[pi] = 0;
                if (pi == c.np-1) pdone = true;
            }
            else
            {
                break;
            }
        }
        if (c.np == 0) pdone = true;
    }

    generator_stop;
}

OperatorSum::OperatorSum(const OperatorSum& other)
: l(other.l), r(other.r)
{
    l->increment();
    r->increment();
}

OperatorSum::OperatorSum(Operator& _l, Operator& _r)
: l(&_l), r(&_r)
{
    l->increment();
    r->increment();
}

OperatorSum& OperatorSum::operator=(const OperatorSum& other)
{
    l->decrement();
    r->decrement();
    l = other.l;
    r = other.r;
    l->increment();
    r->increment();
    return *this;
}

OperatorSum::~OperatorSum()
{
    l->decrement();
    r->decrement();
}

Diagram OperatorSum::resolve(const Manifold& left, const Manifold& right) const
{
    return l->resolve(left, right) + r->resolve(left, right);
}

ManifoldGenerator* OperatorSum::matching(const Manifold& leftMin, const Manifold& leftMax,
                                         const Manifold& rightMin, const Manifold& rightMax) const
{
    return new SumManifoldGenerator(leftMin, leftMax, rightMin, rightMax, *this);
}

}
}
