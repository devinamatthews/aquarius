#include "autocc.hpp"

namespace aquarius
{
namespace autocc
{

void Operator::canonicalize(Term& t)
{
    vector<Line> oldinds(t.indices());
    vector<Line> newinds(oldinds);

    int p = 0;
    int h = 0;
    for (int pass = 0;pass < 2;pass++)
    for (vector<Line>::iterator l = newinds.begin();l != newinds.end();++l)
    {
        int n = 0;
        for (vector<Fragment>::const_iterator f = t.getFragments().begin();f != t.getFragments().end();++f)
        {
            vector<Line> finds(f->indices());
            if (find(finds.begin(), finds.end(), *l) != finds.end()) n++;
        }

        if (pass == 1 && n > 1)
        {
            if (l->isVirtual()) l->toIndex(p++);
            else                l->toIndex(h++);
        }
        else if (pass == 0)
        {
            if (l->isVirtual()) l->toIndex(p++);
            else                l->toIndex(h++);
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

Diagram BasicOperator::resolve(Manifold& left, Manifold& right) const
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

Diagram ComplexOperator::resolve(Manifold& left, Manifold& right) const
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

Diagram ExponentialOperator::resolve(Manifold& left, Manifold& right) const
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

    for (vector<Fragment>::const_iterator i1 = f.begin();i1 != f.end();++i1)
    {
        bool found = false;
        vector<Line> inds1(i1->indices());
        for (vector<Line>::iterator l1 = inds1.begin();l1 != inds1.end();++l1)
        {
            for (vector<Fragment>::const_iterator i2 = f.begin();i2 != f.end();++i2)
            {
                if (i1 == i2) continue;
                vector<Line> inds2(i2->indices());
                found |= find(inds2.begin(), inds2.end(), *l1) != inds2.end();
            }
        }
        if (!found) return false;
    }

    return true;
}

bool OperatorProduct::isOpen(const Term& t)
{
    const vector<Fragment>& f(t.getFragments());

    if (f.size() == 1) return true;

    for (vector<Fragment>::const_iterator i1 = f.begin();i1 != f.end();++i1)
    {
        bool found = true;
        vector<Line> inds1(i1->indices());
        for (vector<Line>::iterator l1 = inds1.begin();l1 != inds1.end();++l1)
        {
            bool found2 = false;
            for (vector<Fragment>::const_iterator i2 = f.begin();i2 != f.end();++i2)
            {
                if (i1 == i2) continue;
                vector<Line> inds2(i2->indices());
                found2 |= find(inds2.begin(), inds2.end(), *l1) != inds2.end();
            }
            found &= found2;
        }
        if (!found) return true;
    }

    return false;
}

Diagram OperatorProduct::resolve(Manifold& left, Manifold& right) const
{
    Diagram res(Diagram::SPINORBITAL);
    ManifoldGenerator *lGen, *rGen;
    Manifold ll, lr, rl, rr, c;

    Manifold _zero, _max;
    _max.np[0] = numeric_limits<int>::max();
    _max.nh[0] = numeric_limits<int>::max();

    /*
     * Find all manifold pairs for l which don't overrun left.
     * The loop over lr may be cut short if it is determined that no higher values can produce
     * a valid contraction.
     */
    lGen = l->matching(_zero, left, _zero, _max);
    while (lGen->next(ll, lr))
    {
        /*
         * For this value of lr, loop over all possible ways to contract.
         */
        for (c.np[0] = max(0,lr.np[0]-right.np[0]);c.np[0] <= lr.np[0];c.np[0]++)
        {
            for (c.nh[0] = max(0,lr.nh[0]-right.nh[0]);c.nh[0] <= lr.nh[0];c.nh[0]++)
            {
                rGen = r->matching(left-ll+c, left-ll+c, right-lr+c, right-lr+c);
                while (rGen->next(rl, rr))
                {
                    if (ll == _zero && lr == c &&
                        rr == _zero && rl == c)
                    {
                        if (!(flags&CLOSED)) continue;
                    }
                    else
                    {
                        if (!(flags&OPEN)) continue;
                    }

                    if (c == _zero && !(ll == _zero && lr == _zero) && !(rl == _zero && rr == _zero))
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

void OperatorProduct::doProduct(Diagram& res, Term l, Term r, Manifold& c) const
{
    vector<Line> which_p_lr(c.np[0]);
    vector<Line> which_h_lr(c.nh[0]);
    vector<Line> which_p_rl(c.np[0]);
    vector<Line> which_h_rl(c.nh[0]);

    /*
     * Increase indices in r so that it does not conflict with l.
     */
    int maxp = -1;
    int maxh = -1;
    vector<Line> linds = l.indices();
    vector<Line> oldrinds(linds);
    for (vector<Line>::iterator i = linds.begin();i != linds.end();++i)
    {
        if (i->isVirtual()) i->toIndex(++maxp);
        else                i->toIndex(++maxh);
    }

    l.translate(oldrinds, linds);

    vector<Line> rinds = r.indices();
    oldrinds = rinds;
    for (vector<Line>::iterator i = rinds.begin();i != rinds.end();++i)
    {
        if (i->isVirtual()) i->toIndex(++maxp);
        else                i->toIndex(++maxh);
    }

    r.translate(oldrinds, rinds);

    vector<Line> c_p(c.np[0]);
    vector<Line> c_h(c.nh[0]);
    for (vector<Line>::iterator i = c_p.begin();i != c_p.end();++i) *i = Line(++maxp, 0, Line::VIRTUAL, Line::ALPHA);
    for (vector<Line>::iterator i = c_h.begin();i != c_h.end();++i) *i = Line(++maxh, 0, Line::OCCUPIED, Line::ALPHA);

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

OperatorProduct::ProductGenerator::ProductGenerator(Manifold& c, const Term& term, const Side side)
: c(c), term(term), side(side)
{
    vector<Fragment> fragments(term.getFragments());

    for (vector<Fragment>::iterator f = fragments.begin();f != fragments.end();++f)
    {
        vector<Line> finds;

        if (side == LEFT)
            pindices.push_back(filter_copy(f->getIndicesOut(), isVirtual()));
        else
            pindices.push_back(filter_copy(f->getIndicesIn(), isVirtual()));

        how_many_p.push_back(0);
    }

    for (vector<Fragment>::iterator f = fragments.begin();f != fragments.end();++f)
    {
        vector<Line> finds;

        if (side == LEFT)
            hindices.push_back(filter_copy(f->getIndicesIn(), isOccupied()));
        else
            hindices.push_back(filter_copy(f->getIndicesOut(), isOccupied()));

        how_many_h.push_back(0);
    }

    which_bin_p.resize(c.np[0], 0);
    which_bin_h.resize(c.nh[0], 0);
}

bool OperatorProduct::ProductGenerator::generator_next(Line which_p[], Line which_h[])
{
    pdone = false;
    while (!pdone)
    {
        how_many_p.assign(how_many_p.size(), 0);
        for (pi = 0;pi < c.np[0];pi++)
            how_many_p[which_bin_p[pi]]++;

        hdone = false;
        while (!hdone)
        {
            how_many_h.assign(how_many_h.size(), 0);
            for (hi = 0;hi < c.nh[0];hi++)
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

            for (hi = 0;hi < c.nh[0];hi++)
            {
                which_bin_h[hi]++;
                if (which_bin_h[hi] >= how_many_h.size())
                {
                    which_bin_h[hi] = 0;
                    if (hi == c.nh[0]-1) hdone = true;
                }
                else
                {
                    break;
                }
            }
            if (c.nh[0] == 0) hdone = true;
        }

        for (pi = 0;pi < c.np[0];pi++)
        {
            which_bin_p[pi]++;
            if (which_bin_p[pi] >= how_many_p.size())
            {
                which_bin_p[pi] = 0;
                if (pi == c.np[0]-1) pdone = true;
            }
            else
            {
                break;
            }
        }
        if (c.np[0] == 0) pdone = true;
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

Diagram OperatorSum::resolve(Manifold& left, Manifold& right) const
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
