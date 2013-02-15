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
#include <iostream>

#include "autocc.hpp"

using namespace std;

namespace aquarius
{
namespace autocc
{

bool operator<(const pair<Manifold,Manifold>& p1, const pair<Manifold,Manifold>& p2)
{
    if (p1.first < p2.first)
    {
        return true;
    }
    else if (p2.first < p1.first)
    {
        return false;
    }
    else
    {
        return p1.second < p1.first;
    }
}

bool BasicManifoldGenerator::generator_next(Manifold& left, Manifold& right)
{
    choices.resize(op.terms.size());

    for (int i = 0;i < op.terms.size();i++)
    {
        choices[i] = op.terms[i].getShape();
    }
    sort(choices.begin(), choices.end());

    for (it = choices.begin();it != choices.end();++it)
    {
        if (it->first.np  >= leftMin.np  && it->first.nh  >= leftMin.nh  &&
            it->first.np  <= leftMax.np  && it->first.nh  <= leftMax.nh  &&
            it->second.np >= rightMin.np && it->second.nh >= rightMin.nh &&
            it->second.np <= rightMax.np && it->second.nh <= rightMax.nh)
            generator_yield(left=it->first, right=it->second);
    }

    choices.clear();

    generator_stop;
}

bool ProductManifoldGenerator::generator_next(Manifold& left, Manifold& right)
{
    _ldone = false;
    for (_l.np = leftMin.np;_l.np <= leftMax.np && !_ldone;_l.np++)
    {
        _lhdone = false;
        for (_l.nh = leftMin.nh;_l.nh <= leftMax.nh && !_ldone && !_lhdone;_l.nh++)
        {
            _ldone = true;

            _rdone = false;
            for (_r.np = rightMin.np;_r.np <= rightMax.np && !_rdone;_r.np++)
            {
                _rhdone = false;
                for (_r.nh = rightMin.nh;_r.nh <= rightMax.nh && !_rdone && !_rhdone;_r.nh++)
                {
                    if ((!op.flags&CLOSED) &&  (_l == 0 && _r == 0)) continue;
                    if ((!op.flags&OPEN)   && !(_l == 0 && _r == 0))
                    {
                        _ldone = true;
                        _rdone = true;
                        continue;
                    }

                    //cout << "shape: " << _l << ' ' << _r << endl;

                    _rdone = true;

                    cdone = false;
                    for (c.np = 0;c.np <= Manifold::MAX_VALUE.np && !cdone;c.np++)
                    {
                        chdone = false;
                        for (c.nh = 0;c.nh <= Manifold::MAX_VALUE.nh && !cdone && !chdone;c.nh++)
                        {
                            if ((!op.flags&CONNECTED) && !(c == 0))
                            {
                                cdone = true;
                                continue;
                            }

                            cdone = true;

                            lGen = op.l->matching(Manifold::MIN_VALUE, _l, c, c+_r);
                            while (lGen->next(ll,lr))
                            {
                                rGen = op.r->matching(_l+c-ll, _l+c-ll, _r+c-lr, _r+c-lr);
                                if (rGen->next(rl,rr))
                                {
                                    _ldone = false;
                                    _rdone = false;
                                    cdone = false;

                                    if ((op.flags&DISCONNECTED) || !(c == 0) ||
                                        (ll == 0 && lr == 0) || (rl == 0 && rr == 0))
                                    {
                                        generator_yield(left=_l, right=_r);
                                    }
                                }
                                delete rGen;
                            }
                            delete lGen;

                            chdone = cdone;
                            if (cdone)
                            {
                                lGen = op.l->matching(Manifold::MIN_VALUE,                  _l,
                                                          c+Manifold(0,1), Manifold::MAX_VALUE);
                                rGen = op.r->matching(    c+Manifold(0,1), Manifold::MAX_VALUE,
                                                      Manifold::MIN_VALUE,                  _r);
                                if (lGen->next(ll,lr) && rGen->next(rl,rr))
                                {
                                    cdone = false;
                                    chdone = false;
                                }
                                else
                                {
                                    delete lGen;
                                    delete rGen;
                                    lGen = op.l->matching(Manifold::MIN_VALUE,                  _l,
                                                           Manifold(c.np+1,0), Manifold::MAX_VALUE);
                                    rGen = op.r->matching( Manifold(c.np+1,0), Manifold::MAX_VALUE,
                                                          Manifold::MIN_VALUE,                  _r);
                                    if (lGen->next(ll,lr) && rGen->next(rl,rr))
                                    {
                                        cdone = false;
                                    }
                                }
                                delete lGen;
                                delete rGen;
                            }
                        }
                    }

                    //cout << "rdone: " << _rdone << ' ';

                    _rhdone = _rdone;
                    if (_rdone)
                    {
                        lGen = op.l->matching(Manifold::MIN_VALUE,                  _l,
                                                 _r+Manifold(0,1), Manifold::MAX_VALUE);
                        if (lGen->next(ll,lr))
                        {
                            _rdone = false;
                            _rhdone = false;
                        }
                        else
                        {
                            delete lGen;
                            lGen = op.l->matching(Manifold::MIN_VALUE,                  _l,
                                                  Manifold(_r.np+1,0), Manifold::MAX_VALUE);
                            if (lGen->next(ll,lr))
                            {
                                _rdone = false;
                            }
                        }
                        delete lGen;
                    }

                    //cout << _rdone << endl;
                }
            }

            //cout << "ldone: " << _ldone << ' ';

            _lhdone = _ldone;
            if (_ldone)
            {
                for (__l.np = 0;__l.np <= _l.np && _ldone && _lhdone;__l.np++)
                {
                    for (__l.nh = 0;__l.nh <= _l.nh && _ldone && _lhdone;__l.nh++)
                    {
                        lGen = op.l->matching(                __l, Manifold::MAX_VALUE,
                                              Manifold::MIN_VALUE, Manifold::MAX_VALUE);
                        if (lGen->next(ll,lr))
                        {
                            //cout << "found1 " << ll << ' ' << lr << endl;

                            if (ll.np > _l.np && ll.nh > _l.nh)
                            {
                                _ldone = false;
                                _lhdone = false;
                            }
                            else
                            {
                                rGen = op.r->matching(Manifold(max(0,_l.np-ll.np),max(0,_l.nh+1-ll.nh)), Manifold::MAX_VALUE,
                                                                                  Manifold::MIN_VALUE, Manifold::MAX_VALUE);
                                if (rGen->next(rl,rr))
                                {
                                    //cout << "found2 " << rl << ' ' << rr << endl;
                                    _ldone = false;
                                    _lhdone = false;
                                }
                                delete rGen;
                            }
                        }
                        delete lGen;

                        lGen = op.l->matching( Manifold(__l.np,0), Manifold::MAX_VALUE,
                                              Manifold::MIN_VALUE, Manifold::MAX_VALUE);
                        if (lGen->next(ll,lr))
                        {
                            //cout << "found3 " << ll << ' ' << lr << endl;
                            if (ll.np > _l.np)
                            {
                                _ldone = false;
                            }
                            else
                            {
                                rGen = op.r->matching(Manifold(max(0,_l.np+1-ll.np),0), Manifold::MAX_VALUE,
                                                                 Manifold::MIN_VALUE, Manifold::MAX_VALUE);
                                if (rGen->next(rl,rr))
                                {
                                    //cout << "found4 " << rl << ' ' << rr << endl;
                                    _ldone = false;
                                }
                                delete rGen;
                            }
                        }
                        delete lGen;
                    }
                }
            }

            //cout << _ldone << endl;
        }
    }

    generator_stop;
}

bool SumManifoldGenerator::generator_next(Manifold& left, Manifold& right)
{
    lGenerator = op.l->matching(leftMin, leftMax, rightMin, rightMax);
    rGenerator = op.r->matching(leftMin, leftMax, rightMin, rightMax);
    lHasMore = lGenerator->next(lNext.first, lNext.second);
    rHasMore = rGenerator->next(rNext.first, rNext.second);

    while (lHasMore || rHasMore)
    {
        if (lHasMore)
        {
            if (rHasMore)
            {
                if (lNext.first < rNext.first || (lNext.first == rNext.first && lNext.second < rNext.second))
                {
                    /*
                     * l has the lower manifold pair
                     */
                    generator_yield(left=lNext.first, right=lNext.second);
                    lHasMore = lGenerator->next(lNext.first, lNext.second);
                }
                else if (!(lNext.first == rNext.first && lNext.second == rNext.second))
                {
                    /*
                     * r has the lower manifold pair
                     */
                    generator_yield(left=rNext.first, right=rNext.second);
                    rHasMore = rGenerator->next(rNext.first, rNext.second);
                }
                else
                {
                    /*
                     * l and r have equal manifold pairs, pick either one and discard the other
                     */
                    generator_yield(left=lNext.first, right=lNext.second);
                    lHasMore = lGenerator->next(lNext.first, lNext.second);
                    rHasMore = rGenerator->next(rNext.first, rNext.second);
                }
            }
            else
            {
                generator_yield(left=lNext.first, right=lNext.second);
                lHasMore = lGenerator->next(lNext.first, lNext.second);
            }
        }
        else
        {
            generator_yield(left=rNext.first, right=rNext.second);
            rHasMore = rGenerator->next(rNext.first, rNext.second);
        }
    }

    delete lGenerator;
    delete rGenerator;

    generator_stop;
}

}
}
