#include "autocc.hpp"

namespace aquarius
{
namespace autocc
{

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
        if (it->first.np[0]  >= leftMin.np[0]  && it->first.nh[0]  >= leftMin.nh[0]  &&
            it->first.np[0]  <= leftMax.np[0]  && it->first.nh[0]  <= leftMax.nh[0]  &&
            it->second.np[0] >= rightMin.np[0] && it->second.nh[0] >= rightMin.nh[0] &&
            it->second.np[0] <= rightMax.np[0] && it->second.nh[0] <= rightMax.nh[0])
            generator_yield(left=it->first, right=it->second);
    }

    choices.clear();

    generator_stop;
}

bool ProductManifoldGenerator::generator_next(Manifold& left, Manifold& right)
{
    _one.nh[0] = 1;
    _max.np[0] = numeric_limits<int>::max()-1;
    _max.nh[0] = numeric_limits<int>::max()-1;

    _ldone = false;
    for (_l.np[0] = leftMin.np[0];_l.np[0] <= leftMax.np[0] && !_ldone;_l.np[0]++)
    {
        _lhdone = false;
        for (_l.nh[0] = leftMin.nh[0];_l.nh[0] <= leftMax.nh[0] && !_ldone && !_lhdone;_l.nh[0]++)
        {
            _ldone = true;

            _rdone = false;
            for (_r.np[0] = rightMin.np[0];_r.np[0] <= rightMax.np[0] && !_rdone;_r.np[0]++)
            {
                _rhdone = false;
                for (_r.nh[0] = rightMin.nh[0];_r.nh[0] <= rightMax.nh[0] && !_rdone && !_rhdone;_r.nh[0]++)
                {
                    if ((!op.flags&CLOSED) &&  (_l.np[0]+_l.nh[0] == 0 && _r.np[0]+_r.nh[0] == 0)) continue;
                    if ((!op.flags&OPEN)   && !(_l.np[0]+_l.nh[0] == 0 && _r.np[0]+_r.nh[0] == 0))
                    {
                        _ldone = true;
                        _rdone = true;
                        continue;
                    }

                    //cout << "shape: " << _l << ' ' << _r << endl;

                    _rdone = true;

                    cdone = false;
                    for (c.np[0] = 0;c.np[0] < numeric_limits<int>::max() && !cdone;c.np[0]++)
                    {
                        chdone = false;
                        for (c.nh[0] = 0;c.nh[0] < numeric_limits<int>::max() && !cdone && !chdone;c.nh[0]++)
                        {
                            if ((!op.flags&CONNECTED) && !(c.np[0]+c.nh[0] == 0))
                            {
                                cdone = true;
                                continue;
                            }

                            cdone = true;

                            lGen = op.l->matching(_zero, _l, c, c+_r);
                            while (lGen->next(ll,lr))
                            {
                                rGen = op.r->matching(_l+c-ll, _l+c-ll, _r+c-lr, _r+c-lr);
                                if (rGen->next(rl,rr))
                                {
                                    _ldone = false;
                                    _rdone = false;
                                    cdone = false;

                                    if ((op.flags&DISCONNECTED) ||
                                       !(c.np[0]+c.nh[0] == 0) ||
                                        (ll.np[0]+ll.nh[0] == 0 && lr.np[0]+lr.nh[0] == 0) ||
                                        (rl.np[0]+rl.nh[0] == 0 && rr.np[0]+rr.nh[0] == 0))
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
                                lGen = op.l->matching(_zero, _l, c+_one, _max);
                                rGen = op.r->matching(c+_one, _max, _zero, _r);
                                if (lGen->next(ll,lr) && rGen->next(rl,rr))
                                {
                                    cdone = false;
                                    chdone = false;
                                }
                                else
                                {
                                    delete lGen;
                                    delete rGen;

                                    _next.np[0] = c.np[0]+1;
                                    _next.nh[0] = 0;

                                    lGen = op.l->matching(_zero, _l, _next, _max);
                                    rGen = op.r->matching(_next, _max, _zero, _r);
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
                        lGen = op.l->matching(_zero, _l, _r+_one, _max);
                        if (lGen->next(ll,lr))
                        {
                            _rdone = false;
                            _rhdone = false;
                        }
                        else
                        {
                            delete lGen;

                            _next.np[0] = _r.np[0]+1;
                            _next.nh[0] = 0;

                            lGen = op.l->matching(_zero, _l, _next, _max);
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
                for (__l.np[0] = 0;__l.np[0] <= _l.np[0] && _ldone && _lhdone;__l.np[0]++)
                {
                    for (__l.nh[0] = 0;__l.nh[0] <= _l.nh[0] && _ldone && _lhdone;__l.nh[0]++)
                    {
                        lGen = op.l->matching(__l, _max, _zero, _max);
                        if (lGen->next(ll,lr))
                        {
                            //cout << "found1 " << ll << ' ' << lr << endl;

                            if (ll.np[0] > _l.np[0] && ll.nh[0] > _l.nh[0])
                            {
                                _ldone = false;
                                _lhdone = false;
                            }
                            else
                            {
                                _next.np[0] = max(0,_l.np[0]-ll.np[0]);
                                _next.nh[0] = max(0,_l.nh[0]+1-ll.nh[0]);

                                rGen = op.r->matching(_next, _max, _zero, _max);
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

                        _next.np[0] = __l.np[0];
                        _next.nh[0] = 0;

                        lGen = op.l->matching( _next, _max, _zero, _max);
                        if (lGen->next(ll,lr))
                        {
                            //cout << "found3 " << ll << ' ' << lr << endl;
                            if (ll.np > _l.np)
                            {
                                _ldone = false;
                            }
                            else
                            {
                                _next.np[0] = max(0,_l.np[0]+1-ll.np[0]);
                                _next.nh[0] = 0;

                                rGen = op.r->matching(_next, _max, _zero, _max);
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
