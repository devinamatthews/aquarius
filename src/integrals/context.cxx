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

#include "context.hpp"
#include "internal.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::integrals;
using namespace aquarius::symmetry;

void Context::setCartesianOrdering(int L, const vector<int>& ordering)
{
    assert(L >= 0);
    assert(ordering.size() == (L+1)*(L+2)/2);

    if (L >= cartesian_ordering.size())
    {
        int Lmin = cartesian_ordering.size();
        cartesian_ordering.resize(L+1);
        for (int l = Lmin;l < L;l++)
        {
            cartesian_ordering[l].resize((l+1)*(l+2)/2);
            for (int c = 0;c < (l+1)*(l+2);c++)
            {
                cartesian_ordering[l][c] = c;
            }
        }
    }

    cartesian_ordering[L] = ordering;
}

const vector<int>& Context::getCartesianOrdering(int L) const
{
    assert(L >= 0);

    if (L >= cartesian_ordering.size())
    {
        int Lmin = cartesian_ordering.size();
        cartesian_ordering.resize(L+1);
        for (int l = Lmin;l <= L;l++)
        {
            cartesian_ordering[l].resize((l+1)*(l+2)/2);
            for (int c = 0;c < (l+1)*(l+2)/2;c++)
            {
                cartesian_ordering[l][c] = c;
            }
        }
    }

    return cartesian_ordering[L];
}

void Context::setSphericalOrdering(int L, const vector<int>& ordering)
{
    assert(L >= 0);
    assert(ordering.size() == 2*L+1);

    if (L >= spherical_ordering.size())
    {
        int Lmin = spherical_ordering.size();
        spherical_ordering.resize(L+1);
        for (int l = Lmin;l < L;l++)
        {
            spherical_ordering[l].resize(2*l+1);
            for (int c = 0;c < 2*l+1;c++)
            {
                spherical_ordering[l][c] = c;
            }
        }
    }

    spherical_ordering[L] = ordering;
}

const vector<int>& Context::getSphericalOrdering(int L) const
{
    assert(L >= 0);

    if (L >= spherical_ordering.size())
    {
        int Lmin = spherical_ordering.size();
        spherical_ordering.resize(L+1);
        for (int l = Lmin;l <= L;l++)
        {
            spherical_ordering[l].resize(2*l+1);
            for (int c = 0;c < 2*l+1;c++)
            {
                spherical_ordering[l][c] = c;
            }
        }
    }

    return spherical_ordering[L];
}
