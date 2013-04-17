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

#include "elemental.hpp"
#include "mpi.h"

#include "slide/slide.hpp"
#include "stl_ext/stl_ext.hpp"

using namespace std;
using namespace elem;
using namespace MPI;
using namespace aquarius::slide;

static bool sortIntsByRS(const idx2_t& i1, const idx2_t& i2)
{
    cout << "compare: " << i1.i << " " << i1.j << " " << i2.i << " " << i2.j << " ";

    if (i1.j < i2.j)
    {
        cout << "T" << endl;
        return true;
    }
    else if (i1.j > i2.j)
    {
        cout << "F" << endl;
        return false;
    }
    if (i1.i < i2.i)
    {
        cout << "T" << endl;
        return true;
    }
    else
    {
        cout << "F" << endl;
        return false;
    }
}

template <typename T>
struct sortBackwards
{
    bool operator()(const T& a, const T& b) const
    {
        return a < b;
    }
};

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    SLIDE::init();
    elem::Initialize(argc, argv);

    //idx2_t a[] = {{1,4},{2,3},{5,1},{3,8},{3,7},{3,4},{8,9},{0,3},{9,1},{4,6}};
    int a[] = {10,2,7,4,8,3,9,1,5,6};
    int b[] = {1,2,3,4,5,6,7,8,9,10};

    for (int i = 0;i < 10;i++)
    {
        //cout << a[i].i << " " << a[i].j << " " << b[i] << endl;
        cout << a[i] << " " << b[i] << endl;
    }
    cout << endl;

    //cosort(a, a+10, b, b+10);
    cosort(a, a+10, b, b+10, sortBackwards<int>());

    for (int i = 0;i < 10;i++)
    {
        //cout << a[i].i << " " << a[i].j << " " << b[i] << endl;
        cout << a[i] << " " << b[i] << endl;
    }

    elem::Finalize();
    SLIDE::finish();
    MPI_Finalize();
}
