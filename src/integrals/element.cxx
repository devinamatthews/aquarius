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

#include "element.hpp"
#include "internal.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::integrals;
using namespace aquarius::input;
using namespace aquarius::symmetry;

Element Element::getElement(const string& symbol, int A)
{
    for (int i = 0;i < 101;i++)
    {
        if (strcasecmp(symbol.c_str(), elements[i].symbol) == 0)
        {
            for (int is = 0;is < 25;is++)
            {
                if (elements[is].nucleon[is] == A)
                {
                    return Element(elements[i].charge, elements[i].Z, elements[i].nucleon[is],
                                   elements[i].spin[is], elements[i].symbol, elements[i].name,
                                   elements[i].mass[is]);
                }
            }
            ERROR("Isotope A=%d not found for element %s", A, symbol.c_str());
        }
    }
    ERROR("Element %s not found", symbol.c_str());
}

Element Element::getElement(const string& symbol)
{
    for (int i = 0;i < 101;i++)
    {
        if (strcasecmp(symbol.c_str(), elements[i].symbol) == 0)
        {
            int is = elements[i].isotope;
            return Element(elements[i].charge, elements[i].Z, elements[i].nucleon[is],
                           elements[i].spin[is], elements[i].symbol, elements[i].name,
                           elements[i].mass[is]);
        }
    }
    ERROR("Element %s not found", symbol.c_str());
}
