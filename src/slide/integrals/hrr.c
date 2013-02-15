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

#include "internal.h"

void hrr(const int la, const int lb,
         const double* posa, const double* posb,
         double* integrals1, double* integrals2)
{
    assert(la >= lb);

    double* ints1 = integrals1;
    double* ints2 = integrals2;

    for (int la1 = la+lb;la1 > la;la1--)
    {
        hrr1(la, la1, la+lb-la1, posa, posb, ints1, ints2);

        double* tmp = ints1;
        ints1 = ints2;
        ints2 = tmp;
    }

    if (ints1 != integrals2)
    {
        size_t sz = ((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2);
        dcopy(sz, ints1, 1, integrals2, 1);
    }
}

void hrr1(const int la0, const int la1, const int lb,
          const double* posa, const double* posb,
          double* integrals1, double* integrals2)
{
    double fac[3];
    fac[0] = posa[0]-posb[0];
    fac[1] = posa[1]-posb[1];
    fac[2] = posa[2]-posb[2];

    size_t inca1 = 1;
    size_t incb1 = inca1*((la1+1)*(la1+2)*(la1+3)/6-la0*(la0+1)*(la0+2)/6);

    size_t inca2 = 1;
    size_t incb2 = inca2*(la1*(la1+1)*(la1+2)/6-la0*(la0+1)*(la0+2)/6);

    double *ap1b, *abp1;
    double *ab = integrals1+incb1*(lb+1)*(lb+2)/2-1-inca1*(la1+1)*(la1+2)/2;
    for (int la = la1-1;la >= la0;la--)
    {
        for (int ax = 0;ax <= la;ax++)
        {
            for (int ay = 0;ay <= la-ax;ay++)
            {
                for (int bx = 0;bx <= lb;bx++)
                {
                    for (int by = 0;by <= lb-bx;by++)
                    {
                        ap1b = ab+inca1*(la+1)*(la+2)/2+inca1*ax;
                        abp1 = ab+incb2*(lb+1)*(lb+2)/2+incb2*bx;
                        *abp1 = *ap1b + fac[2]*(*ab);
                        ab -= incb1;
                    }
                    ap1b = ab+inca1*(la+1)*(la+2)/2+inca1*(ax+1);
                    abp1 = ab+incb2*(lb+1)*(lb+2)/2+incb2*(bx+1);
                    *abp1 = *ap1b + fac[1]*(*ab);
                }
                ap1b = ab+inca1*(la+1)*(la+2)/2+inca1*(la+2);
                abp1 = ab+incb2*(lb+1)*(lb+2)/2+incb2*(lb+2);
                *abp1 = *ap1b + fac[0]*(*ab);
                ab -= inca1-incb1*(lb+1)*(lb+2)/2;
            }
        }
    }
}
