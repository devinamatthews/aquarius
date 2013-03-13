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

#ifndef _AQUARIUS_SYMMETRY_HPP_
#define _AQUARIUS_SYMMETRY_HPP_

#include <stdint.h>
#include <cassert>

#include "util/math.hpp"

namespace aquarius
{
namespace symmetry
{

class Representation;

class PointGroup
{
    friend class Representation;

    public:
        static const PointGroup C1;
        static const PointGroup Cs;
        static const PointGroup Ci;
        static const PointGroup Td;
        static const PointGroup Oh;
        static const PointGroup Ih;

        static const PointGroup C2;
        static const PointGroup C3;
        static const PointGroup C4;
        static const PointGroup C5;
        static const PointGroup C6;

        static const PointGroup C2v;
        static const PointGroup C3v;
        static const PointGroup C4v;
        static const PointGroup C5v;
        static const PointGroup C6v;

        static const PointGroup C2h;
        static const PointGroup C3h;
        static const PointGroup C4h;
        static const PointGroup C5h;
        static const PointGroup C6h;

        static const PointGroup D2;
        static const PointGroup D3;
        static const PointGroup D4;
        static const PointGroup D5;
        static const PointGroup D6;

        static const PointGroup D2h;
        static const PointGroup D3h;
        static const PointGroup D4h;
        static const PointGroup D5h;
        static const PointGroup D6h;

        static const PointGroup D2d;
        static const PointGroup D3d;
        static const PointGroup D4d;
        static const PointGroup D5d;
        static const PointGroup D6d;

        static const PointGroup S4;
        static const PointGroup S6;
        static const PointGroup S8;
        static const PointGroup S10;
        static const PointGroup S12;

    protected:
        const int order;
        const int nirrep;
        const char *name;
        const Representation *dirprd;
        const Representation *irreps;
        const char **irrep_names;
        const double *characters;
        const int *irrep_degen;
        const mat3x3 *ops;
        const char **op_names;

        PointGroup(int order, int nirrep, const char *name, const Representation *dirprd,
                   const Representation *irreps, const char **irrep_names, const double *characters,
                   const int *irrep_degen, const mat3x3 *ops, const char **op_names);

    public:
        bool operator==(const PointGroup& other) const;

        int getOrder() const;

        int getNumIrreps() const;

        const char* getName() const;

        const Representation* getIrreps() const;

        const Representation& getIrrep(int i) const;

        const char * getIrrepName(int i) const;

        const mat3x3* getOps() const;

        const mat3x3& getOp(int i) const;

        const char * getOpName(int i) const;

        const Representation& totallySymmetricIrrep() const;
};

class Representation
{
    protected:
        const PointGroup& group;
        uint_fast32_t bits;

    public:
        Representation(const PointGroup& group, uint_fast32_t bits);

        Representation& operator=(const Representation& other);

        bool isReducible() const;

        bool isTotallySymmetric() const;

        Representation operator*(const Representation& other) const;

        Representation& operator*=(const Representation& other);

        Representation operator+(const Representation& other) const;

        Representation& operator+=(const Representation& other);
};

}
}

#endif
