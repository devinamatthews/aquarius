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
#include <cfloat>
#include <cmath>

#include "util/math_ext.h"
#include "util/stl_ext.hpp"

namespace aquarius
{
namespace symmetry
{

class PointGroup;

class Representation : public std::vector<double>
{
    friend class PointGroup;

    protected:
        const PointGroup& group;

        Representation(const PointGroup& group, const std::vector<double>& rep);

    public:
        Representation(const PointGroup& group);

        Representation(const PointGroup& group, int L, int m);

        Representation(const PointGroup& group, int x, int y, int z);

        const PointGroup& getPointGroup() const;

        Representation& operator=(const Representation& other);

        bool isReducible() const;

        bool isTotallySymmetric() const;

        bool transformsAs(const Representation r) const;

        Representation operator^(int p) const;

        Representation& operator^=(int p);

        Representation operator*(const Representation& other) const;

        Representation& operator*=(const Representation& other);

        Representation operator+(const Representation& other) const;

        Representation& operator+=(const Representation& other);

        Representation operator-() const;

        Representation operator-(const Representation& other) const;

        Representation& operator-=(const Representation& other);
};

class PointGroup
{
    friend class Representation;

    public:
        static const PointGroup& C1();
        static const PointGroup& Cs();
        static const PointGroup& Ci();
        static const PointGroup& Td();
        static const PointGroup& Oh();
        static const PointGroup& Ih();
        static const PointGroup& C2();
        static const PointGroup& C3();
        static const PointGroup& C4();
        static const PointGroup& C5();
        static const PointGroup& C6();
        static const PointGroup& C2v();
        static const PointGroup& C3v();
        static const PointGroup& C4v();
        static const PointGroup& C5v();
        static const PointGroup& C6v();
        static const PointGroup& C2h();
        static const PointGroup& C3h();
        static const PointGroup& C4h();
        static const PointGroup& C5h();
        static const PointGroup& C6h();
        static const PointGroup& D2();
        static const PointGroup& D3();
        static const PointGroup& D4();
        static const PointGroup& D5();
        static const PointGroup& D6();
        static const PointGroup& D2h();
        static const PointGroup& D3h();
        static const PointGroup& D4h();
        static const PointGroup& D5h();
        static const PointGroup& D6h();
        static const PointGroup& D2d();
        static const PointGroup& D3d();
        static const PointGroup& S4();
        static const PointGroup& S6();

    protected:
        const int order;
        const int nirrep;
        const char *name;
        const char **irrep_names;
        const int *generators;
        const double *generator_reps;
        std::vector<std::vector<double> > characters;
        std::vector<std::vector<std::vector<double> > > reps;
        const int *irrep_degen;
        const mat3x3 *ops;
        const char **op_names;
        std::vector<int> op_inverse;
        std::vector<int> op_product;
        std::vector<Representation> irreps;

        PointGroup(int order, int nirrep, int ngenerators, const char *name,
                   const char **irrep_names, const int *generators, const double *generator_reps,
                   const int *irrep_degen, const mat3x3 *ops, const char **op_names);

    public:
        bool operator==(const PointGroup& other) const { return this == &other; }

        std::vector<int> DCR(const std::vector<int>& U, const std::vector<int>& V, int& lambda) const;

        double character(int irrep, int op) const { return characters[irrep][op]; }

        double sphericalParity(int L, int m, int op) const;

        double cartesianParity(int x, int y, int z, int op) const;

        int getOrder() const { return order; }

        int getNumIrreps() const { return nirrep; }

        const char* getName() const { return name; }

        bool areConjugate(int i, int j) const;

        int getIrrepDegeneracy(int i) const { return irrep_degen[i]; }

        /*
         * Get the chacter representation of the given irrep
         */
        const Representation& getIrrep(int i) const;

        /*
         * Get the representation of the given row of the given irrep
         */
        Representation getIrrep(int i, int r) const;

        /*
         * Get the off-diagonal representation of the given given irrep
         */
        Representation getIrrep(int i, int r, int s) const;

        const char * getIrrepName(int i) const { return irrep_names[i]; }

        const mat3x3* getOps() const { return ops; }

        const mat3x3& getOp(int i) const { return ops[i]; }

        const char * getOpName(int i) const { return op_names[i]; }

        int getOpInverse(int i) const { return op_inverse[i]; }

        int getOpProduct(int i, int j) const { return op_product[i*order+j]; }

        const Representation& totallySymmetricIrrep() const { return getIrrep(0); }
};

inline mat3x3 Rotation(vec3 axis, double degrees)
{
    axis.normalize();
    double c = cos(degrees*M_PI/180);
    double s = sin(degrees*M_PI/180);
    double x = axis[0];
    double y = axis[1];
    double z = axis[2];
    return mat3x3(x*x*(1-c)+c  , x*y*(1-c)-z*s, x*z*(1-c)+y*s,
                  x*y*(1-c)+z*s, y*y*(1-c)+c  , y*z*(1-c)-x*s,
                  x*z*(1-c)-y*s, y*z*(1-c)+x*s, z*z*(1-c)+c  );
}

inline mat3x3 Rotation(vec3 a, vec3 b)
{
    a.normalize();
    b.normalize();
    vec3 axis = a^b;
    double cost = a*b;
    return cost-(b|a)+(a|b)+(1-cost)*(axis|axis)/(axis*axis+DBL_MIN);
}

inline mat3x3 Reflection(vec3 axis)
{
    axis.normalize();
    double x = axis[0];
    double y = axis[1];
    double z = axis[2];
    return mat3x3(1-2*x*x,  -2*x*y,  -2*x*z,
                   -2*x*y, 1-2*y*y,  -2*y*z,
                   -2*x*z,  -2*y*z, 1-2*z*z);
}

inline mat3x3 Identity()
{
    return mat3x3(1, 0, 0,
                  0, 1, 0,
                  0, 0, 1);
}

inline mat3x3 Inversion()
{
    return mat3x3(-1,  0,  0,
                   0, -1,  0,
                   0,  0, -1);
}

template <int n>
mat3x3 C(const vec3& axis)
{
    return Rotation(axis, 360.0/n);
}

template <int n>
mat3x3 S(const vec3& axis)
{
    return Rotation(axis, 360.0/n)*Reflection(axis);
}

}
}

#endif
