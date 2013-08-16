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

#include <cmath>

#include "symmetry.hpp"

namespace aquarius
{
namespace symmetry
{

PointGroup::PointGroup(int order, int nirrep, const char *name, const Representation *dirprd,
                       const Representation *irreps, const char **irrep_names, const double *characters,
                       const int *irrep_degen, const mat3x3 *ops, const char **op_names)
: order(order), nirrep(nirrep), name(name), dirprd(dirprd), irreps(irreps), irrep_names(irrep_names),
  characters(characters), irrep_degen(irrep_degen), ops(ops), op_names(op_names) {}

bool PointGroup::operator==(const PointGroup& other) const
{
    /*
     * Comparing pointers is sufficient since only the static instances are allowed
     */
    return this == &other;
}

int PointGroup::getOrder() const
{
    return order;
}

int PointGroup::getNumIrreps() const
{
    return nirrep;
}

const char* PointGroup::getName() const
{
    return name;
}

const Representation* PointGroup::getIrreps() const
{
    return irreps;
}

const Representation& PointGroup::getIrrep(int i) const
{
    return irreps[i];
}

const char * PointGroup::getIrrepName(int i) const
{
    return irrep_names[i];
}

const mat3x3* PointGroup::getOps() const
{
    return ops;
}

const mat3x3& PointGroup::getOp(int i) const
{
    return ops[i];
}

const char * PointGroup::getOpName(int i) const
{
    return op_names[i];
}

const Representation& PointGroup::totallySymmetricIrrep() const
{
    return irreps[0];
}

Representation::Representation(const PointGroup& group, uint_fast32_t bits)
: group(group), bits(bits) {}

Representation& Representation::operator=(const Representation& other)
{
    assert(group == other.group);
    bits = other.bits;
    return *this;
}

bool Representation::isReducible() const
{
    bool found = false;

    uint_fast32_t mask = 1;
    for (int i = 0;i < group.nirrep;i++)
    {
        if (mask&bits)
        {
            if (found) return true;
            found = true;
        }
        mask <<= 1;
    }

    return false;
}

bool Representation::isTotallySymmetric() const
{
    return bits&1;
}

Representation Representation::operator*(const Representation& other) const
{
    return Representation(*this) *= other;
}

Representation& Representation::operator*=(const Representation& other)
{
    assert(group == other.group);

    int n = group.nirrep;

    uint_fast32_t oldbits = bits;
    bits = 0;

    uint_fast32_t maski = 1;
    for (int i = 0;i < n;i++)
    {
        if (maski&oldbits)
        {
            uint_fast32_t maskj = 1;
            for (int j = 0;j < n;j++)
            {
                if (maskj&other.bits)
                {
                    *this += group.dirprd[i*n+j];
                }
                maskj <<= 1;
            }
        }
        maski <<= 1;
    }

    return *this;
}

Representation Representation::operator+(const Representation& other) const
{
    return Representation(*this) += other;
}

Representation& Representation::operator+=(const Representation& other)
{
    assert(group == other.group);
    bits |= other.bits;
    return *this;
}

mat3x3 Rotation(vec3 axis, double degrees)
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

mat3x3 Reflection(vec3 axis)
{
    axis.normalize();
    double x = axis[0];
    double y = axis[1];
    double z = axis[2];
    return mat3x3(1-2*x*x,  -2*x*y,  -2*x*z,
                   -2*x*y, 1-2*y*y,  -2*y*z,
                   -2*x*z,  -2*y*z, 1-2*z*z);
}

mat3x3 Identity()
{
    return mat3x3(1, 0, 0,
                  0, 1, 0,
                  0, 0, 1);
}

mat3x3 Inversion()
{
    return mat3x3(-1,  0,  0,
                   0, -1,  0,
                   0,  0, -1);
}

template <int n>
mat3x3 C(const vec3& axis)
{
    return Rotation(axis, 360/n);
}

template <int n>
mat3x3 S(const vec3& axis)
{
    return Rotation(axis, 360/n)*Reflection(axis);
}

/*
 * C1
 */
#define A Representation(PointGroup::C1, 0x1)

static const char *c1_name = "C1";
static const Representation c1_dirprd[] = {A};
static const Representation c1_irreps[] = {A};
static const char *c1_irrep_names[] = {"A"};
static const double c1_characters[] = {1};
static const int c1_irrep_degen[] = {1};
static const mat3x3 c1_ops[] = {Identity()};
static const char *c1_op_names[] = {"E"};

const PointGroup PointGroup::C1 = PointGroup(1, 1, c1_name, c1_dirprd, c1_irreps,
                                             c1_irrep_names, c1_characters, c1_irrep_degen,
                                             c1_ops, c1_op_names);

#undef A

/*
 * Cs
 */
#define Ap  Representation(PointGroup::Cs, 0x1)
#define App Representation(PointGroup::Cs, 0x2)

static const char *cs_name = "Cs";
static const Representation cs_dirprd[] = { Ap,App,
                                           App, Ap};
static const Representation cs_irreps[] = {Ap,App};
static const char *cs_irrep_names[] = {"Ap","App"};
static const double cs_characters[] = { 1, 1,
                                        1,-1};
static const int cs_irrep_degen[] = {1,1};
static const mat3x3 cs_ops[] = {Identity(),
                                Reflection(vec3(0,0,1))};
static const char *cs_op_names[] = {"E", "sxy"};

const PointGroup PointGroup::Cs = PointGroup(2, 2, cs_name, cs_dirprd, cs_irreps,
                                             cs_irrep_names, cs_characters, cs_irrep_degen,
                                             cs_ops, cs_op_names);

#undef Ap
#undef App

/*
 * Ci
 */
#define Ag Representation(PointGroup::Ci, 0x1)
#define Au Representation(PointGroup::Ci, 0x2)

static const char *ci_name = "Ci";
static const Representation ci_dirprd[] = {Ag,Au,
                                           Au,Ag};
static const Representation ci_irreps[] = {Ag,Au};
static const char *ci_irrep_names[] = {"Ag","Au"};
static const double ci_characters[] = { 1, 1,
                                        1,-1};
static const int ci_irrep_degen[] = {1,1};
static const mat3x3 ci_ops[] = {Identity(),
                                Inversion()};
static const char *ci_op_names[] = {"E", "i"};

const PointGroup PointGroup::Ci = PointGroup(2, 2, ci_name, ci_dirprd, ci_irreps,
                                             ci_irrep_names, ci_characters, ci_irrep_degen,
                                             ci_ops, ci_op_names);

#undef Ag
#undef Au

/*
 * Td
 */
#define A1 Representation(PointGroup::Td, 0x01)
#define A2 Representation(PointGroup::Td, 0x02)
#define E  Representation(PointGroup::Td, 0x04)
#define T1 Representation(PointGroup::Td, 0x08)
#define T2 Representation(PointGroup::Td, 0x10)

static const char *td_name = "Td";
static const Representation td_dirprd[] = {A1, A2,       E,         T1,         T2,
                                           A2, A1,       E,         T2,         T1,
                                            E,  E, A1+A2+E,      T1+T2,      T1+T2,
                                           T1, T2,   T1+T2, A1+E+T1+T2, A2+E+T1+T2,
                                           T2, T1,   T1+T2, A2+E+T1+T2, A1+E+T1+T2};
static const Representation td_irreps[] = {A1,A2,E,T1,T2};
static const char *td_irrep_names[] = {"A1","A2","E","T1","T2"};
static const double td_characters[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                        2,-1,-1,-1,-1,-1,-1,-1,-1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        3, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,
                                        3, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1};
static const int td_irrep_degen[] = {1,1,2,3,3};
static const mat3x3 td_ops[] = {Identity(),
                                C<3>(vec3( 1, 1, 1)),
                                C<3>(vec3(-1,-1,-1)),
                                C<3>(vec3( 1,-1,-1)),
                                C<3>(vec3(-1, 1, 1)),
                                C<3>(vec3(-1, 1,-1)),
                                C<3>(vec3( 1,-1, 1)),
                                C<3>(vec3(-1,-1, 1)),
                                C<3>(vec3( 1, 1,-1)),
                                C<2>(vec3( 1, 0, 0)),
                                C<2>(vec3( 0, 1, 0)),
                                C<2>(vec3( 0, 0, 1)),
                                S<4>(vec3( 1, 0, 0)),
                                S<4>(vec3(-1, 0, 0)),
                                S<4>(vec3( 0, 1, 0)),
                                S<4>(vec3( 0,-1, 0)),
                                S<4>(vec3( 0, 0, 1)),
                                S<4>(vec3( 0, 0,-1)),
                                Reflection(vec3( 1, 1, 0)),
                                Reflection(vec3( 1,-1, 0)),
                                Reflection(vec3( 1, 0, 1)),
                                Reflection(vec3( 1, 0,-1)),
                                Reflection(vec3( 0, 1, 1)),
                                Reflection(vec3( 0, 1,-1))};
static const char *td_op_names[] = {"E",
                                    "C3+++", "C3---", "C3+--", "C3-++", "C3-+-", "C3+-+", "C3--+", "C3++-",
                                    "C2x", "C2y", "C2z",
                                    "S4x", "S4x^3", "S4y", "S4y^3", "S4z", "S4z^3",
                                    "sx+y", "sx-y", "sx+z", "sx-z", "sy+z", "sy-z"};

const PointGroup PointGroup::Td = PointGroup(24, 5, td_name, td_dirprd, td_irreps,
                                             td_irrep_names, td_characters, td_irrep_degen,
                                             td_ops, td_op_names);

#undef A1
#undef A2
#undef E
#undef T1
#undef T2

/*
 * Oh
 */
#define A1g Representation(PointGroup::Oh, 0x001)
#define A2g Representation(PointGroup::Oh, 0x002)
#define Eg  Representation(PointGroup::Oh, 0x004)
#define T1g Representation(PointGroup::Oh, 0x008)
#define T2g Representation(PointGroup::Oh, 0x010)
#define A1u Representation(PointGroup::Oh, 0x020)
#define A2u Representation(PointGroup::Oh, 0x040)
#define Eu  Representation(PointGroup::Oh, 0x080)
#define T1u Representation(PointGroup::Oh, 0x100)
#define T2u Representation(PointGroup::Oh, 0x200)

static const char *oh_name = "Oh";
static const Representation oh_dirprd[] = {A1g, A2g,         Eg,            T1g,            T2g, A1u, A2u,         Eu,            T1u,            T2u,
                                           A2g, A1g,         Eg,            T2g,            T1g, A2u, A1u,         Eu,            T2u,            T1u,
                                            Eg,  Eg, A1g+A2g+Eg,        T1g+T2g,        T1g+T2g,  Eu,  Eu, A1u+A2u+Eu,        T1u+T2u,        T1u+T2u,
                                           T1g, T2g,    T1g+T2g, A1g+Eg+T1g+T2g, A2g+Eg+T1g+T2g, T1u, T2u,    T1u+T2u, A1u+Eu+T1u+T2u, A2u+Eu+T1u+T2u,
                                           T2g, T1g,    T1g+T2g, A2g+Eg+T1g+T2g, A1g+Eg+T1g+T2g, T2u, T1u,    T1u+T2u, A2u+Eu+T1u+T2u, A1u+Eu+T1u+T2u,
                                           A1u, A2u,         Eu,            T1u,            T2u, A1g, A2g,         Eg,            T1g,            T2g,
                                           A2u, A1u,         Eu,            T2u,            T1u, A2g, A1g,         Eg,            T2g,            T1g,
                                            Eu,  Eu, A1u+A2u+Eu,        T1u+T2u,        T1u+T2u,  Eg,  Eg, A1g+A2g+Eg,        T1g+T2g,        T1g+T2g,
                                           T1u, T2u,    T1u+T2u, A1u+Eu+T1u+T2u, A2u+Eu+T1u+T2u, T1g, T2g,    T1g+T2g, A1g+Eg+T1g+T2g, A2g+Eg+T1g+T2g,
                                           T2u, T1u,    T1u+T2u, A2u+Eu+T1u+T2u, A1u+Eu+T1u+T2u, T2g, T1g,    T1g+T2g, A2g+Eg+T1g+T2g, A1g+Eg+T1g+T2g};
static const Representation oh_irreps[] = {A1g,A2g,Eg,T1g,T2g,A1u,A2u,Eu,T1u,T2u};
static const char *oh_irrep_names[] = {"A1g","A2g","Eg","T1g","T2g","A1u","A2u","Eu","T1u","T2u"};
static const double oh_characters[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,
                                        2,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1, 2, 2, 2, 0, 0, 0, 0, 0, 0,
                                        3, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1,-1,-1,-1, 3, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                        3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                        1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1,
                                        2,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,-2,-2,-2, 0, 0, 0, 0, 0, 0,
                                        3, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-3,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1};
static const int oh_irrep_degen[] = {1,1,2,3,3,1,1,2,3,3};
static const mat3x3 oh_ops[] = {Identity(),
                                C<3>(vec3( 1, 1, 1)),
                                C<3>(vec3(-1,-1,-1)),
                                C<3>(vec3( 1,-1,-1)),
                                C<3>(vec3(-1, 1, 1)),
                                C<3>(vec3(-1, 1,-1)),
                                C<3>(vec3( 1,-1, 1)),
                                C<3>(vec3(-1,-1, 1)),
                                C<3>(vec3( 1, 1,-1)),
                                C<2>(vec3( 1, 1, 0)),
                                C<2>(vec3( 1, 0, 1)),
                                C<2>(vec3( 0, 1, 1)),
                                C<2>(vec3( 1,-1, 0)),
                                C<2>(vec3( 1, 0,-1)),
                                C<2>(vec3( 0, 1,-1)),
                                C<4>(vec3( 1, 0, 0)),
                                C<4>(vec3(-1, 0, 0)),
                                C<4>(vec3( 0, 1, 0)),
                                C<4>(vec3( 0,-1, 0)),
                                C<4>(vec3( 0, 0, 1)),
                                C<4>(vec3( 0, 0,-1)),
                                C<2>(vec3( 1, 0, 0)),
                                C<2>(vec3( 0, 1, 0)),
                                C<2>(vec3( 0, 0, 1)),
                                Inversion(),
                                S<4>(vec3( 1, 0, 0)),
                                S<4>(vec3(-1, 0, 0)),
                                S<4>(vec3( 0, 1, 0)),
                                S<4>(vec3( 0,-1, 0)),
                                S<4>(vec3( 0, 0, 1)),
                                S<4>(vec3( 0, 0,-1)),
                                S<6>(vec3( 1, 1, 1)),
                                S<6>(vec3(-1,-1,-1)),
                                S<6>(vec3( 1,-1,-1)),
                                S<6>(vec3(-1, 1, 1)),
                                S<6>(vec3(-1, 1,-1)),
                                S<6>(vec3( 1,-1, 1)),
                                S<6>(vec3(-1,-1, 1)),
                                S<6>(vec3( 1, 1,-1)),
                                Reflection(vec3( 0, 0, 1)),
                                Reflection(vec3( 0, 1, 0)),
                                Reflection(vec3( 1, 0, 0)),
                                Reflection(vec3( 1, 1, 0)),
                                Reflection(vec3( 1,-1, 0)),
                                Reflection(vec3( 1, 0, 1)),
                                Reflection(vec3( 1, 0,-1)),
                                Reflection(vec3( 0, 1, 1)),
                                Reflection(vec3( 0, 1,-1))};
static const char *oh_op_names[] = {"E",
                                    "C3+++", "C3---", "C3+--", "C3-++", "C3-+-", "C3+-+", "C3--+", "C3++-",
                                    "C2x+y", "C2x+z", "C2y+z", "C2x-y", "C2x-z", "C2y-z",
                                    "C4x", "C4x^3", "C4y", "C4y^3", "C4z", "C4z^3",
                                    "C2x", "C2y", "C2z",
                                    "i",
                                    "S4x", "S4x^3", "S4y", "S4y^3", "S4z", "S4z^3",
                                    "S6+++", "S6---", "S6+--", "S6-++", "S6-+-", "S6+-+", "S6--+", "S6++-",
                                    "sxy", "sxz", "syz",
                                    "sx+y", "sx-y", "sx+z", "sx-z", "sy+z", "sy-z"};

const PointGroup PointGroup::Oh = PointGroup(48, 10, oh_name, oh_dirprd, oh_irreps,
                                             oh_irrep_names, oh_characters, oh_irrep_degen,
                                             oh_ops, oh_op_names);

#undef A1g
#undef A2g
#undef Eg
#undef T1g
#undef T2g
#undef A1u
#undef A2u
#undef Eu
#undef T1u
#undef T2u

/*
 * Ih
 */
#define Ag  Representation(PointGroup::Ih, 0x001)
#define T1g Representation(PointGroup::Ih, 0x002)
#define T2g Representation(PointGroup::Ih, 0x004)
#define Gg  Representation(PointGroup::Ih, 0x008)
#define Hg  Representation(PointGroup::Ih, 0x010)
#define Au  Representation(PointGroup::Ih, 0x020)
#define T1u Representation(PointGroup::Ih, 0x040)
#define T2u Representation(PointGroup::Ih, 0x080)
#define Gu  Representation(PointGroup::Ih, 0x100)
#define Hu  Representation(PointGroup::Ih, 0x200)

#define e1 (2*cos(M_PI*0.4)) // also equals 1/phi = phi-1
#define e2 (2*cos(M_PI*0.8)) // also equals -1/e1 = -phi

#define avv 1.1071487177940904 // angle between vertices
#define avf 0.6523581397843685 // angle between vertex and face center
#define aff 0.7297276562269659 // angle between face centers
                               // avv+aff+2*avf = pi

#define a cos( avv)       // vector to vertex in yz plane
#define b sin( avv)       //
#define c cos( avf)       // vector to center of upper face in yz plane
#define d sin(-avf)       //
#define e cos( avf+aff)   // vector to center of upper-mid face in yz plane
#define f sin(-avf-aff)   //
#define g cos( avv/2)     // vector to midpoint of upper edge in yz plane
#define h sin( avv/2)     //
#define i cos( avf+aff/2) // vector to midpoint of upper edge perp. to yz plane
#define j sin(-avf-aff/2) //

static const char *ih_name = "Ih";
static const Representation ih_dirprd[] = { Ag,           T1g,           T2g,               Gg,               Hg,  Au,           T1u,           T2u,               Gu,               Hu,
                                           T1g,     Ag+T1g+Hg,         Gg+Hg,        T2g+Gg+Hg,    T1g+T2g+Gg+Hg, T1u,     Au+T1u+Hu,         Gu+Hu,        T2u+Gu+Hu,    T1u+T2u+Gu+Hu,
                                           T2g,         Gg+Hg,     Ag+T2g+Hg,        T1g+Gg+Hg,    T1g+T2g+Gg+Hg, T2u,         Gu+Hu,     Au+T2u+Hu,        T1u+Gu+Hu,    T1u+T2u+Gu+Hu,
                                            Gg,     T2g+Gg+Hg,     T1g+Gg+Hg, Ag+T1g+T2g+Gg+Hg,    T1g+T2g+Gg+Hg,  Gu,     T2u+Gu+Hu,     T1u+Gu+Hu, Au+T1u+T2u+Gu+Hu,    T1u+T2u+Gu+Hu,
                                            Hg, T1g+T2g+Gg+Hg, T1g+T2g+Gg+Hg,    T1g+T2g+Gg+Hg, Ag+T1g+T2g+Gg+Hg,  Hu, T1u+T2u+Gu+Hu, T1u+T2u+Gu+Hu,    T1u+T2u+Gu+Hu, Au+T1u+T2u+Gu+Hu,
                                            Au,           T1u,           T2u,               Gu,               Hu,  Ag,           T1g,           T2g,               Gg,               Hg,
                                           T1u,     Au+T1u+Hu,         Gu+Hu,        T2u+Gu+Hu,    T1u+T2u+Gu+Hu, T1g,     Ag+T1g+Hg,         Gg+Hg,        T2g+Gg+Hg,    T1g+T2g+Gg+Hg,
                                           T2u,         Gu+Hu,     Au+T2u+Hu,        T1u+Gu+Hu,    T1u+T2u+Gu+Hu, T2g,         Gg+Hg,     Ag+T2g+Hg,        T1g+Gg+Hg,    T1g+T2g+Gg+Hg,
                                            Gu,     T2u+Gu+Hu,     T1u+Gu+Hu, Au+T1u+T2u+Gu+Hu,    T1u+T2u+Gu+Hu,  Gg,     T2g+Gg+Hg,     T1g+Gg+Hg, Ag+T1g+T2g+Gg+Hg,    T1g+T2g+Gg+Hg,
                                            Hu, T1u+T2u+Gu+Hu, T1u+T2u+Gu+Hu,    T1u+T2u+Gu+Hu, Au+T1u+T2u+Gu+Hu,  Hg, T1g+T2g+Gg+Hg, T1g+T2g+Gg+Hg,    T1g+T2g+Gg+Hg, Ag+T1g+T2g+Gg+Hg};
static const Representation ih_irreps[] = {Ag,T1g,T2g,Gg,Hg,Au,T1u,T2u,Gu,Hu};
static const char *ih_irrep_names[] = {"Ag","T1g","T2g","Hg","Gg","Au","T1u","T2u","Hu","Gu"};
static const double ih_characters[] = { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        3,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                        3,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                        4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                        3,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-3, e1, e1, e1, e1, e1, e1, e1, e1, e1, e1, e1, e1, e2, e2, e2, e2, e2, e2, e2, e2, e2, e2, e2, e2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        3,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e1,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2,-e2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-3, e2, e2, e2, e2, e2, e2, e2, e2, e2, e2, e2, e2, e1, e1, e1, e1, e1, e1, e1, e1, e1, e1, e1, e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                        4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                        5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
static const int ih_irrep_degen[] = {1,3,3,4,5,1,3,3,4,5};
static const mat3x3 ih_ops[] = {Identity(),
                                C<5>(vec3( 0, 0, 1)),
                                C<5>(vec3( 0, 0,-1)),
                                C<5>(vec3( 0, b, a)),
                                C<5>(vec3( 0,-b,-a)),
                                C<5>(vec3( 0, b, a)*C<5>(vec3(0,0, 1))),
                                C<5>(vec3( 0,-b,-a)*C<5>(vec3(0,0, 1))),
                                C<5>(vec3( 0, b, a)*C<5>(vec3(0,0,-1))),
                                C<5>(vec3( 0,-b,-a)*C<5>(vec3(0,0,-1))),
                                C<5>(vec3( 0, b, a)*(C<5>(vec3(0,0, 1))^2)),
                                C<5>(vec3( 0,-b,-a)*(C<5>(vec3(0,0, 1))^2)),
                                C<5>(vec3( 0, b, a)*(C<5>(vec3(0,0,-1))^2)),
                                C<5>(vec3( 0,-b,-a)*(C<5>(vec3(0,0,-1))^2)),
                                (C<5>(vec3( 0, 0, 1))^2),
                                (C<5>(vec3( 0, 0,-1))^2),
                                (C<5>(vec3( 0, b, a))^2),
                                (C<5>(vec3( 0,-b,-a))^2),
                                C<5>(vec3( 0, b, a)*C<5>(vec3(0,0, 1)))^2,
                                C<5>(vec3( 0,-b,-a)*C<5>(vec3(0,0, 1)))^2,
                                C<5>(vec3( 0, b, a)*C<5>(vec3(0,0,-1)))^2,
                                C<5>(vec3( 0,-b,-a)*C<5>(vec3(0,0,-1)))^2,
                                C<5>(vec3( 0, b, a)*(C<5>(vec3(0,0, 1))^2))^2,
                                C<5>(vec3( 0,-b,-a)*(C<5>(vec3(0,0, 1))^2))^2,
                                C<5>(vec3( 0, b, a)*(C<5>(vec3(0,0,-1))^2))^2,
                                C<5>(vec3( 0,-b,-a)*(C<5>(vec3(0,0,-1))^2))^2,
                                C<3>(vec3( 0, d, c)),
                                C<3>(vec3( 0,-d,-c)),
                                C<3>(vec3( 0, d, c)*C<5>(vec3(0,0, 1))),
                                C<3>(vec3( 0,-d,-c)*C<5>(vec3(0,0, 1))),
                                C<3>(vec3( 0, d, c)*C<5>(vec3(0,0,-1))),
                                C<3>(vec3( 0,-d,-c)*C<5>(vec3(0,0,-1))),
                                C<3>(vec3( 0, d, c)*(C<5>(vec3(0,0, 1))^2)),
                                C<3>(vec3( 0,-d,-c)*(C<5>(vec3(0,0, 1))^2)),
                                C<3>(vec3( 0, d, c)*(C<5>(vec3(0,0,-1))^2)),
                                C<3>(vec3( 0,-d,-c)*(C<5>(vec3(0,0,-1))^2)),
                                C<3>(vec3( 0, f, e)),
                                C<3>(vec3( 0,-f,-e)),
                                C<3>(vec3( 0, f, e)*C<5>(vec3(0,0, 1))),
                                C<3>(vec3( 0,-f,-e)*C<5>(vec3(0,0, 1))),
                                C<3>(vec3( 0, f, e)*C<5>(vec3(0,0,-1))),
                                C<3>(vec3( 0,-f,-e)*C<5>(vec3(0,0,-1))),
                                C<3>(vec3( 0, f, e)*(C<5>(vec3(0,0, 1))^2)),
                                C<3>(vec3( 0,-f,-e)*(C<5>(vec3(0,0, 1))^2)),
                                C<3>(vec3( 0, f, e)*(C<5>(vec3(0,0,-1))^2)),
                                C<3>(vec3( 0,-f,-e)*(C<5>(vec3(0,0,-1))^2)),
                                C<2>(vec3( 0, h, g)),
                                C<2>(vec3( 0, h, g)*C<5>(vec3(0,0, 1))),
                                C<2>(vec3( 0, h, g)*C<5>(vec3(0,0,-1))),
                                C<2>(vec3( 0, h, g)*(C<5>(vec3(0,0, 1))^2)),
                                C<2>(vec3( 0, h, g)*(C<5>(vec3(0,0,-1))^2)),
                                C<2>(vec3( 0, j, i)),
                                C<2>(vec3( 0, j, i)*C<5>(vec3(0,0, 1))),
                                C<2>(vec3( 0, j, i)*C<5>(vec3(0,0,-1))),
                                C<2>(vec3( 0, j, i)*(C<5>(vec3(0,0, 1))^2)),
                                C<2>(vec3( 0, j, i)*(C<5>(vec3(0,0,-1))^2)),
                                C<2>(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))),
                                C<2>(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))*C<5>(vec3(0,0, 1))),
                                C<2>(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))*C<5>(vec3(0,0,-1))),
                                C<2>(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))*(C<5>(vec3(0,0, 1))^2)),
                                C<2>(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))*(C<5>(vec3(0,0,-1))^2)),
                                Inversion(),
                                S<10>(vec3( 0, 0, 1)),
                                S<10>(vec3( 0, 0,-1)),
                                S<10>(vec3( 0, b, a)),
                                S<10>(vec3( 0,-b,-a)),
                                S<10>(vec3( 0, b, a)*C<5>(vec3(0,0, 1))),
                                S<10>(vec3( 0,-b,-a)*C<5>(vec3(0,0, 1))),
                                S<10>(vec3( 0, b, a)*C<5>(vec3(0,0,-1))),
                                S<10>(vec3( 0,-b,-a)*C<5>(vec3(0,0,-1))),
                                S<10>(vec3( 0, b, a)*(C<5>(vec3(0,0, 1))^2)),
                                S<10>(vec3( 0,-b,-a)*(C<5>(vec3(0,0, 1))^2)),
                                S<10>(vec3( 0, b, a)*(C<5>(vec3(0,0,-1))^2)),
                                S<10>(vec3( 0,-b,-a)*(C<5>(vec3(0,0,-1))^2)),
                                S<10>(vec3( 0, 0, 1))^3,
                                S<10>(vec3( 0, 0,-1))^3,
                                S<10>(vec3( 0, b, a))^3,
                                S<10>(vec3( 0,-b,-a))^3,
                                S<10>(vec3( 0, b, a)*C<5>(vec3(0,0, 1)))^3,
                                S<10>(vec3( 0,-b,-a)*C<5>(vec3(0,0, 1)))^3,
                                S<10>(vec3( 0, b, a)*C<5>(vec3(0,0,-1)))^3,
                                S<10>(vec3( 0,-b,-a)*C<5>(vec3(0,0,-1)))^3,
                                S<10>(vec3( 0, b, a)*(C<5>(vec3(0,0, 1))^2))^3,
                                S<10>(vec3( 0,-b,-a)*(C<5>(vec3(0,0, 1))^2))^3,
                                S<10>(vec3( 0, b, a)*(C<5>(vec3(0,0,-1))^2))^3,
                                S<10>(vec3( 0,-b,-a)*(C<5>(vec3(0,0,-1))^2))^3,
                                S<6>(vec3( 0, d, c)),
                                S<6>(vec3( 0,-d,-c)),
                                S<6>(vec3( 0, d, c)*C<5>(vec3(0,0, 1))),
                                S<6>(vec3( 0,-d,-c)*C<5>(vec3(0,0, 1))),
                                S<6>(vec3( 0, d, c)*C<5>(vec3(0,0,-1))),
                                S<6>(vec3( 0,-d,-c)*C<5>(vec3(0,0,-1))),
                                S<6>(vec3( 0, d, c)*(C<5>(vec3(0,0, 1))^2)),
                                S<6>(vec3( 0,-d,-c)*(C<5>(vec3(0,0, 1))^2)),
                                S<6>(vec3( 0, d, c)*(C<5>(vec3(0,0,-1))^2)),
                                S<6>(vec3( 0,-d,-c)*(C<5>(vec3(0,0,-1))^2)),
                                S<6>(vec3( 0, f, e)),
                                S<6>(vec3( 0,-f,-e)),
                                S<6>(vec3( 0, f, e)*C<5>(vec3(0,0, 1))),
                                S<6>(vec3( 0,-f,-e)*C<5>(vec3(0,0, 1))),
                                S<6>(vec3( 0, f, e)*C<5>(vec3(0,0,-1))),
                                S<6>(vec3( 0,-f,-e)*C<5>(vec3(0,0,-1))),
                                S<6>(vec3( 0, f, e)*(C<5>(vec3(0,0, 1))^2)),
                                S<6>(vec3( 0,-f,-e)*(C<5>(vec3(0,0, 1))^2)),
                                S<6>(vec3( 0, f, e)*(C<5>(vec3(0,0,-1))^2)),
                                S<6>(vec3( 0,-f,-e)*(C<5>(vec3(0,0,-1))^2)),
                                Reflection(vec3( 0, h, g)),
                                Reflection(vec3( 0, h, g)*C<5>(vec3(0,0, 1))),
                                Reflection(vec3( 0, h, g)*C<5>(vec3(0,0,-1))),
                                Reflection(vec3( 0, h, g)*(C<5>(vec3(0,0, 1))^2)),
                                Reflection(vec3( 0, h, g)*(C<5>(vec3(0,0,-1))^2)),
                                Reflection(vec3( 0, j, i)),
                                Reflection(vec3( 0, j, i)*C<5>(vec3(0,0, 1))),
                                Reflection(vec3( 0, j, i)*C<5>(vec3(0,0,-1))),
                                Reflection(vec3( 0, j, i)*(C<5>(vec3(0,0, 1))^2)),
                                Reflection(vec3( 0, j, i)*(C<5>(vec3(0,0,-1))^2)),
                                Reflection(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))),
                                Reflection(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))*C<5>(vec3(0,0, 1))),
                                Reflection(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))*C<5>(vec3(0,0,-1))),
                                Reflection(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))*(C<5>(vec3(0,0, 1))^2)),
                                Reflection(vec3( 0, 1, 0)*C<20>(vec3(0,0,1))*(C<5>(vec3(0,0,-1))^2))};
static const char *ih_op_names[] = {"E",
                                    "C5", "C5", "C5", "C5", "C5", "C5", "C5", "C5", "C5", "C5", "C5", "C5",
                                    "C5^2", "C5^2", "C5^2", "C5^2", "C5^2", "C5^2", "C5^2", "C5^2", "C5^2", "C5^2", "C5^2", "C5^2",
                                    "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3",
                                    "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2",
                                    "i",
                                    "S10", "S10", "S10", "S10", "S10", "S10", "S10", "S10", "S10", "S10", "S10", "S10",
                                    "S10^3", "S10^3", "S10^3", "S10^3", "S10^3", "S10^3", "S10^3", "S10^3", "S10^3", "S10^3", "S10^3", "S10^3",
                                    "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6",
                                    "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg"};

const PointGroup PointGroup::Ih = PointGroup(120, 10, ih_name, ih_dirprd, ih_irreps,
                                             ih_irrep_names, ih_characters, ih_irrep_degen,
                                             ih_ops, ih_op_names);

#undef a
#undef b
#undef c
#undef d
#undef e
#undef f
#undef g
#undef h
#undef i
#undef j

#undef avv
#undef avf
#undef aff

#undef e1
#undef e2

#undef Ag
#undef T1g
#undef T2g
#undef Gg
#undef Hg
#undef Au
#undef T1u
#undef T2u
#undef Gu
#undef Hu

/*
 * C2
 */
#define A1 Representation(PointGroup::C2, 0x1)
#define A2 Representation(PointGroup::C2, 0x2)

static const char *c2_name = "C2";
static const Representation c2_dirprd[] = {A1,A1,
                                           A2,A2};
static const Representation c2_irreps[] = {A1,A2};
static const char *c2_irrep_names[] = {"A1","A2"};
static const double c2_characters[] = { 1, 1,
                                        1,-1};
static const int c2_irrep_degen[] = {1,1};
static const mat3x3 c2_ops[] = {Identity(),
                                C<2>(vec3(0,0,1))};
static const char *c2_op_names[] = {"E", "C2z"};

const PointGroup PointGroup::C2 = PointGroup(2, 2, c2_name, c2_dirprd, c2_irreps,
                                             c2_irrep_names, c2_characters, c2_irrep_degen,
                                             c2_ops, c2_op_names);

#undef A1
#undef A2

/*
 * C3
 */

/*
 * C4
 */

/*
 * C5
 */

/*
 * C6
 */

/*
 * C2v
 */
#define A1 Representation(PointGroup::C2v, 0x1)
#define A2 Representation(PointGroup::C2v, 0x2)
#define B1 Representation(PointGroup::C2v, 0x4)
#define B2 Representation(PointGroup::C2v, 0x8)

static const char *c2v_name = "C2v";
static const Representation c2v_dirprd[] = {A1,A1,B1,B2,
                                            A2,A2,B2,B1,
                                            B1,B2,A1,A2,
                                            B2,B1,A2,A1};
static const Representation c2v_irreps[] = {A1,A2,B1,B2};
static const char *c2v_irrep_names[] = {"A1","A2","B1","B2"};
static const double c2v_characters[] = { 1, 1, 1, 1,
                                         1, 1,-1,-1,
                                         1,-1, 1,-1,
                                         1,-1,-1, 1};
static const int c2v_irrep_degen[] = {1,1,1,1};
static const mat3x3 c2v_ops[] = {Identity(),
                                 C<2>(vec3(0,0,1)),
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(1,0,0))};
static const char *c2v_op_names[] = {"E", "C2z", "sxz", "syz"};

const PointGroup PointGroup::C2v = PointGroup(4, 4, c2v_name, c2v_dirprd, c2v_irreps,
                                              c2v_irrep_names, c2v_characters, c2v_irrep_degen,
                                              c2v_ops, c2v_op_names);

#undef A1
#undef A2
#undef B1
#undef B2

/*
 * C3v
 */

/*
 * C4v
 */

/*
 * C5v
 */

/*
 * C6v
 */

/*
 * C2h
 */
#define Ag Representation(PointGroup::C2h, 0x1)
#define Bg Representation(PointGroup::C2h, 0x2)
#define Au Representation(PointGroup::C2h, 0x4)
#define Bu Representation(PointGroup::C2h, 0x8)

static const char *c2h_name = "C2h";
static const Representation c2h_dirprd[] = {Ag,Bg,Au,Bu,
                                            Bg,Ag,Bu,Au,
                                            Au,Bu,Ag,Bg,
                                            Bu,Au,Bg,Ag};
static const Representation c2h_irreps[] = {Ag,Bg,Au,Bu};
static const char *c2h_irrep_names[] = {"Ag","Bg","Au","Bu"};
static const double c2h_characters[] = { 1, 1, 1, 1,
                                         1,-1, 1,-1,
                                         1, 1,-1,-1,
                                         1,-1,-1, 1};
static const int c2h_irrep_degen[] = {1,1,1,1};
static const mat3x3 c2h_ops[] = {Identity(),
                                 C<2>(vec3(0,0,1)),
                                 Inversion(),
                                 Reflection(vec3(0,0,1))};
static const char *c2h_op_names[] = {"E", "C2z", "i", "sxy"};

const PointGroup PointGroup::C2h = PointGroup(4, 4, c2h_name, c2h_dirprd, c2h_irreps,
                                              c2h_irrep_names, c2h_characters, c2h_irrep_degen,
                                              c2h_ops, c2h_op_names);

#undef Ag
#undef Bg
#undef Au
#undef Bu

/*
 * C3h
 */

/*
 * C4h
 */

/*
 * C5h
 */

/*
 * C6h
 */

/*
 * D2
 */
#define A  Representation(PointGroup::D2, 0x1)
#define B1 Representation(PointGroup::D2, 0x2)
#define B2 Representation(PointGroup::D2, 0x4)
#define B3 Representation(PointGroup::D2, 0x8)

static const char *d2_name = "D2";
static const Representation d2_dirprd[] = {A ,B1,B2,B3,
                                           B1,A ,B3,B2,
                                           B2,B3,A ,B1,
                                           B3,B2,B1,A };
static const Representation d2_irreps[] = {A,B1,B2,B3};
static const char *d2_irrep_names[] = {"A","B1","B2","B3"};
static const double d2_characters[] = { 1, 1, 1, 1,
                                        1, 1,-1,-1,
                                        1,-1, 1,-1,
                                        1,-1,-1, 1};
static const int d2_irrep_degen[] = {1,1,1,1};
static const mat3x3 d2_ops[] = {Identity(),
                                C<2>(vec3(0,0,1)),
                                C<2>(vec3(0,1,0)),
                                C<2>(vec3(1,0,0))};
static const char *d2_op_names[] = {"E", "C2z", "C2y", "C2x"};

const PointGroup PointGroup::D2 = PointGroup(4, 4, d2_name, d2_dirprd, d2_irreps,
                                             d2_irrep_names, d2_characters, d2_irrep_degen,
                                             d2_ops, d2_op_names);

#undef A
#undef B1
#undef B2
#undef B3

/*
 * D3
 */

/*
 * D4
 */

/*
 * D5
 */

/*
 * D6
 */

/*
 * D2h
 */
#define Ag  Representation(PointGroup::D2h, 0x01)
#define B1g Representation(PointGroup::D2h, 0x02)
#define B2g Representation(PointGroup::D2h, 0x04)
#define B3g Representation(PointGroup::D2h, 0x08)
#define Au  Representation(PointGroup::D2h, 0x10)
#define B1u Representation(PointGroup::D2h, 0x20)
#define B2u Representation(PointGroup::D2h, 0x40)
#define B3u Representation(PointGroup::D2h, 0x80)

static const char *d2h_name = "D2h";
static const Representation d2h_dirprd[] = {Ag ,B1g,B2g,B3g,Au ,B1u,B2u,B3u,
                                            B1g,Ag ,B3g,B2g,B1u,Au ,B3u,B2u,
                                            B2g,B3g,Ag ,B1g,B2u,B3u,Au ,B1u,
                                            B3g,B2g,B1g,Ag ,B3u,B2u,B1u,Au ,
                                            Au ,B1u,B2u,B3u,Ag ,B1g,B2g,B3g,
                                            B1u,Au ,B3u,B2u,B1g,Ag ,B3g,B2g,
                                            B2u,B3u,Au ,B1u,B2g,B3g,Ag ,B1g,
                                            B3u,B2u,B1u,Au ,B3g,B2g,B1g,Ag };
static const Representation d2h_irreps[] = {Ag, B1g,B2g,B3g,Au, B1u,B2u,B3u};
static const char *d2h_irrep_names[] = {"Ag","B1g","B2g","B3g","Au","B1u","B2u","B3u"};
static const double d2h_characters[] = { 1, 1, 1, 1, 1, 1, 1, 1,
                                         1, 1,-1,-1, 1, 1,-1,-1,
                                         1,-1, 1,-1, 1,-1, 1,-1,
                                         1,-1,-1, 1, 1,-1,-1, 1,
                                         1, 1, 1, 1,-1,-1,-1,-1,
                                         1, 1,-1,-1,-1,-1, 1, 1,
                                         1,-1, 1,-1,-1, 1,-1, 1,
                                         1,-1,-1, 1,-1, 1, 1,-1};
static const int d2h_irrep_degen[] = {1,1,1,1,1,1,1,1};
static const mat3x3 d2h_ops[] = {Identity(),
                                 C<2>(vec3(0,0,1)),
                                 C<2>(vec3(0,1,0)),
                                 C<2>(vec3(1,0,0)),
                                 Inversion(),
                                 Reflection(vec3(0,0,1)),
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(1,0,0))};
static const char *d2h_op_names[] = {"E", "C2z", "C2y", "C2x", "i", "sxy", "sxz", "syz"};

const PointGroup PointGroup::D2h = PointGroup(8, 8, d2h_name, d2h_dirprd, d2h_irreps,
                                              d2h_irrep_names, d2h_characters, d2h_irrep_degen,
                                              d2h_ops, d2h_op_names);

#undef Ag
#undef B1g
#undef B2g
#undef B3g
#undef Au
#undef B1u
#undef B2u
#undef B3u

/*
 * D3h
 */

/*
 * D4h
 */

/*
 * D5h
 */

/*
 * D6h
 */

/*
 * D2d
 */

/*
 * D3d
 */

/*
 * D4d
 */

/*
 * D5d
 */

/*
 * D6d
 */

/*
 * S4
 */

/*
 * S6
 */

/*
 * S8
 */

/*
 * S10
 */

/*
 * S12
 */

}
}
