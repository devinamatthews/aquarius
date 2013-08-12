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
#include <set>
#include <vector>

#include "util/util.h"
#include "util/blas.h"
#include "symmetry.hpp"
#include "stl_ext/stl_ext.hpp"

using namespace std;

namespace aquarius
{
namespace symmetry
{

PointGroup::PointGroup(int order, int nirrep, int ngenerators, const char *name,
                       const char **irrep_names, const int *generators, const double *generator_reps,
                       const int *irrep_degen, const mat3x3 *ops, const char **op_names)
: order(order), nirrep(nirrep), name(name), irrep_names(irrep_names),
  characters(nirrep, std::vector<double>(order)), reps(nirrep, std::vector<std::vector<double> >(order)),
  irrep_degen(irrep_degen), ops(ops), op_names(op_names),
  op_inverse(order, -1), op_product(order*order, -1)
{
    ADPRINTF("Creating point group %s\n", name);

    /*
     * Fill in representation information for the specified generators
     */
    int idx = 0;
    for (int i = 0;i < ngenerators;i++)
    {
        for (int irrep = 0;irrep < nirrep;irrep++)
        {
            ADPRINTF("Using hard-coded representation of %s(%s):\n", irrep_names[irrep], op_names[generators[i]]);
            int n = irrep_degen[irrep];
            for (int j = 0;j < n;j++)
            {
                for (int k = 0;k < n;k++)
                {
                    ADPRINTFC("% f ", generator_reps[idx+j*n+k]);
                }
                ADPRINTFC("\n");
            }
            reps[irrep][generators[i]].assign(generator_reps+idx, generator_reps+idx+n*n);
            idx += n*n;
        }
    }

    /*
     * For everybody except separably-degenerate groups, D(E) = I
     */
    if (reps[0][0].empty())
    {
        for (int irrep = 0;irrep < nirrep;irrep++)
        {
            ADPRINTF("Setting representation of %s(E) to I\n", irrep_names[irrep]);
            int n = irrep_degen[irrep];
            reps[irrep][0].resize(n*n, 0.0);
            for (int i = 0;i < n;i++) reps[irrep][0][i+i*n] = 1;
        }
    }

    /*
     * Find products and inverses of operators:
     *
     * ops[op_product[i*order+j]] = ops[i]*ops[j]
     *
     * ops[op_inverse[i]]*ops[i] = ops[i]*ops[op_inverse[i]] = ops[0] (E)
     */
    for (int i = 0;i < order;i++)
    {
        for (int j = 0;j < order;j++)
        {
            mat3x3 ij = ops[i]*ops[j];

            for (int k = 0;k < order;k++)
            {
                if (norm(ij-ops[k]) < 1e-10) op_product[i*order+j] = k;
            }
            if (op_product[i*order+j] == -1)
                ERROR("No product found for operations %s and %s of group %s", op_names[i], op_names[j], name);

            if (op_product[i*order+j] == 0) op_inverse[i] = j;
        }

        if (op_inverse[i] == -1) ERROR("No inverse found for operation %s of group %s", op_names[i], name);
    }

    /*
     * Generate remaining representations
     */
    for (int i = 1;i < order;i++)
    {
        for (int j = 1;j < order;j++)
        {
            int k = op_product[i*order+j];

            if (!reps[0][i].empty() && !reps[0][j].empty() && reps[0][k].empty())
            {
                for (int irrep = 0;irrep < nirrep;irrep++)
                {
                    ADPRINTF("Generating representation of %s(%s) from %s*%s:\n",
                            irrep_names[irrep], op_names[k], op_names[i], op_names[j]);
                    int n = irrep_degen[irrep];
                    reps[irrep][k].resize(n*n);
                    gemm('N', 'N', n, n, n, 1.0, reps[irrep][i].data(), n,
                                                 reps[irrep][j].data(), n,
                                            0.0, reps[irrep][k].data(), n);
                    for (int l = 0;l < n;l++)
                    {
                        for (int m = 0;m < n;m++)
                        {
                            ADPRINTFC("% f ", reps[irrep][k][l*n+m]);
                        }
                        ADPRINTFC("\n");
                    }
                }
            }
        }
    }

    ADPRINTF("Character table:\n\n        %9s ", op_names[0]);
    for (int i = 1;i < order;i++)
    {
        if (!areConjugate(i-1,i)) ADPRINTFC("%9s ", op_names[i]);
    }
    ADPRINTFC("\n      + ");
    for (int i = 0;i < nirrep*10-1;i++) ADPRINTFC("-");
    ADPRINTFC("\n");

    /*
     * Form character table
     */
    for (int irrep = 0;irrep < nirrep;irrep++)
    {
        ADPRINTFC("%5s | ", irrep_names[irrep]);
        for (int g = 0;g < order;g++)
        {
            if (reps[irrep][g].empty()) ERROR("No representation for %s(%s)", irrep_names[irrep], op_names[g]);
            double c = 0;
            int n = irrep_degen[irrep];
            for (int j = 0;j < n;j++) c += reps[irrep][g][j+j*n];
            characters[irrep][g] = c;
            if (!areConjugate(g-1,g)) ADPRINTFC("% f ", c);
        }
        ADPRINTFC("\n");
    }

    /*
     * Check normalization of character table
     */
    for (int irrep = 0;irrep < nirrep;irrep++)
    {
        double dp = 0;
        for (int g = 0;g < order;g++)
        {
            dp += characters[irrep][g]*characters[irrep][g];
        }
        if (abs(dp-order) > 1e-10)
            ERROR("Characters for irrep %s of group %s have the wrong norm", irrep_names[irrep], name);
    }

    /*
     * Check orthogonality of the representations
     */
    for (int i1 = 0;i1 < nirrep;i1++)
    {
        int n1 = irrep_degen[i1];
        for (int i2 = i1;i2 < nirrep;i2++)
        {
            int n2 = irrep_degen[i2];
            for (int j1 = 0;j1 < n1*n1;j1++)
            {
                for (int j2 = (i1 == i2 ? j1 : 0);j2 < n2*n2;j2++)
                {
                    double dp = 0;
                    for (int g = 0;g < order;g++)
                    {
                        dp += reps[i1][g][j1]*reps[i2][g][j2];
                    }
                    if (i1 == i2 && j1 == j2)
                    {
                        if (abs(dp-order/n1) > 1e-10)
                            ERROR("Representation %s[%d] of group %s has the wrong norm",
                                  irrep_names[i1], j1, name);
                    }
                    else
                    {
                        if (abs(dp) > 1e-10)
                            ERROR("Representations %s[%d] and %s[%d] of group %s are not orthognal",
                                  irrep_names[i1], j1, irrep_names[i2], j2, name);
                    }
                }
            }
        }
    }

    /*
     * Check that characters of elements are the same within classes
     */
    int nclasses = 1;
    for (int g = 1;g < order;g++)
    {
        if (!areConjugate(g-1, g))
        {
            nclasses++;
        }
        else
        {
            for (int irrep = 0;irrep < nirrep;irrep++)
            {
                if (abs(characters[irrep][g-1]-characters[irrep][g]) > 1e-10)
                    ERROR("Characters of same-class operations %s and %s are not the same for irrep %s of group %s",
                          op_names[g-1], op_names[g], irrep_names[irrep], name);
            }
        }
    }
    if (nclasses != nirrep) ERROR("The number of classes and irreps is not the same for group %s", name);
}

bool PointGroup::areConjugate(int i, int j) const
{
    for (int k = 0;k < order;k++)
    {
        // Check if i = k*j*k^1 for some k in G
        if (op_product[k*order+op_product[j*order+op_inverse[k]]] == i) return true;
    }
    return false;
}

vector<int> PointGroup::DCR(const vector<int>& U, const vector<int>& V, int& lambda) const
{
    //TODO: non-D2h

    // \lambda_R = |U^V|
    lambda = intersection(U,V).size();

    set<int> dcr;

    for (int g = 0;g < order;g++)
    {
        set<int> dcd;

        for (int u = 0;u < U.size();u++)
        {
            for (int v = 0;v < V.size();v++)
            {
                // R(G) = UGV
                // each element will appear \lambda_R times, but we already
                // know how many this is so only keep unique ones
                dcd.insert(op_product[U[u]*order+op_product[g*order+V[v]]]);
            }
        }

        // since R(G) are a (degenerate) disjoint partitioning,
        // R will have one unique operation per unique DCD
        dcr.insert(*dcd.begin());
    }

    return vector<int>(dcr.begin(), dcr.end());
}

double PointGroup::sphericalParity(int L, int m, int op) const
{
    /*
     * E   - sign(1)*(-1)^(0)
     * C2z - sign(1)*(-1)^(m)
     * C2y - sign(m)*(-1)^(l)
     * C2x - sign(m)*(-1)^(l+m)
     * i   - sign(1)*(-1)^(l)
     * sxy - sign(1)*(-1)^(l+m)
     * sxz - sign(m)*(-1)^(0)
     * syz - sign(m)*(-1)^(m)
     */

    int xneg = (ops[op][0][0] < 0 ? 1 : 0);
    int yneg = (ops[op][1][1] < 0 ? 1 : 0);
    int zneg = (ops[op][2][2] < 0 ? 1 : 0);

    return (xneg^yneg && m   <  0 ? -1 : 1)*
           (xneg^zneg && m&1 == 1 ? -1 : 1)*
           (zneg      && L&1 == 1 ? -1 : 1);
}

double PointGroup::cartesianParity(int x, int y, int z, int op) const
{
    /*
     * signs for cartesian functions:;
     *
     * E   - (-1)^(0)
     * C2z - (-1)^(lx+ly)
     * C2y - (-1)^(lx+lz)
     * C2x - (-1)^(ly+lz)
     * i   - (-1)^(lx+ly+lz)
     * sxy - (-1)^(lz)
     * sxz - (-1)^(ly)
     * syz - (-1)^(lx)
     */

    return (x&1 == 1 && ops[op][0][0] < 0 ? -1 : 1)*
           (y&1 == 1 && ops[op][1][1] < 0 ? -1 : 1)*
           (z&1 == 1 && ops[op][2][2] < 0 ? -1 : 1);
}

Representation PointGroup::getIrrep(int i) const
{
    return Representation(*this, characters[i]);
}

Representation PointGroup::getIrrep(int i, int r) const
{
    Representation rep(*this);

    int n = irrep_degen[i];

    for (int g = 0;g < order;g++)
    {
        rep[g] = reps[i][g][r*n+r];
    }

    return rep;
}

Representation PointGroup::getIrrep(int i, int r, int s) const
{
    Representation rep(*this);

    int n = irrep_degen[i];

    for (int g = 0;g < order;g++)
    {
        rep[g] = reps[i][g][r*n+s];
    }

    return rep;
}

Representation::Representation(const PointGroup& group)
: group(group), vector<double>(group.getOrder(), 0.0) {}

Representation::Representation(const PointGroup& group, const vector<double>& rep)
: group(group), vector<double>(rep) {}

Representation::Representation(const PointGroup& group, int L, int m)
: group(group)
{
    assert(0);
}

Representation::Representation(const PointGroup& group, int x, int y, int z)
: group(group)
{
    assert(0);
}

const PointGroup& Representation::getPointGroup() const
{
    return group;
}

Representation& Representation::operator=(const Representation& other)
{
    assert(group == other.group);
    assign(other.begin(), other.end());
    return *this;
}

bool Representation::isReducible() const
{
    int nmatch = 0;
    for (int i = 0;i < group.nirrep;i++)
    {
        if (transformsAs(group.getIrrep(i))) nmatch++;
        if (nmatch > 1) return true;
    }
    if (nmatch == 0) return true;
    return false;
}

bool Representation::isTotallySymmetric() const
{
    return transformsAs(group.totallySymmetricIrrep());
}

bool Representation::transformsAs(const Representation r) const
{
    assert(group == r.group);
    double dp = 0;
    for (int i = 0;i < size();i++)
    {
        dp += (*this)[i]*r[i];
    }
    return abs(dp) > 1e-10;
}

Representation Representation::operator^(int p) const
{
    return Representation(*this) ^= p;
}

Representation& Representation::operator^=(int p)
{
    for (int i = 0;i < size();i++) (*this)[i] = pow((*this)[i],p);
    return *this;
}

Representation Representation::operator*(const Representation& other) const
{
    return Representation(*this) *= other;
}

Representation& Representation::operator*=(const Representation& other)
{
    for (int i = 0;i < size();i++)
    {
        (*this)[i] *= other[i];
    }
    return *this;
}

Representation Representation::operator+(const Representation& other) const
{
    return Representation(*this) += other;
}

Representation& Representation::operator+=(const Representation& other)
{
    for (int i = 0;i < size();i++)
    {
        (*this)[i] += other[i];
    }
    return *this;
}

Representation Representation::operator-() const
{
    Representation neg(*this);

    for (int i = 0;i < size();i++)
    {
        neg[i] = -neg[i];
    }

    return neg;
}

Representation Representation::operator-(const Representation& other) const
{
    return Representation(*this) -= other;
}

Representation& Representation::operator-=(const Representation& other)
{
    for (int i = 0;i < size();i++)
    {
        (*this)[i] -= other[i];
    }
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
    return Rotation(axis, 360.0/n);
}

template <int n>
mat3x3 S(const vec3& axis)
{
    return Rotation(axis, 360.0/n)*Reflection(axis);
}

/*
 * C1
 */
static const char *c1_name = "C1";
static const char *c1_irrep_names[] = {"A"};
static const int c1_irrep_degen[] = {1};
static const mat3x3 c1_ops[] = {Identity()};
static const char *c1_op_names[] = {"E"};

const PointGroup PointGroup::C1 = PointGroup(1, 1, 0, c1_name, c1_irrep_names,
                                             NULL, NULL, c1_irrep_degen,
                                             c1_ops, c1_op_names);

/*
 * Cs
 */
static const char *cs_name = "Cs";
static const char *cs_irrep_names[] = {"A'","A''"};
static const int cs_generators[] = {1};
static const double cs_generator_reps[] = {1,-1};
static const int cs_irrep_degen[] = {1,1};
static const mat3x3 cs_ops[] = {Identity(),
                                Reflection(vec3(0,0,1))};
static const char *cs_op_names[] = {"E", "sxy"};

const PointGroup PointGroup::Cs = PointGroup(2, 2, 1, cs_name, cs_irrep_names,
                                             cs_generators, cs_generator_reps, cs_irrep_degen,
                                             cs_ops, cs_op_names);

/*
 * Ci
 */
static const char *ci_name = "Ci";
static const char *ci_irrep_names[] = {"Ag","Au"};
static const int ci_generators[] = {1};
static const double ci_generator_reps[] = {1,-1};
static const int ci_irrep_degen[] = {1,1};
static const mat3x3 ci_ops[] = {Identity(),
                                Inversion()};
static const char *ci_op_names[] = {"E", "i"};

const PointGroup PointGroup::Ci = PointGroup(2, 2, 1, ci_name, ci_irrep_names,
                                             ci_generators, ci_generator_reps, ci_irrep_degen,
                                             ci_ops, ci_op_names);

/*
 * Td
 */
static const char *td_name = "Td";
static const char *td_irrep_names[] = {"A1","A2","E","T1","T2"};
static const int td_generators[] = {17,1};
static const double td_generator_reps[] = { // S4z
                                            1,       // A1
                                           -1,       // A2
                                            1, 0,    // E
                                            0,-1,
                                            0,-1, 0, // T1
                                            1, 0, 0,
                                            0, 0, 1,
                                            0, 1, 0, // T2
                                           -1, 0, 0,
                                            0, 0,-1,
                                            // C3xyz
                                            1,                    // A1
                                            1,                    // A2
                                               -1.0/2,-sqrt(3)/2, // E
                                            sqrt(3)/2,    -1.0/2,
                                            0, 0, 1,              // T1
                                            1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1,              // T2
                                            1, 0, 0,
                                            0, 1, 0};
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

const PointGroup PointGroup::Td = PointGroup(24, 5, 2, td_name, td_irrep_names,
                                             td_generators, td_generator_reps, td_irrep_degen,
                                             td_ops, td_op_names);

/*
 * Oh
 */
static const char *oh_name = "Oh";
static const char *oh_irrep_names[] = {"A1g","A2g","Eg","T1g","T2g","A1u","A2u","Eu","T1u","T2u"};
static const int oh_generators[] = {19,1,24};
static const double oh_generator_reps[] = { // C4z
                                            1,       // A1g
                                           -1,       // A2g
                                            1, 0,    // Eg
                                            0,-1,
                                            0,-1, 0, // T1g
                                            1, 0, 0,
                                            0, 0, 1,
                                            0, 1, 0, // T2g
                                           -1, 0, 0,
                                            0, 0,-1,
                                            1,       // A1u
                                           -1,       // A2u
                                            1, 0,    // Eu
                                            0,-1,
                                            0,-1, 0, // T1u
                                            1, 0, 0,
                                            0, 0, 1,
                                            0, 1, 0, // T2u
                                           -1, 0, 0,
                                            0, 0,-1,
                                            // C3xyz
                                            1,                    // A1g
                                            1,                    // A2g
                                               -1.0/2,-sqrt(3)/2, // Eg
                                            sqrt(3)/2,    -1.0/2,
                                            0, 0, 1,              // T1g
                                            1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1,              // T2g
                                            1, 0, 0,
                                            0, 1, 0,
                                            1,                    // A1u
                                            1,                    // A2u
                                               -1.0/2,-sqrt(3)/2, // Eu
                                            sqrt(3)/2,    -1.0/2,
                                            0, 0, 1,              // T1u
                                            1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1,              // T2u
                                            1, 0, 0,
                                            0, 1, 0,
                                            // i
                                            1,       // A1g
                                            1,       // A2g
                                            1, 0,    // Eg
                                            0, 1,
                                            1, 0, 0, // T1g
                                            0, 1, 0,
                                            0, 0, 1,
                                            1, 0, 0, // T2g
                                            0, 1, 0,
                                            0, 0, 1,
                                           -1,       // A1u
                                           -1,       // A2u
                                           -1, 0,    // Eu
                                            0,-1,
                                           -1, 0, 0, // T1u
                                            0,-1, 0,
                                            0, 0,-1,
                                           -1, 0, 0, // T2u
                                            0,-1, 0,
                                            0, 0,-1};
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

const PointGroup PointGroup::Oh = PointGroup(48, 10, 3, oh_name, oh_irrep_names,
                                             oh_generators, oh_generator_reps, oh_irrep_degen,
                                             oh_ops, oh_op_names);

/*
 * Ih
 */
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
static const char *ih_irrep_names[] = {"Ag","T1g","T2g","Hg","Gg","Au","T1u","T2u","Hu","Gu"};
static const int ih_generators[] = {};
static const double ih_generator_reps[] = {};
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
                                    "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6", "S6",
                                    "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg", "sg"};

//const PointGroup PointGroup::Ih = PointGroup(120, 10, 0, ih_name, ih_irrep_names,
//                                             ih_generators, ih_generator_reps, ih_irrep_degen,
//                                             ih_ops, ih_op_names);

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

/*
 * C2
 */
static const char *c2_name = "C2";
static const char *c2_irrep_names[] = {"A1","A2"};
static const int c2_generators[] = {1};
static const double c2_generator_reps[] = {1,-1};
static const int c2_irrep_degen[] = {1,1};
static const mat3x3 c2_ops[] = {Identity(),
                                C<2>(vec3(0,0,1))};
static const char *c2_op_names[] = {"E", "C2z"};

const PointGroup PointGroup::C2 = PointGroup(2, 2, 1, c2_name, c2_irrep_names,
                                             c2_generators, c2_generator_reps, c2_irrep_degen,
                                             c2_ops, c2_op_names);

/*
 * C3
 */
static const char *c3_name = "C3";
static const char *c3_irrep_names[] = {"A","E+","E-"};
static const int c3_generators[] = {0,1,2};
static const double c3_generator_reps[] = {1,   sqrt(2),         0,  // E
                                           1,-1/sqrt(2), sqrt(1.5),  // C3
                                           1,-1/sqrt(2),-sqrt(1.5)}; // C3^2
static const int c3_irrep_degen[] = {1,1,1};
static const mat3x3 c3_ops[] = {Identity(),
                                C<3>(vec3( 0, 0, 1)),
                                C<3>(vec3( 0, 0,-1))};
static const char *c3_op_names[] = {"E", "C3z", "C3z^2"};

const PointGroup PointGroup::C3 = PointGroup(3, 3, 3, c3_name, c3_irrep_names,
                                             c3_generators, c3_generator_reps, c3_irrep_degen,
                                             c3_ops, c3_op_names);

/*
 * C4
 */
static const char *c4_name = "C4";
static const char *c4_irrep_names[] = {"A","B","E+","E-"};
static const int c4_generators[] = {0,1,2,3};
static const double c4_generator_reps[] = {1, 1, sqrt(2),       0,  // E
                                           1,-1,       0, sqrt(2),  // C4
                                           1, 1,-sqrt(2),       0,  // C2
                                           1,-1,       0,-sqrt(2)}; // C4^3
static const int c4_irrep_degen[] = {1,1,1,1};
static const mat3x3 c4_ops[] = {Identity(),
                                C<4>(vec3( 0, 0, 1)),
                                C<2>(vec3( 0, 0, 1)),
                                C<4>(vec3( 0, 0,-1))};
static const char *c4_op_names[] = {"E", "C4z", "C2z", "C4z^3"};

const PointGroup PointGroup::C4 = PointGroup(4, 4, 4, c4_name, c4_irrep_names,
                                             c4_generators, c4_generator_reps, c4_irrep_degen,
                                             c4_ops, c4_op_names);

/*
 * C5
 */
#define ce  sqrt(2)*cos(0.4*M_PI)
#define se  sqrt(2)*sin(0.4*M_PI)
#define c2e sqrt(2)*cos(0.8*M_PI)
#define s2e sqrt(2)*sin(0.8*M_PI)

static const char *c5_name = "C5";
static const char *c5_irrep_names[] = {"A","E1+","E1-","E2+","E2-"};
static const int c5_generators[] = {0,1,2,3,4};
static const double c5_generator_reps[] = {1, sqrt(2),   0, sqrt(2),   0,  // E
                                           1,      ce,  se,     c2e, s2e,  // C5
                                           1,     c2e, s2e,      ce, -se,  // C5^2
                                           1,     c2e,-s2e,      ce,  se,  // C5^3
                                           1,      ce, -se,     c2e,-s2e}; // C5^4
static const int c5_irrep_degen[] = {1,1,1,1,1};
static const mat3x3 c5_ops[] = {Identity(),
                                C<5>(vec3( 0, 0, 1)),
                                (C<5>(vec3( 0, 0, 1))^2),
                                (C<5>(vec3( 0, 0,-1))^2),
                                C<5>(vec3( 0, 0,-1))};
static const char *c5_op_names[] = {"E", "C5z", "C5z^2", "C5z^3", "C5z^4"};

const PointGroup PointGroup::C5 = PointGroup(5, 5, 5, c5_name, c5_irrep_names,
                                             c5_generators, c5_generator_reps, c5_irrep_degen,
                                             c5_ops, c5_op_names);

#undef ce
#undef se
#undef c2e
#undef s2e

/*
 * C6
 */
static const char *c6_name = "C6";
static const char *c6_irrep_names[] = {"A","B","E1+","E1-","E2+","E2-"};
static const int c6_generators[] = {0,1,2,3,4,5};
static const double c6_generator_reps[] = {1, 1,   sqrt(2),         0,   sqrt(2),         0,  // E
                                           1,-1, 1/sqrt(2), sqrt(1.5),-1/sqrt(2), sqrt(1.5),  // C6
                                           1, 1,-1/sqrt(2), sqrt(1.5),-1/sqrt(2),-sqrt(1.5),  // C3
                                           1,-1,  -sqrt(2),         0,   sqrt(2),         0,  // C2
                                           1, 1,-1/sqrt(2),-sqrt(1.5),-1/sqrt(2), sqrt(1.5),  // C3^2
                                           1,-1, 1/sqrt(2),-sqrt(1.5),-1/sqrt(2),-sqrt(1.5)}; // C6^5
static const int c6_irrep_degen[] = {1,1,1,1,1,1};
static const mat3x3 c6_ops[] = {Identity(),
                                C<6>(vec3( 0, 0, 1)),
                                C<3>(vec3( 0, 0, 1)),
                                C<2>(vec3( 0, 0, 1)),
                                C<3>(vec3( 0, 0,-1)),
                                C<6>(vec3( 0, 0,-1))};
static const char *c6_op_names[] = {"E", "C6z", "C3z", "C2z", "C3z^2", "C6z^5"};

const PointGroup PointGroup::C6 = PointGroup(6, 6, 6, c6_name, c6_irrep_names,
                                             c6_generators, c6_generator_reps, c6_irrep_degen,
                                             c6_ops, c6_op_names);

/*
 * C2v
 */
static const char *c2v_name = "C2v";
static const char *c2v_irrep_names[] = {"A1","A2","B1","B2"};
static const int c2v_generators[] = {1,2};
static const double c2v_generator_reps[] = {1, 1,-1,-1,  // C2z
                                            1,-1, 1,-1}; // sxz
static const int c2v_irrep_degen[] = {1,1,1,1};
static const mat3x3 c2v_ops[] = {Identity(),
                                 C<2>(vec3(0,0,1)),
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(1,0,0))};
static const char *c2v_op_names[] = {"E", "C2z", "sxz", "syz"};

const PointGroup PointGroup::C2v = PointGroup(4, 4, 2, c2v_name, c2v_irrep_names,
                                              c2v_generators, c2v_generator_reps, c2v_irrep_degen,
                                              c2v_ops, c2v_op_names);

/*
 * C3v
 */
static const char *c3v_name = "C3v";
static const char *c3v_irrep_names[] = {"A1","A2","E"};
static const int c3v_generators[] = {1,3};
static const double c3v_generator_reps[] = {// C3z
                                             1,                       // A1
                                             1,                       // A2
                                                   -0.5, -sqrt(0.75), // E
                                             sqrt(0.75),        -0.5,
                                            // sxz
                                             1,     // A1
                                            -1,     // A2
                                             1, 0,  // E
                                             0,-1};
static const int c3v_irrep_degen[] = {1,1,2};
static const mat3x3 c3v_ops[] = {Identity(),
                                 C<3>(vec3(0,0, 1)),
                                 C<3>(vec3(0,0,-1)),
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(0,1,0)*C<3>(vec3(0,0, 1))),
                                 Reflection(vec3(0,1,0)*C<3>(vec3(0,0,-1)))};
static const char *c3v_op_names[] = {"E", "C3z", "C3z^2", "sxz", "s2", "s3"};

const PointGroup PointGroup::C3v = PointGroup(6, 3, 2, c3v_name, c3v_irrep_names,
                                              c3v_generators, c3v_generator_reps, c3v_irrep_degen,
                                              c3v_ops, c3v_op_names);

/*
 * C4v
 */
static const char *c4v_name = "C4v";
static const char *c4v_irrep_names[] = {"A1","A2","B1","B2","E"};
static const int c4v_generators[] = {1,4};
static const double c4v_generator_reps[] = {// C4z
                                             1,     // A1
                                             1,     // A2
                                            -1,     // B1
                                            -1,     // B2
                                             0,-1,  // E
                                             1, 0,
                                             // sxz
                                             1,     // A1
                                            -1,     // A2
                                             1,     // B1
                                            -1,     // B2
                                             1, 0,  // E
                                             0,-1};
static const int c4v_irrep_degen[] = {1,1,1,1,2};
static const mat3x3 c4v_ops[] = {Identity(),
                                 C<4>(vec3(0,0, 1)),
                                 C<4>(vec3(0,0,-1)),
                                 C<2>(vec3(0,0, 1)),
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(1,0,0)),
                                 Reflection(vec3(1, 1,0)),
                                 Reflection(vec3(1,-1,0))};
static const char *c4v_op_names[] = {"E", "C4z", "C4z^3", "C2z", "sxz", "syz", "sx+y", "sx-y"};

const PointGroup PointGroup::C4v = PointGroup(8, 5, 2, c4v_name, c4v_irrep_names,
                                              c4v_generators, c4v_generator_reps, c4v_irrep_degen,
                                              c4v_ops, c4v_op_names);

/*
 * C5v
 */
static const char *c5v_name = "C5v";
static const char *c5v_irrep_names[] = {"A1","A2","E1","E2"};
static const int c5v_generators[] = {1,5};
static const double c5v_generator_reps[] = {// C5z
                                             1,                            // A1
                                             1,                            // A2
                                             cos(0.4*M_PI),-sin(0.4*M_PI), // E1
                                             sin(0.4*M_PI), cos(0.4*M_PI),
                                             cos(0.8*M_PI),-sin(0.8*M_PI), // E2
                                             sin(0.8*M_PI), cos(0.8*M_PI),
                                            // sxz
                                             1,     // A1
                                            -1,     // A2
                                             1, 0,  // E1
                                             0,-1,
                                             1, 0,  // E2
                                             0,-1};
static const int c5v_irrep_degen[] = {1,1,2,2};
static const mat3x3 c5v_ops[] = {Identity(),
                                 C<5>(vec3(0,0, 1)),
                                 C<5>(vec3(0,0,-1)),
                                 C<5>(vec3(0,0, 1))^2,
                                 C<5>(vec3(0,0,-1))^2,
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(0,1,0)*C<5>(vec3(0,0, 1))),
                                 Reflection(vec3(0,1,0)*C<5>(vec3(0,0,-1))),
                                 Reflection(vec3(0,1,0)*(C<5>(vec3(0,0, 1))^2)),
                                 Reflection(vec3(0,1,0)*(C<5>(vec3(0,0,-1))^2))};
static const char *c5v_op_names[] = {"E", "C5z", "C5z^4", "C5z^2", "C5z^3", "sxz", "s2", "s3", "s4", "s5"};

const PointGroup PointGroup::C5v = PointGroup(10, 4, 2, c5v_name, c5v_irrep_names,
                                              c5v_generators, c5v_generator_reps, c5v_irrep_degen,
                                              c5v_ops, c5v_op_names);

/*
 * C6v
 */
static const char *c6v_name = "C6v";
static const char *c6v_irrep_names[] = {"A1","A2","B1","B2","E1","E2"};
static const int c6v_generators[] = {1,6};
static const double c6v_generator_reps[] = {// C6z
                                             1,                      // A1
                                             1,                      // A2
                                            -1,                      // B1
                                            -1,                      // B2
                                                    0.5,-sqrt(0.75), // E1
                                             sqrt(0.75),        0.5,
                                                   -0.5,-sqrt(0.75), // E2
                                             sqrt(0.75),       -0.5,
                                            // sxz
                                             1,     // A1
                                            -1,     // A2
                                             1,     // B1
                                            -1,     // B2
                                             1, 0,  // E1
                                             0,-1,
                                             1, 0,  // E2
                                             0,-1};
static const int c6v_irrep_degen[] = {1,1,1,1,2,2};
static const mat3x3 c6v_ops[] = {Identity(),
                                 C<6>(vec3(0,0, 1)),
                                 C<6>(vec3(0,0,-1)),
                                 C<3>(vec3(0,0, 1)),
                                 C<3>(vec3(0,0,-1)),
                                 C<2>(vec3(0,0, 1)),
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(0,1,0)*C<6>(vec3(0,0,1))),
                                 Reflection(vec3(0,1,0)*C<3>(vec3(0,0,1))),
                                 Reflection(vec3(0,1,0)*C<12>(vec3(0,0,1))),
                                 Reflection(vec3(1,0,0)),
                                 Reflection(vec3(0,1,0)*(C<12>(vec3(0,0,1))^5))};
static const char *c6v_op_names[] = {"E", "C6z", "C6z^5", "C3z", "C3z^2", "C2z", "sxz", "s2", "s3", "syz", "s5", "s6"};

const PointGroup PointGroup::C6v = PointGroup(12, 6, 2, c6v_name, c6v_irrep_names,
                                              c6v_generators, c6v_generator_reps, c6v_irrep_degen,
                                              c6v_ops, c6v_op_names);

/*
 * C2h
 */
static const char *c2h_name = "C2h";
static const char *c2h_irrep_names[] = {"Ag","Bg","Au","Bu"};
static const int c2h_generators[] = {1,2};
static const double c2h_generator_reps[] = {1,-1, 1,-1,  // C2z
                                            1, 1,-1,-1}; // i
static const int c2h_irrep_degen[] = {1,1,1,1};
static const mat3x3 c2h_ops[] = {Identity(),
                                 C<2>(vec3(0,0,1)),
                                 Inversion(),
                                 Reflection(vec3(0,0,1))};
static const char *c2h_op_names[] = {"E", "C2z", "i", "sxy"};

const PointGroup PointGroup::C2h = PointGroup(4, 4, 2, c2h_name, c2h_irrep_names,
                                              c2h_generators, c2h_generator_reps, c2h_irrep_degen,
                                              c2h_ops, c2h_op_names);

/*
 * C3h
 */
static const char *c3h_name = "C3h";
static const char *c3h_irrep_names[] = {"A'","E'+","E'-","A''","E''+","E''-"};
static const int c3h_generators[] = {0,1,2,3,4,5};
static const double c3h_generator_reps[] = {1,   sqrt(2),         0, 1,   sqrt(2),         0,  // E
                                            1,-1/sqrt(2), sqrt(1.5), 1,-1/sqrt(2), sqrt(1.5),  // C3
                                            1,-1/sqrt(2),-sqrt(1.5), 1,-1/sqrt(2),-sqrt(1.5),  // C3^2
                                            1,   sqrt(2),         0,-1,  -sqrt(2),         0,  // sxy
                                            1,-1/sqrt(2), sqrt(1.5),-1, 1/sqrt(2),-sqrt(1.5),  // S3
                                            1,-1/sqrt(2),-sqrt(1.5),-1, 1/sqrt(2), sqrt(1.5)}; // S3^5
static const int c3h_irrep_degen[] = {1,1,1,1,1,1};
static const mat3x3 c3h_ops[] = {Identity(),
                                 C<3>(vec3( 0, 0, 1)),
                                 C<3>(vec3( 0, 0,-1)),
                                 Reflection(vec3(0,0,1)),
                                 S<3>(vec3( 0, 0, 1)),
                                 S<3>(vec3( 0, 0,-1))};
static const char *c3h_op_names[] = {"E", "C3z", "C3z^2", "sxy", "S3z", "S3z^5"};

const PointGroup PointGroup::C3h = PointGroup(6, 6, 6, c3h_name, c3h_irrep_names,
                                              c3h_generators, c3h_generator_reps, c3h_irrep_degen,
                                              c3h_ops, c3h_op_names);

/*
 * C4h
 */
static const char *c4h_name = "C4h";
static const char *c4h_irrep_names[] = {"Ag","Bg","Eg+","Eg-","Au","Bu","Eu+","Eu-"};
static const int c4h_generators[] = {0,1,2,3,4,5,6,7};
static const double c4h_generator_reps[] = {1, 1, sqrt(2),       0, 1, 1, sqrt(2),       0,  // E
                                            1,-1,       0, sqrt(2), 1,-1,       0, sqrt(2),  // C4
                                            1, 1,-sqrt(2),       0, 1, 1,-sqrt(2),       0,  // C2
                                            1,-1,       0,-sqrt(2), 1,-1,       0,-sqrt(2),  // C4^3
                                            1, 1, sqrt(2),       0,-1,-1,-sqrt(2),       0,  // i
                                            1,-1,       0, sqrt(2),-1, 1,       0,-sqrt(2),  // S4^3
                                            1, 1,-sqrt(2),       0,-1,-1, sqrt(2),       0,  // sxy
                                            1,-1,       0,-sqrt(2),-1, 1,       0, sqrt(2)}; // S4
static const int c4h_irrep_degen[] = {1,1,1,1,1,1,1,1};
static const mat3x3 c4h_ops[] = {Identity(),
                                 C<4>(vec3( 0, 0, 1)),
                                 C<2>(vec3( 0, 0, 1)),
                                 C<4>(vec3( 0, 0,-1)),
                                 Inversion(),
                                 S<4>(vec3( 0, 0,-1)),
                                 Reflection(vec3(0,0,1)),
                                 S<4>(vec3( 0, 0, 1))};
static const char *c4h_op_names[] = {"E", "C4z", "C2z", "C4z^3", "i", "S4z^3", "sxy", "S4z"};

const PointGroup PointGroup::C4h = PointGroup(8, 8, 8, c4h_name, c4h_irrep_names,
                                              c4h_generators, c4h_generator_reps, c4h_irrep_degen,
                                              c4h_ops, c4h_op_names);

/*
 * C5h
 */
#define ce  sqrt(2)*cos(0.4*M_PI)
#define se  sqrt(2)*sin(0.4*M_PI)
#define c2e sqrt(2)*cos(0.8*M_PI)
#define s2e sqrt(2)*sin(0.8*M_PI)

static const char *c5h_name = "C5h";
static const char *c5h_irrep_names[] = {"A'","E1'+","E1'-","E2'+","E2'-","A''","E1''+","E1''-","E2''+","E2''-"};
static const int c5h_generators[] = {0,1,2,3,4,5,6,7,8,9};
static const double c5h_generator_reps[] = {1, sqrt(2),   0, sqrt(2),   0, 1, sqrt(2),   0, sqrt(2),   0,  // E
                                            1,      ce,  se,     c2e, s2e, 1,      ce,  se,     c2e, s2e,  // C5
                                            1,     c2e, s2e,      ce, -se, 1,     c2e, s2e,      ce, -se,  // C5^2
                                            1,     c2e,-s2e,      ce,  se, 1,     c2e,-s2e,      ce,  se,  // C5^3
                                            1,      ce, -se,     c2e,-s2e, 1,      ce, -se,     c2e,-s2e,  // C5^4
                                            1, sqrt(2),   0, sqrt(2),   0,-1,-sqrt(2),   0,-sqrt(2),   0,  // sxy
                                            1,      ce,  se,     c2e, s2e,-1,     -ce, -se,    -c2e,-s2e,  // S5
                                            1,     c2e, s2e,      ce, -se,-1,    -c2e,-s2e,     -ce,  se,  // S5^7
                                            1,     c2e,-s2e,      ce,  se,-1,    -c2e, s2e,     -ce, -se,  // S5^3
                                            1,      ce, -se,     c2e,-s2e,-1,     -ce,  se,    -c2e, s2e}; // S5^9
static const int c5h_irrep_degen[] = {1,1,1,1,1,1,1,1,1,1};
static const mat3x3 c5h_ops[] = {Identity(),
                                 C<5>(vec3( 0, 0, 1)),
                                 (C<5>(vec3( 0, 0, 1))^2),
                                 (C<5>(vec3( 0, 0,-1))^2),
                                 C<5>(vec3( 0, 0,-1)),
                                 Reflection(vec3(0,0,1)),
                                 S<5>(vec3( 0, 0, 1)),
                                 S<5>(vec3( 0, 0,-1))^3,
                                 S<5>(vec3( 0, 0, 1))^3,
                                 S<5>(vec3( 0, 0,-1))};
static const char *c5h_op_names[] = {"E", "C5z", "C5z^2", "C5z^3", "C5z^4", "sxy", "S5z", "S5z^7", "S5z^3", "S5z^9"};

const PointGroup PointGroup::C5h = PointGroup(10, 10, 10, c5h_name, c5h_irrep_names,
                                              c5h_generators, c5h_generator_reps, c5h_irrep_degen,
                                              c5h_ops, c5h_op_names);

#undef ce
#undef se
#undef c2e
#undef s2e

/*
 * C6h
 */
static const char *c6h_name = "C6h";
static const char *c6h_irrep_names[] = {"Ag","Bg","E1g+","E1g-","E2g+","E2g-","Au","Bu","E1u+","E1u-","E2u+","E2u-"};
static const int c6h_generators[] = {0,1,2,3,4,5,6,7,8,9,10,11};
static const double c6h_generator_reps[] = {1, 1,   sqrt(2),         0,   sqrt(2),         0, 1, 1,   sqrt(2),         0,   sqrt(2),         0,  // E
                                            1,-1, 1/sqrt(2), sqrt(1.5),-1/sqrt(2), sqrt(1.5), 1,-1, 1/sqrt(2), sqrt(1.5),-1/sqrt(2), sqrt(1.5),  // C6
                                            1, 1,-1/sqrt(2), sqrt(1.5),-1/sqrt(2),-sqrt(1.5), 1, 1,-1/sqrt(2), sqrt(1.5),-1/sqrt(2),-sqrt(1.5),  // C3
                                            1,-1,  -sqrt(2),         0,   sqrt(2),         0, 1,-1,  -sqrt(2),         0,   sqrt(2),         0,  // C2
                                            1, 1,-1/sqrt(2),-sqrt(1.5),-1/sqrt(2), sqrt(1.5), 1, 1,-1/sqrt(2),-sqrt(1.5),-1/sqrt(2), sqrt(1.5),  // C3^2
                                            1,-1, 1/sqrt(2),-sqrt(1.5),-1/sqrt(2),-sqrt(1.5), 1,-1, 1/sqrt(2),-sqrt(1.5),-1/sqrt(2),-sqrt(1.5),  // C6^5
                                            1, 1,   sqrt(2),         0,   sqrt(2),         0,-1,-1,  -sqrt(2),         0,  -sqrt(2),         0,  // i
                                            1,-1, 1/sqrt(2), sqrt(1.5),-1/sqrt(2), sqrt(1.5),-1, 1,-1/sqrt(2),-sqrt(1.5), 1/sqrt(2),-sqrt(1.5),  // S3^5
                                            1, 1,-1/sqrt(2), sqrt(1.5),-1/sqrt(2),-sqrt(1.5),-1,-1, 1/sqrt(2),-sqrt(1.5), 1/sqrt(2), sqrt(1.5),  // S6^5
                                            1,-1,  -sqrt(2),         0,   sqrt(2),         0,-1, 1,   sqrt(2),         0,  -sqrt(2),         0,  // sxy
                                            1, 1,-1/sqrt(2),-sqrt(1.5),-1/sqrt(2), sqrt(1.5),-1,-1, 1/sqrt(2), sqrt(1.5), 1/sqrt(2),-sqrt(1.5),  // S6
                                            1,-1, 1/sqrt(2),-sqrt(1.5),-1/sqrt(2),-sqrt(1.5),-1, 1,-1/sqrt(2), sqrt(1.5), 1/sqrt(2), sqrt(1.5)}; // S3
static const int c6h_irrep_degen[] = {1,1,1,1,1,1,1,1,1,1,1,1};
static const mat3x3 c6h_ops[] = {Identity(),
                                 C<6>(vec3( 0, 0, 1)),
                                 C<3>(vec3( 0, 0, 1)),
                                 C<2>(vec3( 0, 0, 1)),
                                 C<3>(vec3( 0, 0,-1)),
                                 C<6>(vec3( 0, 0,-1)),
                                 Inversion(),
                                 S<3>(vec3( 0, 0,-1)),
                                 S<6>(vec3( 0, 0,-1)),
                                 Reflection(vec3(0,0,1)),
                                 S<6>(vec3( 0, 0, 1)),
                                 S<3>(vec3( 0, 0, 1))};
static const char *c6h_op_names[] = {"E", "C6z", "C3z", "C2z", "C3z^2", "C6z^5", "i", "S3^5", "S6^5", "sxy", "S6", "S3"};

const PointGroup PointGroup::C6h = PointGroup(12, 12, 12, c6h_name, c6h_irrep_names,
                                              c6h_generators, c6h_generator_reps, c6h_irrep_degen,
                                              c6h_ops, c6h_op_names);

/*
 * D2
 */
static const char *d2_name = "D2";
static const char *d2_irrep_names[] = {"A","B1","B2","B3"};
static const int d2_generators[] = {1,2};
static const double d2_generator_reps[] = {1, 1,-1,-1,  // C2z
                                           1,-1, 1,-1}; // C2y
static const int d2_irrep_degen[] = {1,1,1,1};
static const mat3x3 d2_ops[] = {Identity(),
                                C<2>(vec3(0,0,1)),
                                C<2>(vec3(0,1,0)),
                                C<2>(vec3(1,0,0))};
static const char *d2_op_names[] = {"E", "C2z", "C2y", "C2x"};

const PointGroup PointGroup::D2 = PointGroup(4, 4, 2, d2_name, d2_irrep_names,
                                             d2_generators, d2_generator_reps, d2_irrep_degen,
                                             d2_ops, d2_op_names);

/*
 * D3
 */
static const char *d3_name = "D3";
static const char *d3_irrep_names[] = {"A1","A2","E"};
static const int d3_generators[] = {1,3};
static const double d3_generator_reps[] = {// C3z
                                            1,                       // A1
                                            1,                       // A2
                                                  -0.5, -sqrt(0.75), // E
                                            sqrt(0.75),        -0.5,
                                           // C2x
                                            1,     // A1
                                           -1,     // A2
                                            1, 0,  // E
                                            0,-1};
static const int d3_irrep_degen[] = {1,1,2};
static const mat3x3 d3_ops[] = {Identity(),
                                C<3>(vec3(0,0, 1)),
                                C<3>(vec3(0,0,-1)),
                                C<2>(vec3(1,0,0)),
                                C<2>(vec3(1,0,0)*C<3>(vec3(0,0, 1))),
                                C<2>(vec3(1,0,0)*C<3>(vec3(0,0,-1)))};
static const char *d3_op_names[] = {"E", "C3z", "C3z^2", "C2x", "C22", "C23"};

const PointGroup PointGroup::D3 = PointGroup(6, 3, 2, d3_name, d3_irrep_names,
                                             d3_generators, d3_generator_reps, d3_irrep_degen,
                                             d3_ops, d3_op_names);

/*
 * D4
 */
static const char *d4_name = "D4";
static const char *d4_irrep_names[] = {"A1","A2","B1","B2","E"};
static const int d4_generators[] = {1,4};
static const double d4_generator_reps[] = {// C4z
                                           1,     // A1
                                           1,     // A2
                                          -1,     // B1
                                          -1,     // B2
                                           0,-1,  // E
                                           1, 0,
                                           // C2x
                                           1,     // A1
                                          -1,     // A2
                                           1,     // B1
                                          -1,     // B2
                                           1, 0,  // E
                                           0,-1};
static const int d4_irrep_degen[] = {1,1,1,1,2};
static const mat3x3 d4_ops[] = {Identity(),
                                C<4>(vec3(0,0, 1)),
                                C<4>(vec3(0,0,-1)),
                                C<2>(vec3(0,0, 1)),
                                C<2>(vec3(1,0,0)),
                                C<2>(vec3(0,1,0)),
                                C<2>(vec3(1, 1,0)),
                                C<2>(vec3(1,-1,0))};
static const char *d4_op_names[] = {"E", "C4z", "C4z^3", "C2z", "C2x", "C2y", "C2x+y", "C2x-y"};

const PointGroup PointGroup::D4 = PointGroup(8, 5, 2, d4_name, d4_irrep_names,
                                             d4_generators, d4_generator_reps, d4_irrep_degen,
                                             d4_ops, d4_op_names);

/*
 * D5
 */
static const char *d5_name = "D5";
static const char *d5_irrep_names[] = {"A1","A2","E1","E2"};
static const int d5_generators[] = {1,5};
static const double d5_generator_reps[] = {// C5z
                                            1,                            // A1
                                            1,                            // A2
                                            cos(0.4*M_PI),-sin(0.4*M_PI), // E1
                                            sin(0.4*M_PI), cos(0.4*M_PI),
                                            cos(0.8*M_PI),-sin(0.8*M_PI), // E2
                                            sin(0.8*M_PI), cos(0.8*M_PI),
                                           // C2x
                                            1,     // A1
                                           -1,     // A2
                                            1, 0,  // E1
                                            0,-1,
                                            1, 0,  // E2
                                            0,-1};
static const int d5_irrep_degen[] = {1,1,2,2};
static const mat3x3 d5_ops[] = {Identity(),
                                C<5>(vec3(0,0, 1)),
                                C<5>(vec3(0,0,-1)),
                                C<5>(vec3(0,0, 1))^2,
                                C<5>(vec3(0,0,-1))^2,
                                C<2>(vec3(1,0,0)),
                                C<2>(vec3(1,0,0)*C<5>(vec3(0,0, 1))),
                                C<2>(vec3(1,0,0)*C<5>(vec3(0,0,-1))),
                                C<2>(vec3(1,0,0)*(C<5>(vec3(0,0, 1))^2)),
                                C<2>(vec3(1,0,0)*(C<5>(vec3(0,0,-1))^2))};
static const char *d5_op_names[] = {"E", "C5z", "C5z^4", "C5z^2", "C5z^3", "C2x", "C22", "C23", "C24", "C25"};

const PointGroup PointGroup::D5 = PointGroup(10, 4, 2, d5_name, d5_irrep_names,
                                             d5_generators, d5_generator_reps, d5_irrep_degen,
                                             d5_ops, d5_op_names);

/*
 * D6
 */
static const char *d6_name = "D6";
static const char *d6_irrep_names[] = {"A1","A2","B1","B2","E1","E2"};
static const int d6_generators[] = {1,6};
static const double d6_generator_reps[] = {// C6z
                                            1,                      // A1
                                            1,                      // A2
                                           -1,                      // B1
                                           -1,                      // B2
                                                   0.5,-sqrt(0.75), // E1
                                            sqrt(0.75),        0.5,
                                                  -0.5,-sqrt(0.75), // E2
                                            sqrt(0.75),       -0.5,
                                           // C2x
                                            1,     // A1
                                           -1,     // A2
                                            1,     // B1
                                           -1,     // B2
                                            1, 0,  // E1
                                            0,-1,
                                            1, 0,  // E2
                                            0,-1};
static const int d6_irrep_degen[] = {1,1,1,1,2,2};
static const mat3x3 d6_ops[] = {Identity(),
                                C<6>(vec3(0,0, 1)),
                                C<6>(vec3(0,0,-1)),
                                C<3>(vec3(0,0, 1)),
                                C<3>(vec3(0,0,-1)),
                                C<2>(vec3(0,0, 1)),
                                C<2>(vec3(1,0,0)),
                                C<2>(vec3(1,0,0)*C<6>(vec3(0,0, 1))),
                                C<2>(vec3(1,0,0)*C<6>(vec3(0,0,-1))),
                                C<2>(vec3(0,1,0)),
                                C<2>(vec3(0,1,0)*C<6>(vec3(0,0, 1))),
                                C<2>(vec3(0,1,0)*C<6>(vec3(0,0,-1)))};
static const char *d6_op_names[] = {"E", "C6z", "C6z^5", "C3z", "C3z^2", "C2z", "C2x", "C22", "C23", "C2y", "C25", "C26"};

const PointGroup PointGroup::D6 = PointGroup(12, 6, 2, d6_name, d6_irrep_names,
                                             d6_generators, d6_generator_reps, d6_irrep_degen,
                                             d6_ops, d6_op_names);

/*
 * D2h
 */
static const char *d2h_name = "D2h";
static const char *d2h_irrep_names[] = {"Ag","B1g","B2g","B3g","Au","B1u","B2u","B3u"};
static const int d2h_generators[] = {1,2,4};
static const double d2h_generator_reps[] = {1, 1,-1,-1, 1, 1,-1,-1,  // C2z
                                            1,-1, 1,-1, 1,-1, 1,-1,  // C2y
                                            1, 1, 1, 1,-1,-1,-1,-1}; // i
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

const PointGroup PointGroup::D2h = PointGroup(8, 8, 3, d2h_name, d2h_irrep_names,
                                              d2h_generators, d2h_generator_reps, d2h_irrep_degen,
                                              d2h_ops, d2h_op_names);

/*
 * D3h
 */
static const char *d3h_name = "D3h";
static const char *d3h_irrep_names[] = {"A1'","A2'","E'","A1''","A2''","E''"};
static const int d3h_generators[] = {1,3,6};
static const double d3h_generator_reps[] = {// C3z
                                             1,                       // A1'
                                             1,                       // A2'
                                                   -0.5, -sqrt(0.75), // E'
                                             sqrt(0.75),        -0.5,
                                             1,                       // A1''
                                             1,                       // A2''
                                                   -0.5, -sqrt(0.75), // E''
                                             sqrt(0.75),        -0.5,
                                            // C2x
                                             1,     // A1'
                                            -1,     // A2'
                                             1, 0,  // E'
                                             0,-1,
                                             1,     // A1''
                                            -1,     // A2''
                                             1, 0,  // E''
                                             0,-1,
                                            // sxy
                                             1,     // A1'
                                             1,     // A2'
                                             1, 0,  // E'
                                             0, 1,
                                            -1,     // A1''
                                            -1,     // A2''
                                            -1, 0,  // E''
                                             0,-1};
static const int d3h_irrep_degen[] = {1,1,2,1,1,2};
static const mat3x3 d3h_ops[] = {Identity(),
                                 C<3>(vec3(0,0, 1)),
                                 C<3>(vec3(0,0,-1)),
                                 C<2>(vec3(1,0,0)),
                                 C<2>(vec3(1,0,0)*C<3>(vec3(0,0, 1))),
                                 C<2>(vec3(1,0,0)*C<3>(vec3(0,0,-1))),
                                 Reflection(vec3(0,0,1)),
                                 S<3>(vec3(0,0, 1)),
                                 S<3>(vec3(0,0,-1)),
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(0,1,0)*C<3>(vec3(0,0, 1))),
                                 Reflection(vec3(0,1,0)*C<3>(vec3(0,0,-1)))};
static const char *d3h_op_names[] = {"E", "C3z", "C3z^2", "C2x", "C22", "C23",
                                     "sxy", "S3z", "S3z^5", "sxz", "s2", "s3"};

const PointGroup PointGroup::D3h = PointGroup(12, 6, 3, d3h_name, d3h_irrep_names,
                                              d3h_generators, d3h_generator_reps, d3h_irrep_degen,
                                              d3h_ops, d3h_op_names);

/*
 * D4h
 */
static const char *d4h_name = "D4h";
static const char *d4h_irrep_names[] = {"A1g","A2g","B1g","B2g","Eg","A1u","A2u","B1u","B2u","Eu"};
static const int d4h_generators[] = {1,4,8};
static const double d4h_generator_reps[] = {// C4z
                                             1,     // A1g
                                             1,     // A2g
                                            -1,     // B1g
                                            -1,     // B2g
                                             0,-1,  // Eg
                                             1, 0,
                                             1,     // A1u
                                             1,     // A2u
                                            -1,     // B1u
                                            -1,     // B2u
                                             0,-1,  // Eu
                                             1, 0,
                                            // C2x
                                             1,     // A1g
                                            -1,     // A2g
                                             1,     // B1g
                                            -1,     // B2g
                                             1, 0,  // Eg
                                             0,-1,
                                             1,     // A1u
                                            -1,     // A2u
                                             1,     // B1u
                                            -1,     // B2u
                                             1, 0,  // Eu
                                             0,-1,
                                            // i
                                             1,     // A1g
                                             1,     // A2g
                                             1,     // B1g
                                             1,     // B2g
                                             1, 0,  // Eg
                                             0, 1,
                                            -1,     // A1u
                                            -1,     // A2u
                                            -1,     // B1u
                                            -1,     // B2u
                                            -1, 0,  // Eu
                                             0,-1};
static const int d4h_irrep_degen[] = {1,1,1,1,2,1,1,1,1,2};
static const mat3x3 d4h_ops[] = {Identity(),
                                 C<4>(vec3(0,0, 1)),
                                 C<4>(vec3(0,0,-1)),
                                 C<2>(vec3(0,0, 1)),
                                 C<2>(vec3(1,0,0)),
                                 C<2>(vec3(0,1,0)),
                                 C<2>(vec3(1, 1,0)),
                                 C<2>(vec3(1,-1,0)),
                                 Inversion(),
                                 S<4>(vec3(0,0, 1)),
                                 S<4>(vec3(0,0,-1)),
                                 Reflection(vec3(0,0,1)),
                                 Reflection(vec3(1,0,0)),
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(1, 1,0)),
                                 Reflection(vec3(1,-1,0))};
static const char *d4h_op_names[] = {"E", "C4z", "C4z^3", "C2z", "C2x", "C2y", "C2x+y", "C2x-y",
                                     "i", "S4z", "S4z^7", "sxy", "syz", "sxz", "sx+y", "sx-y"};

const PointGroup PointGroup::D4h = PointGroup(16, 10, 3, d4h_name, d4h_irrep_names,
                                              d4h_generators, d4h_generator_reps, d4h_irrep_degen,
                                              d4h_ops, d4h_op_names);

/*
 * D5h
 */
static const char *d5h_name = "D5h";
static const char *d5h_irrep_names[] = {"A1'","A2'","E1'","E2'","A1''","A2''","E1''","E2''"};
static const int d5h_generators[] = {1,5,10};
static const double d5h_generator_reps[] = {// C5z
                                             1,                            // A1'
                                             1,                            // A2'
                                             cos(0.4*M_PI),-sin(0.4*M_PI), // E1'
                                             sin(0.4*M_PI), cos(0.4*M_PI),
                                             cos(0.8*M_PI),-sin(0.8*M_PI), // E2'
                                             sin(0.8*M_PI), cos(0.8*M_PI),
                                             1,                            // A1''
                                             1,                            // A2''
                                             cos(0.4*M_PI),-sin(0.4*M_PI), // E1''
                                             sin(0.4*M_PI), cos(0.4*M_PI),
                                             cos(0.8*M_PI),-sin(0.8*M_PI), // E2''
                                             sin(0.8*M_PI), cos(0.8*M_PI),
                                            // C2x
                                             1,     // A1'
                                            -1,     // A2'
                                             1, 0,  // E1'
                                             0,-1,
                                             1, 0,  // E2'
                                             0,-1,
                                             1,     // A1''
                                            -1,     // A2''
                                             1, 0,  // E1''
                                             0,-1,
                                             1, 0,  // E2''
                                             0,-1,
                                            // sxy
                                             1,     // A1'
                                             1,     // A2'
                                             1, 0,  // E1'
                                             0, 1,
                                             1, 0,  // E2'
                                             0, 1,
                                            -1,     // A1''
                                            -1,     // A2''
                                            -1, 0,  // E1''
                                             0,-1,
                                            -1, 0,  // E2''
                                             0,-1};
static const int d5h_irrep_degen[] = {1,1,2,2,1,1,2,2};
static const mat3x3 d5h_ops[] = {Identity(),
                                 C<5>(vec3(0,0, 1)),
                                 C<5>(vec3(0,0,-1)),
                                 C<5>(vec3(0,0, 1))^2,
                                 C<5>(vec3(0,0,-1))^2,
                                 C<2>(vec3(1,0,0)),
                                 C<2>(vec3(1,0,0)*C<5>(vec3(0,0, 1))),
                                 C<2>(vec3(1,0,0)*C<5>(vec3(0,0,-1))),
                                 C<2>(vec3(1,0,0)*(C<5>(vec3(0,0, 1))^2)),
                                 C<2>(vec3(1,0,0)*(C<5>(vec3(0,0,-1))^2)),
                                 Reflection(vec3(0,0,1)),
                                 S<5>(vec3(0,0, 1)),
                                 S<5>(vec3(0,0,-1)),
                                 S<5>(vec3(0,0, 1))^3,
                                 S<5>(vec3(0,0,-1))^3,
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(0,1,0)*C<5>(vec3(0,0, 1))),
                                 Reflection(vec3(0,1,0)*C<5>(vec3(0,0,-1))),
                                 Reflection(vec3(0,1,0)*(C<5>(vec3(0,0, 1))^2)),
                                 Reflection(vec3(0,1,0)*(C<5>(vec3(0,0,-1))^2))};
static const char *d5h_op_names[] = {"E", "C5z", "C5z^4", "C5z^2", "C5z^3", "C2x", "C22", "C23", "C24", "C25",
                                     "sxy", "S5z", "S5z^9", "S5z^3", "S5z^7", "sxz", "s2", "s3", "s4", "s5"};

const PointGroup PointGroup::D5h = PointGroup(20, 8, 3, d5h_name, d5h_irrep_names,
                                              d5h_generators, d5h_generator_reps, d5h_irrep_degen,
                                              d5h_ops, d5h_op_names);

/*
 * D6h
 */
static const char *d6h_name = "D6h";
static const char *d6h_irrep_names[] = {"A1g","A2g","B1g","B2g","E1g","E2g","A1u","A2u","B1u","B2u","E1u","E2u"};
static const int d6h_generators[] = {1,6,12};
static const double d6h_generator_reps[] = {// C6z
                                             1,                      // A1g
                                             1,                      // A2g
                                            -1,                      // B1g
                                            -1,                      // B2g
                                                    0.5,-sqrt(0.75), // E1g
                                             sqrt(0.75),        0.5,
                                                   -0.5,-sqrt(0.75), // E2g
                                             sqrt(0.75),       -0.5,
                                             1,                      // A1u
                                             1,                      // A2u
                                            -1,                      // B1u
                                            -1,                      // B2u
                                                    0.5,-sqrt(0.75), // E1u
                                             sqrt(0.75),        0.5,
                                                   -0.5,-sqrt(0.75), // E2u
                                             sqrt(0.75),       -0.5,
                                            // C2x
                                             1,     // A1g
                                            -1,     // A2g
                                             1,     // B1g
                                            -1,     // B2g
                                             1, 0,  // E1g
                                             0,-1,
                                             1, 0,  // E2g
                                             0,-1,
                                             1,     // A1u
                                            -1,     // A2u
                                             1,     // B1u
                                            -1,     // B2u
                                             1, 0,  // E1u
                                             0,-1,
                                             1, 0,  // E2u
                                             0,-1,
                                            // i
                                             1,     // A1g
                                             1,     // A2g
                                             1,     // B1g
                                             1,     // B2g
                                             1, 0,  // E1g
                                             0, 1,
                                             1, 0,  // E2g
                                             0, 1,
                                            -1,     // A1u
                                            -1,     // A2u
                                            -1,     // B1u
                                            -1,     // B2u
                                            -1, 0,  // E1u
                                             0,-1,
                                            -1, 0,  // E2u
                                             0,-1};
static const int d6h_irrep_degen[] = {1,1,1,1,2,2,1,1,1,1,2,2};
static const mat3x3 d6h_ops[] = {Identity(),
                                 C<6>(vec3(0,0, 1)),
                                 C<6>(vec3(0,0,-1)),
                                 C<3>(vec3(0,0, 1)),
                                 C<3>(vec3(0,0,-1)),
                                 C<2>(vec3(0,0, 1)),
                                 C<2>(vec3(1,0,0)),
                                 C<2>(vec3(1,0,0)*C<6>(vec3(0,0,1))),
                                 C<2>(vec3(1,0,0)*C<3>(vec3(0,0,1))),
                                 C<2>(vec3(1,0,0)*C<12>(vec3(0,0,1))),
                                 C<2>(vec3(0,1,0)),
                                 C<2>(vec3(1,0,0)*(C<12>(vec3(0,0,1))^5)),
                                 Inversion(),
                                 S<3>(vec3(0,0, 1)),
                                 S<3>(vec3(0,0,-1)),
                                 S<6>(vec3(0,0, 1)),
                                 S<6>(vec3(0,0,-1)),
                                 Reflection(vec3(0,0,1)),
                                 Reflection(vec3(1,0,0)),
                                 Reflection(vec3(1,0,0)*C<6>(vec3(0,0, 1))),
                                 Reflection(vec3(1,0,0)*C<6>(vec3(0,0,-1))),
                                 Reflection(vec3(0,1,0)),
                                 Reflection(vec3(0,1,0)*C<6>(vec3(0,0, 1))),
                                 Reflection(vec3(0,1,0)*C<6>(vec3(0,0,-1)))};
static const char *d6h_op_names[] = {"E", "C6z", "C6z^5", "C3z", "C3z^2", "C2z", "C2x", "C22", "C23", "C2y", "C25", "C26",
                                     "i", "S3z", "S3z^5", "S6z", "S6z^5", "sxy", "syz", "s2", "s3", "sxz", "s5", "s6"};

const PointGroup PointGroup::D6h = PointGroup(24, 12, 3, d6h_name, d6h_irrep_names,
                                              d6h_generators, d6h_generator_reps, d6h_irrep_degen,
                                              d6h_ops, d6h_op_names);

/*
 * D2d
 */
static const char *d2d_name = "D2d";
static const char *d2d_irrep_names[] = {"A1","A2","B1","B2","E"};
static const int d2d_generators[] = {1,4};
static const double d2d_generator_reps[] = {// S4z
                                            1,     // A1
                                            1,     // A2
                                           -1,     // B1
                                           -1,     // B2
                                            0,-1,  // E
                                            1, 0,
                                            // C2x
                                            1,     // A1
                                           -1,     // A2
                                            1,     // B1
                                           -1,     // B2
                                            1, 0,  // E
                                            0,-1};
static const int d2d_irrep_degen[] = {1,1,1,1,2};
static const mat3x3 d2d_ops[] = {Identity(),
                                S<4>(vec3(0,0, 1)),
                                S<4>(vec3(0,0,-1)),
                                C<2>(vec3(0,0, 1)),
                                C<2>(vec3(1,0,0)),
                                C<2>(vec3(0,1,0)),
                                Reflection(vec3(1, 1,0)),
                                Reflection(vec3(1,-1,0))};
static const char *d2d_op_names[] = {"E", "S4z", "S4z^3", "C2z", "C2x", "C2y", "sx+y", "sx-y"};

const PointGroup PointGroup::D2d = PointGroup(8, 5, 2, d2d_name, d2d_irrep_names,
                                              d2d_generators, d2d_generator_reps, d2d_irrep_degen,
                                              d2d_ops, d2d_op_names);

/*
 * D3d
 */
static const char *d3d_name = "D3d";
static const char *d3d_irrep_names[] = {"A1g","A2g","Eg","A1u","A2u","Eu"};
static const int d3d_generators[] = {1,6};
static const double d3d_generator_reps[] = {// S6z
                                             1,                      // A1g
                                             1,                      // A2g
                                                   -0.5,-sqrt(0.75), // Eg
                                             sqrt(0.75),       -0.5,
                                            -1,                      // A1u
                                            -1,                      // A2u
                                                    0.5,-sqrt(0.75), // Eu
                                             sqrt(0.75),        0.5,
                                            // C2x
                                             1,     // A1g
                                            -1,     // A2g
                                             1, 0,  // Eg
                                             0,-1,
                                             1,     // A1u
                                            -1,     // A2u
                                             1, 0,  // Eu
                                             0,-1};
static const int d3d_irrep_degen[] = {1,1,2,1,1,2};
static const mat3x3 d3d_ops[] = {Identity(),
                                 S<6>(vec3(0,0, 1)),
                                 S<6>(vec3(0,0,-1)),
                                 C<3>(vec3(0,0, 1)),
                                 C<3>(vec3(0,0,-1)),
                                 Inversion(),
                                 C<2>(vec3(1,0,0)),
                                 C<2>(vec3(1,0,0)*C<6>(vec3(0,0, 1))),
                                 C<2>(vec3(1,0,0)*C<6>(vec3(0,0,-1))),
                                 Reflection(vec3(1,0,0)),
                                 Reflection(vec3(1,0,0)*C<6>(vec3(0,0, 1))),
                                 Reflection(vec3(1,0,0)*C<6>(vec3(0,0,-1)))};
static const char *d3d_op_names[] = {"E", "S6z", "S6z^5", "C3z", "C3z^2", "i", "C2x", "C22", "C23", "syz", "s2", "s3"};

const PointGroup PointGroup::D3d = PointGroup(12, 6, 2, d3d_name, d3d_irrep_names,
                                              d3d_generators, d3d_generator_reps, d3d_irrep_degen,
                                              d3d_ops, d3d_op_names);

/*
 * S4
 */
static const char *s4_name = "S4";
static const char *s4_irrep_names[] = {"A","B","E+","E-"};
static const int s4_generators[] = {0,1,2,3};
static const double s4_generator_reps[] = {1, 1, sqrt(2),       0,  // E
                                           1,-1,       0, sqrt(2),  // S4
                                           1, 1,-sqrt(2),       0,  // C2
                                           1,-1,       0,-sqrt(2)}; // S4^3
static const int s4_irrep_degen[] = {1,1,1,1};
static const mat3x3 s4_ops[] = {Identity(),
                                S<4>(vec3( 0, 0, 1)),
                                C<2>(vec3( 0, 0, 1)),
                                S<4>(vec3( 0, 0,-1))};
static const char *s4_op_names[] = {"E", "S4z", "C2z", "S4z^3"};

const PointGroup PointGroup::S4 = PointGroup(4, 4, 4, s4_name, s4_irrep_names,
                                             s4_generators, s4_generator_reps, s4_irrep_degen,
                                             s4_ops, s4_op_names);

/*
 * S6
 */
static const char *s6_name = "S6";
static const char *s6_irrep_names[] = {"Ag","Eg+","Eg-","Au","Eu+","Eu-"};
static const int s6_generators[] = {0,1,2,3,4,5};
static const double s6_generator_reps[] = {1,   sqrt(2),         0, 1,   sqrt(2),         0,  // E
                                           1,-1/sqrt(2), sqrt(1.5),-1, 1/sqrt(2), sqrt(1.5),  // S6
                                           1,-1/sqrt(2),-sqrt(1.5), 1,-1/sqrt(2), sqrt(1.5),  // C3
                                           1,   sqrt(2),         0,-1,  -sqrt(2),         0,  // i
                                           1,-1/sqrt(2), sqrt(1.5), 1,-1/sqrt(2),-sqrt(1.5),  // C3^2
                                           1,-1/sqrt(2),-sqrt(1.5),-1, 1/sqrt(2),-sqrt(1.5)}; // S6^5
static const int s6_irrep_degen[] = {1,1,1,1,1,1};
static const mat3x3 s6_ops[] = {Identity(),
                                S<6>(vec3( 0, 0, 1)),
                                C<3>(vec3( 0, 0, 1)),
                                Inversion(),
                                C<3>(vec3( 0, 0,-1)),
                                S<6>(vec3( 0, 0,-1))};
static const char *s6_op_names[] = {"E", "S6z", "C3z", "i", "C3z^2", "S6z^5"};

const PointGroup PointGroup::S6 = PointGroup(6, 6, 6, s6_name, s6_irrep_names,
                                             s6_generators, s6_generator_reps, s6_irrep_degen,
                                             s6_ops, s6_op_names);

}
}
