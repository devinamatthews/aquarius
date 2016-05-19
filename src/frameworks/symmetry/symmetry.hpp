#ifndef _AQUARIUS_SYMMETRY_HPP_
#define _AQUARIUS_SYMMETRY_HPP_

#include "util/global.hpp"

namespace aquarius
{
namespace symmetry
{

class PointGroup;

class Representation : public vector<double>
{
    friend class PointGroup;

    protected:
        const PointGroup& group;

        Representation(const PointGroup& group, const vector<double>& rep);

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
        vector<vector<double>> characters;
        vector<vector<vector<double>>> reps;
        const int *irrep_degen;
        const mat3x3 *ops;
        const char **op_names;
        vector<int> op_inverse;
        vector<int> op_product;
        vector<Representation> irreps;

        PointGroup(int order, int nirrep, int ngenerators, const char *name,
                   const char **irrep_names, const int *generators, const double *generator_reps,
                   const int *irrep_degen, const mat3x3 *ops, const char **op_names);

    public:
        bool operator==(const PointGroup& other) const { return this == &other; }

        vector<int> DCR(const vector<int>& U, const vector<int>& V, int& lambda) const;

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

mat3x3 Rotation(vec3 axis, double degrees);

mat3x3 Rotation(vec3 a, vec3 b);

mat3x3 Reflection(vec3 axis);

mat3x3 Identity();

mat3x3 Inversion();

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
