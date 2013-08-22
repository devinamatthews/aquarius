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

#ifndef _AQUARIUS_UTIL_MATH_EXT_H_
#define _AQUARIUS_UTIL_MATH_EXT_H_

#include <stdint.h>

#ifdef __cplusplus
#include <cmath>
#else
#include <math.h>
#endif

#include "util.h"

#ifdef __cplusplus
extern "C"
{
#endif

int roundup(int x, int y);

double dist2(const double* a, const double* b);

double dist(const double* a, const double* b);

int binom(int a, int b);

int64_t fact(int n);

int64_t dfact(int n);

#ifdef __cplusplus
}

#include <ostream>

namespace aquarius
{

template <typename T>
void transpose(const size_t m, const size_t n, const T alpha, const T* A, const size_t lda,
                                               const T  beta,       T* B, const size_t ldb);

class vec3;
class mat3x3;

class vec3
{
    private:
        double v[3];

    public:
        vec3();

        vec3(double pos[3]);

        vec3(double x, double y, double z);

        vec3(const vec3& other);

        vec3& operator+=(const vec3& other);

        vec3& operator-=(const vec3& other);

        vec3 operator+(const vec3& other) const;

        vec3 operator-(const vec3& other) const;

        vec3 operator^(const vec3& other) const;

        double& operator[](int i);

        const double& operator[](int i) const;

        bool operator==(const vec3& other) const;

        operator double*();

        operator const double*() const;

        double operator*(const vec3& other) const;

        vec3 operator*(const mat3x3& other) const;

        vec3& operator*=(const double a);

        vec3& operator/=(const double a);

        vec3 operator*(const double a) const;

        vec3 operator/(const double a) const;

        vec3 operator/(const vec3& other) const;

        double norm() const;

        vec3& normalize();

        vec3& orthogonalize(const vec3& other);
};

std::ostream& operator<<(std::ostream& os, const vec3& v);

vec3 operator*(const double a, const vec3& v);

vec3 operator/(const double a, const vec3& v);

vec3 unit(const vec3& v);

double norm(const vec3& v);

class mat3x3
{
    protected:
        double m[3][3];

        class column
        {
            friend class mat3x3;

            protected:
                double &r1, &r2, &r3;
                column(double& r1, double& r2, double& r3) : r1(r1), r2(r2), r3(r3) {};

            public:
                double& operator[](int i);
        };

        class const_column
        {
            friend class mat3x3;

            protected:
                const double &r1, &r2, &r3;
                const_column(const double& r1, const double& r2, const double& r3) : r1(r1), r2(r2), r3(r3) {};

            public:
                const double& operator[](int i) const;
        };

    public:
        mat3x3();

        mat3x3(double m00, double m01, double m02,
               double m10, double m11, double m12,
               double m20, double m21, double m22);

        const mat3x3 operator^(int p) const;

        column operator[](int i);

        const_column operator[](int i) const;

        bool operator==(const mat3x3& other) const;

        mat3x3 operator*(const mat3x3& other) const;

        mat3x3 operator+(const mat3x3& other) const;

        mat3x3 operator-(const mat3x3& other) const;

        vec3 operator*(const vec3& other) const;

        operator double*();

        operator const double*() const;

        void diagonalize(vec3& eigenvalues, mat3x3& eigenvectors) const;

        static mat3x3 identity();

        double norm() const;
};

double norm(const mat3x3& m);

std::ostream& operator<<(std::ostream& os, const mat3x3& v);

}

#endif

#endif
