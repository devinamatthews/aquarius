#ifndef _AQUARIUS_UTIL_MATH_EXT_H_
#define _AQUARIUS_UTIL_MATH_EXT_H_

#include "stl_ext.hpp"
#include "lawrap.hpp"
#include "marray.hpp"

namespace aquarius
{

template <typename T>
typename enable_if<is_integral<T>::value,T>::type
roundup(T x, T y)
{
    return ((x+y-1)/y)*y;
}

template <typename T>
typename enable_if<is_integral<T>::value,T>::type
binom(T a, T b)
{
    T i, j;

    if (b < 0 || b > a) return 0;

    j = 1;
    for (i = 1;i <= min(b,a-b);i++)
    {
        j = (j*(a-i+1))/i;
    }

    return j;
}

template <typename T>
typename enable_if<is_integral<T>::value,T>::type
binomial(T a, T b)
{
    return binom(a, b);
}

template <typename T>
typename enable_if<is_integral<T>::value,T>::type
fact(T n)
{
    T i, j;

    j = 1;

    for (i = n;i > 1;i--)
    {
        j = j*i;
    }

    return j;
}

template <typename T>
typename enable_if<is_integral<T>::value,T>::type
factorial(T n)
{
    return fact(n);
}

template <typename T>
typename enable_if<is_integral<T>::value,T>::type
dfact(T n)
{
    T i, j;

    j = 1;

    for (i = n;i > 1;i -= 2)
    {
        j = j*i;
    }

    return j;
}

template <typename T>
void transpose(size_t data_, size_t n, const T& alpha, const T* restrict A, size_t lda,
                                   const T&  beta,       T* restrict B, size_t ldb)
{
    size_t i, j;

    if (alpha == 1.0)
    {
        if (beta == 0.0)
        {
            for (i = 0;i < data_;i++)
            {
                for (j = 0;j < n;j++)
                {
                    B[i*ldb + j] = A[j*lda + i];
                }
            }
        }
        else if (beta == 1.0)
        {
            for (i = 0;i < data_;i++)
            {
                for (j = 0;j < n;j++)
                {
                    B[i*ldb + j] += A[j*lda + i];
                }
            }
        }
        else
        {
            for (i = 0;i < data_;i++)
            {
                for (j = 0;j < n;j++)
                {
                    B[i*ldb + j] = beta*B[i*ldb + j] + A[j*lda + i];
                }
            }
        }
    }
    else
    {
        if (beta == 0.0)
        {
            for (i = 0;i < data_;i++)
            {
                for (j = 0;j < n;j++)
                {
                    B[i*ldb + j] = alpha*A[j*lda + i];
                }
            }
        }
        else if (beta == 1.0)
        {
            for (i = 0;i < data_;i++)
            {
                for (j = 0;j < n;j++)
                {
                    B[i*ldb + j] += alpha*A[j*lda + i];
                }
            }
        }
        else
        {
            for (i = 0;i < data_;i++)
            {
                for (j = 0;j < n;j++)
                {
                    B[i*ldb + j] = beta*B[i*ldb + j] + alpha*A[j*lda + i];
                }
            }
        }
    }
}

template <typename T>
void zero(size_t n, T* a, size_t inca)
{
    if (inca == 1)
    {
        fill(a, a+n, T());
    }
    else
    {
        size_t ia = 0;
        for (size_t i = 0;i < n;i++)
        {
            a[ia] = 0.0;
            ia += inca;
        }
    }
}

template <typename T>
void xypz(size_t n, const T& alpha, const T* restrict a, size_t inca,
                                    const T* restrict b, size_t incb,
                    const T&  beta,       T* restrict c, size_t incc)
{
    if (inca == 1 && incb == 1 && incc == 1)
    {
        if (alpha == 1.0)
        {
            if (beta == 0.0)
            {
                for (size_t i = 0;i < n;i++)
                {
                    c[i] = a[i]*b[i];
                }
            }
            else
            {
                for (size_t i = 0;i < n;i++)
                {
                    c[i] = beta*c[i] + a[i]*b[i];
                }
            }
        }
        else
        {
            if (beta == 0.0)
            {
                for (size_t i = 0;i < n;i++)
                {
                    c[i] = alpha*a[i]*b[i];
                }
            }
            else
            {
                for (size_t i = 0;i < n;i++)
                {
                    c[i] = beta*c[i] + alpha*a[i]*b[i];
                }
            }
        }
    }
    else
    {
        size_t ia = 0;
        size_t ib = 0;
        size_t ic = 0;
        for (size_t i = 0;i < n;i++)
        {
            c[ic] = beta*c[ic] + alpha*a[ia]*b[ib];
            ia += inca;
            ib += incb;
            ic += incc;
        }
    }
}

class vec3;
class mat3x3;

class vec3 : public marray<double, 1>
{
    friend class mat3x3;

    public:
        vec3() : marray<double, 1>(3) {}

        vec3(vec3&& other) : marray<double, 1>(other) {}

        vec3(const vec3& other) : marray<double, 1>(other, construct_copy) {}

        vec3(double pos[3]) : marray<double, 1>(3)
        {
            data_[0] = pos[0];
            data_[1] = pos[1];
            data_[2] = pos[2];
        }

        vec3(double x, double y, double z) : marray<double, 1>(3)
        {
            data_[0] = x;
            data_[1] = y;
            data_[2] = z;
        }

        bool operator==(const vec3& other) const
        {
            return data_[0] == other[0] &&
                   data_[1] == other[1] &&
                   data_[2] == other[2];
        }

        vec3& operator=(double other)
        {
            data_[0] = other;
            data_[1] = other;
            data_[2] = other;
            return *this;
        }

        vec3& operator=(const vec3& other)
        {
            data_[0] = other.data_[0];
            data_[1] = other.data_[1];
            data_[2] = other.data_[2];
            return *this;
        }

        vec3& operator+=(const vec3& other)
        {
            data_[0] += other[0];
            data_[1] += other[1];
            data_[2] += other[2];
            return *this;
        }

        vec3& operator-=(const vec3& other)
        {
            data_[0] -= other[0];
            data_[1] -= other[1];
            data_[2] -= other[2];
            return *this;
        }

        vec3 operator+(const vec3& other) const
        {
            return vec3(*this) += other;
        }

        vec3 operator-(const vec3& other) const
        {
            return vec3(*this) -= other;
        }

        vec3 operator^(const vec3& other) const
        {
            return vec3(data_[1]*other[2] - data_[2]*other[1],
                        data_[2]*other[0] - data_[0]*other[2],
                        data_[0]*other[1] - data_[1]*other[0]);
        }

        vec3 operator%(const vec3& other) const
        {
            return vec3(*this).orthogonalize(other);
        }

        vec3 operator-() const
        {
            return vec3(-data_[0], -data_[1], -data_[2]);
        }

        double operator*(const vec3& other) const
        {
            return data_[0]*other[0] + data_[1]*other[1] + data_[2]*other[2];
        }

        vec3 operator*(const mat3x3& other) const;

        vec3& operator*=(double a)
        {
            data_[0] *= a;
            data_[1] *= a;
            data_[2] *= a;
            return *this;
        }

        vec3& operator/=(double a)
        {
            data_[0] /= a;
            data_[1] /= a;
            data_[2] /= a;
            return *this;
        }

        vec3 operator*(double a) const
        {
            return vec3(*this) *= a;
        }

        vec3 operator/(double a) const
        {
            return vec3(*this) /= a;
        }

        friend vec3 operator*(double a, const vec3& data_)
        {
            return vec3(data_) *= a;
        }

        vec3 operator/(const vec3& other) const
        {
            return vec3(*this).orthogonalize(other);
        }

        mat3x3 operator|(const vec3& other) const;

        friend vec3 unit(const vec3& data_)
        {
            return vec3(data_).normalize();
        }

        friend double norm(const vec3& data_)
        {
            return data_.norm();
        }

        double norm() const
        {
            return sqrt(norm2());
        }

        friend double norm2(const vec3& data_)
        {
            return data_.norm2();
        }

        double norm2() const
        {
            return data_[0]*data_[0] + data_[1]*data_[1] + data_[2]*data_[2];
        }

        vec3& normalize()
        {
            *this /= norm();
            return *this;
        }

        vec3& orthogonalize(const vec3& other)
        {
            double n = other.norm();
            *this -= other*(*this*other)/n;
            return *this;
        }

        friend ostream& operator<<(ostream& os, const vec3& data_)
        {
            os << "(" << data_[0] << ", " << data_[1] << ", " << data_[2] << ")";
            return os;
        }
};

class mat3x3 : public marray<double, 2>
{
    friend class vec3;

    public:
        mat3x3() : marray<double, 2>(3,3) {}

        mat3x3(mat3x3&& other) : marray<double, 2>(other) {}

        mat3x3(const mat3x3& other) : marray<double, 2>(other, construct_copy) {}

        mat3x3(double m00, double m01, double m02,
               double m10, double m11, double m12,
               double m20, double m21, double m22) : marray<double, 2>(3,3)
        {
            (*this)[0][0] = m00; (*this)[0][1] = m01; (*this)[0][2] = m02;
            (*this)[1][0] = m10; (*this)[1][1] = m11; (*this)[1][2] = m12;
            (*this)[2][0] = m20; (*this)[2][1] = m21; (*this)[2][2] = m22;
        }

        const mat3x3 operator^(int p) const
        {
            assert(p >= 0);

            mat3x3 ret(1,0,0,
                       0,1,0,
                       0,0,1);

            for (;p > 0;p--) ret = ret*(*this);

            return ret;
        }

        bool operator==(const mat3x3& other) const
        {
            return (*this)[0][0] == other[0][0] && (*this)[0][1] == other[0][1] && (*this)[0][2] == other[0][2] &&
                   (*this)[1][0] == other[1][0] && (*this)[1][1] == other[1][1] && (*this)[1][2] == other[1][2] &&
                   (*this)[2][0] == other[2][0] && (*this)[2][1] == other[2][1] && (*this)[2][2] == other[2][2];
        }

        mat3x3& operator=(const mat3x3& other)
        {
            (*this)[0][0] = other[0][0]; (*this)[0][1] = other[0][1]; (*this)[0][2] = other[0][2];
            (*this)[1][0] = other[1][0]; (*this)[1][1] = other[1][1]; (*this)[1][2] = other[1][2];
            (*this)[2][0] = other[2][0]; (*this)[2][1] = other[2][1]; (*this)[2][2] = other[2][2];
            return *this;
        }

        mat3x3& operator=(double other)
        {
            (*this)[0][0] = other;
            (*this)[1][0] = other;
            (*this)[2][0] = other;
            (*this)[0][1] = other;
            (*this)[1][1] = other;
            (*this)[2][1] = other;
            (*this)[0][2] = other;
            (*this)[1][2] = other;
            (*this)[2][2] = other;
            return *this;
        }

        mat3x3 operator*(const mat3x3& other) const
        {
            mat3x3 r;
            r[0][0] = (*this)[0][0]*other[0][0] + (*this)[0][1]*other[1][0] + (*this)[0][2]*other[2][0];
            r[0][1] = (*this)[0][0]*other[0][1] + (*this)[0][1]*other[1][1] + (*this)[0][2]*other[2][1];
            r[0][2] = (*this)[0][0]*other[0][2] + (*this)[0][1]*other[1][2] + (*this)[0][2]*other[2][2];
            r[1][0] = (*this)[1][0]*other[0][0] + (*this)[1][1]*other[1][0] + (*this)[1][2]*other[2][0];
            r[1][1] = (*this)[1][0]*other[0][1] + (*this)[1][1]*other[1][1] + (*this)[1][2]*other[2][1];
            r[1][2] = (*this)[1][0]*other[0][2] + (*this)[1][1]*other[1][2] + (*this)[1][2]*other[2][2];
            r[2][0] = (*this)[2][0]*other[0][0] + (*this)[2][1]*other[1][0] + (*this)[2][2]*other[2][0];
            r[2][1] = (*this)[2][0]*other[0][1] + (*this)[2][1]*other[1][1] + (*this)[2][2]*other[2][1];
            r[2][2] = (*this)[2][0]*other[0][2] + (*this)[2][1]*other[1][2] + (*this)[2][2]*other[2][2];
            return r;
        }

        mat3x3 operator+(const mat3x3& other) const
        {
            mat3x3 r;
            r[0][0] = (*this)[0][0] + other[0][0];
            r[0][1] = (*this)[0][1] + other[0][1];
            r[0][2] = (*this)[0][2] + other[0][2];
            r[1][0] = (*this)[1][0] + other[1][0];
            r[1][1] = (*this)[1][1] + other[1][1];
            r[1][2] = (*this)[1][2] + other[1][2];
            r[2][0] = (*this)[2][0] + other[2][0];
            r[2][1] = (*this)[2][1] + other[2][1];
            r[2][2] = (*this)[2][2] + other[2][2];
            return r;
        }

        mat3x3 operator-(const mat3x3& other) const
        {
            mat3x3 r;
            r[0][0] = (*this)[0][0] - other[0][0];
            r[0][1] = (*this)[0][1] - other[0][1];
            r[0][2] = (*this)[0][2] - other[0][2];
            r[1][0] = (*this)[1][0] - other[1][0];
            r[1][1] = (*this)[1][1] - other[1][1];
            r[1][2] = (*this)[1][2] - other[1][2];
            r[2][0] = (*this)[2][0] - other[2][0];
            r[2][1] = (*this)[2][1] - other[2][1];
            r[2][2] = (*this)[2][2] - other[2][2];
            return r;
        }

        mat3x3& operator+=(const mat3x3& other)
        {
            (*this)[0][0] += other[0][0];
            (*this)[0][1] += other[0][1];
            (*this)[0][2] += other[0][2];
            (*this)[1][0] += other[1][0];
            (*this)[1][1] += other[1][1];
            (*this)[1][2] += other[1][2];
            (*this)[2][0] += other[2][0];
            (*this)[2][1] += other[2][1];
            (*this)[2][2] += other[2][2];
            return *this;
        }

        mat3x3& operator-=(const mat3x3& other)
        {
            (*this)[0][0] -= other[0][0];
            (*this)[0][1] -= other[0][1];
            (*this)[0][2] -= other[0][2];
            (*this)[1][0] -= other[1][0];
            (*this)[1][1] -= other[1][1];
            (*this)[1][2] -= other[1][2];
            (*this)[2][0] -= other[2][0];
            (*this)[2][1] -= other[2][1];
            (*this)[2][2] -= other[2][2];
            return *this;
        }

        friend mat3x3 operator+(double other, const mat3x3& m)
        {
            return m+other;
        }

        friend mat3x3 operator-(double other, const mat3x3& m)
        {
            return -(m-other);
        }

        friend mat3x3 operator*(double other, const mat3x3& m)
        {
            return m*other;
        }

        mat3x3 operator+(double other) const
        {
            mat3x3 r;
            r[0][0] = (*this)[0][0] + other;
            r[0][1] = (*this)[0][1];
            r[0][2] = (*this)[0][2];
            r[1][0] = (*this)[1][0];
            r[1][1] = (*this)[1][1] + other;
            r[1][2] = (*this)[1][2];
            r[2][0] = (*this)[2][0];
            r[2][1] = (*this)[2][1];
            r[2][2] = (*this)[2][2] + other;
            return r;
        }

        mat3x3 operator-(double other) const
        {
            mat3x3 r;
            r[0][0] = (*this)[0][0] - other;
            r[0][1] = (*this)[0][1];
            r[0][2] = (*this)[0][2];
            r[1][0] = (*this)[1][0];
            r[1][1] = (*this)[1][1] - other;
            r[1][2] = (*this)[1][2];
            r[2][0] = (*this)[2][0];
            r[2][1] = (*this)[2][1];
            r[2][2] = (*this)[2][2] - other;
            return r;
        }

        mat3x3 operator*(double other) const
        {
            mat3x3 r;
            r[0][0] = (*this)[0][0]*other;
            r[0][1] = (*this)[0][1]*other;
            r[0][2] = (*this)[0][2]*other;
            r[1][0] = (*this)[1][0]*other;
            r[1][1] = (*this)[1][1]*other;
            r[1][2] = (*this)[1][2]*other;
            r[2][0] = (*this)[2][0]*other;
            r[2][1] = (*this)[2][1]*other;
            r[2][2] = (*this)[2][2]*other;
            return r;
        }

        mat3x3 operator/(double other) const
        {
            mat3x3 r;
            r[0][0] = (*this)[0][0]/other;
            r[0][1] = (*this)[0][1]/other;
            r[0][2] = (*this)[0][2]/other;
            r[1][0] = (*this)[1][0]/other;
            r[1][1] = (*this)[1][1]/other;
            r[1][2] = (*this)[1][2]/other;
            r[2][0] = (*this)[2][0]/other;
            r[2][1] = (*this)[2][1]/other;
            r[2][2] = (*this)[2][2]/other;
            return r;
        }

        mat3x3 operator-() const
        {
            mat3x3 r;
            r[0][0] = -(*this)[0][0];
            r[0][1] = -(*this)[0][1];
            r[0][2] = -(*this)[0][2];
            r[1][0] = -(*this)[1][0];
            r[1][1] = -(*this)[1][1];
            r[1][2] = -(*this)[1][2];
            r[2][0] = -(*this)[2][0];
            r[2][1] = -(*this)[2][1];
            r[2][2] = -(*this)[2][2];
            return r;
        }

        vec3 operator*(const vec3& other) const
        {
            vec3 r;
            r[0] = (*this)[0][0]*other[0] + (*this)[0][1]*other[1] + (*this)[0][2]*other[2];
            r[1] = (*this)[1][0]*other[0] + (*this)[1][1]*other[1] + (*this)[1][2]*other[2];
            r[2] = (*this)[2][0]*other[0] + (*this)[2][1]*other[1] + (*this)[2][2]*other[2];
            return r;
        }

        void diagonalize(vec3& eigenvalues, mat3x3& eigenvectors) const
        {
            eigenvectors = *this;
            heev('V', 'U', 3, eigenvectors.data(), 3, eigenvalues.data());
        }

        mat3x3 identity()
        {
            mat3x3 i;
            i[0][0] = 1;
            i[1][1] = 1;
            i[2][2] = 1;
            return i;
        }

        double norm() const
        {
            return sqrt((*this)[0][0]*(*this)[0][0] + (*this)[0][1]*(*this)[0][1] + (*this)[0][2]*(*this)[0][2] +
                        (*this)[1][0]*(*this)[1][0] + (*this)[1][1]*(*this)[1][1] + (*this)[1][2]*(*this)[1][2] +
                        (*this)[2][0]*(*this)[2][0] + (*this)[2][1]*(*this)[2][1] + (*this)[2][2]*(*this)[2][2]);
        }

        friend double norm(const mat3x3& m)
        {
            return m.norm();
        }

        friend ostream& operator<<(ostream& os, const mat3x3& m)
        {
            auto oldfmt = os.flags();
            os << scientific << showpos << setprecision(12) << setw(20);
            os <<  "/" << m[0][0] << " " << m[0][1] << " " << m[0][2] << "\\\n";
            os <<  "|" << m[1][0] << " " << m[1][1] << " " << m[1][2] << "|\n";
            os << "\\" << m[2][0] << " " << m[2][1] << " " << m[2][2] << "/\n";
            os.flags(oldfmt);
            return os;
        }
};

inline vec3 vec3::operator*(const mat3x3& other) const
{
    vec3 r;
    r[0] = data_[0]*other[0][0] + data_[1]*other[1][0] + data_[2]*other[2][0];
    r[1] = data_[0]*other[0][1] + data_[1]*other[1][1] + data_[2]*other[2][1];
    r[2] = data_[0]*other[0][2] + data_[1]*other[1][2] + data_[2]*other[2][2];
    return r;
}

inline mat3x3 vec3::operator|(const vec3& other) const
{
    return mat3x3(data_[0]*other[0],data_[0]*other[1],data_[0]*other[2],
                  data_[1]*other[0],data_[1]*other[1],data_[1]*other[2],
                  data_[2]*other[0],data_[2]*other[1],data_[2]*other[2]);
}

}

#endif
