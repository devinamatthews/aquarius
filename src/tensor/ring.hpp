#ifndef _AQUARIUS_TENSOR_RING_HPP_
#define _AQUARIUS_TENSOR_RING_HPP_

#include "util/global.hpp"

namespace aquarius
{
namespace tensor
{

typedef complex<      float>  scomplex;
typedef complex<     double>  dcomplex;
typedef complex<long double> ldcomplex;

struct Ring
{
    const int type;
    const size_t size;

    Ring(int type, size_t size) : type(type), size(size) {}

    bool operator==(const Ring& other) const { return type == other.type; }
};

namespace detail
{
    static size_t field_sizes[] = {sizeof(      float),
                                   sizeof(     double),
                                   sizeof(long double),
                                   sizeof(   scomplex),
                                   sizeof(   dcomplex),
                                   sizeof(  ldcomplex)};
}

struct Field : Ring
{
    enum field : int {SINGLE, DOUBLE, LDOUBLE, SCOMPLEX, DCOMPLEX, LDCOMPLEX};

    Field(field F) : Ring(F, detail::field_sizes[F]) {}

    Field(const       float& val) : Ring(   SINGLE, sizeof(      float)) {}
    Field(const      double& val) : Ring(   DOUBLE, sizeof(     double)) {}
    Field(const long double& val) : Ring(  LDOUBLE, sizeof(long double)) {}
    Field(const    scomplex& val) : Ring( SCOMPLEX, sizeof(   scomplex)) {}
    Field(const    dcomplex& val) : Ring( DCOMPLEX, sizeof(   dcomplex)) {}
    Field(const   ldcomplex& val) : Ring(LDCOMPLEX, sizeof(  ldcomplex)) {}

    static Datatype MPI_TYPE(const Field& F)
    {
        return MPI_TYPE((field)F.type);
    }

    static Datatype MPI_TYPE(field F)
    {
        switch (F)
        {
            case    SINGLE: return MPI_TYPE_<      float>::value();
            case    DOUBLE: return MPI_TYPE_<     double>::value();
            case   LDOUBLE: return MPI_TYPE_<long double>::value();
            case  SCOMPLEX: return MPI_TYPE_<   scomplex>::value();
            case  DCOMPLEX: return MPI_TYPE_<   dcomplex>::value();
            case LDCOMPLEX: return MPI_TYPE_<  ldcomplex>::value();
        }
        return MPI_TYPE_<double>::value();
    }
};

template <typename T>
struct is_field : is_floating_point<typename real_type<T>::type> {};

class Scalar
{
    protected:
        Field F;
        union
        {
            float fval;
            double dval;
            long double ldval;
            scomplex fcval;
            dcomplex dcval;
            ldcomplex ldcval;
        };

        static Field::field resultType(int f1, int f2)
        {
            /*
             * If the other argument is from a real or complex field, up-convert to
             * complex if either argument is and set the floating-point type to the
             * widest of the operands. Otherwise, keep the field's type.
             */
            return Field::field((f1 >= Field::SCOMPLEX || f2 >= Field::SCOMPLEX ? Field::SCOMPLEX : 0)+
                                max(f1 >= Field::SCOMPLEX ? f1-Field::SCOMPLEX : f1,
                                         f2 >= Field::SCOMPLEX ? f2-Field::SCOMPLEX : f2));
        }

        template <typename T>
        enable_if_t<is_field<T>::value,Field::field> resultType(T other) const
        {
            return resultType(F.type, Field(other).type);
        }

        template <typename T>
        enable_if_t<!is_field<T>::value,Field::field> resultType(T other) const
        {
            return (Field::field)F.type;
        }

    public:
        Scalar(Field::field type = Field::DOUBLE) : F(type)
        {
            *this = 0.0;
        }

        template <typename T>
        Scalar(Field::field type, T val) : F(type)
        {
            *this = val;
        }

        template <typename T>
        Scalar(T val, enable_if_t<is_arithmetic<real_type_t<T>>::value>* foo = 0) : F(is_field<T>::value ? val : double())
        {
            *this = val;
        }

        friend Scalar abs(const Scalar& s)
        {
            switch (s.F.type)
            {
                case Field::SINGLE:    return Scalar(aquarius::abs(  s.fval)); break;
                case Field::DOUBLE:    return Scalar(aquarius::abs(  s.dval)); break;
                case Field::LDOUBLE:   return Scalar(aquarius::abs( s.ldval)); break;
                case Field::SCOMPLEX:  return Scalar(aquarius::abs( s.fcval)); break;
                case Field::DCOMPLEX:  return Scalar(aquarius::abs( s.dcval)); break;
                case Field::LDCOMPLEX: return Scalar(aquarius::abs(s.ldcval)); break;
            }
            return Scalar(0.0);
        }

        friend Scalar conj(Scalar s)
        {
            switch (s.F.type)
            {
                case Field::SCOMPLEX:   s.fcval = conj( s.fcval); break;
                case Field::DCOMPLEX:   s.dcval = conj( s.dcval); break;
                case Field::LDCOMPLEX: s.ldcval = conj(s.ldcval); break;
            }
            return s;
        }

        friend Scalar sqrt(Scalar s)
        {
            switch (s.F.type)
            {
                case Field::SINGLE:      s.fval = sqrt(  s.fval); break;
                case Field::DOUBLE:      s.dval = sqrt(  s.dval); break;
                case Field::LDOUBLE:    s.ldval = sqrt( s.ldval); break;
                case Field::SCOMPLEX:   s.fcval = sqrt( s.fcval); break;
                case Field::DCOMPLEX:   s.dcval = sqrt( s.dcval); break;
                case Field::LDCOMPLEX: s.ldcval = sqrt(s.ldcval); break;
            }
            return s;
        }

        friend Scalar min(const Scalar& s1, const Scalar& s2)
        {
            switch (s1.F.type)
            {
                case Field::SINGLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(min((      float)  s1.fval , (      float)  s2.fval) ); break;
                        case Field::DOUBLE:    return Scalar(min((     double)  s1.fval , (     double)  s2.dval) ); break;
                        case Field::LDOUBLE:   return Scalar(min((long double)  s1.fval , (long double) s2.ldval) ); break;
                        case Field::SCOMPLEX:  return Scalar(min(    scomplex(  s1.fval),     scomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(min(    dcomplex(  s1.fval),     dcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(min(   ldcomplex(  s1.fval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::DOUBLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(min((     double)  s1.dval , (     double)  s2.fval) ); break;
                        case Field::DOUBLE:    return Scalar(min((     double)  s1.dval , (     double)  s2.dval) ); break;
                        case Field::LDOUBLE:   return Scalar(min((long double)  s1.dval , (long double) s2.ldval) ); break;
                        case Field::SCOMPLEX:  return Scalar(min(    dcomplex(  s1.dval),     dcomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(min(    dcomplex(  s1.dval),     dcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(min(   ldcomplex(  s1.dval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::LDOUBLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(min((long double) s1.ldval , (long double)  s2.fval) ); break;
                        case Field::DOUBLE:    return Scalar(min((long double) s1.ldval , (long double)  s2.dval) ); break;
                        case Field::LDOUBLE:   return Scalar(min((long double) s1.ldval , (long double) s2.ldval) ); break;
                        case Field::SCOMPLEX:  return Scalar(min(   ldcomplex( s1.ldval),    ldcomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(min(   ldcomplex( s1.ldval),    ldcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(min(   ldcomplex( s1.ldval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::SCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(min(    scomplex( s1.fcval),     scomplex(  s2.fval))); break;
                        case Field::DOUBLE:    return Scalar(min(    dcomplex( s1.fcval),     dcomplex(  s2.dval))); break;
                        case Field::LDOUBLE:   return Scalar(min(   ldcomplex( s1.fcval),    ldcomplex( s2.ldval))); break;
                        case Field::SCOMPLEX:  return Scalar(min(    scomplex( s1.fcval),     scomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(min(    dcomplex( s1.fcval),     dcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(min(   ldcomplex( s1.fcval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::DCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(min(    dcomplex( s1.dcval),     dcomplex(  s2.fval))); break;
                        case Field::DOUBLE:    return Scalar(min(    dcomplex( s1.dcval),     dcomplex(  s2.dval))); break;
                        case Field::LDOUBLE:   return Scalar(min(   ldcomplex( s1.dcval),    ldcomplex( s2.ldval))); break;
                        case Field::SCOMPLEX:  return Scalar(min(    dcomplex( s1.dcval),     dcomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(min(    dcomplex( s1.dcval),     dcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(min(   ldcomplex( s1.dcval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::LDCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(min(   ldcomplex(s1.ldcval),    ldcomplex(  s2.fval))); break;
                        case Field::DOUBLE:    return Scalar(min(   ldcomplex(s1.ldcval),    ldcomplex(  s2.dval))); break;
                        case Field::LDOUBLE:   return Scalar(min(   ldcomplex(s1.ldcval),    ldcomplex( s2.ldval))); break;
                        case Field::SCOMPLEX:  return Scalar(min(   ldcomplex(s1.ldcval),    ldcomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(min(   ldcomplex(s1.ldcval),    ldcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(min(   ldcomplex(s1.ldcval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
            }

            return Scalar();
        }

        friend Scalar max(const Scalar& s1, const Scalar& s2)
        {
            switch (s1.F.type)
            {
                case Field::SINGLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(max((      float)  s1.fval , (      float)  s2.fval) ); break;
                        case Field::DOUBLE:    return Scalar(max((     double)  s1.fval , (     double)  s2.dval) ); break;
                        case Field::LDOUBLE:   return Scalar(max((long double)  s1.fval , (long double) s2.ldval) ); break;
                        case Field::SCOMPLEX:  return Scalar(max(    scomplex(  s1.fval),     scomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(max(    dcomplex(  s1.fval),     dcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(max(   ldcomplex(  s1.fval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::DOUBLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(max((     double)  s1.dval , (     double)  s2.fval) ); break;
                        case Field::DOUBLE:    return Scalar(max((     double)  s1.dval , (     double)  s2.dval) ); break;
                        case Field::LDOUBLE:   return Scalar(max((long double)  s1.dval , (long double) s2.ldval) ); break;
                        case Field::SCOMPLEX:  return Scalar(max(    dcomplex(  s1.dval),     dcomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(max(    dcomplex(  s1.dval),     dcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(max(   ldcomplex(  s1.dval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::LDOUBLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(max((long double) s1.ldval , (long double)  s2.fval) ); break;
                        case Field::DOUBLE:    return Scalar(max((long double) s1.ldval , (long double)  s2.dval) ); break;
                        case Field::LDOUBLE:   return Scalar(max((long double) s1.ldval , (long double) s2.ldval) ); break;
                        case Field::SCOMPLEX:  return Scalar(max(   ldcomplex( s1.ldval),    ldcomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(max(   ldcomplex( s1.ldval),    ldcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(max(   ldcomplex( s1.ldval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::SCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(max(    scomplex( s1.fcval),     scomplex(  s2.fval))); break;
                        case Field::DOUBLE:    return Scalar(max(    dcomplex( s1.fcval),     dcomplex(  s2.dval))); break;
                        case Field::LDOUBLE:   return Scalar(max(   ldcomplex( s1.fcval),    ldcomplex( s2.ldval))); break;
                        case Field::SCOMPLEX:  return Scalar(max(    scomplex( s1.fcval),     scomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(max(    dcomplex( s1.fcval),     dcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(max(   ldcomplex( s1.fcval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::DCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(max(    dcomplex( s1.dcval),     dcomplex(  s2.fval))); break;
                        case Field::DOUBLE:    return Scalar(max(    dcomplex( s1.dcval),     dcomplex(  s2.dval))); break;
                        case Field::LDOUBLE:   return Scalar(max(   ldcomplex( s1.dcval),    ldcomplex( s2.ldval))); break;
                        case Field::SCOMPLEX:  return Scalar(max(    dcomplex( s1.dcval),     dcomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(max(    dcomplex( s1.dcval),     dcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(max(   ldcomplex( s1.dcval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
                case Field::LDCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return Scalar(max(   ldcomplex(s1.ldcval),    ldcomplex(  s2.fval))); break;
                        case Field::DOUBLE:    return Scalar(max(   ldcomplex(s1.ldcval),    ldcomplex(  s2.dval))); break;
                        case Field::LDOUBLE:   return Scalar(max(   ldcomplex(s1.ldcval),    ldcomplex( s2.ldval))); break;
                        case Field::SCOMPLEX:  return Scalar(max(   ldcomplex(s1.ldcval),    ldcomplex( s2.fcval))); break;
                        case Field::DCOMPLEX:  return Scalar(max(   ldcomplex(s1.ldcval),    ldcomplex( s2.dcval))); break;
                        case Field::LDCOMPLEX: return Scalar(max(   ldcomplex(s1.ldcval),    ldcomplex(s2.ldcval))); break;
                    }
                    break;
            }

            return Scalar();
        }

        friend bool operator==(const Scalar& s1, const Scalar& s2)
        {
            switch (s1.F.type)
            {
                case Field::SINGLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return (      float)  s1.fval  == (      float)  s2.fval ; break;
                        case Field::DOUBLE:    return (     double)  s1.fval  == (     double)  s2.dval ; break;
                        case Field::LDOUBLE:   return (long double)  s1.fval  == (long double) s2.ldval ; break;
                        case Field::SCOMPLEX:  return     scomplex(  s1.fval) ==     scomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return     dcomplex(  s1.fval) ==     dcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex(  s1.fval) ==    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::DOUBLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return (     double)  s1.dval  == (     double)  s2.fval ; break;
                        case Field::DOUBLE:    return (     double)  s1.dval  == (     double)  s2.dval ; break;
                        case Field::LDOUBLE:   return (long double)  s1.dval  == (long double) s2.ldval ; break;
                        case Field::SCOMPLEX:  return     dcomplex(  s1.dval) ==     dcomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return     dcomplex(  s1.dval) ==     dcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex(  s1.dval) ==    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::LDOUBLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return (long double) s1.ldval  == (long double)  s2.fval ; break;
                        case Field::DOUBLE:    return (long double) s1.ldval  == (long double)  s2.dval ; break;
                        case Field::LDOUBLE:   return (long double) s1.ldval  == (long double) s2.ldval ; break;
                        case Field::SCOMPLEX:  return    ldcomplex( s1.ldval) ==    ldcomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return    ldcomplex( s1.ldval) ==    ldcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex( s1.ldval) ==    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::SCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return     scomplex( s1.fcval) ==     scomplex(  s2.fval); break;
                        case Field::DOUBLE:    return     dcomplex( s1.fcval) ==     dcomplex(  s2.dval); break;
                        case Field::LDOUBLE:   return    ldcomplex( s1.fcval) ==    ldcomplex( s2.ldval); break;
                        case Field::SCOMPLEX:  return     scomplex( s1.fcval) ==     scomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return     dcomplex( s1.fcval) ==     dcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex( s1.fcval) ==    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::DCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return     dcomplex( s1.dcval) ==     dcomplex(  s2.fval); break;
                        case Field::DOUBLE:    return     dcomplex( s1.dcval) ==     dcomplex(  s2.dval); break;
                        case Field::LDOUBLE:   return    ldcomplex( s1.dcval) ==    ldcomplex( s2.ldval); break;
                        case Field::SCOMPLEX:  return     dcomplex( s1.dcval) ==     dcomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return     dcomplex( s1.dcval) ==     dcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex( s1.dcval) ==    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::LDCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return    ldcomplex(s1.ldcval) ==    ldcomplex(  s2.fval); break;
                        case Field::DOUBLE:    return    ldcomplex(s1.ldcval) ==    ldcomplex(  s2.dval); break;
                        case Field::LDOUBLE:   return    ldcomplex(s1.ldcval) ==    ldcomplex( s2.ldval); break;
                        case Field::SCOMPLEX:  return    ldcomplex(s1.ldcval) ==    ldcomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return    ldcomplex(s1.ldcval) ==    ldcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex(s1.ldcval) ==    ldcomplex(s2.ldcval); break;
                    }
                    break;
            }

            return false;
        }

        friend bool operator!=(const Scalar& s1, const Scalar& s2)
        {
            return !(s1 == s2);
        }

        friend bool operator<(const Scalar& s1, const Scalar& s2)
        {
            switch (s1.F.type)
            {
                case Field::SINGLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return (      float)  s1.fval  < (      float)  s2.fval ; break;
                        case Field::DOUBLE:    return (     double)  s1.fval  < (     double)  s2.dval ; break;
                        case Field::LDOUBLE:   return (long double)  s1.fval  < (long double) s2.ldval ; break;
                        case Field::SCOMPLEX:  return     scomplex(  s1.fval) <     scomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return     dcomplex(  s1.fval) <     dcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex(  s1.fval) <    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::DOUBLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return (     double)  s1.dval  < (     double)  s2.fval ; break;
                        case Field::DOUBLE:    return (     double)  s1.dval  < (     double)  s2.dval ; break;
                        case Field::LDOUBLE:   return (long double)  s1.dval  < (long double) s2.ldval ; break;
                        case Field::SCOMPLEX:  return     dcomplex(  s1.dval) <     dcomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return     dcomplex(  s1.dval) <     dcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex(  s1.dval) <    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::LDOUBLE:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return (long double) s1.ldval  < (long double)  s2.fval ; break;
                        case Field::DOUBLE:    return (long double) s1.ldval  < (long double)  s2.dval ; break;
                        case Field::LDOUBLE:   return (long double) s1.ldval  < (long double) s2.ldval ; break;
                        case Field::SCOMPLEX:  return    ldcomplex( s1.ldval) <    ldcomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return    ldcomplex( s1.ldval) <    ldcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex( s1.ldval) <    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::SCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return     scomplex( s1.fcval) <     scomplex(  s2.fval); break;
                        case Field::DOUBLE:    return     dcomplex( s1.fcval) <     dcomplex(  s2.dval); break;
                        case Field::LDOUBLE:   return    ldcomplex( s1.fcval) <    ldcomplex( s2.ldval); break;
                        case Field::SCOMPLEX:  return     scomplex( s1.fcval) <     scomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return     dcomplex( s1.fcval) <     dcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex( s1.fcval) <    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::DCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return     dcomplex( s1.dcval) <     dcomplex(  s2.fval); break;
                        case Field::DOUBLE:    return     dcomplex( s1.dcval) <     dcomplex(  s2.dval); break;
                        case Field::LDOUBLE:   return    ldcomplex( s1.dcval) <    ldcomplex( s2.ldval); break;
                        case Field::SCOMPLEX:  return     dcomplex( s1.dcval) <     dcomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return     dcomplex( s1.dcval) <     dcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex( s1.dcval) <    ldcomplex(s2.ldcval); break;
                    }
                    break;
                case Field::LDCOMPLEX:
                    switch (s2.F.type)
                    {
                        case Field::SINGLE:    return    ldcomplex(s1.ldcval) <    ldcomplex(  s2.fval); break;
                        case Field::DOUBLE:    return    ldcomplex(s1.ldcval) <    ldcomplex(  s2.dval); break;
                        case Field::LDOUBLE:   return    ldcomplex(s1.ldcval) <    ldcomplex( s2.ldval); break;
                        case Field::SCOMPLEX:  return    ldcomplex(s1.ldcval) <    ldcomplex( s2.fcval); break;
                        case Field::DCOMPLEX:  return    ldcomplex(s1.ldcval) <    ldcomplex( s2.dcval); break;
                        case Field::LDCOMPLEX: return    ldcomplex(s1.ldcval) <    ldcomplex(s2.ldcval); break;
                    }
                    break;
            }

            return false;
        }

        friend bool operator<=(const Scalar& s1, const Scalar& s2)
        {
            return !(s2 < s1);
        }

        friend bool operator>(const Scalar& s1, const Scalar& s2)
        {
            return s2 < s1;
        }

        friend bool operator>=(const Scalar& s1, const Scalar& s2)
        {
            return !(s1 < s2);
        }

        Field field() const { return F; }

        /*
        template <typename T> operator T() const
        {
            static_assert(is_field<T>::value, "");
            return to<T>();
        }
        */

        void* data() { return (void*)&fval; }

        const void* data() const { return (void*)&fval; }

        template <typename T> typename enable_if<is_complex<T>::value,T>::type to() const
        {
            switch (F.type)
            {
                case Field::SINGLE:    return T(  fval); break;
                case Field::DOUBLE:    return T(  dval); break;
                case Field::LDOUBLE:   return T( ldval); break;
                case Field::SCOMPLEX:  return T( fcval); break;
                case Field::DCOMPLEX:  return T( dcval); break;
                case Field::LDCOMPLEX: return T(ldcval); break;
            }
        }

        template <typename T> typename enable_if<!is_complex<T>::value,T>::type to() const
        {
            switch (F.type)
            {
                case Field::SINGLE:    return T(       fval ); break;
                case Field::DOUBLE:    return T(       dval ); break;
                case Field::LDOUBLE:   return T(      ldval ); break;
                case Field::SCOMPLEX:  return T(real( fcval)); break;
                case Field::DCOMPLEX:  return T(real( dcval)); break;
                case Field::LDCOMPLEX: return T(real(ldcval)); break;
            }
        }

        template <typename T>
        typename enable_if<is_complex<T>::value,Scalar&>::type operator=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval = real(other); break;
                case Field::DOUBLE:      dval = real(other); break;
                case Field::LDOUBLE:    ldval = real(other); break;
                case Field::SCOMPLEX:   fcval =      other ; break;
                case Field::DCOMPLEX:   dcval =      other ; break;
                case Field::LDCOMPLEX: ldcval =      other ; break;
            }

            return *this;
        }

        template <typename T>
        typename enable_if<!is_complex<T>::value,Scalar&>::type operator=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval = other; break;
                case Field::DOUBLE:      dval = other; break;
                case Field::LDOUBLE:    ldval = other; break;
                case Field::SCOMPLEX:   fcval = other; break;
                case Field::DCOMPLEX:   dcval = other; break;
                case Field::LDCOMPLEX: ldcval = other; break;
            }

            return *this;
        }

        Scalar& operator=(const Scalar& other)
        {
            switch (other.F.type)
            {
                case Field::SINGLE:    *this =   other.fval; break;
                case Field::DOUBLE:    *this =   other.dval; break;
                case Field::LDOUBLE:   *this =  other.ldval; break;
                case Field::SCOMPLEX:  *this =  other.fcval; break;
                case Field::DCOMPLEX:  *this =  other.dcval; break;
                case Field::LDCOMPLEX: *this = other.ldcval; break;
            }

            return *this;
        }

        template <typename T>
        typename enable_if<is_complex<T>::value,Scalar&>::type operator+=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval += abs(other); break;
                case Field::DOUBLE:      dval += abs(other); break;
                case Field::LDOUBLE:    ldval += abs(other); break;
                case Field::SCOMPLEX:   fcval +=          other ; break;
                case Field::DCOMPLEX:   dcval +=          other ; break;
                case Field::LDCOMPLEX: ldcval +=          other ; break;
            }

            return *this;
        }

        template <typename T>
        typename enable_if<!is_complex<T>::value,Scalar&>::type operator+=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval += other; break;
                case Field::DOUBLE:      dval += other; break;
                case Field::LDOUBLE:    ldval += other; break;
                case Field::SCOMPLEX:   fcval += other; break;
                case Field::DCOMPLEX:   dcval += other; break;
                case Field::LDCOMPLEX: ldcval += other; break;
            }

            return *this;
        }

        Scalar& operator+=(const Scalar& other)
        {
            switch (other.F.type)
            {
                case Field::SINGLE:    *this +=   other.fval; break;
                case Field::DOUBLE:    *this +=   other.dval; break;
                case Field::LDOUBLE:   *this +=  other.ldval; break;
                case Field::SCOMPLEX:  *this +=  other.fcval; break;
                case Field::DCOMPLEX:  *this +=  other.dcval; break;
                case Field::LDCOMPLEX: *this += other.ldcval; break;
            }

            return *this;
        }

        template <typename T>
        typename enable_if<is_complex<T>::value,Scalar&>::type operator-=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval -= abs(other); break;
                case Field::DOUBLE:      dval -= abs(other); break;
                case Field::LDOUBLE:    ldval -= abs(other); break;
                case Field::SCOMPLEX:   fcval -=          other ; break;
                case Field::DCOMPLEX:   dcval -=          other ; break;
                case Field::LDCOMPLEX: ldcval -=          other ; break;
            }

            return *this;
        }

        template <typename T>
        typename enable_if<!is_complex<T>::value,Scalar&>::type operator-=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval -= other; break;
                case Field::DOUBLE:      dval -= other; break;
                case Field::LDOUBLE:    ldval -= other; break;
                case Field::SCOMPLEX:   fcval -= other; break;
                case Field::DCOMPLEX:   dcval -= other; break;
                case Field::LDCOMPLEX: ldcval -= other; break;
            }

            return *this;
        }

        Scalar& operator-=(const Scalar& other)
        {
            switch (other.F.type)
            {
                case Field::SINGLE:    *this -=   other.fval; break;
                case Field::DOUBLE:    *this -=   other.dval; break;
                case Field::LDOUBLE:   *this -=  other.ldval; break;
                case Field::SCOMPLEX:  *this -=  other.fcval; break;
                case Field::DCOMPLEX:  *this -=  other.dcval; break;
                case Field::LDCOMPLEX: *this -= other.ldcval; break;
            }

            return *this;
        }

        template <typename T>
        typename enable_if<is_complex<T>::value,Scalar&>::type operator*=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval *= abs(other); break;
                case Field::DOUBLE:      dval *= abs(other); break;
                case Field::LDOUBLE:    ldval *= abs(other); break;
                case Field::SCOMPLEX:   fcval *=          other ; break;
                case Field::DCOMPLEX:   dcval *=          other ; break;
                case Field::LDCOMPLEX: ldcval *=          other ; break;
            }

            return *this;
        }

        template <typename T>
        typename enable_if<!is_complex<T>::value,Scalar&>::type operator*=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval *= other; break;
                case Field::DOUBLE:      dval *= other; break;
                case Field::LDOUBLE:    ldval *= other; break;
                case Field::SCOMPLEX:   fcval *= other; break;
                case Field::DCOMPLEX:   dcval *= other; break;
                case Field::LDCOMPLEX: ldcval *= other; break;
            }

            return *this;
        }

        Scalar& operator*=(const Scalar& other)
        {
            switch (other.F.type)
            {
                case Field::SINGLE:    *this *=   other.fval; break;
                case Field::DOUBLE:    *this *=   other.dval; break;
                case Field::LDOUBLE:   *this *=  other.ldval; break;
                case Field::SCOMPLEX:  *this *=  other.fcval; break;
                case Field::DCOMPLEX:  *this *=  other.dcval; break;
                case Field::LDCOMPLEX: *this *= other.ldcval; break;
            }

            return *this;
        }

        template <typename T>
        typename enable_if<is_complex<T>::value,Scalar&>::type operator/=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval /= abs(other); break;
                case Field::DOUBLE:      dval /= abs(other); break;
                case Field::LDOUBLE:    ldval /= abs(other); break;
                case Field::SCOMPLEX:   fcval /=          other ; break;
                case Field::DCOMPLEX:   dcval /=          other ; break;
                case Field::LDCOMPLEX: ldcval /=          other ; break;
            }

            return *this;
        }

        template <typename T>
        typename enable_if<!is_complex<T>::value,Scalar&>::type operator/=(T other)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            switch (F.type)
            {
                case Field::SINGLE:      fval /= other; break;
                case Field::DOUBLE:      dval /= other; break;
                case Field::LDOUBLE:    ldval /= other; break;
                case Field::SCOMPLEX:   fcval /= other; break;
                case Field::DCOMPLEX:   dcval /= other; break;
                case Field::LDCOMPLEX: ldcval /= other; break;
            }

            return *this;
        }

        Scalar& operator/=(const Scalar& other)
        {
            switch (other.F.type)
            {
                case Field::SINGLE:    *this /=   other.fval; break;
                case Field::DOUBLE:    *this /=   other.dval; break;
                case Field::LDOUBLE:   *this /=  other.ldval; break;
                case Field::SCOMPLEX:  *this /=  other.fcval; break;
                case Field::DCOMPLEX:  *this /=  other.dcval; break;
                case Field::LDCOMPLEX: *this /= other.ldcval; break;
            }

            return *this;
        }

        Scalar operator-() const
        {
            Scalar n(*this);

            switch (n.F.type)
            {
                case Field::SINGLE:      n.fval =   -n.fval; break;
                case Field::DOUBLE:      n.dval =   -n.dval; break;
                case Field::LDOUBLE:    n.ldval =  -n.ldval; break;
                case Field::SCOMPLEX:   n.fcval =  -n.fcval; break;
                case Field::DCOMPLEX:   n.dcval =  -n.dcval; break;
                case Field::LDCOMPLEX: n.ldcval = -n.ldcval; break;
            }

            return n;
        }

        template <typename T>
        Scalar operator+(T other) const
        {
            return other+(*this);
        }

        template <typename T> friend
        Scalar operator+(T other, const Scalar& s)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            Field::field new_type = s.resultType(other);

            switch (s.F.type)
            {
                case Field::SINGLE:    return Scalar(new_type, other+  s.fval); break;
                case Field::DOUBLE:    return Scalar(new_type, other+  s.dval); break;
                case Field::LDOUBLE:   return Scalar(new_type, other+ s.ldval); break;
                case Field::SCOMPLEX:  return Scalar(new_type, other+ s.fcval); break;
                case Field::DCOMPLEX:  return Scalar(new_type, other+ s.dcval); break;
                case Field::LDCOMPLEX: return Scalar(new_type, other+s.ldcval); break;
            }

            return Scalar(0.0);
        }

        Scalar operator+(const Scalar& other) const
        {
            switch (other.F.type)
            {
                case Field::SINGLE:    return *this +   other.fval; break;
                case Field::DOUBLE:    return *this +   other.dval; break;
                case Field::LDOUBLE:   return *this +  other.ldval; break;
                case Field::SCOMPLEX:  return *this +  other.fcval; break;
                case Field::DCOMPLEX:  return *this +  other.dcval; break;
                case Field::LDCOMPLEX: return *this + other.ldcval; break;
            }

            return Scalar(0.0);
        }

        template <typename T>
        Scalar operator-(T other) const
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            Field::field new_type = resultType(other);

            switch (F.type)
            {
                case Field::SINGLE:    return Scalar(new_type,   fval-other); break;
                case Field::DOUBLE:    return Scalar(new_type,   dval-other); break;
                case Field::LDOUBLE:   return Scalar(new_type,  ldval-other); break;
                case Field::SCOMPLEX:  return Scalar(new_type,  fcval-other); break;
                case Field::DCOMPLEX:  return Scalar(new_type,  dcval-other); break;
                case Field::LDCOMPLEX: return Scalar(new_type, ldcval-other); break;
            }

            return Scalar(0.0);
        }

        template <typename T> friend
        Scalar operator-(T other, const Scalar& s)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            Field::field new_type = s.resultType(other);

            switch (s.F.type)
            {
                case Field::SINGLE:    return Scalar(new_type, other-  s.fval); break;
                case Field::DOUBLE:    return Scalar(new_type, other-  s.dval); break;
                case Field::LDOUBLE:   return Scalar(new_type, other- s.ldval); break;
                case Field::SCOMPLEX:  return Scalar(new_type, other- s.fcval); break;
                case Field::DCOMPLEX:  return Scalar(new_type, other- s.dcval); break;
                case Field::LDCOMPLEX: return Scalar(new_type, other-s.ldcval); break;
            }

            return Scalar(0.0);
        }

        Scalar operator-(const Scalar& other) const
        {
            switch (other.F.type)
            {
                case Field::SINGLE:    return *this -   other.fval; break;
                case Field::DOUBLE:    return *this -   other.dval; break;
                case Field::LDOUBLE:   return *this -  other.ldval; break;
                case Field::SCOMPLEX:  return *this -  other.fcval; break;
                case Field::DCOMPLEX:  return *this -  other.dcval; break;
                case Field::LDCOMPLEX: return *this - other.ldcval; break;
            }

            return Scalar(0.0);
        }

        template <typename T>
        Scalar operator*(T other) const
        {
            return other*(*this);
        }

        template <typename T> friend
        Scalar operator*(T other, const Scalar& s)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            Field::field new_type = s.resultType(other);

            switch (s.F.type)
            {
                case Field::SINGLE:    return Scalar(new_type, other*  s.fval); break;
                case Field::DOUBLE:    return Scalar(new_type, other*  s.dval); break;
                case Field::LDOUBLE:   return Scalar(new_type, other* s.ldval); break;
                case Field::SCOMPLEX:  return Scalar(new_type, other* s.fcval); break;
                case Field::DCOMPLEX:  return Scalar(new_type, other* s.dcval); break;
                case Field::LDCOMPLEX: return Scalar(new_type, other*s.ldcval); break;
            }

            return Scalar(0.0);
        }

        Scalar operator*(const Scalar& other) const
        {
            switch (other.F.type)
            {
                case Field::SINGLE:    return *this *   other.fval; break;
                case Field::DOUBLE:    return *this *   other.dval; break;
                case Field::LDOUBLE:   return *this *  other.ldval; break;
                case Field::SCOMPLEX:  return *this *  other.fcval; break;
                case Field::DCOMPLEX:  return *this *  other.dcval; break;
                case Field::LDCOMPLEX: return *this * other.ldcval; break;
            }

            return Scalar(0.0);
        }

        template <typename T>
        Scalar operator/(T other) const
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            Field::field new_type = resultType(other);

            switch (F.type)
            {
                case Field::SINGLE:    return Scalar(new_type,   fval/other); break;
                case Field::DOUBLE:    return Scalar(new_type,   dval/other); break;
                case Field::LDOUBLE:   return Scalar(new_type,  ldval/other); break;
                case Field::SCOMPLEX:  return Scalar(new_type,  fcval/other); break;
                case Field::DCOMPLEX:  return Scalar(new_type,  dcval/other); break;
                case Field::LDCOMPLEX: return Scalar(new_type, ldcval/other); break;
            }

            return Scalar(0.0);
        }

        template <typename T> friend
        Scalar operator/(T other, const Scalar& s)
        {
            static_assert(is_arithmetic<typename real_type<T>::type>::value, "");

            Field::field new_type = s.resultType(other);

            switch (s.F.type)
            {
                case Field::SINGLE:    return Scalar(new_type, other/  s.fval); break;
                case Field::DOUBLE:    return Scalar(new_type, other/  s.dval); break;
                case Field::LDOUBLE:   return Scalar(new_type, other/ s.ldval); break;
                case Field::SCOMPLEX:  return Scalar(new_type, other/ s.fcval); break;
                case Field::DCOMPLEX:  return Scalar(new_type, other/ s.dcval); break;
                case Field::LDCOMPLEX: return Scalar(new_type, other/s.ldcval); break;
            }

            return Scalar(0.0);
        }

        Scalar operator/(const Scalar& other) const
        {
            switch (other.F.type)
            {
                case Field::SINGLE:    return *this /   other.fval; break;
                case Field::DOUBLE:    return *this /   other.dval; break;
                case Field::LDOUBLE:   return *this /  other.ldval; break;
                case Field::SCOMPLEX:  return *this /  other.fcval; break;
                case Field::DCOMPLEX:  return *this /  other.dcval; break;
                case Field::LDCOMPLEX: return *this / other.ldcval; break;
            }

            return Scalar(0.0);
        }

        friend ostream& operator<<(ostream& os, const Scalar& s)
        {
            switch (s.F.type)
            {
                case Field::SINGLE:    return os << s.fval;
                case Field::DOUBLE:    return os << s.dval;
                case Field::LDOUBLE:   return os << s.ldval;
                case Field::SCOMPLEX:  return os << s.fcval;
                case Field::DCOMPLEX:  return os << s.dcval;
                case Field::LDCOMPLEX: return os << s.ldcval;
            }
            return os;
        }
};

}
}

#endif
