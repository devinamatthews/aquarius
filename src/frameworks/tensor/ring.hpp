#ifndef _AQUARIUS_TENSOR_RING_HPP_
#define _AQUARIUS_TENSOR_RING_HPP_

#include "frameworks/util.hpp"

namespace aquarius
{

class Scalar;
namespace tensor { class KeyValueVector; }

typedef complex< float> scomplex;
typedef complex<double> dcomplex;

class Ring
{
    public:
        int type() const { return type_; }

        int size() const { return size_; }

        bool operator==(const Ring& other) const { return type_ == other.type_; }

        bool operator!=(const Ring& other) const { return type_ != other.type_; }

    protected:
        Ring(int type, size_t size) : type_(type), size_(size) {}

        int type_;
        size_t size_;
};

namespace detail
{
    static size_t field_sizes[] = {sizeof(   float),
                                   sizeof(  double),
                                   sizeof(scomplex),
                                   sizeof(dcomplex)};
}

class Field : public Ring
{
    friend Scalar;
    friend tensor::KeyValueVector;

    public:
        enum field : int {SINGLE=0, DOUBLE=1, SCOMPLEX=2, DCOMPLEX=3};

        Field(field F) : Ring(F, detail::field_sizes[F]) {}

        Field(const    float& val) : Ring(  SINGLE, sizeof(   float)) {}
        Field(const   double& val) : Ring(  DOUBLE, sizeof(  double)) {}
        Field(const scomplex& val) : Ring(SCOMPLEX, sizeof(scomplex)) {}
        Field(const dcomplex& val) : Ring(DCOMPLEX, sizeof(dcomplex)) {}

        field type() const { return (field)type_; }

        static Datatype MPI_TYPE(const Field& F)
        {
            return MPI_TYPE(F.type());
        }

        static Datatype MPI_TYPE(field F)
        {
            switch (F)
            {
                case   SINGLE: return MPI_TYPE_<   float>::value();
                case   DOUBLE: return MPI_TYPE_<  double>::value();
                case SCOMPLEX: return MPI_TYPE_<scomplex>::value();
                case DCOMPLEX: return MPI_TYPE_<dcomplex>::value();
            }
            return MPI_TYPE_<double>::value();
        }
};

template <typename T> struct field_type           { static constexpr Field::field val =   Field::DOUBLE; };
template <>           struct field_type<   float> { static constexpr Field::field val =   Field::SINGLE; };
template <>           struct field_type<  double> { static constexpr Field::field val =   Field::DOUBLE; };
template <>           struct field_type<scomplex> { static constexpr Field::field val = Field::SCOMPLEX; };
template <>           struct field_type<dcomplex> { static constexpr Field::field val = Field::DCOMPLEX; };

template <typename T> struct is_field           : false_type {};
template <>           struct is_field<   float> :  true_type {};
template <>           struct is_field<  double> :  true_type {};
template <>           struct is_field<scomplex> :  true_type {};
template <>           struct is_field<dcomplex> :  true_type {};

template <typename T, typename U=void>
using enable_if_field_t = enable_if_t<is_field<T>::value,U>;

class Scalar
{
    protected:
        Field F = Field::DOUBLE;
        union
        {
            float fval;
            double dval;
            scomplex fcval;
            dcomplex dcval;
        };

        constexpr static Field::field resultType(int f1, int f2)
        {
            return Field::field(f1|f2);
        }

        template <typename T>
        enable_if_t<is_field<T>::value,Field::field> resultType(T other) const
        {
            return resultType(F.type(), field_type<T>::val);
        }

        template <typename T>
        enable_if_t<!is_field<T>::value,Field::field> resultType(T other) const
        {
            return Field::field(F.type());
        }

    public:
        explicit Scalar(Field::field type = Field::DOUBLE)
        {
            reset(type);
        }

        template <typename T>
        Scalar(Field::field type, T val)
        {
            reset(type, val);
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        explicit Scalar(T val)
        {
            reset(val);
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        Scalar(T r, T i)
        {
            reset(r,i);
        }

        void reset(Field::field type = Field::DOUBLE)
        {
            F = type;
            *this = 0.0;
        }

        template <typename T>
        void reset(Field::field type, T val)
        {
            F = type;
            *this = val;
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        void reset(T val)
        {
            F = field_type<T>::val;
            *this = val;
        }

        template <typename T, typename=enable_if_arithmetic_t<T>>
        void reset(T r, T i)
        {
            F = field_type<std::complex<T>>::val;
            *this = std::complex<T>(r,i);
        }

        friend Scalar abs(const Scalar& s)
        {
            switch (s.F.type())
            {
                case Field::SINGLE:   return Scalar(std::abs( s.fval)); break;
                case Field::DOUBLE:   return Scalar(std::abs( s.dval)); break;
                case Field::SCOMPLEX: return Scalar(std::abs(s.fcval)); break;
                case Field::DCOMPLEX: return Scalar(std::abs(s.dcval)); break;
            }
            return Scalar(0.0);
        }

        friend Scalar conj(Scalar s)
        {
            switch (s.F.type())
            {
                case Field::SINGLE:                            break;
                case Field::DOUBLE:                            break;
                case Field::SCOMPLEX: s.fcval = conj(s.fcval); break;
                case Field::DCOMPLEX: s.dcval = conj(s.dcval); break;
            }
            return s;
        }

        friend Scalar sqrt(Scalar s)
        {
            switch (s.F.type())
            {
                case Field::SINGLE:    s.fval = sqrt( s.fval); break;
                case Field::DOUBLE:    s.dval = sqrt( s.dval); break;
                case Field::SCOMPLEX: s.fcval = sqrt(s.fcval); break;
                case Field::DCOMPLEX: s.dcval = sqrt(s.dcval); break;
            }
            return s;
        }

        friend Scalar real(const Scalar& s)
        {
            Scalar r(s.F.type() & Field::DOUBLE);
            switch (s.F.type())
            {
                case Field::SINGLE:   r.fval = real( s.fval); break;
                case Field::DOUBLE:   r.dval = real( s.dval); break;
                case Field::SCOMPLEX: r.fval = real(s.fcval); break;
                case Field::DCOMPLEX: r.dval = real(s.dcval); break;
            }
            return r;
        }

        friend Scalar imag(const Scalar& s)
        {
            Scalar i(s.F.type() & Field::DOUBLE);
            switch (s.F.type())
            {
                case Field::SINGLE:   i.fval = imag( s.fval); break;
                case Field::DOUBLE:   i.dval = imag( s.dval); break;
                case Field::SCOMPLEX: i.fval = imag(s.fcval); break;
                case Field::DCOMPLEX: i.dval = imag(s.dcval); break;
            }
            return i;
        }

        friend Scalar min(const Scalar& s1, const Scalar& s2)
        {
            switch (s1.F.type())
            {
                case Field::SINGLE:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return Scalar(min(   float(s1.fval),    float( s2.fval))); break;
                        case Field::DOUBLE:   return Scalar(min(  double(s1.fval),   double( s2.dval))); break;
                        case Field::SCOMPLEX: return Scalar(min(scomplex(s1.fval), scomplex(s2.fcval))); break;
                        case Field::DCOMPLEX: return Scalar(min(dcomplex(s1.fval), dcomplex(s2.dcval))); break;
                    }
                    break;
                case Field::DOUBLE:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return Scalar(min(  double(s1.dval),   double( s2.fval))); break;
                        case Field::DOUBLE:   return Scalar(min(  double(s1.dval),   double( s2.dval))); break;
                        case Field::SCOMPLEX: return Scalar(min(dcomplex(s1.dval), dcomplex(s2.fcval))); break;
                        case Field::DCOMPLEX: return Scalar(min(dcomplex(s1.dval), dcomplex(s2.dcval))); break;
                    }
                    break;
                case Field::SCOMPLEX:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return Scalar(min(scomplex(s1.fcval), scomplex( s2.fval))); break;
                        case Field::DOUBLE:   return Scalar(min(dcomplex(s1.fcval), dcomplex( s2.dval))); break;
                        case Field::SCOMPLEX: return Scalar(min(scomplex(s1.fcval), scomplex(s2.fcval))); break;
                        case Field::DCOMPLEX: return Scalar(min(dcomplex(s1.fcval), dcomplex(s2.dcval))); break;
                    }
                    break;
                case Field::DCOMPLEX:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return Scalar(min(dcomplex(s1.dcval), dcomplex( s2.fval))); break;
                        case Field::DOUBLE:   return Scalar(min(dcomplex(s1.dcval), dcomplex( s2.dval))); break;
                        case Field::SCOMPLEX: return Scalar(min(dcomplex(s1.dcval), dcomplex(s2.fcval))); break;
                        case Field::DCOMPLEX: return Scalar(min(dcomplex(s1.dcval), dcomplex(s2.dcval))); break;
                    }
                    break;
            }

            return Scalar();
        }

        friend Scalar max(const Scalar& s1, const Scalar& s2)
        {
            switch (s1.F.type())
            {
                case Field::SINGLE:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return Scalar(max(   float(s1.fval),    float( s2.fval))); break;
                        case Field::DOUBLE:   return Scalar(max(  double(s1.fval),   double( s2.dval))); break;
                        case Field::SCOMPLEX: return Scalar(max(scomplex(s1.fval), scomplex(s2.fcval))); break;
                        case Field::DCOMPLEX: return Scalar(max(dcomplex(s1.fval), dcomplex(s2.dcval))); break;
                    }
                    break;
                case Field::DOUBLE:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return Scalar(max(  double(s1.dval),   double( s2.fval))); break;
                        case Field::DOUBLE:   return Scalar(max(  double(s1.dval),   double( s2.dval))); break;
                        case Field::SCOMPLEX: return Scalar(max(dcomplex(s1.dval), dcomplex(s2.fcval))); break;
                        case Field::DCOMPLEX: return Scalar(max(dcomplex(s1.dval), dcomplex(s2.dcval))); break;
                    }
                    break;
                case Field::SCOMPLEX:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return Scalar(max(scomplex(s1.fcval), scomplex( s2.fval))); break;
                        case Field::DOUBLE:   return Scalar(max(dcomplex(s1.fcval), dcomplex( s2.dval))); break;
                        case Field::SCOMPLEX: return Scalar(max(scomplex(s1.fcval), scomplex(s2.fcval))); break;
                        case Field::DCOMPLEX: return Scalar(max(dcomplex(s1.fcval), dcomplex(s2.dcval))); break;
                    }
                    break;
                case Field::DCOMPLEX:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return Scalar(max(dcomplex(s1.dcval), dcomplex( s2.fval))); break;
                        case Field::DOUBLE:   return Scalar(max(dcomplex(s1.dcval), dcomplex( s2.dval))); break;
                        case Field::SCOMPLEX: return Scalar(max(dcomplex(s1.dcval), dcomplex(s2.fcval))); break;
                        case Field::DCOMPLEX: return Scalar(max(dcomplex(s1.dcval), dcomplex(s2.dcval))); break;
                    }
                    break;
            }

            return Scalar();
        }

        friend bool operator==(const Scalar& s1, const Scalar& s2)
        {
            switch (s1.F.type())
            {
                case Field::SINGLE:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return    float(s1.fval) ==    float( s2.fval); break;
                        case Field::DOUBLE:   return   double(s1.fval) ==   double( s2.dval); break;
                        case Field::SCOMPLEX: return scomplex(s1.fval) == scomplex(s2.fcval); break;
                        case Field::DCOMPLEX: return dcomplex(s1.fval) == dcomplex(s2.dcval); break;
                    }
                    break;
                case Field::DOUBLE:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return   double(s1.dval) ==   double( s2.fval); break;
                        case Field::DOUBLE:   return   double(s1.dval) ==   double( s2.dval); break;
                        case Field::SCOMPLEX: return dcomplex(s1.dval) == dcomplex(s2.fcval); break;
                        case Field::DCOMPLEX: return dcomplex(s1.dval) == dcomplex(s2.dcval); break;
                    }
                    break;
                case Field::SCOMPLEX:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return scomplex(s1.fcval) == scomplex( s2.fval); break;
                        case Field::DOUBLE:   return dcomplex(s1.fcval) == dcomplex( s2.dval); break;
                        case Field::SCOMPLEX: return scomplex(s1.fcval) == scomplex(s2.fcval); break;
                        case Field::DCOMPLEX: return dcomplex(s1.fcval) == dcomplex(s2.dcval); break;
                    }
                    break;
                case Field::DCOMPLEX:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return dcomplex(s1.dcval) == dcomplex( s2.fval); break;
                        case Field::DOUBLE:   return dcomplex(s1.dcval) == dcomplex( s2.dval); break;
                        case Field::SCOMPLEX: return dcomplex(s1.dcval) == dcomplex(s2.fcval); break;
                        case Field::DCOMPLEX: return dcomplex(s1.dcval) == dcomplex(s2.dcval); break;
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
            switch (s1.F.type())
            {
                case Field::SINGLE:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return    float(s1.fval) <    float( s2.fval); break;
                        case Field::DOUBLE:   return   double(s1.fval) <   double( s2.dval); break;
                        case Field::SCOMPLEX: return scomplex(s1.fval) < scomplex(s2.fcval); break;
                        case Field::DCOMPLEX: return dcomplex(s1.fval) < dcomplex(s2.dcval); break;
                    }
                    break;
                case Field::DOUBLE:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return   double(s1.dval) <   double( s2.fval); break;
                        case Field::DOUBLE:   return   double(s1.dval) <   double( s2.dval); break;
                        case Field::SCOMPLEX: return dcomplex(s1.dval) < dcomplex(s2.fcval); break;
                        case Field::DCOMPLEX: return dcomplex(s1.dval) < dcomplex(s2.dcval); break;
                    }
                    break;
                case Field::SCOMPLEX:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return scomplex(s1.fcval) < scomplex( s2.fval); break;
                        case Field::DOUBLE:   return dcomplex(s1.fcval) < dcomplex( s2.dval); break;
                        case Field::SCOMPLEX: return scomplex(s1.fcval) < scomplex(s2.fcval); break;
                        case Field::DCOMPLEX: return dcomplex(s1.fcval) < dcomplex(s2.dcval); break;
                    }
                    break;
                case Field::DCOMPLEX:
                    switch (s2.F.type())
                    {
                        case Field::SINGLE:   return dcomplex(s1.dcval) < dcomplex( s2.fval); break;
                        case Field::DOUBLE:   return dcomplex(s1.dcval) < dcomplex( s2.dval); break;
                        case Field::SCOMPLEX: return dcomplex(s1.dcval) < dcomplex(s2.fcval); break;
                        case Field::DCOMPLEX: return dcomplex(s1.dcval) < dcomplex(s2.dcval); break;
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

        void* data() { return (void*)&fval; }

        const void* data() const { return (void*)&fval; }

        template <typename T> enable_if_t<is_complex<T>::value,T> to() const
        {
            switch (F.type())
            {
                case Field::SINGLE:   return T( fval); break;
                case Field::DOUBLE:   return T( dval); break;
                case Field::SCOMPLEX: return T(fcval); break;
                case Field::DCOMPLEX: return T(dcval); break;
            }
            return T();
        }

        template <typename T> enable_if_t<!is_complex<T>::value,T> to() const
        {
            switch (F.type())
            {
                case Field::SINGLE:   return T(      fval ); break;
                case Field::DOUBLE:   return T(      dval ); break;
                case Field::SCOMPLEX: return T(real(fcval)); break;
                case Field::DCOMPLEX: return T(real(dcval)); break;
            }
            return T();
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        explicit operator T() const
        {
            return to<T>();
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        Scalar& operator=(T other)
        {
            switch (F.type())
            {
                case Field::SINGLE:    fval = real(other); break;
                case Field::DOUBLE:    dval = real(other); break;
                case Field::SCOMPLEX: fcval =      other ; break;
                case Field::DCOMPLEX: dcval =      other ; break;
            }

            return *this;
        }

        Scalar& operator=(const Scalar& other)
        {
            switch (other.F.type())
            {
                case Field::SINGLE:   *this =  other.fval; break;
                case Field::DOUBLE:   *this =  other.dval; break;
                case Field::SCOMPLEX: *this = other.fcval; break;
                case Field::DCOMPLEX: *this = other.dcval; break;
            }

            return *this;
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        Scalar& operator+=(T other)
        {
            switch (F.type())
            {
                case Field::SINGLE:    fval += real(other); break;
                case Field::DOUBLE:    dval += real(other); break;
                case Field::SCOMPLEX: fcval +=      other ; break;
                case Field::DCOMPLEX: dcval +=      other ; break;
            }

            return *this;
        }

        Scalar& operator+=(const Scalar& other)
        {
            switch (other.F.type())
            {
                case Field::SINGLE:   *this +=  other.fval; break;
                case Field::DOUBLE:   *this +=  other.dval; break;
                case Field::SCOMPLEX: *this += other.fcval; break;
                case Field::DCOMPLEX: *this += other.dcval; break;
            }

            return *this;
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        Scalar& operator-=(T other)
        {
            switch (F.type())
            {
                case Field::SINGLE:    fval -= real(other); break;
                case Field::DOUBLE:    dval -= real(other); break;
                case Field::SCOMPLEX: fcval -=      other ; break;
                case Field::DCOMPLEX: dcval -=      other ; break;
            }

            return *this;
        }

        Scalar& operator-=(const Scalar& other)
        {
            switch (other.F.type())
            {
                case Field::SINGLE:   *this -=  other.fval; break;
                case Field::DOUBLE:   *this -=  other.dval; break;
                case Field::SCOMPLEX: *this -= other.fcval; break;
                case Field::DCOMPLEX: *this -= other.dcval; break;
            }

            return *this;
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        enable_if_complex_t<T,Scalar&>
        operator*=(T other)
        {
            switch (F.type())
            {
                case Field::SINGLE:    fval *= std::abs(other); break;
                case Field::DOUBLE:    dval *= std::abs(other); break;
                case Field::SCOMPLEX: fcval *=          other ; break;
                case Field::DCOMPLEX: dcval *=          other ; break;
            }

            return *this;
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        enable_if_not_complex_t<T,Scalar&>
        operator*=(T other)
        {
            switch (F.type())
            {
                case Field::SINGLE:    fval *= other; break;
                case Field::DOUBLE:    dval *= other; break;
                case Field::SCOMPLEX: fcval *= other; break;
                case Field::DCOMPLEX: dcval *= other; break;
            }

            return *this;
        }

        Scalar& operator*=(const Scalar& other)
        {
            switch (other.F.type())
            {
                case Field::SINGLE:   *this *=  other.fval; break;
                case Field::DOUBLE:   *this *=  other.dval; break;
                case Field::SCOMPLEX: *this *= other.fcval; break;
                case Field::DCOMPLEX: *this *= other.dcval; break;
            }

            return *this;
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        enable_if_complex_t<T,Scalar&>
        operator/=(T other)
        {
            switch (F.type())
            {
                case Field::SINGLE:    fval /= std::abs(other); break;
                case Field::DOUBLE:    dval /= std::abs(other); break;
                case Field::SCOMPLEX: fcval /=          other ; break;
                case Field::DCOMPLEX: dcval /=          other ; break;
            }

            return *this;
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        enable_if_not_complex_t<T,Scalar&>
        operator/=(T other)
        {
            switch (F.type())
            {
                case Field::SINGLE:    fval /= other; break;
                case Field::DOUBLE:    dval /= other; break;
                case Field::SCOMPLEX: fcval /= other; break;
                case Field::DCOMPLEX: dcval /= other; break;
            }

            return *this;
        }

        Scalar& operator/=(const Scalar& other)
        {
            switch (other.F.type())
            {
                case Field::SINGLE:   *this /=  other.fval; break;
                case Field::DOUBLE:   *this /=  other.dval; break;
                case Field::SCOMPLEX: *this /= other.fcval; break;
                case Field::DCOMPLEX: *this /= other.dcval; break;
            }

            return *this;
        }

        Scalar operator-() const
        {
            Scalar n(*this);

            switch (n.F.type())
            {
                case Field::SINGLE:    n.fval =  -n.fval; break;
                case Field::DOUBLE:    n.dval =  -n.dval; break;
                case Field::SCOMPLEX: n.fcval = -n.fcval; break;
                case Field::DCOMPLEX: n.dcval = -n.dcval; break;
            }

            return n;
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        Scalar operator+(T other) const
        {
            return other+(*this);
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        friend Scalar operator+(T other, const Scalar& s)
        {
            Field::field new_type = s.resultType(other);

            switch (s.F.type())
            {
                case Field::SINGLE:   return Scalar(new_type, other+ s.fval); break;
                case Field::DOUBLE:   return Scalar(new_type, other+ s.dval); break;
                case Field::SCOMPLEX: return Scalar(new_type, other+s.fcval); break;
                case Field::DCOMPLEX: return Scalar(new_type, other+s.dcval); break;
            }

            return Scalar(0.0);
        }

        Scalar operator+(const Scalar& other) const
        {
            switch (other.F.type())
            {
                case Field::SINGLE:   return *this +  other.fval; break;
                case Field::DOUBLE:   return *this +  other.dval; break;
                case Field::SCOMPLEX: return *this + other.fcval; break;
                case Field::DCOMPLEX: return *this + other.dcval; break;
            }

            return Scalar(0.0);
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        Scalar operator-(T other) const
        {
            Field::field new_type = resultType(other);

            switch (F.type())
            {
                case Field::SINGLE:   return Scalar(new_type,  fval-other); break;
                case Field::DOUBLE:   return Scalar(new_type,  dval-other); break;
                case Field::SCOMPLEX: return Scalar(new_type, fcval-other); break;
                case Field::DCOMPLEX: return Scalar(new_type, dcval-other); break;
            }

            return Scalar(0.0);
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        friend Scalar operator-(T other, const Scalar& s)
        {
            Field::field new_type = s.resultType(other);

            switch (s.F.type())
            {
                case Field::SINGLE:   return Scalar(new_type, other- s.fval); break;
                case Field::DOUBLE:   return Scalar(new_type, other- s.dval); break;
                case Field::SCOMPLEX: return Scalar(new_type, other-s.fcval); break;
                case Field::DCOMPLEX: return Scalar(new_type, other-s.dcval); break;
            }

            return Scalar(0.0);
        }

        Scalar operator-(const Scalar& other) const
        {
            switch (other.F.type())
            {
                case Field::SINGLE:   return *this -  other.fval; break;
                case Field::DOUBLE:   return *this -  other.dval; break;
                case Field::SCOMPLEX: return *this - other.fcval; break;
                case Field::DCOMPLEX: return *this - other.dcval; break;
            }

            return Scalar(0.0);
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        Scalar operator*(T other) const
        {
            return other*(*this);
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        friend Scalar operator*(T other, const Scalar& s)
        {
            Field::field new_type = s.resultType(other);

            switch (s.F.type())
            {
                case Field::SINGLE:   return Scalar(new_type, other* s.fval); break;
                case Field::DOUBLE:   return Scalar(new_type, other* s.dval); break;
                case Field::SCOMPLEX: return Scalar(new_type, other*s.fcval); break;
                case Field::DCOMPLEX: return Scalar(new_type, other*s.dcval); break;
            }

            return Scalar(0.0);
        }

        Scalar operator*(const Scalar& other) const
        {
            switch (other.F.type())
            {
                case Field::SINGLE:   return *this *  other.fval; break;
                case Field::DOUBLE:   return *this *  other.dval; break;
                case Field::SCOMPLEX: return *this * other.fcval; break;
                case Field::DCOMPLEX: return *this * other.dcval; break;
            }

            return Scalar(0.0);
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        Scalar operator/(T other) const
        {
            Field::field new_type = resultType(other);

            switch (F.type())
            {
                case Field::SINGLE:   return Scalar(new_type,  fval/other); break;
                case Field::DOUBLE:   return Scalar(new_type,  dval/other); break;
                case Field::SCOMPLEX: return Scalar(new_type, fcval/other); break;
                case Field::DCOMPLEX: return Scalar(new_type, dcval/other); break;
            }

            return Scalar(0.0);
        }

        template <typename T, typename=enable_if_arithmetic_t<real_type_t<T>>>
        friend Scalar operator/(T other, const Scalar& s)
        {
            Field::field new_type = s.resultType(other);

            switch (s.F.type())
            {
                case Field::SINGLE:   return Scalar(new_type, other/ s.fval); break;
                case Field::DOUBLE:   return Scalar(new_type, other/ s.dval); break;
                case Field::SCOMPLEX: return Scalar(new_type, other/s.fcval); break;
                case Field::DCOMPLEX: return Scalar(new_type, other/s.dcval); break;
            }

            return Scalar(0.0);
        }

        Scalar operator/(const Scalar& other) const
        {
            switch (other.F.type())
            {
                case Field::SINGLE:   return *this / other.fval; break;
                case Field::DOUBLE:   return *this / other.dval; break;
                case Field::SCOMPLEX: return *this /other.fcval; break;
                case Field::DCOMPLEX: return *this /other.dcval; break;
            }

            return Scalar(0.0);
        }

        friend ostream& operator<<(ostream& os, const Scalar& s)
        {
            switch (s.F.type())
            {
                case Field::SINGLE:   return os << s.fval;
                case Field::DOUBLE:   return os << s.dval;
                case Field::SCOMPLEX: return os << s.fcval;
                case Field::DCOMPLEX: return os << s.dcval;
            }
            return os;
        }
};

}

#endif
