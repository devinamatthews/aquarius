#ifndef _MPIWRAP_DATATYPE_HPP_
#define _MPIWRAP_DATATYPE_HPP_

#include "mpiwrap_common.hpp"

namespace MPIWrap
{

namespace internal
{
    template <typename T, typename=void> struct is_complex : std::false_type {};
    template <typename T> struct is_complex<std::complex<T>, typename std::enable_if<std::is_arithmetic<T>::value>::type> : std::true_type {};

    template <typename Derived> class Comm;
}

template <typename... T> struct MPI_TYPE_;

class Datatype
{
    friend class Intercomm;
    friend class Intracomm;
    template <typename Derived> friend class internal::Comm;
    template <typename... T> friend struct MPI_TYPE_;

    protected:
        MPI_Datatype type;
        bool predefined;

        static std::vector<MPI_Aint> displacements(const MPI_Int* counts,
                                                   const Datatype* types, MPI_Int nelem)
        {
            std::vector<MPI_Aint> bdispls(nelem);
            for (MPI_Int i = 1;i < nelem;i++)
            {
                bdispls[i] = bdispls[i-1]+counts[i-1]*types[i-1].extent();
            }
            return bdispls;
        }

        static std::vector<MPI_Aint> displacements(const std::vector<MPI_Int>& counts,
                                                   const std::vector<Datatype>& types)
        {
            MPIWRAP_ASSERT(counts.size() == types.size(),
                           "counts and types must have the same size.");
            return displacements(&counts[0], &types[0], counts.size());
        }

        static std::vector<MPI_Aint> displacements(const Datatype* types, MPI_Int nelem)
        {
            std::vector<MPI_Aint> bdispls(nelem);
            for (MPI_Int i = 1;i < nelem;i++)
            {
                bdispls[i] = bdispls[i-1]+types[i-1].extent();
            }
            return bdispls;
        }

        static std::vector<MPI_Aint> displacements(const std::vector<Datatype>& types)
        {
            return displacements(&types[0], types.size());
        }

        explicit Datatype(const MPI_Datatype& type, bool predefined = false)
        : type(type), predefined(predefined)
        {
            if (type != MPI_DATATYPE_NULL && !predefined)
            {
                MPIWRAP_CALL(MPI_Type_commit(&this->type));
            }
        }

        template <typename T>
        explicit Datatype(const T& x, typename std::enable_if<std::is_arithmetic<T>::value ||
                                                              internal::is_complex<T>::value>::type* foo = NULL)
        : type(MPI_TYPE_<T>::value().type), predefined(true) {}

    public:
#if MPIWRAP_CXX11

        Datatype(Datatype&& other)
        : type(other.type), predefined(other.predefined)
        {
            other.type = MPI_DATATYPE_NULL;
        }

#else

        Datatype(const Datatype& other)
        : type(other.type), predefined(other.predefined)
        {
        	if (!predefined)
        	{
                MPIWRAP_CALL(MPI_Type_dup(other, *this));
        	}
        }

#endif

        ~Datatype()
        {
            if (type != MPI_DATATYPE_NULL && !predefined)
            {
                MPIWRAP_CALL(MPI_Type_free(&type));
            }
        }

        operator MPI_Datatype&() { return type; }

        operator const MPI_Datatype&() const { return type; }

        operator MPI_Datatype*() { return &type; }

        operator const MPI_Datatype*() const { return &type; }

        friend Datatype operator*(const Datatype& old, MPI_Int n)
        {
            return old.contiguous(n);
        }

        friend Datatype operator*(MPI_Int n, const Datatype& old)
        {
            return old.contiguous(n);
        }

        Datatype contiguous(MPI_Int n) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_contiguous(n, type, &t));
            return Datatype(t);
        }

        Datatype vector(MPI_Int m, MPI_Int n, MPI_Int stride) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_vector(m, n, stride, type, &t));
            return Datatype(t);
        }

        Datatype hvector(MPI_Int m, MPI_Int n, MPI_Aint bstride) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_create_hvector(m, n, bstride, type, &t));
            return Datatype(t);
        }

        Datatype indexed(const MPI_Int* counts, const MPI_Int* displs, MPI_Int nblock) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_indexed(nblock, nc(counts), nc(displs), type, &t));
            return Datatype(t);
        }

        Datatype indexed(const std::vector<MPI_Int>& counts, const std::vector<MPI_Int>& displs) const
        {
            MPIWRAP_ASSERT(counts.size() == displs.size(),
                           "counts and displs must have the same size.");
            return indexed(&counts[0], &displs[0], counts.size());
        }

        Datatype hindexed(const MPI_Int* counts, const MPI_Aint* bdispls, MPI_Int nblock) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_create_hindexed(nblock, nc(counts), nc(bdispls), type, &t));
            return Datatype(t);
        }

        Datatype hindexed(const std::vector<MPI_Int>& counts, const std::vector<MPI_Aint>& bdispls) const
        {
            MPIWRAP_ASSERT(counts.size() == bdispls.size(),
                           "counts and bdispls must have the same size.");
            return hindexed(&counts[0], &bdispls[0], counts.size());
        }

        Datatype indexed(MPI_Int count, const MPI_Int* displs, MPI_Int nblock) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_create_indexed_block(nblock, count, nc(displs), type, &t));
            return Datatype(t);
        }

        Datatype indexed(MPI_Int count, const std::vector<MPI_Int>& displs) const
        {
            return indexed(count, &displs[0], displs.size());
        }

#if MPIWRAP_HAVE_MPI_TYPE_CREATE_HINDEXED_BLOCK

        Datatype hindexed(MPI_Int count, const MPI_Aint* bdispls, MPI_Int nblock) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_create_hindexed_block(nblock, count, nc(bdispls), type, &t));
            return Datatype(t);
        }

        Datatype hindexed(MPI_Int count, const std::vector<MPI_Aint>& bdispls) const
        {
            return hindexed(count, &bdispls[0], bdispls.size());
        }

#endif

        Datatype indexed(const MPI_Int* displs, MPI_Int nblock) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_create_indexed_block(nblock, 1, nc(displs), type, &t));
            return Datatype(t);
        }

        Datatype indexed(const std::vector<MPI_Int>& displs) const
        {
            return indexed(&displs[0], displs.size());
        }

#if MPIWRAP_HAVE_MPI_TYPE_CREATE_HINDEXED_BLOCK

        Datatype hindexed(const MPI_Aint* bdispls, MPI_Int nblock) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_create_hindexed_block(nblock, 1, nc(bdispls), type, &t));
            return Datatype(t);
        }

        Datatype hindexed(const std::vector<MPI_Aint>& bdispls) const
        {
            return hindexed(&bdispls[0], bdispls.size());
        }

#endif

        static Datatype structure(const MPI_Int* counts, const MPI_Aint* bdispls,
                                  const Datatype* types, MPI_Int nelem)
        {
            std::vector<MPI_Datatype> mpi_types(nelem);
            for (MPI_Int i = 0;i < nelem;i++) mpi_types[i] = types[i].type;
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_create_struct(nelem, nc(counts), nc(bdispls), nc(&mpi_types[0]), &t));
            return Datatype(t);
        }

        static Datatype structure(const std::vector<MPI_Int>& counts, const std::vector<MPI_Aint>& bdispls,
                                  const std::vector<Datatype>& types)
        {
            MPIWRAP_ASSERT(counts.size() == bdispls.size() &&
                           counts.size() == types.size(),
                           "counts, bdispls, and types must have the same size.");
            return structure(&counts[0], &bdispls[0], &types[0], counts.size());
        }

        static Datatype structure(const MPI_Aint* bdispls,
                                  const Datatype* types, MPI_Int nelem)
        {
            std::vector<MPI_Int> counts(nelem, 1);
            return structure(&counts[0], bdispls, types, nelem);
        }

        static Datatype structure(const std::vector<MPI_Aint>& bdispls,
                                  const std::vector<Datatype>& types)
        {
            return structure(std::vector<MPI_Int>(types.size(), 1), bdispls, types);
        }

        static Datatype structure(const MPI_Int* counts, const Datatype* types, MPI_Int nelem)
        {
            std::vector<MPI_Aint> bdispls = displacements(counts, types, nelem);
            return structure(counts, &bdispls[0], types, nelem);
        }

        static Datatype structure(const std::vector<MPI_Int>& counts, const std::vector<Datatype>& types)
        {
            return structure(counts, displacements(counts, types), types);
        }

        static Datatype structure(const Datatype* types, MPI_Int nelem)
        {
            std::vector<MPI_Int> counts(nelem, 1);
            std::vector<MPI_Aint> bdispls = displacements(types, nelem);
            return structure(&counts[0], &bdispls[0], types, nelem);
        }

        static Datatype structure(const std::vector<Datatype>& types)
        {
            return structure(std::vector<MPI_Int>(types.size(), 1), displacements(types), types);
        }

        Datatype subarray(const MPI_Int* sizes, const MPI_Int* subsizes,
                          const MPI_Int* starts, MPI_Int ndim, MPI_Int order = MPI_ORDER_C) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_create_subarray(ndim, nc(sizes), nc(subsizes), nc(starts), order, type, &t));
            return Datatype(t);
        }

        Datatype subarray(const std::vector<MPI_Int>& sizes, const std::vector<MPI_Int>& subsizes,
                          const std::vector<MPI_Int>& starts, MPI_Int order = MPI_ORDER_C) const
        {
            MPIWRAP_ASSERT(sizes.size() == subsizes.size() &&
                           starts.size() == sizes.size(),
                           "sizes, subsizes, and starts must have the same size.");
            return subarray(&sizes[0], &subsizes[0], &starts[0], sizes.size(), order);
        }

#if MPIWRAP_HAVE_MPI_COUNT

        MPI_Count size() const
        {
            MPI_Count c;
            MPIWRAP_CALL(MPI_Type_size_x(type, &c));
            return c;
        }

#else

        MPI_Int size() const
        {
            MPI_Int c;
            MPIWRAP_CALL(MPI_Type_size(type, &c));
            return c;
        }

#endif

        void extent(MPI_Aint& lb, MPI_Aint& extent) const
        {
            MPIWRAP_CALL(MPI_Type_get_extent(type, &lb, &extent));
        }

#if MPIWRAP_HAVE_MPI_COUNT

        void extent(MPI_Count& lb, MPI_Count& extent) const
        {
            MPIWRAP_CALL(MPI_Type_get_extent_x(type, &lb, &extent));
        }

        MPI_Count extent() const
        {
        	MPI_Count lb, extent;
            this->extent(lb, extent);
            return extent;
        }

        MPI_Count lb() const
        {
        	MPI_Count lb, extent;
            this->extent(lb, extent);
            return lb;
        }

        MPI_Count ub() const
        {
        	MPI_Count lb, extent;
            this->extent(lb, extent);
            return lb+extent;
        }

#else

        MPI_Aint extent() const
        {
            MPI_Aint lb, extent;
            this->extent(lb, extent);
            return extent;
        }

        MPI_Aint lb() const
        {
            MPI_Aint lb, extent;
            this->extent(lb, extent);
            return lb;
        }

        MPI_Aint ub() const
        {
            MPI_Aint lb, extent;
            this->extent(lb, extent);
            return lb+extent;
        }

#endif

        Datatype resize(MPI_Aint lb, MPI_Aint extent) const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_create_resized(type, lb, extent, &t));
            return Datatype(t);
        }

        void true_extent(MPI_Aint& lb, MPI_Aint& extent) const
        {
            MPIWRAP_CALL(MPI_Type_get_true_extent(type, &lb, &extent));
        }

#if MPIWRAP_HAVE_MPI_COUNT

        void true_extent(MPI_Count& lb, MPI_Count& extent) const
        {
            MPIWRAP_CALL(MPI_Type_get_true_extent_x(type, &lb, &extent));
        }

        MPI_Count true_extent() const
        {
        	MPI_Count lb, extent;
            this->true_extent(lb, extent);
            return extent;
        }

        MPI_Count true_lb() const
        {
        	MPI_Count lb, extent;
            this->true_extent(lb, extent);
            return lb;
        }

        MPI_Count true_ub() const
        {
        	MPI_Count lb, extent;
            this->true_extent(lb, extent);
            return lb+extent;
        }

#else

        MPI_Aint true_extent() const
        {
            MPI_Aint lb, extent;
            this->true_extent(lb, extent);
            return extent;
        }

        MPI_Aint true_lb() const
        {
            MPI_Aint lb, extent;
            this->true_extent(lb, extent);
            return lb;
        }

        MPI_Aint true_ub() const
        {
            MPI_Aint lb, extent;
            this->true_extent(lb, extent);
            return lb+extent;
        }

#endif

        Datatype duplicate() const
        {
            MPI_Datatype t;
            MPIWRAP_CALL(MPI_Type_dup(type, &t));
            return Datatype(t);
        }
};

#if MPIWRAP_CXX11

template <typename... T> struct MPI_TYPE_;

#else

template <typename T> struct MPI_TYPE_;

#endif

namespace internal
{
    #if MPIWRAP_CXX11

    template <typename T, typename... Ts>
    struct get_types
    {
        get_types(std::vector<Datatype>& types)
        {
            types.push_back(MPI_TYPE_<T>::value());
            get_types<Ts...> getter(types);
        }
    };

    #else

    template <typename T, typename U=void>
    struct get_types
    {
        get_types(std::vector<Datatype>& types)
        {
            types.push_back(MPI_TYPE_<T>::value());
            types.push_back(MPI_TYPE_<U>::value());
        }
    };

    #endif

    template <typename T> struct get_types<T>
    {
        get_types(std::vector<Datatype>& types)
        {
            types.push_back(MPI_TYPE_<T>::value());
        }
    };
}

template <typename T, typename U>
struct MPI_TYPE_<std::pair<T,U>>
{
    static Datatype value()
    {
        std::vector<Datatype> types;
        internal::get_types<T,U> getter(types);
        return Datatype::structure(types);
    }
};

#if MPIWRAP_CXX11

template <typename... T>
struct MPI_TYPE_
{
    static Datatype value()
    {
        std::vector<Datatype> types;
        internal::get_types<T...> getter(types);
        return Datatype::structure(types);
    }
};

template <typename... T>
struct MPI_TYPE_<std::tuple<T...>>
{
    static Datatype value()
    {
        std::vector<Datatype> types;
        internal::get_types<T...> getter(types);
        return Datatype::structure(types);
    }
};

#endif

template <>
struct MPI_TYPE_<float>
{
    static Datatype value() { return Datatype(MPI_FLOAT, true); }
};

template <>
struct MPI_TYPE_<double>
{
    static Datatype value() { return Datatype(MPI_DOUBLE, true); }
};

template <>
struct MPI_TYPE_<long double>
{
    static Datatype value() { return Datatype(MPI_LONG_DOUBLE, true); }
};

/*
 * Cheat and use the C99 complex types. The binary layout should be identical.
 */

template <>
struct MPI_TYPE_< std::complex<float> >
{
    static Datatype value() { return Datatype(MPI_C_FLOAT_COMPLEX, true); }
};

template <>
struct MPI_TYPE_< std::complex<double> >
{
    static Datatype value() { return Datatype(MPI_C_DOUBLE_COMPLEX, true); }
};

template <>
struct MPI_TYPE_< std::complex<long double> >
{
    static Datatype value() { return Datatype(MPI_C_LONG_DOUBLE_COMPLEX, true); }
};

template <>
struct MPI_TYPE_<void>
{
    static Datatype value() { return Datatype(MPI_BYTE, true); }
};

template <>
struct MPI_TYPE_<char>
{
    static Datatype value() { return Datatype(MPI_CHAR, true); }
};

template <>
struct MPI_TYPE_<signed char>
{
    static Datatype value() { return Datatype(MPI_SIGNED_CHAR, true); }
};

template <>
struct MPI_TYPE_<unsigned char>
{
    static Datatype value() { return Datatype(MPI_UNSIGNED_CHAR, true); }
};

template <>
struct MPI_TYPE_<wchar_t>
{
    static Datatype value() { return Datatype(MPI_WCHAR, true); }
};

template <>
struct MPI_TYPE_<short>
{
    static Datatype value() { return Datatype(MPI_SHORT, true); }
};

template <>
struct MPI_TYPE_<unsigned short>
{
    static Datatype value() { return Datatype(MPI_UNSIGNED_SHORT, true); }
};

template <>
struct MPI_TYPE_<int>
{
    static Datatype value() { return Datatype(MPI_INT, true); }
};

template <>
struct MPI_TYPE_<unsigned int>
{
    static Datatype value() { return Datatype(MPI_UNSIGNED, true); }
};

template <>
struct MPI_TYPE_<long>
{
    static Datatype value() { return Datatype(MPI_LONG, true); }
};

template <>
struct MPI_TYPE_<unsigned long>
{
    static Datatype value() { return Datatype(MPI_UNSIGNED_LONG, true); }
};

template <>
struct MPI_TYPE_<long long>
{
    static Datatype value() { return Datatype(MPI_LONG_LONG, true); }
};

template <>
struct MPI_TYPE_<unsigned long long>
{
    static Datatype value() { return Datatype(MPI_UNSIGNED_LONG_LONG, true); }
};

template <typename T>
void Reduce(const T* in, T* inout, MPI_Int count, MPI_Op op)
{
    Reduce(in, inout, count, op, MPI_TYPE_<T>::value());
}

template <typename T>
void Reduce(const std::vector<T>& in, std::vector<T>& inout, MPI_Op op)
{
    Reduce(in, inout, op, MPI_TYPE_<T>::value());
}

template <typename T>
void Reduce(const T* in, T* inout, MPI_Int count, MPI_Op op, const Datatype& type)
{
    MPIWRAP_CALL(MPI_Reduce_local(nc(in), inout, count, type, op));
}

template <typename T>
void Reduce(const std::vector<T>& in, std::vector<T>& inout, MPI_Op op, const Datatype& type)
{
    MPIWRAP_ASSERT(in.size() == inout.size(), "Buffers must have the same size.");
    Reduce(&in[0], &inout[0], in.size(), op, type);
}

}

#endif
