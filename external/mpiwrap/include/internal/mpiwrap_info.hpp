#ifndef _MPIWRAP_INFO_HPP_
#define _MPIWRAP_INFO_HPP_

#include "mpiwrap_common.hpp"

namespace MPIWrap
{

class Info
{
    friend class Window;

    protected:
        MPI_Info info;

        Info(const MPI_Info& info) : info(info) {}

    public:
        class ConstInfoRef
        {
            friend class Info;
            template <typename T> friend class iterator_;

            protected:
                MPI_Info info;
                std::string key;

                ConstInfoRef(const MPI_Info& info, const std::string& key)
                : info(info), key(key) {}

            public:
                operator std::string() const
                {
                    std::string val;
                    MPI_Int len, flag;
                    MPIWRAP_CALL(MPI_Info_get_valuelen(info, nc(key.c_str()), &len, &flag));
                    if (!flag) return val;

                    val.resize(len);
                    MPIWRAP_CALL(MPI_Info_get(info, nc(key.c_str()), len, &val[0], &flag));
                    return val;
                }
        };

        class InfoRef : public ConstInfoRef
        {
            friend class Info;
            template <typename T> friend class iterator_;

            protected:
                InfoRef(const MPI_Info& info, const std::string& key)
                : ConstInfoRef(info, key) {}

            public:
                using ConstInfoRef::operator std::string;

                InfoRef& operator=(const std::string& val)
                {
                    MPIWRAP_CALL(MPI_Info_set(info, nc(key.c_str()), nc(val.c_str())));
                    return *this;
                }
        };

        Info()
        : info(MPI_INFO_NULL)
        {
            MPIWRAP_CALL(MPI_Info_create(&info));
        }

#if MPIWRAP_CXX11

        Info(Info&& other)
        : info(other.info)
        {
            other.info = MPI_INFO_NULL;
        }

#endif

        ~Info()
        {
            if (info != MPI_INFO_NULL)
            {
                MPIWRAP_CALL(MPI_Info_free(&info));
            }
        }

        operator MPI_Info&() { return info; }

        operator const MPI_Info&() const { return info; }

        operator MPI_Info*() { return &info; }

        operator const MPI_Info*() const { return &info; }

        static Info null()
        {
            return Info(MPI_INFO_NULL);
        }

        InfoRef operator[](const std::string& key)
        {
            return InfoRef(info, key);
        }

        ConstInfoRef operator[](const std::string& key) const
        {
            MPI_Int len, flag;
            MPIWRAP_CALL(MPI_Info_get_valuelen(info, nc(key.c_str()), &len, &flag));
            if (!flag) throw std::out_of_range(std::string("no such key: ") + key);
            return ConstInfoRef(info, key);
        }
};

}

#endif
