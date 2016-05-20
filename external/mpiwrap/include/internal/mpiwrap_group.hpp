#ifndef _MPIWRAP_GROUP_HPP_
#define _MPIWRAP_GROUP_HPP_

#include "mpiwrap_common.hpp"

namespace MPIWrap
{

namespace internal
{
    template <typename Derived> class Comm;
}

class Group
{
    friend class Intercomm;
    friend class Intracomm;
    template <typename Derived> friend class internal::Comm;

    protected:
        MPI_Group group;

        static int getRank(const MPI_Group& group)
        {
            MPI_Int rank;
            MPIWRAP_CALL(MPI_Group_rank(group, &rank));
            return rank;
        }

        static int getSize(const MPI_Group& group)
        {
            MPI_Int size;
            MPIWRAP_CALL(MPI_Group_size(group, &size));
            return size;
        }

        explicit Group(const MPI_Group& group)
        : group(group), rank(getRank(group)), size(getSize(group)) {}

    public:
        struct Range
        {
            MPI_Int from, to, stride;
            Range(MPI_Int from, MPI_Int to, MPI_Int stride)
            : from(from), to(to), stride(stride) {}
        };

        const MPI_Int rank;
        const MPI_Int size;

#if MPIWRAP_CXX11

        Group(Group&& other) : group(other.group), rank(other.rank), size(other.size)
        {
            other.group = MPI_GROUP_NULL;
        }

#endif

        ~Group()
        {
            if (group != MPI_GROUP_NULL)
            {
                MPIWRAP_CALL(MPI_Group_free(&group));
            }
        }

        operator MPI_Group&() { return group; }

        operator const MPI_Group&() const { return group; }

        operator MPI_Group*() { return &group; }

        operator const MPI_Group*() const { return &group; }

        std::vector<MPI_Int> translateRanks(const MPI_Int* ranks, MPI_Int count, const Group& other) const
        {
            std::vector<MPI_Int> other_ranks(count);
            MPIWRAP_CALL(MPI_Group_translate_ranks(group, count, nc(ranks), other, &other_ranks.front()));
            return other_ranks;
        }

        MPI_Int compare(const Group& other) const
        {
            MPI_Int result;
            MPIWRAP_CALL(MPI_Group_compare(group, other, &result));
            return result;
        }

        Group together(const Group& other) const
        {
            MPI_Group g;
            MPIWRAP_CALL(MPI_Group_union(group, other, &g));
            return Group(g);
        }

        Group intersection(const Group& other) const
        {
            MPI_Group g;
            MPIWRAP_CALL(MPI_Group_intersection(group, other, &g));
            return Group(g);
        }

        Group difference(const Group& other) const
        {
            MPI_Group g;
            MPIWRAP_CALL(MPI_Group_difference(group, other, &g));
            return Group(g);
        }

        Group including(const MPI_Int* ranks, MPI_Int count)
        {
            MPI_Group g;
            MPIWRAP_CALL(MPI_Group_incl(group, count, nc(ranks), &g));
            return Group(g);
        }

        Group including(const std::vector<MPI_Int>& ranks)
        {
            return including(&ranks.front(), ranks.size());
        }

        Group excluding(const MPI_Int* ranks, MPI_Int count)
        {
            MPI_Group g;
            MPIWRAP_CALL(MPI_Group_excl(group, count, nc(ranks), &g));
            return Group(g);
        }

        Group excluding(const std::vector<MPI_Int>& ranks)
        {
            return excluding(&ranks.front(), ranks.size());
        }

        Group including(const Range* ranges, MPI_Int count)
        {
            MPI_Group g;
            MPIWRAP_CALL(MPI_Group_incl(group, count, reinterpret_cast<MPI_Int*>(nc(ranges)), &g));
            return Group(g);
        }

        Group including(const std::vector<Range>& ranges)
        {
            return including(&ranges.front(), ranges.size());
        }

        Group excluding(const Range* ranges, MPI_Int count)
        {
            MPI_Group g;
            MPIWRAP_CALL(MPI_Group_excl(group, count, reinterpret_cast<MPI_Int*>(nc(ranges)), &g));
            return Group(g);
        }

        Group excluding(const std::vector<Range>& ranges)
        {
            return excluding(&ranges.front(), ranges.size());
        }
};

}

#endif
