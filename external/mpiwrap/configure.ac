AC_INIT(mpiwrap, 1.0, dmatthews@utexas.edu, mpiwrap, http://www.github.com/devinamatthews/mpiwrap)

AC_CONFIG_MACRO_DIR([m4])

AX_PROG_CXX_MPI([true], [], [AC_MSG_FAILURE([No MPI C++ compiler found.])])
AC_LANG([C++])

AC_ARG_ENABLE([c++11],
              [AS_HELP_STRING([--disable-c++11],
                              [Disable the use of C++11 features.
                               These features are enable by default
                               if the compiler supports them.])])
                               
AS_IF([test "x$enable_c__11" != "xno"],
[
    AX_CXX_COMPILE_STDCXX_11([noext], [optional])
], [HAVE_CXX11=0])

AC_CHECK_TYPE([MPI_Count], [HAVE_MPI_COUNT=1], [HAVE_MPI_COUNT=0], [[#include "mpi.h"]])
AC_SUBST([HAVE_MPI_COUNT])

HAVE_MPI_WIN=1
AC_CHECK_FUNCS([MPI_Compare_and_swap \
                MPI_Fetch_and_op \
                MPI_Get_accumulate \
                MPI_Raccumulate \
                MPI_Rget \
                MPI_Rget_accumulate \
                MPI_Rput \
                MPI_Win_allocate \
                MPI_Win_allocate_shared \
                MPI_Win_attach \
                MPI_Win_create_dynamic \
                MPI_Win_detach \
                MPI_Win_flush \
                MPI_Win_flush_all \
                MPI_Win_flush_local \
                MPI_Win_flush_local_all \
                MPI_Win_get_info \
                MPI_Win_lock_all \
                MPI_Win_set_info \
                MPI_Win_shared_query \
                MPI_Win_sync \
                MPI_Win_unlock_all], [], [HAVE_MPI_WIN=0])
AC_SUBST([HAVE_MPI_WIN])

HAVE_MPI_ICOLLECTIVES=1
AC_CHECK_FUNCS([MPI_Iallgather \
                MPI_Iallgatherv \
                MPI_Iallreduce \
                MPI_Ialltoall \
                MPI_Ialltoallv \
                MPI_Ialltoallw \
                MPI_Ibarrier \
                MPI_Iexscan \
                MPI_Igather \
                MPI_Igatherv \
                MPI_Ireduce \
                MPI_Ireduce_scatter \
                MPI_Ireduce_scatter_block \
                MPI_Iscatter \
                MPI_Iscatterv], [], [HAVE_MPI_ICOLLECTIVES=0])
AC_SUBST([HAVE_MPI_ICOLLECTIVES])

HAVE_MPI_MESSAGE=1
AC_CHECK_FUNCS([MPI_Mprobe \
                MPI_Improbe \
                MPI_Mrecv \
                MPI_Imrecv], [], [HAVE_MPI_MESSAGE=0])
AC_SUBST([HAVE_MPI_MESSAGE])

HAVE_MPI_TYPE_CREATE_HINDEXED_BLOCK=1
AC_CHECK_FUNCS([MPI_Type_create_hindexed_block], [], [HAVE_MPI_TYPE_CREATE_HINDEXED_BLOCK=0])
AC_SUBST([HAVE_MPI_TYPE_CREATE_HINDEXED_BLOCK])

AC_CONFIG_FILES([
include/mpiwrap.hpp
include/internal/mpiwrap_buffer.hpp:include/internal/mpiwrap_buffer.hpp
include/internal/mpiwrap_comm.hpp:include/internal/mpiwrap_comm.hpp
include/internal/mpiwrap_common.hpp:include/internal/mpiwrap_common.hpp
include/internal/mpiwrap_datatype.hpp:include/internal/mpiwrap_datatype.hpp
include/internal/mpiwrap_group.hpp:include/internal/mpiwrap_group.hpp
include/internal/mpiwrap_info.hpp:include/internal/mpiwrap_info.hpp
include/internal/mpiwrap_intercomm.hpp:include/internal/mpiwrap_intercomm.hpp
include/internal/mpiwrap_intracomm.hpp:include/internal/mpiwrap_intracomm.hpp
include/internal/mpiwrap_message.hpp:include/internal/mpiwrap_message.hpp
include/internal/mpiwrap_request.hpp:include/internal/mpiwrap_request.hpp
include/internal/mpiwrap_status.hpp:include/internal/mpiwrap_status.hpp
include/internal/mpiwrap_window.hpp:include/internal/mpiwrap_window.hpp
])
AC_OUTPUT
