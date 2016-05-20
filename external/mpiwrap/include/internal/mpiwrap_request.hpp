#ifndef _MPIWRAP_REQUEST_HPP_
#define _MPIWRAP_REQUEST_HPP_

#include "mpiwrap_common.hpp"
#include "mpiwrap_status.hpp"

namespace MPIWrap
{

namespace internal
{
    template <typename Derived> class Comm;
}

class Request
{
    friend class Intracomm;
    friend class Intercomm;
    friend class Message;
    template <typename Derived> friend class internal::Comm;

    protected:
        MPI_Request req;

        explicit Request(const MPI_Request& req) : req(req) {}

    public:
#if MPIWRAP_CXX11

        Request(Request&& other) : req(other.req)
        {
            other.req = MPI_REQUEST_NULL;
        }

#endif

        ~Request()
        {
            if (req != MPI_REQUEST_NULL)
            {
                MPIWRAP_CALL(MPI_Request_free(&req));
            }
        }

        operator MPI_Request&() { return req; }

        operator const MPI_Request&() const { return req; }

        operator MPI_Request*() { return &req; }

        operator const MPI_Request*() const { return &req; }

        void start()
        {
            MPIWRAP_CALL(MPI_Start(&req));
        }

        void cancel()
        {
            MPIWRAP_CALL(MPI_Cancel(&req));
        }

        void wait()
        {
#if MPIWRAP_VERSION_AT_LEAST(2,0)
            MPIWRAP_CALL(MPI_Wait(&req, MPI_STATUS_IGNORE));
#else
            Status status;
            wait(status);
#endif
        }

        void wait(Status& status)
        {
            MPIWRAP_CALL(MPI_Wait(&req, status));
        }

        bool test()
        {
#if MPIWRAP_VERSION_AT_LEAST(2,0)
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Test(&req, &flag, MPI_STATUS_IGNORE));
            return flag;
#else
            Status status;
            return test(status);
#endif
        }

        bool test(Status& status)
        {
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Test(&req, &flag, status));
            return flag;
        }

        friend Request& waitAny(Request* reqs, MPI_Int count)
        {
#if MPIWRAP_VERSION_AT_LEAST(2,0)
            MPI_Int i;
            MPIWRAP_CALL(MPI_Waitany(count, reinterpret_cast<MPI_Request*>(reqs), &i, MPI_STATUS_IGNORE));
            return reqs[i];
#else
            Status status;
            return waitAny(reqs, count, status);
#endif
        }

        friend Request& waitAny(std::vector<Request>& reqs)
        {
            return waitAny(&reqs.front(), reqs.size());
        }

        friend Request& waitAny(Request* reqs, MPI_Int count, Status& status)
        {
            MPI_Int i;
            MPIWRAP_CALL(MPI_Waitany(count, reinterpret_cast<MPI_Request*>(reqs), &i, status));
            return reqs[i];
        }

        friend Request& waitAny(std::vector<Request>& reqs, Status& status)
        {
            return waitAny(&reqs.front(), reqs.size(), status);
        }

        friend Request* testAny(Request* reqs, MPI_Int count)
        {
#if MPIWRAP_VERSION_AT_LEAST(2,0)
            MPI_Int i, flag;
            MPIWRAP_CALL(MPI_Testany(count, reinterpret_cast<MPI_Request*>(reqs), &i, &flag, MPI_STATUS_IGNORE));
            return (flag ? &reqs[i] : NULL);
#else
            Status status;
            return testAny(reqs, count, status);
#endif
        }

        friend Request* testAny(std::vector<Request>& reqs)
        {
            return testAny(&reqs.front(), reqs.size());
        }

        friend Request* testAny(Request* reqs, MPI_Int count, Status& status)
        {
            MPI_Int i, flag;
            MPIWRAP_CALL(MPI_Testany(count, reinterpret_cast<MPI_Request*>(reqs), &i, &flag, status));
            return (flag ? &reqs[i] : NULL);
        }

        friend Request* testAny(std::vector<Request>& reqs, Status& status)
        {
            return testAny(&reqs.front(), reqs.size(), status);
        }

        friend void waitAll(Request* reqs, MPI_Int count)
        {
#if MPIWRAP_VERSION_AT_LEAST(2,0)
            MPIWRAP_CALL(MPI_Waitall(count, reinterpret_cast<MPI_Request*>(reqs), MPI_STATUSES_IGNORE));
#else
            std::vector<Status> stats(count);
            waitAll(reqs, &stats.front(), count);
#endif
        }

        friend void waitAll(std::vector<Request>& reqs)
        {
            waitAll(&reqs.front(), reqs.size());
        }

        friend void waitAll(Request* reqs, Status* stats, MPI_Int count)
        {
            MPIWRAP_CALL(MPI_Waitall(count, reinterpret_cast<MPI_Request*>(reqs), reinterpret_cast<MPI_Status*>(stats)));
        }

        friend void waitAll(std::vector<Request>& reqs, std::vector<Status>& stats)
        {
            MPIWRAP_ASSERT(reqs.size() == stats.size(),
                           "There must be exactly one status for each request.");
            waitAll(&reqs.front(), &stats.front(), reqs.size());
        }

        friend bool testAll(Request* reqs, MPI_Int count)
        {
#if MPIWRAP_VERSION_AT_LEAST(2,0)
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Testall(count, reinterpret_cast<MPI_Request*>(reqs), &flag, MPI_STATUSES_IGNORE));
            return flag;
#else
            std::vector<Status> stats(count);
            return testAll(reqs, &stats.front(), count);
#endif
        }

        friend bool testAll(std::vector<Request>& reqs)
        {
            return testAll(&reqs.front(), reqs.size());
        }

        friend bool testAll(Request* reqs, Status* stats, MPI_Int count)
        {
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Testall(count, reinterpret_cast<MPI_Request*>(reqs), &flag, reinterpret_cast<MPI_Status*>(stats)));
            return flag;
        }

        friend bool testAll(std::vector<Request>& reqs, std::vector<Status>& stats)
        {
            MPIWRAP_ASSERT(reqs.size() == stats.size(),
                           "There must be exactly one status for each request.");
            return testAll(&reqs.front(), &stats.front(), reqs.size());
        }

        friend std::vector<Request*> waitSome(Request* reqs, MPI_Int count)
        {
#if MPIWRAP_VERSION_AT_LEAST(2,0)
            std::vector<MPI_Int> i(count);
            MPI_Int n;
            MPIWRAP_CALL(MPI_Waitsome(count, reinterpret_cast<MPI_Request*>(reqs), &n, &i.front(), MPI_STATUSES_IGNORE))
            std::vector<Request*> ret(n);
            for (MPI_Int j = 0;j < n;j++)
            {
                ret[j] = reqs+i[j];
            }
            return ret;
#else
            std::vector<Status> stats(count);
            return waitSome(reqs, &stats.front(), count);
#endif
        }

        friend std::vector<Request*> waitSome(std::vector<Request>& reqs)
        {
            return waitSome(&reqs.front(), reqs.size());
        }

        friend std::vector<Request*> waitSome(Request* reqs, Status* stats, MPI_Int count)
        {
            std::vector<MPI_Int> i(count);
            MPI_Int n;
            MPIWRAP_CALL(MPI_Waitsome(count, reinterpret_cast<MPI_Request*>(reqs), &n, &i.front(), reinterpret_cast<MPI_Status*>(stats)))
            std::vector<Request*> ret(n);
            for (MPI_Int j = 0;j < n;j++)
            {
                ret[j] = reqs+i[j];
            }
            return ret;
        }

        friend std::vector<Request*> waitSome(std::vector<Request>& reqs, std::vector<Status>& stats)
        {
            MPIWRAP_ASSERT(reqs.size() == stats.size(),
                           "There must be exactly one status for each request.");
            return waitSome(&reqs.front(), &stats.front(), reqs.size());
        }

        friend std::vector<Request*> testSome(Request* reqs, MPI_Int count)
        {
#if MPIWRAP_VERSION_AT_LEAST(2,0)
            std::vector<MPI_Int> i(count);
            MPI_Int n;
            MPIWRAP_CALL(MPI_Testsome(count, reinterpret_cast<MPI_Request*>(reqs), &n, &i.front(), MPI_STATUSES_IGNORE))
            std::vector<Request*> ret(n);
            for (MPI_Int j = 0;j < n;j++)
            {
                ret[j] = reqs+i[j];
            }
            return ret;
#else
            std::vector<Status> stats(count);
            return testSome(reqs, &stats.front(), count);
#endif
        }

        friend std::vector<Request*> testSome(std::vector<Request>& reqs)
        {
            return testSome(&reqs.front(), reqs.size());
        }

        friend std::vector<Request*> testSome(Request* reqs, Status* stats, MPI_Int count)
        {
            std::vector<MPI_Int> i(count);
            MPI_Int n;
            MPIWRAP_CALL(MPI_Testsome(count, reinterpret_cast<MPI_Request*>(reqs), &n, &i.front(), reinterpret_cast<MPI_Status*>(stats)))
            std::vector<Request*> ret(n);
            for (MPI_Int j = 0;j < n;j++)
            {
                ret[j] = reqs+i[j];
            }
            return ret;
        }

        friend std::vector<Request*> testSome(std::vector<Request>& reqs, std::vector<Status>& stats)
        {
            MPIWRAP_ASSERT(reqs.size() == stats.size(),
                           "There must be exactly one status for each request.");
            return testSome(&reqs.front(), &stats.front(), reqs.size());
        }

#if MPIWRAP_VERSION_AT_LEAST(2,0)
        bool getStatus() const
        {
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Request_get_status(*this, &flag, MPI_STATUS_IGNORE));
            return flag;
        }

        bool getStatus(Status& status) const
        {
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Request_get_status(*this, &flag, status));
            return flag;
        }
#endif
};

}

#endif
