#ifndef _AQUARIUS_FRAMEWORKS_LOGGING_LOGGING_HPP_
#define _AQUARIUS_FRAMEWORKS_LOGGING_LOGGING_HPP_

#include "frameworks/util.hpp"

namespace aquarius
{
namespace logging
{

class Logger
{
    protected:
        class LogToStreamBuffer : public streambuf
        {
            protected:
                ostream& os;

                int sync();

                streamsize xsputn (const char* s, streamsize n);

                int overflow (int c = EOF);

            public:
                LogToStreamBuffer(ostream& os);
        };

        class SelfDestructBuffer : public LogToStreamBuffer
        {
            protected:
                int sync();

            public:
                SelfDestructBuffer();
        };

        class NullBuffer : public streambuf
        {
            protected:
                streamsize xsputn (const char* s, streamsize n);

                int overflow (int c = EOF);
        };

        class NullStream : public ostream
        {
            public:
                NullStream();

                ~NullStream();
        };

        class SelfDestructStream : public ostream
        {
            public:
                SelfDestructStream();
        };

        static string dateTime();

        static NullStream nullstream;
        static SelfDestructStream sddstream;

    public:
        static ostream& log(const Arena& arena);

        static ostream& warn(const Arena& arena);

        static ostream& error(const Arena& arena);
};

class Printer
{
    protected:
        ostream* out;

    public:
        Printer() : out(NULL) {}

        Printer(ostream& out) : out(&out) {}

        ostream& getOutputStream() { return *out; }

        void setOutputStream(ostream& os) { out = &os; }

        template <class T>
        Printer& operator<<(const T& t)
        {
            if (out != NULL) out << t;
            return *this;
        }
};

}
}

#endif
