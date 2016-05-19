#include "logging.hpp"

namespace aquarius
{
namespace logging
{

Logger::NullStream Logger::nullstream;
Logger::SelfDestructStream Logger::sddstream;

int Logger::LogToStreamBuffer::sync()
{
    if (os.flush()) return 0;
    else return -1;
}

streamsize Logger::LogToStreamBuffer::xsputn (const char* s, streamsize n)
{
    if (os.write(s, n)) return n;
    else return 0;
}

int Logger::LogToStreamBuffer::overflow (int c)
{
    if (c != EOF) os << traits_type::to_char_type(c);
    if (c == EOF || !os) return EOF;
    else return traits_type::to_int_type(c);
}

Logger::LogToStreamBuffer::LogToStreamBuffer(ostream& os)
: os(os) {}

int Logger::SelfDestructBuffer::sync()
{
    LogToStreamBuffer::sync();
    abort();
    return -1;
}

Logger::SelfDestructBuffer::SelfDestructBuffer()
: LogToStreamBuffer(cerr) {}

streamsize Logger::NullBuffer::xsputn (const char* s, streamsize n)
{
    return n;
}

int Logger::NullBuffer::overflow (int c)
{
    if (c == EOF) return EOF;
    else return traits_type::to_int_type(c);
}

Logger::NullStream::NullStream()
: ostream(new NullBuffer()) {}

Logger::NullStream::~NullStream()
{
    delete rdbuf();
}

Logger::SelfDestructStream::SelfDestructStream()
: ostream(new SelfDestructBuffer()) {}

string Logger::dateTime()
{
    char buf[256];
    time_t t = ::time(NULL);
    tm* timeptr = localtime(&t);
    strftime(buf, 256, "%c", timeptr);
    return string(buf);
}

ostream& Logger::log(const Arena& arena)
{
    if (arena.rank == 0)
    {
        cout << dateTime() << ": ";
        return cout;
    }
    else
    {
        return nullstream;
    }
}

ostream& Logger::warn(const Arena& arena)
{
    if (arena.rank == 0)
    {
        cerr << dateTime() << ": warning: ";
        return cerr;
    }
    else
    {
        return nullstream;
    }
}

ostream& Logger::error(const Arena& arena)
{
    arena.comm().Barrier();
    if (arena.rank == 0)
    {
        sddstream << dateTime() << ": error: ";
        return sddstream;
        //cout << dateTime() << ": error: ";
        //return cout;
    }
    else
    {
        return nullstream;
    }
}

}
}
