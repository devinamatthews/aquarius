#ifndef _AQUARIUS_TIME_HPP_
#define _AQUARIUS_TIME_HPP_

#include "util/global.hpp"

#ifdef PROFILE

#define PROFILE_SECTION(name) \
{ \
static time::Timer __timer = time::Timer::(name); \
__timer.start();

#define PROFILE_FUNCTION PROFILE_SECTION(__func__)

#define PROFILE_STOP \
__timer.stop(); \
}

#define PROFILE_RETURN \
__timer.stop(); \
return;

#define PROFILE_FLOPS(n) time::do_flops(n)

#else

#define PROFILE_FUNCTION

#define PROFILE_SECTION(name)

#define PROFILE_STOP

#define PROFILE_RETURN return

#define PROFILE_FLOPS(n) time::do_flops(n)

#endif

namespace aquarius
{
namespace time
{

class Timer;

class Interval
{
    friend class Timer;
    friend void do_flops(int64_t flops);
    friend Interval toc();
    friend Interval cputoc();

    protected:
        double dt;
        int64_t flops;

        Interval(double start, int64_t flops) : dt(start), flops(flops) {}

    public:
        Interval() : dt(0), flops(0) {}

        static Interval time();

        static Interval cputime();

        bool operator<(const Interval& other) const;

        Interval& operator+=(const Interval& other);

        Interval& operator-=(const Interval& other);

        Interval& operator*=(int m);

        Interval& operator/=(int m);

        Interval operator+(const Interval& other) const;

        Interval operator-(const Interval& other) const;

        Interval operator*(int m) const;

        Interval operator/(int m) const;

        double seconds() const;

        double gflops() const;

        double seconds(const Arena& arena) const;

        double gflops(const Arena& arena) const;
};

void tic();

Interval toc();

void cputic();

Interval cputoc();

void do_flops(int64_t flops);

class Timer : public Interval
{
    protected:
        static vector<Timer> timers;

        string name;
        int64_t count;
        CTF_Flop_Counter ctfflops;

        Timer(const string& name) : name(name), count(0) {}

    public:
        Timer() : name(), count(0) {}

        static Timer& timer(const string& name);

        void start()
        {
            tic();
            ctfflops.zero();
        }

        void stop()
        {
            #ifdef _OPENMP
            if (!omp_in_parallel())
            {
                *this += toc();//*omp_get_max_threads();
                count++;
                flops += ctfflops.count();
            }
            else
            #pragma omp critical
            #endif
            {
                *this += toc();
                count++;
                flops += ctfflops.count();
            }
        }

        static void printTimers(const Arena& arena);

        static void clearTimers(const Arena& arena);
};

}
}

#endif
