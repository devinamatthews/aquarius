/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_TIME_HPP_
#define _AQUARIUS_TIME_HPP_

#include <string>
#include <stdint.h>
#include <ctime>
#include <vector>

#include "util/stl_ext.hpp"
#include "util/distributed.hpp"

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
        static std::vector<Timer> timers;

        std::string name;
        int64_t count;
        CTF_Flop_Counter ctfflops;

        Timer(const std::string& name) : name(name), count(0) {}

    public:
        Timer() : name(), count(0) {}

        static Timer& timer(const std::string& name);

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
