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

#include <utility>
#include <cstdio>
#include <cassert>
#include <algorithm>

#ifdef __MACH__
#include <mach/mach_time.h>
#endif

#include "util/util.h"
#include "task/task.hpp"

#include "time.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::task;

extern "C"
{

void c_do_flops(int64_t flops);

void c_do_flops(int64_t flops)
{
    aquarius::time::do_flops(flops);
}

}

namespace aquarius
{
namespace time
{

static std::vector<Interval> *tics[128];

Interval Interval::time()
{
    #ifdef __MACH__
    static double conv = -1.0;
    if (conv < 0)
    {
        mach_timebase_info_data_t timebase;
        mach_timebase_info(&timebase);
        conv = (double)timebase.numer / (double)timebase.denom;
    }
    uint64_t nsec = mach_absolute_time();
    return Interval(conv*(double)nsec/1e9, 0);
    #else
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return Interval((double)ts.tv_sec+(double)ts.tv_nsec/1e9, 0);
    #endif
}

Interval Interval::cputime()
{
    #ifdef __MACH__
    static double conv = -1.0;
    if (conv < 0)
    {
        mach_timebase_info_data_t timebase;
        mach_timebase_info(&timebase);
        conv = (double)timebase.numer / (double)timebase.denom;
    }
    uint64_t nsec = mach_absolute_time();
    return Interval(conv*(double)nsec/1e9, 0);
    #else
    timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return Interval((double)ts.tv_sec+(double)ts.tv_nsec/1e9, 0);
    #endif
}

bool Interval::operator<(const Interval& other) const
{
    return dt < other.dt;
}

Interval& Interval::operator+=(const Interval& other)
{
    dt += other.dt;
    flops += other.flops;
    return *this;
}

Interval& Interval::operator-=(const Interval& other)
{
    dt -= other.dt;
    flops -= other.flops;
    return *this;
}

Interval& Interval::operator*=(int m)
{
    dt *= m;
    return *this;
}

Interval& Interval::operator/=(int m)
{
    dt /= m;
    return *this;
}

Interval Interval::operator+(const Interval& other) const
{
    Interval r(*this);
    r += other;
    return r;
}

Interval Interval::operator-(const Interval& other) const
{
    Interval r(*this);
    r -= other;
    return r;
}

Interval Interval::operator*(int m) const
{
    Interval r(*this);
    r *= m;
    return r;
}

Interval Interval::operator/(int m) const
{
    Interval r(*this);
    r /= m;
    return r;
}

double Interval::seconds() const
{
    return dt;
}

double Interval::gflops() const
{
    return (double)flops/1e9/seconds();
}

double Interval::seconds(const Arena& arena) const
{
    double dtmax = seconds();
    arena.Allreduce(&dtmax, 1, MPI::MAX);
    return dtmax;
}

double Interval::gflops(const Arena& arena) const
{
    int64_t fl = flops;
    arena.Allreduce(&fl, 1, MPI::SUM);
    return (double)fl/1e9/seconds(arena);
}

std::vector<Timer> Timer::timers;

Timer& Timer::timer(const std::string& name)
{
    for (vector<Timer>::iterator it = timers.begin();it != timers.end();++it)
    {
        if (it->name == name) return *it;
    }

    timers.push_back(Timer(name));
    return timers.back();
}

void Timer::printTimers(const Arena& arena)
{
    int max_len = 0;
    for (vector<Timer>::iterator it = timers.begin();it != timers.end();++it)
    {
        max_len = max(max_len, (int)it->name.size());
    }

    sort(timers.begin(), timers.end());

    for (vector<Timer>::iterator it = timers.begin();it != timers.end();++it)
    {
        double tot = it->seconds(arena);
        int64_t count = it->count;
        arena.Allreduce(&count, 1, MPI::SUM);
        double gflops = it->gflops(arena);
        Logger::log(arena) << strprintf("%s:%*s %13.6f s %10ld x %11.6f gflops/sec\n", it->name.c_str(), (int)(max_len-it->name.size()), "", tot, count, gflops) << endl;
    }
}

void Timer::clearTimers(const Arena& arena)
{
    for (vector<Timer>::iterator it = timers.begin();it != timers.end();++it)
    {
        it->dt = 0;
        it->flops = 0;
        it->count = 0;
    }
}

void tic()
{
    #ifdef _OPENMP
    int tid = omp_get_thread_num();
    int ntd = (omp_in_parallel() ? 1 : omp_get_max_threads());
    #else
    int tid = 0;
    int ntd = 1;
    #endif
    for (int td = tid;td < tid+ntd;td++)
    {
        if (tics[td] == NULL) tics[td] = new vector<Interval>();
        tics[td]->push_back(td == tid ? Interval::time() : Interval());
    }
}

Interval toc()
{
    #ifdef _OPENMP
    int tid = omp_get_thread_num();
    int ntd = (omp_in_parallel() ? 1 : omp_get_max_threads());
    #else
    int tid = 0;
    int ntd = 1;
    #endif
    Interval dt = Interval::time();
    dt -= tics[tid]->back();
    tics[tid]->pop_back();
    for (int td = tid+1;td < tid+ntd;td++)
    {
        dt -= tics[td]->back();
        tics[td]->pop_back();
    }
    if (!tics[tid]->empty()) tics[tid]->back().flops -= dt.flops;
    return dt;
}

void cputic()
{
    #ifdef _OPENMP
    int tid = omp_get_thread_num();
    int ntd = (omp_in_parallel() ? 1 : omp_get_max_threads());
    #else
    int tid = 0;
    int ntd = 1;
    #endif
    for (int td = tid;td < tid+ntd;td++)
    {
        if (tics[td] == NULL) tics[td] = new vector<Interval>();
        tics[td]->push_back(td == tid ? Interval::cputime() : Interval());
    }
}

Interval cputoc()
{
    #ifdef _OPENMP
    int tid = omp_get_thread_num();
    int ntd = (omp_in_parallel() ? 1 : omp_get_max_threads());
    #else
    int tid = 0;
    int ntd = 1;
    #endif
    Interval dt = Interval::cputime();
    dt -= tics[tid]->back();
    tics[tid]->pop_back();
    for (int td = tid+1;td < tid+ntd;td++)
    {
        dt -= tics[td]->back();
        tics[td]->pop_back();
    }
    if (!tics[tid]->empty()) tics[tid]->back().flops -= dt.flops;
    return dt;
}

void do_flops(int64_t flops)
{
    #ifdef _OPENMP
    int tid = omp_get_thread_num();
    #else
    int tid = 0;
    #endif
    tics[tid]->back().flops -= flops;
}


}
}
