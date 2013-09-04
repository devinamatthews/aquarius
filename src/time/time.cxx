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

#include <vector>
#include <utility>
#include <string>
#include <cstdio>
#include <cassert>
#include <algorithm>

#ifdef __MACH__
#include <mach/mach_time.h>
#endif

#ifdef _OPENMP
#include "omp.h"
#endif

#include "mpi.h"

#include "util/util.h"

#include "time.hpp"

using namespace std;

/*
bool operator<(const timespec& t1, const timespec& t2)
{
    return t1.tv_sec < t2.tv_sec || (!(t2.tv_sec < t1.tv_sec) && t1.tv_nsec < t2.tv_nsec);
}

timespec& operator+=(timespec& t1, const timespec& t2)
{
    t1.tv_sec += t2.tv_sec;
    t1.tv_nsec += t2.tv_nsec;
    if (t1.tv_nsec >= 1000000000)
    {
        t1.tv_nsec -= 1000000000;
        t1.tv_sec++;
    }
    return t1;
}

timespec& operator-=(timespec& t1, const timespec& t2)
{
    if (t1.tv_sec < t2.tv_sec || (t1.tv_sec == t2.tv_sec && t1.tv_nsec < t2.tv_nsec))
    {
        printf("%ld %ld %ld %ld\n", t1.tv_sec, t1.tv_nsec, t2.tv_sec, t2.tv_nsec);
        assert(0);
        t1.tv_sec = 0;
        t1.tv_nsec = 0;
        return t1;
    }

    t1.tv_sec -= t2.tv_sec;
    if (t1.tv_nsec < t2.tv_nsec)
    {
        t1.tv_nsec += 1000000000;
        t1.tv_sec--;
    }
    t1.tv_nsec -= t2.tv_nsec;
    return t1;
}

timespec& operator*=(timespec& t, const int m)
{
    t.tv_sec *= m;
    t.tv_nsec *= m;
    t.tv_sec += t.tv_nsec/1000000000;
    t.tv_nsec = t.tv_nsec%1000000000;
    return t;
}

timespec& operator/=(timespec& t, const int m)
{
    long tmp = t.tv_sec%m;
    t.tv_sec /= m;
    t.tv_nsec = (t.tv_nsec+tmp*1000000000)/m;
    t.tv_sec += t.tv_nsec/1000000000;
    t.tv_nsec = t.tv_nsec%1000000000;
    return t;
}

timespec operator+(const timespec& t1, const timespec& t2)
{
    timespec r = t1;
    r += t2;
    return r;
}

timespec operator-(const timespec& t1, const timespec& t2)
{
    timespec r = t1;
    r -= t2;
    return r;
}

timespec operator*(const timespec& t, const int m)
{
    timespec r = t;
    r *= m;
    return r;
}

timespec operator/(const timespec& t, const int m)
{
    timespec r = t;
    r /= m;
    return r;
}
*/

namespace aquarius
{
namespace time
{

//static timespec delta = {0,-1};
static vector< pair<time_type,time_type> > *tics[128];
static vector< pair<string, pair<time_type,long> > > timers;

double todouble(const time_type& t)
{
    return t;
}

/*
double todouble(const time_type& t)
{
#ifdef __MACH__
    static double conv = -1.0;
    if (conv < 0)
    {
        mach_timebase_info_data_t timebase;
        mach_timebase_info(&timebase);
        conv = (double)timebase.numer / (double)timebase.denom;
    }
    return conv*(double(t.tv_sec)+double(t.tv_nsec)/1e9);
#else
    return double(t.tv_sec)+double(t.tv_nsec)/1e9;
#endif
}
*/

int register_timer(const string& name)
{
    int timer = 0;
    #pragma omp critical
    {
        int i;
        for (i = 0;i < timers.size();i++)
        {
            if (timers[i].first == name) break;
        }
        timer = i;
        //if (timer == timers.size()) timers.push_back(make_pair(name, make_pair((timespec){0,0}, 0l)));
        if (timer == timers.size()) timers.push_back(make_pair(name, make_pair(0.0, 0l)));
    }
    return timer;
}

void inc_timer(const int timer, const time_type& time)
{
    #ifdef _OPENMP
    if (!omp_in_parallel())
    {
        int nt = omp_get_max_threads();
        timers[timer].second.first += time*nt;
        timers[timer].second.second++;
    }
    else
    #endif
    #pragma omp critical
    {
        timers[timer].second.first += time;
        timers[timer].second.second++;
    }
}

void print_timers()
{
    int max_len = 0;
    for (vector< pair<string, pair<time_type,long> > >::iterator it = timers.begin();it != timers.end();++it)
    {
        max_len = max(max_len, (int)it->first.size());
    }

    sort(timers.begin(), timers.end());

    for (vector< pair<string, pair<time_type,long> > >::iterator it = timers.begin();it != timers.end();++it)
    {
        double tot = todouble(it->second.first);
        long count = it->second.second;
        MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &tot, 1, MPI::DOUBLE, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &count, 1, MPI::LONG, MPI::SUM);
        if (MPI::COMM_WORLD.Get_rank() == 0)
            printf("%s:%*s %13.6f s %10ld x\n", it->first.c_str(), (int)(max_len-it->first.size()), "", tot, count);
    }
}

void clear_timers()
{
    for (vector< pair<string, pair<time_type,long> > >::iterator it = timers.begin();it != timers.end();++it)
    {
        it->second.first = 0;
        //it->second.first.tv_sec = 0;
        //it->second.first.tv_nsec = 0;
        it->second.second = 0;
    }
}

void tic()
{
    #ifdef _OPENMP
    int tid = omp_get_thread_num();
    #else
    int tid = 0;
    #endif

    /*
    if (delta.tv_nsec == -1)
    {
        delta.tv_nsec = 0;

        timespec dt = {0,0};
        for (int i = 0;i < 10;i++)
        {
            tic();
            dt += toc();
        }

        delta = dt/10;

        PRINT("timer delta: %es\n", todouble(delta));
    }
    */

    if (tics[tid] == NULL) tics[tid] = new vector< pair<time_type,time_type> >();

    //PRINT("tic\n");
    tics[tid]->push_back(make_pair(time_type(),time_type()));
    tics[tid]->back().first = MPI::Wtime();
/*
#ifdef __MACH__
    uint64_t nsec = mach_absolute_time();
    tics[tid]->back().first.tv_sec = nsec/1000000000;
    tics[tid]->back().first.tv_nsec = nsec%1000000000;
#else
    int ret = clock_gettime(CLOCK_MONOTONIC, &(tics[tid]->back().first));
    assert(ret == 0);
#endif
*/
    //ALLPRINT("tic: %ld %ld\n", tics[tid]->back().first.tv_sec, tics[tid]->back().first.tv_nsec);
}

time_type toc()
{
    #ifdef _OPENMP
    int tid = omp_get_thread_num();
    #else
    int tid = 0;
    #endif

    //PRINT("toc\n");
    tics[tid]->back().second = MPI::Wtime();
/*
#ifdef __MACH__
    uint64_t nsec = mach_absolute_time();
    tics[tid]->back().second.tv_sec = nsec/1000000000;
    tics[tid]->back().second.tv_nsec = nsec%1000000000;
#else
    int ret = clock_gettime(CLOCK_MONOTONIC, &(tics[tid]->back().second));
    assert(ret == 0);
#endif
*/
    //ALLPRINT("toc: %ld %ld\n", tics[tid]->back().second.tv_sec, tics[tid]->back().second.tv_nsec);
    //timespec dt = (tics->back().second-tics[tid]->back().first)-delta;
    time_type dt = tics[tid]->back().second-tics[tid]->back().first;
    tics[tid]->pop_back();
    return dt;
}

}
}
