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

#include "util/util.h"

#include "time.hpp"

using namespace std;

namespace aquarius
{
namespace time
{

static timespec delta = {0,-1};
static vector< pair<timespec,timespec> > tics;
static vector< pair<string, pair<timespec,long> > > timers;

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

double todouble(const timespec& t)
{
	return double(t.tv_sec)+double(t.tv_nsec)/1e9;
}

int register_timer(const string& name)
{
    timers.push_back(make_pair(name, make_pair(timespec(), 0l)));
    return timers.size()-1;
}

void inc_timer(const int timer, const timespec& time)
{
    timers[timer].second.first += time;
    timers[timer].second.second++;
}

void print_timers()
{
    int max_len = 0;
    for (vector< pair<string, pair<timespec,long> > >::iterator it = timers.begin();it != timers.end();++it)
    {
        max_len = max(max_len, (int)it->first.size());
    }

    for (vector< pair<string, pair<timespec,long> > >::iterator it = timers.begin();it != timers.end();++it)
    {
        double tot = allsum(todouble(it->second.first));
        long count = allsum(it->second.second);
        PRINT("%s:%*s %13.6fs %10ld\n", it->first.c_str(), (int)(max_len-it->first.size()), "", tot, count);
    }
}

void tic()
{
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

    tics.push_back(make_pair(timespec(),timespec()));
    clock_gettime(CLOCK_REALTIME, &(tics.back().first));
}

timespec toc()
{
    clock_gettime(CLOCK_REALTIME, &(tics.back().second));
    //timespec dt = (tics.back().second-tics.back().first)-delta;
    timespec dt = tics.back().second-tics.back().first;
    tics.pop_back();
    return dt;
}

}
}
