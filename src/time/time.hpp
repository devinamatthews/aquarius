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
#include <ctime>

typedef double time_type;

#ifdef PROFILE

#define PROFILE_FUNCTION \
{ \
static int __timer = time::register_timer(__func__); \
time::tic();

#define PROFILE_SECTION(name) \
{ \
static int __timer = time::register_timer(#name); \
time::tic();

#define PROFILE_STOP \
time_type __time = time::toc(); \
time::inc_timer(__timer, __time); \
}

#define PROFILE_RETURN \
{ \
time_type __time = time::toc(); \
time::inc_timer(__timer, __time); \
return; \
}

#else

#define PROFILE_FUNCTION

#define PROFILE_SECTION(name)

#define PROFILE_STOP

#define PROFILE_RETURN return

#endif

namespace aquarius
{
namespace time
{

double todouble(const time_type& t);

int register_timer(const std::string& name);

void inc_timer(const int timer, const time_type& time);

void print_timers();

void clear_timers();

void tic();

time_type toc();

}
}

#endif
