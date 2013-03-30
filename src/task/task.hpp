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

#ifndef _AQUARIUS_TASK_HPP_
#define _AQUARIUS_TASK_HPP_

#include <map>
#include <utility>
#include <string>

#include "input/config.hpp"
#include "stl_ext/stl_ext.hpp"

namespace aquarius
{
namespace task
{

class Resource
{
    public:
        virtual ~Resource() {};

        virtual bool operator==(const Resource& other) = 0;

        virtual Task* defaultProvider() = 0;
};

class Task
{
    protected:
        std::map<std::string,Resource*> requirements;
        std::map<std::string,Resource*> products;
        bool hasRun_;
        bool wasSuccessful_;

    public:
        Task() : hasRun_(false), wasSuccessful_(false) {}

        virtual ~Task()
        {
            for (std::map<std::string,Resource*>::iterator i = requirements.begin();i != requirements.end();++i)
            {
                delete i->second;
            }

            for (std::map<std::string,Resource*>::iterator i = products.begin();i != products.end();++i)
            {
                delete i->second;
            }
        }

        virtual void configure(const input::Config& config) = 0;

        virtual void run() = 0;

        bool hasRun() const { return hasRun_; }

        bool wasSuccessful() const { return wasSuccessful_; }

        const std::map<std::string,Resource*>& requires() const { return requirements; }

        const std::map<std::string,Resource*>& provides() const { return products; }
};

}
}

#endif
