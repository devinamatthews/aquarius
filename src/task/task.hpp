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
#include <stdexcept>
#include <ostream>
#include <streambuf>
#include <ios>
#include <iostream>
#include <ctime>
#include <vector>
#include <stdexcept>
#include <unistd.h>
#include <iomanip>

#include "input/config.hpp"
#include "util/stl_ext.hpp"
#include "util/util.h"
#include "util/distributed.hpp"
#include "time/time.hpp"

namespace aquarius
{
namespace task
{

class Logger
{
    protected:
        class LogToStreamBuffer : public std::streambuf
        {
            protected:
                std::ostream& os;

                int sync();

                std::streamsize xsputn (const char* s, std::streamsize n);

                int overflow (int c = EOF);

            public:
                LogToStreamBuffer(std::ostream& os);
        };

        class SelfDestructBuffer : public LogToStreamBuffer
        {
            protected:
                int sync();

            public:
                SelfDestructBuffer();
        };

        class NullBuffer : public std::streambuf
        {
            protected:
                std::streamsize xsputn (const char* s, std::streamsize n);

                int overflow (int c = EOF);
        };

        class NullStream : public std::ostream
        {
            public:
                NullStream();

                ~NullStream();
        };

        class SelfDestructStream : public std::ostream
        {
            public:
                SelfDestructStream();
        };

        static std::string dateTime();

        static NullStream nullstream;
        static SelfDestructStream sddstream;

    public:
        static std::ostream& log(const Arena& arena);

        static std::ostream& warn(const Arena& arena);

        static std::ostream& error(const Arena& arena);
};

class Printer
{
    protected:
        std::ostream* out;

    public:
        Printer() : out(NULL) {}

        Printer(std::ostream& out) : out(&out) {}

        std::ostream& getOutputStream() { return *out; }

        void setOutputStream(std::ostream& os) { out = &os; }

        template <class T>
        Printer& operator<<(const T& t)
        {
            if (out != NULL) out << t;
            return *this;
        }
};

class Destructible
{
    public:
        virtual ~Destructible() {}
};

template <typename T>
class Resource : public Destructible
{
    private:
        Resource(const Resource& other);

        Resource& operator=(const Resource& other);

    public:
        T* const data;

        Resource(T* data) : data(data) {}

        ~Resource() { delete data; }
};

class Product;

class Requirement
{
    protected:
        std::string type;
        std::string name;
        std::global_ptr<Product> product;

    public:
        Requirement(const std::string& type, const std::string& name);

        const std::string& getName() const { return name; }

        const std::string getType() const { return type; }

        bool exists() const;

        void fulfil(const Product& product);

        bool isFulfilled() const { return product; }

        Product& get();
};

class Task;

class Product
{
    friend class Requirement;
    friend class Task;

    protected:
        std::string type;
        std::string name;
        std::global_ptr<Destructible> data;
        std::shared_ptr<std::vector<Requirement> > requirements;
        std::shared_ptr<bool> used;

    public:
        Product(const std::string& type, const std::string& name);

        Product(const std::string& type, const std::string& name, const std::vector<Requirement>& reqs);

        const std::string& getName() const { return name; }

        const std::string& getType() const { return type; }

        bool isUsed() const { return *used; }

        template <typename T> T& put(T* resource)
        {
            data.reset(new Resource<T>(resource));
            return *resource;
        }

        template <typename T> T& get()
        {
            if (!data) throw std::logic_error("Product " + name + " does not exist.");
            return *static_cast<Resource<T>&>(*data).data;
        }

        bool exists() const { return data; }

        void addRequirement(const Requirement& req);

        void addRequirements(const std::vector<Requirement>& reqs);

        std::vector<Requirement>& getRequirements() { return *requirements; }
};

class TaskDAG;

class Task
{
    protected:
        typedef Task* (*factory_func)(const std::string&, const input::Config&);

        std::string type;
        std::string name;
        std::vector<Product> products;
        std::vector<Product> temporaries;

        static std::map<std::string,factory_func>& tasks();

        void addProduct(Product&& product)
        {
            products.push_back(std::forward<Product>(product));
        }

        template <typename... Args>
        void addProduct(Args&&... args)
        {
            products.emplace_back(std::forward<Args>(args)...);
        }

        template <typename T> T& puttmp(const std::string& name, T* resource)
        {
            for (std::vector<Product>::iterator i = temporaries.begin();i != temporaries.end();++i)
            {
                if (i->getName() == name)
                {
                    return i->put(resource);
                }
            }

            temporaries.push_back(Product("temporary", name));
            return temporaries.back().put(resource);
        }

        template <typename T> T& gettmp(const std::string& name)
        {
            for (std::vector<Product>::iterator i = temporaries.begin();i != temporaries.end();++i)
            {
                if (i->getName() == name)
                {
                    return i->get<T>();
                }
            }

            throw std::logic_error("Temporary " + name + " not found on task " + this->name);
        }

        std::ostream& log(const Arena& arena);

        std::ostream& warn(const Arena& arena);

        std::ostream& error(const Arena& arena);

    public:
        Task(const std::string& type, const std::string& name);

        virtual ~Task()
        {
            for (std::vector<Product>::iterator i = products.begin();i != products.end();++i)
            {
                i->requirements->clear();
            }
        }

        const std::string& getType() const { return type; }

        const std::string& getName() const { return name; }

        bool isUsed(const std::string& name) const { return getProduct(name).isUsed(); }

        template <typename T> static std::string type_string();

        static bool registerTask(const std::string& name, factory_func create);

        static const input::Schema& getSchema(const std::string& name);

        template <typename T> T& put(const std::string& name, T* resource)
        {
            return getProduct(name).put(resource);
        }

        template <typename T> T& get(const std::string& name)
        {
            std::vector<Product*> to_search;

            for (std::vector<Product>::iterator i = products.begin();i != products.end();++i)
            {
                if (i->getName() == name) return i->get<T>();

                for (std::vector<Requirement>::iterator j = i->getRequirements().begin();
                     j != i->getRequirements().end();j++)
                {
                    if (j->getName() == name) return j->get().get<T>();
                }
            }

            throw std::logic_error("Product " + name + " not found on task " + this->name + " or its dependencies");
        }

        Product& getProduct(const std::string& name);

        const Product& getProduct (const std::string& name) const;

        std::vector<Product>& getProducts() { return products; }

        const std::vector<Product>& getProducts() const { return products; }

        virtual void run(TaskDAG& dag, const Arena& arena) = 0;

        static Task* createTask(const std::string& type, const std::string& name, const input::Config& config);
};

template <class T>
class TaskFactory
{
    friend class Task;

    protected:
        static bool initialized;

        static Task* create(const std::string& name, const input::Config& config)
        {
            return new T(name, config);
        }
};

#define REGISTER_TASK(task,name) \
template <> bool TaskFactory<task >::initialized = \
    Task::registerTask(name,TaskFactory<task >::create);

class TaskDAG
{
    NON_COPYABLE(TaskDAG);
    NON_ASSIGNABLE(TaskDAG);

    protected:
        std::vector<std::pair<Task*,input::Config> > tasks;

        void parseTasks(const std::string& context, input::Config& config);

        void satisfyRemainingRequirements(const Arena& world);

        void satisfyExplicitRequirements(const Arena& world);

    public:
        TaskDAG() {}

        TaskDAG(const std::string& file);

        ~TaskDAG();

        void addTask(Task* task, const input::Config& config);

        void execute(Arena& world);
};

class CompareScalars : public Task
{
    protected:
        double tolerance;

    public:
        CompareScalars(const std::string& name, const input::Config& config);

        void run(TaskDAG& dag, const Arena& arena);
};

}
}

#endif
