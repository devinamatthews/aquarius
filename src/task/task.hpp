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

#include "input/config.hpp"
#include "stl_ext/stl_ext.hpp"
#include "util/util.h"
#include "util/distributed.hpp"

namespace aquarius
{
namespace task
{

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

class Resource : public Distributed
{
    public:
        Resource(const Arena& arena)
        : Distributed(arena) {}

        virtual ~Resource() {}

        virtual void print(Printer& p) const = 0;
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

class Product
{
    friend class Requirement;

    protected:
        std::string type;
        std::string name;
        std::global_ptr<Resource> data;
        std::vector<Requirement> requirements;
        std::shared_ptr<bool> used;

    public:
        Product(const std::string& type, const std::string& name);

        Product(const std::string& type, const std::string& name, const std::vector<Requirement>& reqs);

        const std::string& getName() const { return name; }

        const std::string& getType() const { return type; }

        bool isUsed() const { return *used; }

        template <typename T> void put(T* resource)
        {
            data.reset(static_cast<Resource*>(resource));
        }

        template <typename T> T& get()
        {
            if (!data) throw std::logic_error("Product " + name + " does not exist.");
            return static_cast<T&>(*data);
        }

        bool exists() const { return data; }

        void addRequirement(const Requirement& req);

        void addRequirements(const std::vector<Requirement>& reqs);

        std::vector<Requirement>& getRequirements() { return requirements; }
};

class TaskDAG;

class Task
{
    protected:
        typedef Task* (*factory_func)(const std::string&, const input::Config&);

        std::string type;
        std::string name;
        std::vector<Product> products;

        static std::map<std::string,factory_func>& tasks();

        void addProduct(const Product& product);

    public:
        Task(const std::string& type, const std::string& name);

        virtual ~Task() {}

        const std::string& getType() const { return type; }

        const std::string& getName() const { return name; }

        static bool registerTask(const std::string& name, factory_func create);

        template <typename T> void put(const std::string& name, T* resource)
        {
            getProduct(name).put(resource);
        }

        template <typename T> T& get(const std::string& name)
        {
            std::vector<Product*> to_search;

            for (std::vector<Product>::iterator i = products.begin();i != products.end();++i) to_search += &(*i);

            while (!to_search.empty())
            {
                std::vector<Product*> new_to_search;

                for (std::vector<Product*>::iterator i = to_search.begin();i != to_search.end();i++)
                {
                    if ((*i)->getName() == name) return (*i)->get<T>();

                    for (std::vector<Requirement>::iterator j = (*i)->getRequirements().begin();
                         j != (*i)->getRequirements().end();j++) new_to_search += &j->get();
                }

                to_search = new_to_search;
            }

            throw std::logic_error("Product " + name + " not found on task " + this->name + " or its dependencies");
        }

        Product& getProduct(const std::string& name);

        std::vector<Product>& getProducts() { return products; }

        virtual void run(TaskDAG& dag, Arena& arena) = 0;

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
        std::vector<Task*> tasks;

    public:
        TaskDAG() {}

        ~TaskDAG();

        void addTask(Task* task);

        void execute(Arena& world);
};

}
}

#endif
