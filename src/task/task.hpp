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

class Product;

class Requirement
{
    protected:
        std::string type;
        std::string name;
        Product* product;

    public:
        Requirement(const std::string& type, const std::string& name);

        ~Requirement();

        const std::string& getName() const { return name; }

        const std::string getType() const { return type; }

        void fulfil(const Product& product);

        bool isFulfilled() const { return product != NULL; }

        Product& get();
};

class Product
{
    protected:
        std::string type;
        std::string name;
        std::shared_ptr<void*> data;
        std::vector<Requirement> requirements;

    public:
        Product(const std::string& type, const std::string& name);

        const std::string& getName() const { return name; }

        const std::string& getType() const { return type; }

        template <typename T> void put(T* resource);

        template <typename T> T& get();

        bool exists() const { return data; }

        void addRequirement(const Requirement& req);

        std::vector<Requirement> getRequirements() { return requirements; }
};

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

        static bool registerTask(const std::string& name, factory_func create);

        template <typename T> void put(const std::string& name, T* resource);

        template <typename T> T& get(const std::string& name);

        Product& getProduct(const std::string& name);

        std::vector<Product>& getProducts() { return products; }

        virtual void run() = 0;

        static Task* createTask(const std::string& type, const std::string& name, const input::Config& config);
};

template <class T>
class TaskFactory
{
    friend class Task;

    protected:
        static bool initialized;

        Task* create(const std::string& name, const input::Config& config)
        {
            return new T(name, config);
        }
};

#define REGISTER_TASK(task,type) \
template <> int TaskFactory<task>::initialized = \
    Task::registerTask(name,TaskFactory<task>::create);

}
}

#endif
