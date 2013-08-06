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

#include "task.hpp"

#include <set>

using namespace std;
using namespace aquarius;
using namespace aquarius::task;

Requirement::Requirement(const string& type, const string& name)
: type(type), name(name) {}

Requirement::~Requirement()
{
    if (product != NULL) delete product;
}

void Requirement::fulfil(const Product& product)
{
    if (this->product != NULL) delete this->product;
    this->product = new Product(product);
}

Product& Requirement::get()
{
    if (product == NULL) throw logic_error("Requirement " + name + " not fulfilled");
    return *product;
}

Product::Product(const string& type, const string& name)
: type(type), name(name), data(new void*(NULL)) {}

template <typename T>
void Product::put(T* resource)
{
    *data = static_cast<T*>(resource);
}

template <typename T>
T& Product::get()
{
    if (!*data) throw logic_error("Product " + name + " does not exist.");
    return *static_cast<T*>(*data);
}

Task::Task(const string& type, const string& name)
: type(type), name(name) {}

map<string,Task::factory_func>& Task::tasks()
{
    static std::map<std::string,factory_func> tasks_;
    return tasks_;
}

void Task::addProduct(const Product& product)
{
    products.push_back(product);
}

bool Task::registerTask(const string& name, factory_func create)
{
    tasks()[name] = create;
    return true;
}

template <typename T>
void Task::put(const string& name, T* resource)
{
    getProduct(name).put(resource);
}

template <typename T>
T& Task::get(const string& name)
{
    vector<Product*> to_search;

    for (vector<Product>::iterator i = products.begin();i != products.end();i++) to_search += &(*i);

    while (!to_search.empty())
    {
        vector<Product*> new_to_search;

        for (vector<Product*>::iterator i = to_search.begin();i != to_search.end();i++)
        {
            if ((*i)->getName() == name) return (*i)->get<T>();

            for (vector<Requirement>::iterator j = (*i)->getRequirements().begin();
                 j != (*i)->getRequirements().end();j++) new_to_search += &j->get();
        }

        to_search = new_to_search;
    }

    throw logic_error("Product " + name " + not found on task " + this->name + " or its dependencies");
}

Product& Task::getProduct(const string& name)
{
    for (vector<Product>::iterator i = products.begin();i != products.end();i++)
    {
        if (i->getName() == name) return *i;
    }
    throw logic_error("Product " + name + " not found on task " + this->name);
}

Task* Task::createTask(const string& type, const string& name, const input::Config& config)
{
    map<string,factory_func>::iterator i = tasks().find(type);

    if (i == tasks().end()) throw logic_error("Task type " + type + " not found");

    return i->second(name, config);
}
