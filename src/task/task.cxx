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

using namespace std;
using namespace aquarius;
using namespace aquarius::task;

template <typename T>
Resource::Resource(const std::string& type, T& data)
: type(type), data(static_cast<void*>(&data)), revision(0), needsDelete(false) {}

template <typename T>
Resource::Resource(const std::string& type, T* data)
: type(type), data(static_cast<void*>(data)), revision(0), needsDelete(true) {}

Resource::~Resource()
{
    if (needsDelete) delete data;
}

template <typename T> T& Resource::get()
{
    return *static_cast<T*>(data);
}

template <typename T> const T& Resource::get() const
{
    return *static_cast<T*>(data);
}

Requirement::Requirement(const std::string& name, const std::string& type)
: name(name), type(type) {}

Requirement::Requirement(const std::string& name, const std::string& type, const std::string& defaultTask, const std::string& defaultName)
: name(name), type(type), defaultProvider(std::make_pair(defaultTask, defaultName)) {}

void Requirement::fulfil(Product& product)
{
    this->product.reset(&product);
}

Resource& Requirement::getResource()
{
    if (!product) throw std::logic_error("No resource has been provided for requirement " + name + ".");
    return product->getResource();
}

void Product::run_check()
{
    for (int i = 0;i < requirements.size();i++)
    {
        if (!requirements[i].isFulfilled())
            throw logic_error("Requirement " + requirements[i].getName() + " of product " + name + " is not fulfilled.");

        Resource& res = requirements[i].getResource();

        if (requirements[i].type != res.getType())
            throw logic_error("Product of type " + res.getType() + " provides the wrong resource for requirement " + requirements[i].getName() + ".");

        revisions[i] = res.getRevision();
    }

    if (parent == NULL) throw std::logic_error("No parent task has been assigned to product " + name + ".");

    parent.run();
}

Product::Product(Task& parent, const std::string& name, const std::shared_ptr<Resource>& resource)
: resource(resource), parent(parent), name(name) {}

Product::Product(Task& parent, const std::string& name, const std::shared_ptr<Resource>& resource, const Requirement& requirement)
: resource(resource), parent(parent), name(name), requirements(1, requirement), revisions(1, -1) {}

Product::Product(Task& parent, const std::string& name, const std::shared_ptr<Resource>& resource, const std::vector<Requirement>& requirements)
: resource(resource), parent(parent), name(name), requirements(requirements), revisions(requirements.size(), -1) {}

Resource& Product::getResource()
{
    if (!resource) run_check();

    for (int i = 0;i < requirements.size();i++)
    {
        if (revisions[i] != requirements[i].getResource().getRevision()) run_check();
    }

    return *resource;
}

Task::Task(const std::string& type, const std::string& name)
: type(type), name(name) {}

Task::Task(const std::string& type, const std::string& name, const Product& product)
: type(type), name(name), products(1, product)
{
    products.back().parent = this;
}

Task::Task(const std::string& type, const std::string& name, const std::vector<Product>& products)
: type(type), name(name), products(products)
{
    for (std::vector<Product>::iterator i = products.begin();i != products.end();++i) i->parent = this;
}

Product& Task::getProduct(const std::string& name)
{
    for (std::vector<Product>::iterator i = products.begin();i != products.end();++i)
    {
        if (i->name == name) return *i;
    }
    throw std::logic_error("No product " + name + " found in task " + this->name + ".");
}

const Product& Task::getProduct(const std::string& name) const
{
    for (std::vector<Product>::const_iterator i = products.begin();i != products.end();++i)
    {
        if (i->name == name) return *i;
    }
    throw std::logic_error("No product " + name + " found in task " + this->name + ".");
}
