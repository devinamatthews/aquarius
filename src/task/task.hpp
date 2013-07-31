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
    protected:
        std::string type;
        void* data;
        int revision;
        bool needsDelete;

    public:
        template <typename T>
        Resource(const std::string& type, T& data);

        template <typename T>
        Resource(const std::string& type, T* data);

        ~Resource();

        const std::string& getType() const { return type; }

        template <typename T> T& get();

        template <typename T> const T& get() const;

        void update() { revision++; }

        int getRevision() const { return revision; }
};

class Product;

class Requirement
{
    protected:
        std::string name;
        std::string type;
        std::pair<std::string,std::string> defaultProvider;
        std::shared_ptr<Product> product;

    public:
        Requirement(const std::string& name, const std::string& type);

        Requirement(const std::string& name, const std::string& type, const std::string& defaultTask, const std::string& defaultName);

        const std::string& getName() const { return name; }

        const std::string getType() const { return type; }

        void fulfil(Product& product);

        bool isFulfilled() const { return product; }

        const std::pair<std::string,std::string>& getDefaultProvider() const { return defaultProvider; }

        Resource& getResource();
};

class Task;

class Product
{
    friend class Task;

    protected:
        std::vector<Requirement> requirements;
        std::vector<int> revisions;
        std::shared_ptr<Resource> resource;
        Task& parent;
        std::string name;

        void run_check();

    public:
        Product(Task& parent, const std::string& name, const std::shared_ptr<Resource>& resource);

        Product(Task& parent, const std::string& name, const std::shared_ptr<Resource>& resource, const Requirement& requirement);

        Product(Task& parent, const std::string& name, const std::shared_ptr<Resource>& resource, const std::vector<Requirement>& requirements);

        const std::string& getName() const { return name; }

        const std::vector<Requirement> getRequirements() const { return requirements; }

        Resource& getResource();
};

class Task
{
    protected:
        std::string type;
        std::string name;
        std::vector<Product> products;

        void addProduct(const Product& product);

    public:
        Task(const std::string& type, const std::string& name);

        Task(const std::string& type, const std::string& name, const Product& product);

        Task(const std::string& type, const std::string& name, const std::vector<Product>& products);

        virtual ~Task() {}

        Product& getProduct(const std::string& name);

        const Product& getProduct(const std::string& name) const;

        std::vector<Product>& getProducts() { return products; }

        const std::vector<Product>& getProducts() const { return products; }

        virtual void run() = 0;
};

}
}

#endif
