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
using namespace aquarius::time;
using namespace aquarius::input;

namespace aquarius
{
namespace task
{

Logger::NullStream Logger::nullstream;
Logger::SelfDestructStream Logger::sddstream;

int Logger::LogToStreamBuffer::sync()
{
    if (os.flush()) return 0;
    else return -1;
}

std::streamsize Logger::LogToStreamBuffer::xsputn (const char* s, std::streamsize n)
{
    if (os.write(s, n)) return n;
    else return 0;
}

int Logger::LogToStreamBuffer::overflow (int c)
{
    if (c != EOF) os << traits_type::to_char_type(c);
    if (c == EOF || !os) return EOF;
    else return traits_type::to_int_type(c);
}

Logger::LogToStreamBuffer::LogToStreamBuffer(std::ostream& os)
: os(os) {}

int Logger::SelfDestructBuffer::sync()
{
    LogToStreamBuffer::sync();
    abort();
    return -1;
}

Logger::SelfDestructBuffer::SelfDestructBuffer()
: LogToStreamBuffer(std::cerr) {}

std::streamsize Logger::NullBuffer::xsputn (const char* s, std::streamsize n)
{
    return n;
}

int Logger::NullBuffer::overflow (int c)
{
    if (c == EOF) return EOF;
    else return traits_type::to_int_type(c);
}

Logger::NullStream::NullStream()
: std::ostream(new NullBuffer()) {}

Logger::NullStream::~NullStream()
{
    delete rdbuf();
}

Logger::SelfDestructStream::SelfDestructStream()
: std::ostream(new SelfDestructBuffer()) {}

std::string Logger::dateTime()
{
    char buf[256];
    time_t t = ::time(NULL);
    tm* timeptr = localtime(&t);
    strftime(buf, 256, "%c", timeptr);
    return std::string(buf);
}

std::ostream& Logger::log(const Arena& arena)
{
    if (arena.rank == 0)
    {
        std::cout << dateTime() << ": ";
        return std::cout;
    }
    else
    {
        return nullstream;
    }
}

std::ostream& Logger::warn(const Arena& arena)
{
    if (arena.rank == 0)
    {
        std::cerr << dateTime() << ": warning: ";
        return std::cerr;
    }
    else
    {
        return nullstream;
    }
}

std::ostream& Logger::error(const Arena& arena)
{
    if (arena.rank == 0)
    {
        sddstream << dateTime() << ": error: ";
        return sddstream;
    }
    else
    {
        pause();
        return nullstream;
    }
}

Requirement::Requirement(const string& type, const string& name)
: type(type), name(name) {}

void Requirement::fulfil(const Product& product)
{
    this->product.reset(new Product(product));
    *this->product->used = true;
}

bool Requirement::exists() const
{
    return product->exists();
}

Product& Requirement::get()
{
    if (!product) throw logic_error("Requirement " + name + " not fulfilled");
    return *product;
}

Product::Product(const string& type, const string& name)
: type(type), name(name), requirements(new vector<Requirement>()), used(new bool(false)) {}

Product::Product(const string& type, const string& name, const vector<Requirement>& reqs)
: type(type), name(name), requirements(new vector<Requirement>(reqs)), used(new bool(false)) {}

void Product::addRequirement(const Requirement& req)
{
    requirements->push_back(req);
}

void Product::addRequirements(const std::vector<Requirement>& reqs)
{
    requirements->insert(requirements->end(), reqs.begin(), reqs.end());
}

template <> string Task::type_string<float>() { return "<float>"; }
template <> string Task::type_string<double>() { return ""; }
template <> string Task::type_string<complex<float> >() { return "<scomplex>"; }
template <> string Task::type_string<complex<double> >() { return "<dcomplex>"; }

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

const Schema& Task::getSchema(const string& name)
{
    static map<string,Schema> schemas;

    if (schemas.empty())
    {
        Config master(TOPDIR "/desc/schemas");

        vector<pair<string,Config> > configs = master.find<Config>("*");
        for (vector<pair<string,Config> >::iterator i = configs.begin();i != configs.end();++i)
        {
            Config schema = master.get(i->first);
            schemas[i->first] = Schema(schema);
        }
    }

    map<string,Schema>::iterator i = schemas.find(name);
    if (i == schemas.end()) throw logic_error("Cannot find schema for task " + name);
    return i->second;
}

Product& Task::getProduct(const string& name)
{
    for (vector<Product>::iterator i = products.begin();i != products.end();i++)
    {
        if (i->getName() == name) return *i;
    }
    throw logic_error("Product " + name + " not found on task " + this->name);
}

const Product& Task::getProduct(const string& name) const
{
    return const_cast<Task&>(*this).getProduct(name);
}

Task* Task::createTask(const string& type, const string& name, const input::Config& config)
{
    map<string,factory_func>::iterator i = tasks().find(type);

    if (i == tasks().end()) throw logic_error("Task type " + type + " not found");

    return i->second(name, config);
}

TaskDAG::TaskDAG(const std::string& file)
{
    Config input(file);

    vector<pair<string,Config> > taskconfs = input.find<Config>("*");

    for (vector<pair<string,Config> >::iterator i = taskconfs.begin();i != taskconfs.end();++i)
    {
        int num = 0;
        for (vector<Task*>::iterator t = tasks.begin();t != tasks.end();++t)
        {
            if ((*t)->getType() == i->first) num++;
        }

        Task::getSchema(i->first).apply(i->second);
        tasks.push_back(Task::createTask(i->first, i->first+str(num), i->second));
    }
}

TaskDAG::~TaskDAG()
{
    for (vector<Task*>::iterator i = tasks.begin();i != tasks.end();++i) delete *i;
    tasks.clear();
}

void TaskDAG::addTask(Task* task)
{
    tasks.push_back(task);
}

void TaskDAG::execute(Arena& world)
{
    /*
     * Attempy to satisfy task requirements greedily. If we are not careful, this could produce cycles.
     */
    for (vector<Task*>::iterator t1 = tasks.begin();t1 != tasks.end();++t1)
    {
        for (vector<Product>::iterator p1 = (*t1)->getProducts().begin();p1 != (*t1)->getProducts().end();++p1)
        {
            for (vector<Requirement>::iterator r1 = p1->getRequirements().begin();r1 != p1->getRequirements().end();++r1)
            {
                if (r1->isFulfilled()) continue;

                for (vector<Task*>::iterator t2 = tasks.begin();t2 != tasks.end();++t2)
                {
                    for (vector<Product>::iterator p2 = (*t2)->getProducts().begin();p2 != (*t2)->getProducts().end();++p2)
                    {
                        if (r1->isFulfilled()) continue;
                        if (r1->getType() == p2->getType())
                        {
                            r1->fulfil(*p2);
                        }

                        for (vector<Requirement>::iterator r2 = p2->getRequirements().begin();r2 != p2->getRequirements().end();++r2)
                        {
                            if (r1->isFulfilled()) continue;
                            if (!r2->isFulfilled()) continue;
                            if (r1->getType() == r2->getType())
                            {
                                r1->fulfil(r2->get());
                            }
                        }
                    }
                }

                if (!r1->isFulfilled())
                    Logger::error(world) << "Could not fulfil requirement " << r1->getName() <<
                                            " of task " << (*t1)->getName() << endl;
            }
        }
    }

    //TODO: check for cycles

    /*
     * Successively search for executable tasks
     */
    while (true)
    {
        set<Task*> to_execute;

        for (vector<Task*>::iterator t = tasks.begin();t != tasks.end();++t)
        {
            bool can_execute = true;

            for (vector<Product>::iterator p = (*t)->getProducts().begin();p != (*t)->getProducts().end();++p)
            {
                for (vector<Requirement>::iterator r = p->getRequirements().begin();r != p->getRequirements().end();++r)
                {
                    if (!r->exists())
                    {
                        can_execute = false;
                        break;
                    }
                }
                if (!can_execute) break;
            }

            if (can_execute) to_execute.insert(*t);
        }

        if (to_execute.empty()) break;

        for (set<Task*>::iterator t = to_execute.begin();t != to_execute.end();++t)
        {
            Logger::log(world) << "Starting task: " << (*t)->getName() << endl;
            tic();

            try
            {
                (*t)->run(*this, world);
            }
            catch (runtime_error& e)
            {
                Logger::error(world) << e.what() << endl;
            }

            double dt = todouble(toc());
            Logger::log(world) << "Finished task: " << (*t)->getName() << " in " << dt << " s" << endl;

            for (vector<Product>::iterator p = (*t)->getProducts().begin();p != (*t)->getProducts().end();++p)
            {
                if (p->isUsed() && !p->exists())
                    Logger::error(world) << "Product " << p->getName() <<
                                            " of task " << (*t)->getName() <<
                                            " was not successfully produced" << endl;
            }

            tasks.erase(std::find(tasks.begin(), tasks.end(), *t));
            delete *t;
        }
    }

    if (!tasks.empty())
    {
        Logger::error(world) << "Some tasks were not executed" << endl;
    }
}

}
}
