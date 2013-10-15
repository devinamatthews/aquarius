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

#ifndef _AQUARIUS_INPUT_CONFIG_HPP_
#define _AQUARIUS_INPUT_CONFIG_HPP_

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <iterator>
#include <utility>
#include <stdexcept>
#include <exception>
#include <cstring>
#include <strings.h>

#include "util/stl_ext.hpp"

namespace aquarius
{
namespace input
{

class Config;
class Schema;

class EntryNotFoundError : public std::runtime_error
{
    public:
        EntryNotFoundError(const std::string& what_arg) : runtime_error(what_arg) {}
};

class NoValueError : public std::runtime_error
{
    public:
        NoValueError(const std::string& what_arg) : runtime_error(what_arg) {}
};

class BadValueError : public std::runtime_error
{
    public:
        BadValueError() : runtime_error("bad value") {}

        BadValueError(const std::string& what_arg) : runtime_error(what_arg) {}
};

class SchemaValidationError : public std::runtime_error
{
    public:
        SchemaValidationError(const std::string& what_arg) : runtime_error(what_arg) {}
};

class FormatError : public std::runtime_error
{
    private:
        static std::string buildString(const std::string& what_arg, const int lineno)
        {
            std::string s;
            std::ostringstream os(s);
            os << what_arg << ": line " << lineno;
            return s;
        }

    public:
        FormatError(const std::string& what_arg, const int lineno)
        : runtime_error(buildString(what_arg, lineno)) {}
};

class Config
{
    friend class Schema;

    public:
        enum {ALL=-1};

    protected:
        struct node_t
        {
            std::string data;
            node_t *parent;
            node_t *children;
            node_t *next;
            node_t *prev;
            int ref_count;
            node_t() : parent(NULL), children(NULL), next(NULL), prev(NULL), ref_count(0) {}
            node_t(const std::string data, node_t* parent, node_t* children, node_t* next, node_t* prev, const int ref_count)
            : data(data), parent(parent), children(children), next(next), prev(prev), ref_count(ref_count) {}
        };
        node_t* root;

        Config(node_t* node);

        node_t* cloneNode(node_t* node) const;

        void attachNode(node_t* node);

        void detachNode(node_t* node);

        void deleteNode(node_t* node);

        static node_t* resolve(node_t* start, const std::string& path, const bool create = false);

        void writeNode(std::ostream& os, const node_t* n, const int level) const;

        std::string readEntry(std::istream& is, std::string& line, int& lineno);

        node_t* addNode(node_t* parent, std::string& data);

        void addNode(node_t* parent, node_t* child);

        template<typename T>
        void emit(std::ostream& os, T& x) const;

        template<typename T>
        void emit(std::ostream& os, std::vector<T>& x) const;

        template<typename T>
        class Parser
        {
            public:
                static T parse(std::istream& is);
                static T parse(std::string& s);
        };

        template<typename T>
        class Extractor
        {
            public:
                static T extract(node_t* node, int which = 0);
        };

        bool matchesWithGlobs(const char* str, const char* pat) const;

        template<typename T>
        void find(node_t* node, const std::string name, const std::string pattern, std::vector< std::pair<std::string,T> >& v) const;

        static std::string path(node_t* node);

    public:
        Config();

        Config(std::istream& is);

        Config(const std::string& file);

        Config(const Config& copy);

        ~Config();

        Config clone() const;

        Config& operator=(const Config& copy);

        bool exists(const std::string& path) const { return resolve(root, path) != NULL; }

        void remove(const std::string& path)
        {
            node_t* node = resolve(root,path);
            if (node == NULL) throw EntryNotFoundError(path);
            deleteNode(node);
        }

        template<typename T>
        T get(const std::string& path, const int which = 0) const;

        template<typename T>
        void set(const std::string& path, const T& data, const int which = 0, const bool create = true);

        template<typename T>
        std::vector< std::pair<std::string,T> > find(const std::string& pattern) const;

        void read(std::istream& is);

        void read(const std::string& file);

        void write(std::ostream& os = std::cout) const;

        void write(const std::string& file) const;

        Config get(const std::string& path, int which = 0);

        const Config get(const std::string& path, int which = 0) const;

        std::vector< std::pair<std::string,Config> > find(const std::string& pattern);

        std::vector< std::pair<std::string,const Config> > find(const std::string& pattern) const;
};

class Schema : public Config
{
    private:
        bool isPrimitive(const node_t* schema) const;

        bool isInt(const std::string& data) const;

        bool isDouble(const std::string& data) const;

        bool isBool(const std::string& data) const;

        bool hasDefault(const node_t* schema) const;

        void apply(Config& config, const node_t* schema, node_t* root) const;

    public:
        Schema() : Config() {}

        Schema(std::istream& is) : Config(is) {}

        Schema(const std::string& file) : Config(file) {}

        Schema(Config& copy) : Config(copy) {}

        void apply(Config& config) const;
};

template<typename S>
class Config::Parser< std::vector<S> >
{
    public:
    static std::vector<S> parse(std::istream& is)
    {
        std::vector<S> x;
        char delim;

        if (!(is >> std::skipws >> delim)) throw BadValueError();
        if (delim != '[') throw BadValueError();

        while (true)
        {
            x.push_back(Parser<S>::parse(is));
            if (!(is >> std::skipws >> delim)) throw BadValueError();
            if (delim == ']') break;
            if (delim != ',') throw BadValueError();
        }

        return x;
    }
};

template<>
class Config::Parser<bool>
{
    public:
    static bool parse(std::istream& is);
};

template<typename T>
T Config::Parser<T>::parse(std::istream& is)
{
    T x;
    if (!(is >> x)) throw BadValueError();
    return x;
}

template<typename T>
T Config::Parser<T>::parse(std::string& s)
{
    std::istringstream iss(s);
    return parse(iss);
}

template<typename S>
class Config::Extractor< std::vector<S> >
{
    public:
    static std::vector<S> extract(node_t* node, int which = 0)
    {
        if (which != ALL)
        {
            // attempt to follow the linked list to the 'which'th element
            for (node = node->children;node != NULL && which > 0;node = node->next, which--);
            if (node == NULL) throw NoValueError(path(node));

            std::istringstream iss(node->data);
            return Parser< std::vector<S> >::parse(iss);
        }
        else
        {
            std::vector<S> v;
            int i = 0;
            for (node = node->children;node != NULL;node = node->next)
            {
                v.push_back(Extractor<S>::extract(node->parent, i++));
            }
            return v;
        }
    }
};

template<>
class Config::Extractor<std::string>
{
    public:
    static std::string extract(node_t* node, int which = 0);
};

template<>
class Config::Extractor<Config>
{
    public:
    static Config extract(node_t* node, int which = 0);
};

template<typename T>
T Config::Extractor<T>::extract(node_t* node, int which)
{
    // attempt to follow the linked list to the 'which'th element
    for (node = node->children;node != NULL && which > 0;node = node->next, which--);
    if (node == NULL) throw NoValueError(path(node));

    std::istringstream iss(node->data);
    return Parser<T>::parse(iss);
}

template<typename T>
void Config::emit(std::ostream& os, std::vector<T>& x) const
{
    typename std::vector<T>::iterator i;

    if (!(os << '[')) throw BadValueError();
    for (i = x.begin();i != x.end();)
    {
        emit(os, *i);
        if (++i != x.end())
        {
            if (!(os << ',')) throw BadValueError();
        }
    }
    if (!(os << ']')) throw BadValueError();
}

template<typename T>
void Config::emit(std::ostream& os, T& x) const
{
    if (!(os << x)) throw BadValueError();
}

template<typename T>
T Config::get(const std::string& path, const int which) const
{
    node_t* n = resolve(root, path);
    if (n == NULL) throw EntryNotFoundError(path);
    return Extractor<T>::extract(n, which);
}

template<>
void Config::set<std::string>(const std::string& path, const std::string& data, const int which, const bool create);

template<typename T>
void Config::set(const std::string& path, const T& data, const int which, const bool create)
{
    std::string s;
    emit(std::ostringstream(s), data);
    set(path, s, create);
}

template<typename T>
void Config::find(node_t* node, const std::string name, const std::string pattern, std::vector< std::pair<std::string,T> >& v) const
{
    node_t* c;
    std::string toMatch;
    std::string remainder;
    size_t pos;

    pos = pattern.find('.');
    if (pos == std::string::npos)
    {
        toMatch = pattern;
    }
    else
    {
        toMatch = pattern.substr(0, pos);
        remainder = pattern.substr(pos+1);
    }

    for (c = node->children;c;c = c->next)
    {
        if (matchesWithGlobs(c->data.c_str(), toMatch.c_str()))
        {
            if (pos == std::string::npos)
            {
                v.push_back(std::make_pair(name + c->data, Extractor<T>::extract(c)));
            }
            else
            {
                if (name == "")
                {
                    find(c, c->data + ".", remainder, v);
                }
                else
                {
                    find(c, name + c->data + ".", remainder, v);
                }
            }
        }
    }
}

template<typename T>
std::vector< std::pair<std::string,T> > Config::find(const std::string& pattern) const
{
    std::vector< std::pair<std::string,T> > v;
    find(root, "", pattern, v);
    return v;
}

}
}

#endif
