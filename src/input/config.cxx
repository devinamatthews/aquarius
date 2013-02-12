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
 * ARE DISCLAIMED. IN NO EVENT SHALL EDGAR SOLOMONIK BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "config.hpp"

#include <cstring>

using namespace std;

namespace aquarius
{
namespace input
{

Config::Config(node_t* node)
{
    root = node;
    attachNode(root);
}

Config::Config()
{
    root = new node_t;
    root->data = "root";
    root->children = NULL;
    root->parent = NULL;
    root->next = NULL;
    root->prev = NULL;
    root->ref_count = 1;
}

Config::Config(istream& is)
{
    root = new node_t;
    root->children = NULL;
    read(is);
}

Config::Config(const string& file)
{
    root = new node_t;
    root->children = NULL;
    read(file);
}

Config::Config(const Config& copy)
{
    root = copy.root;
    attachNode(root);
}

Config::~Config()
{
    deleteNode(root);
}

Config& Config::operator=(const Config& copy)
{
    deleteNode(root);
    root = copy.root;
    attachNode(root);

    return *this;
}

void Config::attachNode(node_t* node)
{
    node_t* child;

    node->ref_count++;

    for (child = node->children;child != NULL;child = child->next)
    {
        attachNode(child);
    }
}

void Config::deleteNode(node_t* node)
{
    node_t *child, *next;

    node->ref_count--;

    /*
     * decrement count for each child, they may self-destruct and remove themselves from node->children
     */
    for (child = node->children;child != NULL;child = next)
    {
        next = child->next;
        deleteNode(child);
    }

    if (node->ref_count == 0)
    {
        /*
         * for each remaining child, make it an orphan
         */
        for (child = node->children;child != NULL;child = child->next)
        {
            child->parent = NULL;
            child->next = NULL;
            child->prev = NULL;
        }

        /*
         * if we have are the main child, update our parent (if we have one)
         */
        if (node->prev == NULL && node->parent != NULL) node->parent->children = node->next;

        /*
         * and update siblings if necessary
         */
        if (node->prev != NULL) node->prev->next = node->next;
        if (node->next != NULL) node->next->prev = node->prev;

        delete node;
    }
}

void Config::writeNode(ostream& os, const node_t* n, const int level) const
{
    int i;
    node_t* c;

    for (i = 0;i < level;i++) os << '\t';
    if (n->data.find(' ') != string::npos)
    {
        os << "\"" << n->data << "\"\n";
    }
    else
    {
        os << n->data << "\n";
    }

    for (c = n->children;c != NULL;c = c->next)
    {
        writeNode(os, c, level + 1);
    }
}

Config Config::get(const string& path, const int which)
{
    node_t* n = resolve(root, path);

    if (n == NULL) throw EntryNotFoundError(path);

    return Config(n);
}

const Config Config::get(const string& path, const int which) const
{
    node_t* n = resolve(root, path);

    if (n == NULL) throw EntryNotFoundError(path);

    return Config(n);
}

string Config::path(node_t* node)
{
    if (node->parent == NULL)
    {
        return node->data;
    }
    else
    {
        return path(node->parent) + "." + node->data;
    }
}

Config::node_t* Config::resolve(node_t* node, const string& path, const bool create)
{
    node_t* child;
    string name;
    size_t pos, pos2;
    int which = 0;

    pos = path.find('.');
    name = path.substr(0, pos);

    if (name[name.length()-1] == ']')
    {
        pos2 = name.find('[');
        if (pos2 == string::npos) throw EntryNotFoundError(path);
        if (!(istringstream(name.substr(pos2+1,name.length()-1)) >> which)) throw EntryNotFoundError(path);
        name = name.substr(0, pos2);
    }

    for (child = node->children;child != NULL;child = child->next)
    {
        if (name == child->data)
        {
            if (which == 0) break;
            which--;
        }
    }

    if (child == NULL && create)
    {
        which++;

        if (node->children == NULL)
        {
            node->children = new node_t;
            *node->children = node_t(name, node, NULL, NULL, NULL, node->ref_count);
            which--;
        }

        for (child = node->children;which > 0;child = child->next, which--)
        {
            child->next = new node_t;
            *child->next = node_t(name, node, NULL, NULL, child, node->ref_count);
        }
    }

    if (pos == string::npos || child == NULL)
    {
        return child;
    }
    else
    {
        return resolve(child, path.substr(pos+1));
    }
}

string Config::readEntry(istream& is, string& line, int& lineno)
{
    int pos;
    string entry;

    if (line[0] == '"')
    {
        line = line.substr(1);
        entry = "";

        while ((pos = line.find('"', 1)) == string::npos)
        {
            entry += line;
            lineno++;
            if (!getline(is, line)) throw FormatError("unbalanced quotation marks", lineno);
        }

        entry += line.substr(pos);
    }
    else
    {
        pos = line.find(' ');
        entry = line.substr(0, pos);
        if (pos == string::npos)
        {
            line = "";
        }
        else
        {
            pos = line.find_first_not_of(' ', pos);
            if (pos == string::npos)
            {
                line = "";
            }
            else
            {
                line = line.substr(pos);
            }
        }
    }

    return entry;
}

void Config::addNode(node_t* parent, string& data)
{
    node_t* child = new node_t;

    child->children = NULL;
    child->data = data;
    child->ref_count = parent->ref_count;

    addNode(parent, child);
}

void Config::addNode(node_t* parent, node_t* child)
{
    node_t* sib;

    if (parent->children == NULL)
    {
        child->prev = NULL;
        parent->children = child;
    }
    else
    {
        for (sib = parent->children;sib->next != NULL;sib = sib->next);
        sib->next = child;
        child->prev = sib;
    }

    child->parent = parent;
    child->next = NULL;
}

void Config::read(const string& file)
{
    ifstream ifs(file.c_str());
    read(ifs);
}

void Config::read(istream& is)
{
    string line;
    node_t* current;
    int oldindent = -1;
    int newindent;
    int lineno;
    int i;
    string entry;

    deleteNode(root);
    root = new node_t;
    root->data = "root";
    root->children = NULL;
    root->parent = NULL;
    root->next = NULL;
    root->prev = NULL;
    root->ref_count = 1;

    current = root;

    for (lineno = 1;getline(is, line);lineno++)
    {
        newindent = line.find_first_not_of("\t");

        if (newindent == string::npos || line[0] == '#') continue;
        if (newindent > oldindent + 1) throw FormatError("too much indentation", lineno);

        line = line.substr(newindent);

        for (i = 0;i < (oldindent+1) - newindent;i++) current = current->parent;

        entry = readEntry(is, line, lineno);
        addNode(current, entry);
        for (current = current->children;current->next != NULL;current = current->next);

        while ((entry = readEntry(is, line, lineno)).length() > 0)
        {
            addNode(current, entry);
        }

        oldindent = newindent;
    }
}

void Config::write(const string& file) const
{
    ofstream ofs(file.c_str());
    write(ofs);
}

void Config::write(ostream& os) const
{
    node_t* child;

    for (child = root->children;child != NULL;child = child->next)
    {
        writeNode(os, child, 0);
    }
}

bool Config::matchesWithGlobs(const char* str, const char* pat) const
{
    char* s = (char*)str;
    char* p = (char*)pat;

    while (true)
    {
        if (*s == '\0')
        {
            if (*p == '\0')
            {
                return true;
            }
            else
            {
                /*
                 * if all remaining chars in pat are *, str matches
                 * otherwise, nothing can be done
                 */
                for (;*p == '*';p++);
                return (*p == '\0');
            }
        }
        else
        {
            if (*p == '\0') return false;

            if (*p == '*')
            {
                /*
                 * find next char in pat that is not '*'
                 * if none exists, then the rest of str automatically matches
                 */
                for (;*p == '*';p++);
                if (*p == '\0') return true;

                /*
                 * try each character in str that matches this next non-'*' character
                 * if none work, then no match
                 */
                do
                {
                    for (;*s != *p;s++);
                    if (*s == '\0') return false;
                }
                while (!matchesWithGlobs(++s, p+1));

                return true;
            }
            else
            {
                /*
                 * if the next char in pat is non-'*', then it must match
                 */
                if (*s++ != *p++) return false;
            }
        }
    }

    return true;
}

//template<>
bool Config::Parser<bool>::parse(istream& is)
{
    char buf[6];
    string bad = " \t\n,]";
    is.get(buf[0]);
    if (bad.find(buf[0]) != string::npos) throw BadValueError();
    if (strncasecmp(buf, "1", 1) == 0) return true;
    if (strncasecmp(buf, "0", 1) == 0) return false;
    is.get(buf[1]);
    if (bad.find(buf[1]) != string::npos) throw BadValueError();
    if (strncasecmp(buf, "on", 2) == 0) return true;
    if (strncasecmp(buf, "no", 2) == 0) return false;
    is.get(buf[2]);
    if (bad.find(buf[2]) != string::npos) throw BadValueError();
    if (strncasecmp(buf, "yes", 3) == 0) return true;
    if (strncasecmp(buf, "off", 3) == 0) return false;
    is.get(buf[3]);
    if (bad.find(buf[3]) != string::npos) throw BadValueError();
    if (strncasecmp(buf, "true", 4) == 0) return true;
    if (strncasecmp(buf, "keep", 4) == 0) return true;
    is.get(buf[4]);
    if (bad.find(buf[4]) != string::npos) throw BadValueError();
    if (strncasecmp(buf, "false", 5) == 0) return false;
    is.get(buf[5]);
    if (bad.find(buf[5]) != string::npos) throw BadValueError();
    if (strncasecmp(buf, "delete", 6) == 0) return false;
    throw BadValueError();
}

//template<>
string Config::Extractor<string>::extract(node_t* node, int which)
{
    // attempt to follow the linked list to the 'which'th element
    for (node = node->children;node != NULL && which > 0;node = node->next, which--);
    if (node == NULL) throw NoValueError(path(node));

    return node->data;
}

template<>
void Config::set<string>(const string& path, const string& data, const int which, const bool create)
{
    node_t* n = resolve(root, path, create);
    node_t* c;
    int k = which;

    if (n == NULL) throw EntryNotFoundError(path);

    // attempt to follow the linked list to the 'which'th element
    for (c = n->children;c != NULL && k > 0;c = c->next, k--);

    if (c == NULL)
    {
        if (!create) throw NoValueError(path);

        k++;

        if (n->children == NULL)
        {
            n->children = new node_t;
            *n->children = node_t(data, n, NULL, NULL, NULL, n->ref_count);
            k--;
        }

        for (c = n->children;k > 0;c = c->next, k--)
        {
            c->next = new node_t;
            *c->next = node_t(data, n, NULL, NULL, c, n->ref_count);
        }
    }

    n->data = data;
}

void Schema::apply(Config& config) const
{
    apply(config, root, config.root);
}

void Schema::apply(Config& config, const node_t* schema, node_t* root) const
{
    node_t* r;
    node_t* s;
    bool found;
    int slen;
    int minwild, maxwild, nwild;
    string sname;

    if (!isPrimitive(schema))
    {
        minwild = 0;
        maxwild = 0;
        for (s = schema->children;s != NULL;s = s->next)
        {
            slen = s->data.length();
            if (s->data[slen - 1] == '*')
            {
                if (s->data[0] == '*')
                {
                    maxwild = -1;
                    continue;
                }

                for (r = root->children;r != NULL;r = r->next)
                {
                    if (s->data.compare(0, slen-1, r->data, 0, slen-1) == 0)
                    {
                        apply(config, s, r);
                    }
                }
            }
            else if (s->data[slen - 1] == '+')
            {
                if (s->data[0] == '*')
                {
                    minwild++;
                    maxwild = -1;
                    continue;
                }

                found = false;
                for (r = root->children;r != NULL;r = r->next)
                {
                    if (s->data.compare(0, slen-1, r->data, 0, slen-1) == 0)
                    {
                        found = true;
                        apply(config, s, r);
                    }
                }
                if (!found)
                    throw SchemaValidationError("required node not found: " + s->data);
            }
            else if (s->data[slen - 1] == '?')
            {
                if (s->data[0] == '*')
                {
                    if (maxwild != -1) maxwild++;
                    continue;
                }

                found = false;
                for (r = root->children;r != NULL;r = r->next)
                {
                    if (s->data.compare(0, slen-1, r->data, 0, slen-1) == 0)
                    {
                        if (found)
                            throw SchemaValidationError("multiple copies of one-time node found: " + s->data);
                        found = true;
                        apply(config, s, r);
                    }
                }

                if (!found)
                {
                    sname = s->data.substr(0, slen-1);
                    config.addNode(root, sname);
                    for (r = root->children;r->next != NULL;r = r->next);
                    apply(config, s, r);
                }
            }
            else
            {
                if (s->data[0] == '*')
                {
                    minwild++;
                    if (maxwild != -1) maxwild++;
                    continue;
                }

                found = false;
                for (r = root->children;r != NULL;r = r->next)
                {
                    if (s->data.compare(0, slen-1, r->data, 0, slen-1) == 0)
                    {
                        if (found)
                            throw SchemaValidationError("multiple copies of one-time node found: " + s->data);
                        found = true;
                        apply(config, s, r);
                    }
                }
                if (!found)
                    throw SchemaValidationError("required node not found: " + s->data);
            }
        }

        nwild = 0;
        for (r = root->children;r != NULL;r = r->next)
        {
            found = false;
            for (s = schema->children;s != NULL;s = s->next)
            {
                slen = s->data.length();
                if (s->data[slen-1] == '*' || s->data[slen-1] == '+' || s->data[slen-1] == '?') slen--;
                if (s->data.compare(0, slen, r->data, 0, slen) == 0)
                {
                    found = true;
                    break;
                }
            }
            if (!found) nwild++;
        }

        if (nwild < minwild || (maxwild != -1 && nwild > maxwild))
        {
            throw SchemaValidationError("too few or too many wildcard nodes: " + root->data);
        }
    }
    else //primitive
    {
        schema = schema->children;

        if (root->children == NULL) //value not specified
        {
            if (hasDefault(schema))
            {
                config.addNode(root, schema->children->data);
            }
            else //no default
            {
                slen = schema->parent->data.length();
                if (schema->parent->data[slen-1] == '?')
                {
                    //ensure deletion
                    root->ref_count = 1;
                    config.deleteNode(root);
                }
                else
                {
                    throw SchemaValidationError("no value specified and no default: " + root->data);
                }
            }
        }
        else //value specified
        {
            if (root->children->next != NULL)
                throw SchemaValidationError("multiple values for primitive: " + root->data);

            if (schema->data == "int")
            {
                if (!isInt(root->children->data))
                    throw SchemaValidationError("data is not an integer: " + root->data);
            }
            else if (schema->data == "double")
            {
                if (!isDouble(root->children->data))
                    throw SchemaValidationError("data is not a double: " + root->data);
            }
            else if (schema->data == "bool")
            {
                if (!isBool(root->children->data))
                    throw SchemaValidationError("data is not a boolean: " + root->data);
            }
            else if (schema->data == "enum")
            {
                found = false;
                for (s = schema->children;s != NULL;s = s->next)
                {
                    if (root->children->data == s->data)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    throw SchemaValidationError("value is not in enum: " + root->data);
            }
            else if (schema->data == "string")
            {
                //anything goes
            }
        }
    }
}

bool Schema::isInt(const string& data) const
{
    int pos = 0;
    if (data[0] == '-') pos = 1;
    if (data.find_first_not_of("0123456789", pos) == string::npos) return true;
    return false;
}

bool Schema::isDouble(const string& data) const
{
    double x;
    if (!(istringstream(data) >> x)) return false;
    return true;
}

bool Schema::isBool(const string& data) const
{
    if (strcasecmp(data.c_str(), "yes") == 0 ||
        strcasecmp(data.c_str(), "true") == 0 ||
        strcasecmp(data.c_str(), "keep") == 0 ||
        strcasecmp(data.c_str(), "1") == 0 ||
        strcasecmp(data.c_str(), "on") == 0 ||
        strcasecmp(data.c_str(), "no") == 0 ||
        strcasecmp(data.c_str(), "false") == 0 ||
        strcasecmp(data.c_str(), "delete") == 0 ||
        strcasecmp(data.c_str(), "off") == 0 ||
        strcasecmp(data.c_str(), "0") == 0) return true;

    return false;
}

bool Schema::isPrimitive(const node_t* schema) const
{
    if (schema->children == NULL) return false;

    if (schema->children->data == "int" ||
        schema->children->data == "double" ||
        schema->children->data == "bool" ||
        schema->children->data == "enum" ||
        schema->children->data == "string") return true;

    return false;
}

bool Schema::hasDefault(const node_t* schema) const
{
    if (schema->children == NULL) return false;

    return true;
}

}
}
