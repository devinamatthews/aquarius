#ifndef _AQUARIUS_INPUT_CONFIG_HPP_
#define _AQUARIUS_INPUT_CONFIG_HPP_

#include "util/global.hpp"

namespace aquarius
{
namespace input
{

class Config;
class Schema;

class EntryNotFoundError : public runtime_error
{
    public:
        EntryNotFoundError(const string& what_arg) : runtime_error(what_arg) {}
};

class NoValueError : public runtime_error
{
    public:
        NoValueError(const string& what_arg) : runtime_error(what_arg) {}
};

class BadValueError : public runtime_error
{
    public:
        BadValueError() : runtime_error("bad value") {}

        BadValueError(const string& what_arg) : runtime_error(what_arg) {}
};

class SchemaValidationError : public runtime_error
{
    public:
        SchemaValidationError(const string& what_arg) : runtime_error(what_arg) {}
};

class FormatError : public runtime_error
{
    private:
        static string buildString(const string& what_arg, const int lineno)
        {
            string s;
            ostringstream os(s);
            os << what_arg << ": line " << lineno;
            return s;
        }

    public:
        FormatError(const string& what_arg, const int lineno)
        : runtime_error(buildString(what_arg, lineno)) {}
};

class Config
{
    friend class Schema;

    public:
        enum {ALL=-1};

    protected:
        struct Node
        {
            string data;
            Node *parent;
            shared_list<Node> children;

            Node() : parent(NULL) {}

            string path() const;

            string fullName() const;

            shared_ptr<Node> clone() const;

            void write(ostream& os, int level = 0) const;

            shared_list<Node>::iterator getChild(int which);

            shared_list<Node>::iterator addChild(const string& data = string());

            void removeChild(const Node& child);
        };

        friend ostream& operator<<(ostream& os, const Node& node);

        shared_ptr<Node> root;

        Config(const shared_ptr<Node>& node) : root(node) {}

        Config(shared_ptr<Node>&& node) : root(move(node)) {}

        static Node* resolve(Node& start, const string& path, bool create = false);

        string readEntry(istream& is, string& line, int& lineno);

        template<typename T>
        void emit(ostream& os, T& x) const;

        template<typename T>
        void emit(ostream& os, vector<T>& x) const;

        template<typename T>
        class Parser
        {
            public:
                static T parse(istream& is);
                static T parse(string& s);
        };

        template<typename T>
        class Extractor
        {
            public:
                static T extract(Node& node, int which = 0);
        };

        bool matchesWithGlobs(const char* str, const char* pat) const;

        template<typename T>
        void find(Node& node, const string& name, const string& pattern, vector<pair<string,T>>& v) const;

    public:

        Config() : root(new Node()) {}

        Config(const Config& config) = default;

        Config(Config&& config) = default;

        Config(istream& is);

        Config(const string& s);

        Config clone() const;

        Config& operator=(const Config& copy) = default;

        Config& operator=(Config&& copy) = default;

        bool exists(const string& path) const
        {
            return resolve(*root, path) != NULL;
        }

        void remove(const string& path)
        {
            Node* node = resolve(*root, path);
            if (node == NULL) throw EntryNotFoundError(path);
            node->parent->removeChild(*node);
        }

        template<typename T>
        T get(const string& path, int which = 0) const;

        template<typename T>
        void set(const string& path, const T& data, const int which = 0, const bool create = true);

        template<typename T>
        vector<pair<string,T>> find(const string& pattern) const;

        void read(const string& cwd, istream& is);

        void read(const string& file);

        void write(ostream& os = cout) const;

        void write(const string& file) const;

        Config get(const string& path);

        vector<pair<string,Config>> find(const string& pattern);
};

class Schema : public Config
{
    private:
        bool isPrimitive(const Node& schema) const;

        bool isInt(const string& data) const;

        bool isDouble(const string& data) const;

        bool isBool(const string& data) const;

        bool hasDefault(const Node& schema) const;

        void apply(const Node& schema, Node& root) const;

    public:
        Schema() : Config() {}

        Schema(istream& is) : Config(is) {}

        Schema(const string& s) : Config(s) {}

        Schema(Config& copy) : Config(copy) {}

        void apply(Config& config) const;
};

template<typename S>
class Config::Parser<vector<S>>
{
    public:
    static vector<S> parse(istream& is)
    {
        vector<S> x;
        char delim;

        if (!(is >> skipws >> delim)) throw BadValueError();
        if (delim != '[') throw BadValueError();

        while (true)
        {
            x.push_back(Parser<S>::parse(is));
            if (!(is >> skipws >> delim)) throw BadValueError();
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
    static bool parse(istream& is);
};

template<typename T>
T Config::Parser<T>::parse(istream& is)
{
    T x;
    if (!(is >> x)) throw BadValueError();
    return x;
}

template<typename T>
T Config::Parser<T>::parse(string& s)
{
    istringstream iss(s);
    return parse(iss);
}

template<typename S>
class Config::Extractor<vector<S>>
{
    public:
    static vector<S> extract(Node& node, int which = 0)
    {
        if (which != ALL)
        {
            auto i = node.getChild(which);
            if (i == node.children.end()) throw NoValueError(node.fullName());

            istringstream iss(i->data);
            return Parser<vector<S>>::parse(iss);
        }
        else
        {
            vector<S> v;
            for (int i = 0;i < node.children.size();i++)
            {
                v.push_back(Extractor<S>::extract(node, i));
            }
            return v;
        }
    }
};

template<>
class Config::Extractor<string>
{
    public:
    static string extract(Node& node, int which = 0);
};

template<>
class Config::Extractor<Config>
{
    public:
    static Config extract(Node& node, int which = 0);
};

template<typename T>
T Config::Extractor<T>::extract(Node& node, int which)
{
    auto i = node.getChild(which);
    if (i == node.children.end()) throw NoValueError(node.fullName());

    istringstream iss(i->data);
    return Parser<T>::parse(iss);
}

template<typename T>
void Config::emit(ostream& os, vector<T>& x) const
{
    typename vector<T>::iterator i;

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
void Config::emit(ostream& os, T& x) const
{
    if (!(os << x)) throw BadValueError();
}

template<typename T>
T Config::get(const string& path, int which) const
{
    Node* n = resolve(*root, path);
    if (n == NULL) throw EntryNotFoundError(path);
    return Extractor<T>::extract(*n, which);
}

template<>
void Config::set<string>(const string& path, const string& data, int which, bool create);

template<typename T>
void Config::set(const string& path, const T& data, int which, bool create)
{
    string s;
    emit(ostringstream(s), data);
    set(path, s, create);
}

template<typename T>
void Config::find(Node& node, const string& name, const string& pattern, vector<pair<string,T>>& v) const
{
    string toMatch;
    string remainder;

    size_t pos = pattern.find('.');
    if (pos == string::npos)
    {
        toMatch = pattern;
    }
    else
    {
        toMatch = pattern.substr(0, pos);
        remainder = pattern.substr(pos+1);
    }

    for (Node& c : node.children)
    {
        if (matchesWithGlobs(c.data.c_str(), toMatch.c_str()))
        {
            if (pos == string::npos)
            {
                v.push_back(make_pair(name + c.data, Extractor<T>::extract(c)));
            }
            else
            {
                find(c, name + c.data + ".", remainder, v);
            }
        }
    }
}

template<typename T>
vector<pair<string,T>> Config::find(const string& pattern) const
{
    vector<pair<string,T>> v;
    find(*root, "", pattern, v);
    return v;
}

}
}

#endif
