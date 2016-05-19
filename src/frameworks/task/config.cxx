#include "config.hpp"

namespace aquarius
{
namespace input
{

class Tokenizer
{
    protected:
        int lineno;
        istream& is;
        string line;
        istringstream iss;

        bool nextLine()
        {
            while (true)
            {
                if (!is) return false;
                lineno++;
                getline(is, line);
                iss.clear();
                iss.str(line);
                char c;
                while (!(!iss.get(c)))
                {
                    if (c != ' ' && c != '\t' && c != '\r')
                    {
                        if (c == '#')
                        {
                            break;
                        }
                        else
                        {
                            iss.unget();
                            return true;
                        }
                    }
                }
            }
            return false;
        }

    public:
        Tokenizer(istream& is)
        : lineno(0), is(is) {}

        int getLine() const { return lineno; }

        bool next(string& str)
        {
            ostringstream oss;
            string token;

            iss >> token;
            while (token.size() == 0)
            {
                if (!nextLine()) return false;
                iss >> token;
            }

            if (token[0] == '"')
            {
                oss << token.substr(1);
                char c;
                while (!(!iss.get(c)))
                {
                    if (c == '"') break;
                    oss << c;
                }
                if (!iss) throw FormatError("End of line reached while looking for closing \" character", lineno);
            }
            else
            {
                oss << token;
            }

            str = oss.str();

            return true;
        }
};

string Config::Node::path() const
{
    if (parent) return parent->fullName();
    else return string();
}

string Config::Node::fullName() const
{
    string p = path();
    if (p.empty()) return data;
    else return p + '.' + data;
}

shared_ptr<Config::Node> Config::Node::clone() const
{
    Node *node = new Node();
    node->data = data;
    for (const Node& child : children) node->children.push_back(child.clone());
    for (Node& child : node->children) child.parent = node;
    return shared_ptr<Node>(node);
}

void Config::Node::write(ostream& os, int level) const
{
    for (int i = 0;i < level;i++) os << '\t';

    if (data.find(' ') != string::npos)
    {
        os << "\"" << data << "\"\n";
    }
    else
    {
        os << data << "\n";
    }

    for (const Node& c : children)
    {
        c.write(os, level+1);
    }
}

shared_list<Config::Node>::iterator Config::Node::getChild(int which)
{
    auto i = children.begin();
    for (;i != children.end() && which > 0;++i, --which);
    return i;
}

shared_list<Config::Node>::iterator Config::Node::addChild(const string& data)
{
    children.emplace_back();
    children.back().parent = this;
    children.back().data = data;
    return --children.end();
}

void Config::Node::removeChild(const Node& child)
{
    children.remove_if([&child](const Node& x) { return &x == &child; });
}

ostream& operator<<(ostream& os, const Config::Node& n)
{
    n.write(os);
    return os;
}

Config::Config(istream& is)
{
    read(".", is);
}

Config::Config(const string& s)
{
    istringstream iss(s);
    read(".", iss);
}

Config Config::clone() const
{
    return Config(root->clone());
}

Config Config::get(const string& path)
{
    Node* n = resolve(*root, path);
    if (n == NULL) throw EntryNotFoundError(path);
    auto i = n->parent->children.pbegin();
    while (i->get() != n) ++i;
    return Config(*i);
}

Config::Node* Config::resolve(Node& node, const string& path, bool create)
{
    size_t pos = path.find('.');
    string name = path.substr(0, pos);

    /*
    if (name[name.length()-1] == ']')
    {
        size_t pos2 = name.find('[');
        if (pos2 == string::npos) throw EntryNotFoundError(path);
        if (!(istringstream(name.substr(pos2+1,name.length()-1)) >> which)) throw EntryNotFoundError(path);
        name = name.substr(0, pos2);
    }
    */

    Node *child = NULL;
    for (Node& c : node.children)
    {
        if (name == c.data)
        {
            child = &c;
            break;
        }
    }

    if (child == NULL && create)
    {
        child = &*node.addChild();
    }

    if (pos == string::npos || child == NULL)
    {
        return child;
    }
    else
    {
        return resolve(*child, path.substr(pos+1));
    }
}

void Config::read(const string& file)
{
    string cwd;
    if (file.find('/') == string::npos)
    {
        cwd = ".";
    }
    else
    {
        cwd = file.substr(0, file.rfind('/'));
    }

    ifstream ifs(file.c_str());
    read(cwd, ifs);
}

void Config::read(const string& cwd, istream& is)
{
    root.reset(new Node());

    ptr_vector<Node> current = {root.get(), root.get()};

    Tokenizer t(is);
    string token;

    while (t.next(token))
    {
        for (string::size_type i = 0;i < token.size();)
        {
            if (token[i] == '{')
            {
                current.push_back(current.pback());
                i++;
            }
            else if (token[i] == '}')
            {
                if (current.size() < 2) throw FormatError("Too many }'s", t.getLine());
                current.pop_back();
                i++;
            }
            else if (token[i] == ',')
            {
                assert(current.size() >= 2);
                current.ptr(current.size()-1) = current.ptr(current.size()-2);
                i++;
            }
            else
            {
                string::size_type pos = token.find_first_of("{},", i);
                string::size_type len = (pos == string::npos ? token.size() : pos)-i;

                string node = token.substr(i, len);

                if (node == "include")
                {
                    if (i+len != token.size())
                    {
                        throw FormatError("\"include\" must be immediately followed by a filename", t.getLine());
                    }

                    t.next(token);
                    pos = token.find_first_of("{},", i);
                    len = (pos == string::npos ? token.size() : pos);

                    if (pos == 0)
                    {
                        throw FormatError("\"include\" must be immediately followed by a filename", t.getLine());
                    }

                    string fname = token.substr(0, len);
                    if (fname[0] != '/') fname = cwd+'/'+fname;

                    Config leaf; leaf.read(fname);
                    current.back().children.splice(current.back().children.end(), leaf.root->children);
                    for (Node& c : current.back().children)
                    {
                        c.parent = current.pback();
                    }
                }
                else
                {
                    current.pback() = &*current.back().addChild(node);
                }

                i += len;
            }
        }
    }

    if (current.size() != 2) throw FormatError("Too few }'s", t.getLine());
}

void Config::write(const string& file) const
{
    ofstream ofs(file.c_str());
    write(ofs);
}

void Config::write(ostream& os) const
{
    os << *root;
}

vector<pair<string,Config>> Config::find(const string& pattern)
{
    vector<pair<string,Config>> v;
    find(*root, "", pattern, v);
    return v;
}

bool Config::matchesWithGlobs(const char* str, const char* pat) const
{
    const char* s = str;
    const char* p = pat;

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

string Config::Extractor<string>::extract(Node& node, int which)
{
    auto i = node.getChild(which);
    if (i == node.children.end()) throw NoValueError(node.fullName());
    return i->data;
}

Config Config::Extractor<Config>::extract(Node& node, int which)
{
    auto i = node.parent->children.pbegin();
    while (i->get() != &node) ++i;
    return Config(*i);
}

template<>
void Config::set<string>(const string& path, const string& data, int which, bool create)
{
    Node* n = resolve(*root, path, create);
    if (n == NULL) throw EntryNotFoundError(path);

    auto i = n->getChild(which);
    if (i == n->children.end())
    {
        if (!create) throw NoValueError(path);
        i = n->addChild();
    }

    i->data = data;
}

void Schema::apply(Config& config) const
{
    apply(*root, *config.root);
}

void Schema::apply(const Node& schema, Node& root) const
{
    if (!isPrimitive(schema))
    {
        int minwild = 0;
        int maxwild = 0;
        for (const Node& s : schema.children)
        {
            int slen = s.data.length();
            if (s.data[slen-1] == '*')
            {
                if (s.data == "*")
                {
                    maxwild = -1;
                    continue;
                }

                for (Node& r : root.children)
                {
                    if (s.data.compare(0, slen-1, r.data, 0, slen-1) == 0)
                    {
                        apply(s, r);
                    }
                }
            }
            else if (s.data[slen-1] == '+')
            {
                if (s.data == "*+")
                {
                    minwild++;
                    maxwild = -1;
                    continue;
                }

                bool found = false;
                for (Node& r : root.children)
                {
                    if (s.data.compare(0, slen-1, r.data, 0, slen-1) == 0)
                    {
                        found = true;
                        apply(s, r);
                    }
                }
                if (!found)
                    throw SchemaValidationError("required node not found: " + s.data);
            }
            else if (s.data[slen-1] == '?')
            {
                if (s.data == "*?")
                {
                    if (maxwild != -1) maxwild++;
                    continue;
                }

                bool found = false;
                for (Node& r : root.children)
                {
                    if (s.data.compare(0, slen-1, r.data, 0, slen-1) == 0)
                    {
                        if (found)
                            throw SchemaValidationError("multiple copies of one-time node found: " + s.data);
                        found = true;
                        apply(s, r);
                    }
                }

                if (!found)
                {
                    Node& r = *root.addChild(s.data.substr(0, slen-1));
                    apply(s, r);
                }
            }
            else
            {
                bool found = false;
                for (Node& r : root.children)
                {
                    if (s.data == r.data)
                    {
                        if (found)
                            throw SchemaValidationError("multiple copies of one-time node found: " + s.data);
                        found = true;
                        apply(s, r);
                    }
                }
                if (!found)
                    throw SchemaValidationError("required node not found: " + s.data);
            }
        }

        int nwild = 0;
        for (Node& r : root.children)
        {
            bool found = false;
            for (const Node& s : schema.children)
            {
                int slen = s.data.length();
                if (s.data == "*" || s.data == "*+" || s.data == "*?") continue;
                if (s.data[slen-1] == '*' || s.data[slen-1] == '+' || s.data[slen-1] == '?') slen--;
                if (s.data.compare(0, slen, r.data, 0, slen) == 0)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                nwild++;

                if (maxwild != -1 && nwild > maxwild)
                {
                    throw SchemaValidationError("The node " + r.data + " is not a valid child of " + root.fullName());
                }
            }
        }

        if (nwild < minwild)
        {
            throw SchemaValidationError("The are too few wildcard nodes on " + root.fullName());
        }
    }
    else //primitive
    {
        const Node& s = schema.children.front();

        if (root.children.empty()) //value not specified
        {
            if (hasDefault(s))
            {
                root.addChild(s.children.front().data);
            }
            else //no default
            {
                int slen = schema.data.length();
                if (schema.data[slen-1] == '?')
                {
                    root.parent->removeChild(root);
                }
                else
                {
                    throw SchemaValidationError("no value specified and no default: " + root.data);
                }
            }
        }
        else //value specified
        {
            if (root.children.size() != 1)
                throw SchemaValidationError("multiple values for primitive: " + root.data);

            Node& r = root.children.front();

            if (s.data == "int")
            {
                if (!isInt(r.data))
                    throw SchemaValidationError("data is not an integer: " + root.data);
            }
            else if (s.data == "double")
            {
                if (!isDouble(r.data))
                    throw SchemaValidationError("data is not a double: " + root.data);
            }
            else if (s.data == "bool")
            {
                if (!isBool(r.data))
                    throw SchemaValidationError("data is not a boolean: " + root.data);
            }
            else if (s.data == "enum")
            {
                bool found = false;
                for (const Node& e : s.children)
                {
                    if (r.data == e.data)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    throw SchemaValidationError("value is not in enum: " + root.data);
            }
            else if (s.data == "string")
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

bool Schema::isPrimitive(const Node& schema) const
{
    if (schema.children.size() != 1) return false;

    if (schema.children.front().data == "int" ||
        schema.children.front().data == "double" ||
        schema.children.front().data == "bool" ||
        schema.children.front().data == "enum" ||
        schema.children.front().data == "string") return true;

    return false;
}

bool Schema::hasDefault(const Node& schema) const
{
    if (schema.children.empty()) return false;
    return true;
}

}
}
