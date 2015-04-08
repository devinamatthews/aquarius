#ifndef _AQUARIUS_TASK_HPP_
#define _AQUARIUS_TASK_HPP_

#include "util/global.hpp"

#include "input/config.hpp"
#include "time/time.hpp"

namespace aquarius
{
namespace task
{

class Logger
{
    protected:
        class LogToStreamBuffer : public streambuf
        {
            protected:
                ostream& os;

                int sync();

                streamsize xsputn (const char* s, streamsize n);

                int overflow (int c = EOF);

            public:
                LogToStreamBuffer(ostream& os);
        };

        class SelfDestructBuffer : public LogToStreamBuffer
        {
            protected:
                int sync();

            public:
                SelfDestructBuffer();
        };

        class NullBuffer : public streambuf
        {
            protected:
                streamsize xsputn (const char* s, streamsize n);

                int overflow (int c = EOF);
        };

        class NullStream : public ostream
        {
            public:
                NullStream();

                ~NullStream();
        };

        class SelfDestructStream : public ostream
        {
            public:
                SelfDestructStream();
        };

        static string dateTime();

        static NullStream nullstream;
        static SelfDestructStream sddstream;

    public:
        static ostream& log(const Arena& arena);

        static ostream& warn(const Arena& arena);

        static ostream& error(const Arena& arena);
};

class Printer
{
    protected:
        ostream* out;

    public:
        Printer() : out(NULL) {}

        Printer(ostream& out) : out(&out) {}

        ostream& getOutputStream() { return *out; }

        void setOutputStream(ostream& os) { out = &os; }

        template <class T>
        Printer& operator<<(const T& t)
        {
            if (out != NULL) out << t;
            return *this;
        }
};

class Destructible
{
    public:
        virtual ~Destructible() {}
};

template <typename T>
class Resource : public Destructible
{
    private:
        Resource(const Resource& other);

        Resource& operator=(const Resource& other);

    public:
        T* const data;

        Resource(T* data) : data(data) {}

        ~Resource() { delete data; }
};

class Product;

class Requirement
{
    protected:
        string type;
        string name;
        global_ptr<Product> product;

    public:
        Requirement(const string& type, const string& name);

        const string& getName() const { return name; }

        const string getType() const { return type; }

        bool exists() const;

        void fulfil(const Product& product);

        bool isFulfilled() const { return product; }

        Product& get();
};

class Task;

class Product
{
    friend class Requirement;
    friend class Task;

    protected:
        string type;
        string name;
        global_ptr<Destructible> data;
        shared_ptr<vector<Requirement>> requirements;
        shared_ptr<bool> used;

    public:
        Product(const string& type, const string& name);

        Product(const string& type, const string& name, const vector<Requirement>& reqs);

        const string& getName() const { return name; }

        const string& getType() const { return type; }

        bool isUsed() const { return *used; }

        template <typename T> T& put(T* resource)
        {
            data.reset(new Resource<T>(resource));
            return *resource;
        }

        template <typename T> T& get()
        {
            if (!data) throw logic_error("Product " + name + " does not exist.");
            return *static_cast<Resource<T>&>(*data).data;
        }

        bool exists() const { return data; }

        void addRequirement(Requirement&& req);

        template <typename... Args>
        void addRequirement(Args&&... args)
        {
            requirements->emplace_back(forward<Args>(args)...);
        }

        void addRequirements(const vector<Requirement>& reqs);

        void addRequirements(vector<Requirement>&& reqs);

        vector<Requirement>& getRequirements() { return *requirements; }
};

class TaskDAG;

class Task
{
    protected:
        typedef unique_ptr<Task> (*factory_func)(const string&, input::Config&);

        string type;
        string name;
        vector<Product> products;
        vector<Product> temporaries;
        input::Config config;

        static map<string,tuple<input::Schema,factory_func>>& tasks();

        void addProduct(const Product& product)
        {
            products.push_back(product);
        }

        void addProduct(Product&& product)
        {
            products.push_back(move(product));
        }

        template <typename... Args>
        void addProduct(Args&&... args)
        {
            products.emplace_back(forward<Args>(args)...);
        }

        template <typename T> T& puttmp(const string& name, T* resource)
        {
            for (vector<Product>::iterator i = temporaries.begin();i != temporaries.end();++i)
            {
                if (i->getName() == name)
                {
                    return i->put(resource);
                }
            }

            temporaries.push_back(Product("temporary", name));
            return temporaries.back().put(resource);
        }

        template <typename T> T& gettmp(const string& name)
        {
            for (vector<Product>::iterator i = temporaries.begin();i != temporaries.end();++i)
            {
                if (i->getName() == name)
                {
                    return i->get<T>();
                }
            }

            throw logic_error("Temporary " + name + " not found on task " + this->name);
        }

        ostream& log(const Arena& arena);

        ostream& warn(const Arena& arena);

        ostream& error(const Arena& arena);

    public:
        Task(const string& name, input::Config& config);

        virtual ~Task()
        {
            for (Product& p : products)
            {
                p.requirements->clear();
            }
        }

        const string& getType() const { return type; }

        const string& getName() const { return name; }

        bool isUsed(const string& name) const { return getProduct(name).isUsed(); }

        template <typename T> static string type_string();

        static bool registerTask(const string& name, input::Schema&& schema, factory_func create);

        static const input::Schema& getSchema(const string& name);

        input::Config& getConfig()
        {
            return config;
        }

        const input::Config& getConfig() const
        {
            return config;
        }

        template <typename T> T& put(const string& name, T* resource)
        {
            return getProduct(name).put(resource);
        }

        template <typename T> T& get(const string& name)
        {
            for (Product& p : products)
            {
                if (p.getName() == name) return p.get<T>();

                for (Requirement& r : p.getRequirements())
                {
                    if (r.getName() == name) return r.get().get<T>();
                }
            }

            throw logic_error("Product " + name + " not found on task " + this->name + " or its dependencies");
        }

        Product& getProduct(const string& name);

        const Product& getProduct (const string& name) const;

        vector<Product>& getProducts() { return products; }

        const vector<Product>& getProducts() const { return products; }

        virtual bool run(TaskDAG& dag, const Arena& arena) = 0;

        static unique_ptr<Task> createTask(const string& type, const string& name, input::Config& config);
};

template <class T>
class TaskFactory
{
    friend class Task;

    protected:
        static bool initialized;

        static unique_ptr<Task> create(const string& name, input::Config& config)
        {
            return unique_ptr<Task>(new T(name, config));
        }
};

#define REGISTER_TASK(type,name,...) \
template <> bool aquarius::task::TaskFactory<type>::initialized = \
    aquarius::task::Task::registerTask(name, aquarius::input::Schema(__VA_ARGS__), aquarius::task::TaskFactory<type>::create)

class TaskDAG
{
    protected:
        unique_vector<Task> tasks;
        vector<tuple<string,string,input::Config>> usings;

        void parseTasks(const string& context, input::Config& config);

        void satisfyRemainingRequirements(const Arena& world);

        void satisfyExplicitRequirements(const Arena& world);

    public:
        TaskDAG() {}

        TaskDAG(const string& file);

        void schedule(unique_ptr<Task>&& task);

        void execute(Arena& world);
};

class CompareScalars : public Task
{
    protected:
        double tolerance;

    public:
        CompareScalars(const string& name, input::Config& config);

        bool run(TaskDAG& dag, const Arena& arena);
};

}
}

#endif
