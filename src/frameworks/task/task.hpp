#ifndef _AQUARIUS_FRAMEWORKS_TASK_TASK_HPP_
#define _AQUARIUS_FRAMEWORKS_TASK_TASK_HPP_

#include "frameworks/util.hpp"
#include "frameworks/time.hpp"

#include "config.hpp"

#define REGISTER_TASK(type,name,...) \
template <> bool aquarius::task::TaskFactory<type>::initialized = \
    aquarius::task::Task::registerTask(name, aquarius::task::Schema(__VA_ARGS__), \
                                       aquarius::task::TaskFactory<type>::create)

namespace aquarius
{
namespace task
{

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
        shared_ptr<any> data;
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
            data = make_shared<any>(resource);
            return *resource;
        }

        template <typename T> decay_t<T>& put(T&& resource)
        {
            data = make_shared<any>(forward(resource));
            return resource;
        }

        template <typename T> T& get()
        {
            if (!exists()) throw logic_error("Product " + name + " does not exist.");
            return data->get<T>();
        }

        bool exists() const { return data; }

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
        typedef unique_ptr<Task> (*factory_func)(const string&, Config&);

        string type;
        string name;
        vector<Product> products;
        vector<Product> temporaries;
        Config config;

        static map<string,tuple<Schema,factory_func>>& tasks();

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
            for (auto& t : temporaries)
            {
                if (t.getName() == name)
                {
                    return t.put(resource);
                }
            }

            temporaries.emplace_back("temporary", name);
            return temporaries.back().put(resource);
        }

        template <typename T> T& gettmp(const string& name)
        {
            for (auto& t : temporaries)
            {
                if (t.getName() == name)
                {
                    return t.get<T>();
                }
            }

            throw logic_error("Temporary " + name + " not found on task " + this->name);
        }

        ostream& log(const Arena& arena);

        ostream& warn(const Arena& arena);

        ostream& error(const Arena& arena);

    public:
        Task(const string& name, Config& config);

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

        static bool registerTask(const string& name, Schema&& schema, factory_func create);

        static const Schema& getSchema(const string& name);

        Config& getConfig()
        {
            return config;
        }

        const Config& getConfig() const
        {
            return config;
        }

        template <typename T> T& put(const string& name, T* resource)
        {
            return getProduct(name).put(resource);
        }

        template <typename T> decay_t<T>& put(const string& name, T&& resource)
        {
            return put(name, new decay_t<T>(forward<T>(resource)));
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

        static unique_ptr<Task> createTask(const string& type, const string& name, Config& config);
};

template <class T>
class TaskFactory
{
    friend class Task;

    protected:
        static bool initialized;

        static unique_ptr<Task> create(const string& name, Config& config)
        {
            return unique_ptr<Task>(new T(name, config));
        }
};

class TaskDAG
{
    protected:
        unique_vector<Task> tasks;
        vector<tuple<string,string,Config>> usings;

        void parseTasks(const string& context, Config& config);

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
        CompareScalars(const string& name, Config& config);

        bool run(TaskDAG& dag, const Arena& arena);
};

}
}

#endif
