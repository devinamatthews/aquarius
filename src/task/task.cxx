#include "task.hpp"

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

streamsize Logger::LogToStreamBuffer::xsputn (const char* s, streamsize n)
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

Logger::LogToStreamBuffer::LogToStreamBuffer(ostream& os)
: os(os) {}

int Logger::SelfDestructBuffer::sync()
{
    LogToStreamBuffer::sync();
    abort();
    return -1;
}

Logger::SelfDestructBuffer::SelfDestructBuffer()
: LogToStreamBuffer(cerr) {}

streamsize Logger::NullBuffer::xsputn (const char* s, streamsize n)
{
    return n;
}

int Logger::NullBuffer::overflow (int c)
{
    if (c == EOF) return EOF;
    else return traits_type::to_int_type(c);
}

Logger::NullStream::NullStream()
: ostream(new NullBuffer()) {}

Logger::NullStream::~NullStream()
{
    delete rdbuf();
}

Logger::SelfDestructStream::SelfDestructStream()
: ostream(new SelfDestructBuffer()) {}

string Logger::dateTime()
{
    char buf[256];
    time_t t = ::time(NULL);
    tm* timeptr = localtime(&t);
    strftime(buf, 256, "%c", timeptr);
    return string(buf);
}

ostream& Logger::log(const Arena& arena)
{
    if (arena.rank == 0)
    {
        cout << dateTime() << ": ";
        return cout;
    }
    else
    {
        return nullstream;
    }
}

ostream& Logger::warn(const Arena& arena)
{
    if (arena.rank == 0)
    {
        cerr << dateTime() << ": warning: ";
        return cerr;
    }
    else
    {
        return nullstream;
    }
}

ostream& Logger::error(const Arena& arena)
{
    arena.comm().Barrier();
    if (arena.rank == 0)
    {
        sddstream << dateTime() << ": error: ";
        return sddstream;
        //cout << dateTime() << ": error: ";
        //return cout;
    }
    else
    {
        return nullstream;
    }
}

Requirement::Requirement(const string& type, const string& name)
: type(type), name(name) {}

void Requirement::fulfil(const Product& product)
{
    this->product.set(new Product(product));
    *this->product->used = true;
}

bool Requirement::exists() const
{
    if (!product) throw logic_error("Requirement " + name + " not fulfilled");
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

void Product::addRequirement(Requirement&& req)
{
    requirements->push_back(forward<Requirement>(req));
}

void Product::addRequirements(const vector<Requirement>& reqs)
{
    requirements->insert(requirements->end(), reqs.begin(), reqs.end());
}

void Product::addRequirements(vector<Requirement>&& reqs)
{
    requirements->reserve(requirements->size()+reqs.size());
    for (Requirement& req : reqs) requirements->push_back(move(req));
}

template <> string Task::type_string<float>() { return "<float>"; }
template <> string Task::type_string<double>() { return ""; }
template <> string Task::type_string<complex<float>>() { return "<scomplex>"; }
template <> string Task::type_string<complex<double>>() { return "<dcomplex>"; }

Task::Task(const string& name, Config& config)
: name(name), config(config.clone()) {}

map<string,tuple<Schema,Task::factory_func>>& Task::tasks()
{
    static map<string,tuple<Schema,Task::factory_func>> tasks_;
    return tasks_;
}

ostream& Task::log(const Arena& arena)
{
    return Logger::log(arena) << name << ": ";
}

ostream& Task::warn(const Arena& arena)
{
    return Logger::warn(arena) << name << ": ";
}

ostream& Task::error(const Arena& arena)
{
    return Logger::error(arena) << name << ": ";
}

bool Task::registerTask(const string& name, Schema&& schema, factory_func create)
{
    tasks()[name] = make_tuple(move(schema), create);
    return true;
}

const Schema& Task::getSchema(const string& name)
{
    auto i = tasks().find(name);
    if (i == tasks().end()) throw logic_error("Cannot find schema for task " + name);
    return aquarius::get<0>(i->second);
}

Product& Task::getProduct(const string& name)
{
    for (auto& p : products)
    {
        if (p.getName() == name) return p;
    }
    throw logic_error("Product " + name + " not found on task " + this->name);
}

const Product& Task::getProduct(const string& name) const
{
    return const_cast<Task&>(*this).getProduct(name);
}

unique_ptr<Task> Task::createTask(const string& type, const string& name, input::Config& config)
{
    auto i = tasks().find(type);
    if (i == tasks().end()) throw logic_error("Task type " + type + " not found");
    aquarius::get<0>(i->second).apply(config);
    unique_ptr<Task> t = aquarius::get<1>(i->second)(name, config);
    t->type = type;
    return t;
}

void TaskDAG::parseTasks(const string& context, Config& input)
{
    for (auto& i : input.find<string>("section"))
    {
        Config section = input.get("section." + i.second);
        parseTasks(context+i.second+".", section);
        input.remove("section");
    }

    for (auto& i : input.find("*"))
    {
        string type = i.first;
        Config config = i.second.clone();
        Config uconfig = i.second.clone();

        while (config.exists("using")) config.remove("using");

        if (config.exists("name"))
        {
            string name = config.get<string>("name");
            if (name.find_first_of(".:") != string::npos)
                Logger::error(Arena()) << "Task names may not contain any of '.:' (" << name << ")" << endl;
            name = context+name;
            config.remove("name");

            addTask(world(), Task::createTask(type, name, config));
        }
        else
        {
            addTask(world(), type, context, config);
        }

        while (uconfig.exists("using"))
        {
            string u = uconfig.get<string>("using");
            Config c = uconfig.get("using");
            usings.emplace_back(tasks.back().getName(), u, c);
        }
    }
}

TaskDAG::TaskDAG(const string& file)
{
    ifstream ifs(file);
    Config input(ifs);
    parseTasks("", input);
}

Task& TaskDAG::addTask(const Arena& arena, const string& type, const string& context, Config& config)
{
    int num = 0;
    for (Task& t : tasks)
    {
        if (t.getName().compare(0, context.size(), context) == 0 &&
            t.getType() == type) num++;
    }

    string name = context + (!context.empty() && context.back() != '.' ? "." : "") +
        type + (num == 0 ? "" : str(num));

    return addTask(arena, Task::createTask(type, name, config));
}

Task& TaskDAG::addTask(const Arena& arena, unique_ptr<Task>&& task)
{
    for (Task& t : tasks)
    {
        if (t.getName() == task->getName())
            Logger::error(arena) << "More than one task with name " << task->getName() << endl;
    }

    cout << "Adding task " << task->getName() << endl;

    Task& t1 = *task;

    string context1;
    size_t sep = t1.getName().find_last_of(".");
    if (sep!=string::npos)
    {
        context1 = t1.getName().substr(0, sep+1);
    }

    for (Task& t2 : tasks)
    {
        string context2;
        size_t sep = t2.getName().find_last_of(".");
        if (sep != string::npos)
        {
            context2 = t2.getName().substr(0, sep+1);
        }

        bool t1_parent_of_t2 = (context2.compare(0, context1.size(), context1) == 0);
        bool t2_parent_of_t1 = (context1.compare(0, context2.size(), context2) == 0);

        for (Product& p1 : t1.getProducts())
        {
            for (Product& p2 : t2.getProducts())
            {
                for (Requirement& r2 : p2.getRequirements())
                {
                    if (t1_parent_of_t2)
                    {
                        cout << "Checking requirement " << r2.getName() << " of task " << t2.getName() <<
                            " against product " << p1.getName() << " of task " << t1.getName() << endl;
                        if (!r2.isFulfilled() && r2.getType() == p1.getType())
                        {
                            cout << "yes" << endl;
                            r2.fulfil(p1);
                        }
                    }
                }
            }

            for (Requirement& r1 : p1.getRequirements())
            {
                for (Product& p2 : t2.getProducts())
                {
                    if (t2_parent_of_t1)
                    {
                        cout << "Checking requirement " << r1.getName() << " of task " << t1.getName() <<
                            " against product " << p2.getName() << " of task " << t2.getName() << endl;
                        if (!r1.isFulfilled() && r1.getType() == p2.getType())
                        {
                            cout << "yes" << endl;
                            r1.fulfil(p2);
                        }
                    }

                    for (Requirement& r2 : p2.getRequirements())
                    {
                        if (t2_parent_of_t1)
                        {
                            cout << "Checking requirement " << r1.getName() << " of task " << t1.getName() <<
                                " against requiremet " << r2.getName() << " of task " << t2.getName() << endl;
                            if (!r1.isFulfilled() && r2.isFulfilled() && r1.getType() == r2.getType())
                            {
                                cout << "yes" << endl;
                                r1.fulfil(r2.get());
                            }
                        }

                        if (t1_parent_of_t2)
                        {
                            cout << "Checking requirement " << r2.getName() << " of task " << t2.getName() <<
                                " against requirement " << r1.getName() << " of task " << t1.getName() << endl;
                            if (!r2.isFulfilled() && r1.isFulfilled() && r2.getType() == r1.getType())
                            {
                                cout << "yes" << endl;
                                r2.fulfil(r1.get());
                            }
                        }
                    }
                }
            }
        }
    }

    tasks.push_back(move(task));
    return tasks.back();
}

void TaskDAG::satisfyExplicitRequirements(const Arena& world)
{
    /*
     * Hook up explicitly fulfilled requirements.
     */
    for (auto& u : usings)
    {
        string context;
        const string& name = get<1>(u);
        Config& config = get<2>(u);

        Task *t1 = NULL;
        for (Task& t : tasks) if (t.getName() == get<0>(u)) t1 = &t;
        assert(t1);

        size_t sep = t1->getName().find_last_of(".");
        if (sep != string::npos)
        {
            context = t1->getName().substr(0, sep+1);
        }

        ptr_vector<Requirement> reqs;

        for (Product& p : t1->getProducts())
        {
            for (Requirement& r : p.getRequirements())
            {
                if (r.getName() == name)
                {
                    if (!reqs.empty() && reqs.back().getType() != r.getType())
                        Logger::error(world) << "Multiple requirements named " << name << " with different types" << endl;
                    reqs.push_back(&r);
                }
            }
        }

        if (reqs.empty())
            Logger::error(world) << "No requirement " << name << " found on task " << t1->getName() << endl;

        Product p("double", name);
        Product* fulfiller = NULL;

        if (config.exists(name+".="))
        {
            if (reqs.back().getType() != "double")
                Logger::error(world) << "Attempting to specify a non-scalar requirement by value" << endl;
            p.put(new double(config.get<double>(name+".=")));
            fulfiller = &p;
        }
        else
        {
            string from = config.get<string>(name+".from");
            string task, req;

            size_t sep = from.find(':');
            if (sep == string::npos)
            {
                task = from;
                req = name;
            }
            else
            {
                task = from.substr(0, sep);
                req = from.substr(sep+1);
            }

            sep = task.find('.');
            if (sep == string::npos)
            {
                task = context+task;
            }

            if (task[0] == '.')
                task = task.substr(1);

            for (Task& t2 : tasks)
            {
                if (t2.getName() == task)
                {
                    for (Product& p2 : t2.getProducts())
                    {
                        if (p2.getName() == req)
                        {
                            if (p2.getType() != reqs.back().getType())
                                Logger::error(world) << "Product " << task << "." << req <<
                                    " is wrong type for requirement " << t1->getName() << "." << req << endl;
                            fulfiller = &p2;
                        }
                    }

                    if (fulfiller == NULL)
                        Logger::error(world) << "Product " << req << " not found on task " << task << endl;
                }
            }

            if (fulfiller == NULL)
                Logger::error(world) << "Task " << task << " not found" << endl;
        }

        for (Requirement& r : reqs)
        {
            r.fulfil(*fulfiller);
        }
    }
}

void TaskDAG::execute(const Arena& world)
{
    satisfyExplicitRequirements(world);

    //TODO: check for cycles

    /*
     * Successively search for executable tasks
     */
    while (!tasks.empty())
    {
        bool ran_something = false;
        for (auto i = tasks.pbegin();i != tasks.pend();)
        {
            bool can_execute = true;
            Task& t = **i;

            for (Product& p : t.getProducts())
            {
                for (Requirement& r : p.getRequirements())
                {
                    if (!r.exists())
                    {
                        can_execute = false;
                        break;
                    }
                }
                if (!can_execute) break;
            }

            if (can_execute)
            {
                ran_something = true;

                Logger::log(world) << "Starting task: " << t.getName() << endl;
                Timer timer;

                bool success = true;
                bool done = false;
                string error;

                timer.start();
                //try
                //{
                    done = t.run(*this, world);
                //}
                //catch (runtime_error& e)
                //{
                //    success = false;
                //    error = e.what();
                //}
                timer.stop();

                double dt = timer.seconds(world);
                double gflops = timer.gflops(world);
                Logger::log(world) << "Finished task: " << t.getName() <<
                           " in " << fixed << setprecision(3) << dt << " s" << endl;
                Logger::log(world) << "Task: " << t.getName() <<
                           " achieved " << fixed << setprecision(3) << gflops << " Gflops/sec" << endl;

                if (!success)
                {
                    throw runtime_error(error);
                }

                if (done)
                {
                    for (Product& p : t.getProducts())
                    {
                        if (p.isUsed() && !p.exists())
                            Logger::error(world) << "Product " << p.getName() <<
                                                    " of task " << t.getName() <<
                                                    " was not successfully produced" << endl;
                    }

                    i = tasks.perase(i);
                }
            }
            else
            {
                ++i;
            }
        }

        if (!ran_something)
        {
            Logger::error(world) << "Some tasks were not executed due to missing dependencies" << endl;
        }
    }
}

CompareScalars::CompareScalars(const string& name, Config& config)
: Task(name, config)
{
    tolerance = config.get<double>("tolerance");

    vector<Requirement> reqs;
    reqs.push_back(Requirement("double", "val1"));
    reqs.push_back(Requirement("double", "val2"));
    addProduct(Product("bool", "match", reqs));
}

bool CompareScalars::run(TaskDAG& dag, const Arena& arena)
{
    double val1 = get<double>("val1");
    double val2 = get<double>("val2");

    bool match = aquarius::abs(val1-val2) < tolerance;

    if (match)
    {
        log(arena) << "passed" << endl;
    }
    else
    {
        error(arena) << "failed: " << fixed <<
                setprecision((int)(0.5-log10(tolerance))) << val1 << " vs " << val2 << endl;
    }

    put("match", new bool(match));

    return true;
}

REGISTER_TASK(CompareScalars,"compare","tolerance double");

}
}
