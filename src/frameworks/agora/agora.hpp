#ifndef _AQUARIUS_FRAMEWORKS_AGORA_AGORA_HPP_
#define _AQUARIUS_FRAMEWORKS_AGORA_AGORA_HPP_

#include "frameworks/util.hpp"

#define REGISTER_MARKETPLACE(marketplace,...) \
template <> bool agora::StaticInitializer<marketplace>::initialized = \
    agora::Agora::registerMarketplace<marketplace>(__VA_ARGS__)

#define REGISTER_VENDOR(marketplace,vendor,type,...) \
template <> bool agora::StaticInitializer<type>::initialized = \
    agora::Marketplace<marketplace>::registerVendor<vendor,type>(__VA_ARGS__)

namespace aquarius
{
namespace agora
{

template <class T>
class StaticInitializer
{
    protected:
        static bool initialized;
};

class Vendor
{
    public:
        Vendor(const string& name) : name(name) {}

        virtual ~Vendor() {}

        const string& getName() const { return name; }

    protected:
        string name;
};

template <class M>
class Marketplace
{
    public:
        virtual ~Marketplace() {}

        template <typename V, typename X, typename... Args>
        static bool registerVendor(Args&&... args)
        {
            vendors().emplace_back(new X(forward<Args>(args)...));
            return true;
        }

        template <typename V, typename... Args>
        bool hasVendor(Args&&... args) const
        {
            for (auto& v : vendors())
            {
                if (v.type() != typeid(V)) continue;
                V& vc = v.get<V>();
                if (vc.check(forward<Args>(args)...)) return true;
            }

            return false;
        }

        template <typename V, typename... Args>
        V& getVendor(Args&&... args) const
        {
            ptr_vector<V> matches;

            for (auto& v : vendors())
            {
                if (v.type() != typeid(V)) continue;
                V& vc = v.get<V>();
                if (vc.check(forward<Args>(args)...)) matches.push_back(&vc);
            }

            if (matches.empty())
            {
                throw logic_error("No vendor matching requirements");
            }
            else
            {
                return static_cast<const M&>(*this).getBestVendor(matches);
            }
        }

    protected:
        static list<any>& vendors()
        {
            static list<any> vendors_;
            return vendors_;
        }
};

class Agora
{
    public:
        template <typename M, typename... Args>
        static bool registerMarketplace(const string& name, Args&&... args)
        {
            marketplaces()[name] = make_any<M>(forward<Args>(args)...);
            return true;
        }

        template <typename M>
        const M& getMarketplace(const string& name) const
        {
            return marketplaces().at(name).get<M>();
        }

    protected:
        static map<string,any>& marketplaces()
        {
            static map<string,any> marketplaces_;
            return marketplaces_;
        }
};

}
}

#endif
