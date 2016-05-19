#ifndef _AQUARIUS_FRAMEWORKS_AGORA_AGORA_HPP_
#define _AQUARIUS_FRAMEWORKS_AGORA_AGORA_HPP_

#include "frameworks/util.hpp"

#define REGISTER_MARKETPLACE(type,name) \
template <> bool aquarius::agora::StaticInitializer<type>::initialized = \
    aquarius::agora::Agora::registerMarketplace(name, type())

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

class Marketplace
{
    public:
        virtual ~Marketplace() {}
};

class Agora
{
    public:
        template <typename T>
        static bool registerMarketplace(const string& name, const T& marketplace)
        {
            marketplaces()[name] = make_any<T>(marketplace);
            return true;
        }

        template <typename T>
        const T& getMarketplace(const string& name) const
        {
            return marketplaces().at(name).get<T>();
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
