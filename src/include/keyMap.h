#pragma once
#include <typeinfo>
#include <typeindex>
#include <tuple>
#include <iostream>

namespace namespaceKeyMap
{

    struct KeyMap
    {
        std::type_index typeFrom = typeid(void);
        std::type_index typeTo = typeid(void);
    };

    inline bool operator==(const KeyMap &lhs, const KeyMap &rhs)
    {
        auto lhsTuple = std::make_tuple(lhs.typeFrom, lhs.typeTo);
        auto rhsTuple = std::make_tuple(rhs.typeFrom, rhs.typeTo);
        return lhsTuple == rhsTuple;
    }

    struct Hash_fn
    {
        std::size_t operator()(const KeyMap &obj) const
        {
            constexpr size_t shift = 7;
            std::size_t h1 = std::hash<std::type_index>()(obj.typeFrom);
            std::size_t h2 = std::hash<std::type_index>()(obj.typeTo);
            return h1 ^ (h2 >> shift);
        }
    };
}