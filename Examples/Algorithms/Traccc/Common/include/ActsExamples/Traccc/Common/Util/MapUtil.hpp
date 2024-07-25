// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Traccc/Common/Util/Functional.hpp"

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <type_traits>

namespace ActsExamples::Traccc::Common::Util {

template <typename K, typename V, typename hash = std::hash<K>, typename equal_to = std::equal_to<K>>
struct ReferenceMap{

    template<typename T1, typename T2>
    ReferenceMap(T1& keys, T2& values){
        for (std::size_t idx = 0; idx < keys.size(); idx++){
            map.emplace(std::piecewise_construct, std::forward_as_tuple(keys[idx]), std::forward_as_tuple(&values[idx]));
        }
    }
    
    V& at(const K& k) const {
        return *map.at(k);
    }

    std::unordered_map<K, V*, hash, equal_to> map;

};

template <typename K, typename V, typename hash = std::hash<K>, typename equal_to = std::equal_to<K>, typename T, typename fn_t>
auto convert(T& keys, fn_t fn){
    std::vector<V> values;
    std::transform(keys.begin(), keys.end(), std::back_inserter(values), fn);
    return std::make_tuple(std::move(values), ReferenceMap<K, V, Hasher<K>, Equals<K>>(keys, values));
}

}