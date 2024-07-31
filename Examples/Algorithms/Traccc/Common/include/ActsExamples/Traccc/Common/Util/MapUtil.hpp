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

template <typename input_container_t, typename output_container_t, typename hash = std::hash<typename input_container_t::value_type>, typename equal_to = std::equal_to<typename input_container_t::value_type>>
struct ConversionData{
    using K = typename input_container_t::value_type;
    using V = typename output_container_t::value_type;

    input_container_t* inputContainer;
    std::unordered_map<K, std::size_t, hash, equal_to> keyToIdxMap;
    output_container_t* outputContainer;

    std::size_t valueToIndex(K& inputValue) const{
        return keyToIdxMap.at(inputValue);
    }

    std::size_t indexToIndex(std::size_t index) const{
        return valueToIndex(inputContainer->at(index));
    }

    auto& indexToValue(std::size_t index) const{
        return outputContainer->at(indexToIndex(index));
    }

    auto& valueToValue(K& inputValue) const{
        return outputContainer->at(valueToIndex(inputValue));
    }
};

template <typename input_container_t, typename hash = std::hash<typename input_container_t::value_type>, typename equal_to = std::equal_to<typename input_container_t::value_type>>
auto create1To1(input_container_t& inputs){
    using K = typename input_container_t::value_type;
    std::unordered_map<K, std::size_t, Hasher<K>, Equals<K>> keyToIdxMap;
    for (std::size_t idx = 0; idx < inputs.size(); idx++){
        keyToIdxMap.emplace(std::piecewise_construct, std::forward_as_tuple(inputs[idx]), std::forward_as_tuple(idx));
    }
    return keyToIdxMap;
}

template <typename input_container_t, typename output_container_t, typename hash = std::hash<typename input_container_t::value_type>, typename equal_to = std::equal_to<typename input_container_t::value_type>, template <typename, typename> class index_map_t>
auto createFromIndexMap(input_container_t& inputs, output_container_t& outputs, index_map_t<std::size_t, std::size_t>& indexMap){
    using K = typename input_container_t::value_type;
    std::unordered_map<K, std::size_t, Hasher<K>, Equals<K>> keyToIdxMap;
    for (std::size_t inputIdx = 0; inputIdx < inputs.size(); inputIdx++){
        auto outputIdx = indexMap.at(inputIdx);
        keyToIdxMap.emplace(std::piecewise_construct, std::forward_as_tuple(inputs[inputIdx]), std::forward_as_tuple(outputIdx));
    }
    return ConversionData{&inputs, std::move(keyToIdxMap), &outputs};
}

template <typename V, typename input_container_t, typename output_container_t, typename hash = std::hash<typename input_container_t::value_type>, typename equal_to = std::equal_to<typename input_container_t::value_type>, typename fn_t>
auto convert(input_container_t& inputs, fn_t fn, output_container_t& outputs){
    using K = typename input_container_t::value_type;
    std::transform(inputs.begin(), inputs.end(), std::back_inserter(outputs), fn);

    ConversionData conv{
        &inputs,
        create1To1<input_container_t, hash, equal_to>(inputs),
        &outputs
    };

    return conv;
}

}