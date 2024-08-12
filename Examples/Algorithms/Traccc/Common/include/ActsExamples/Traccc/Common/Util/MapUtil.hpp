// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

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


/// @brief Class representing a mapping between two collections.
/// Functions for mapping from/to value and indexes are provided.
/// The class is useful as it allows you to obtain a reference to
/// the item mapped in the output collection.
template <typename hash_t, typename equals_t, typename input_container_t, typename output_container_t>
struct ConversionData{
    using K = typename input_container_t::value_type;
    using V = typename output_container_t::value_type;

    input_container_t* inputContainer;
    std::unordered_map<K, std::size_t, hash_t, equals_t> keyToIdxMap;
    output_container_t* outputContainer;

    std::size_t valueToIndex(const K& inputValue) const{
        return keyToIdxMap.at(inputValue);
    }

    std::size_t indexToIndex(std::size_t index) const{
        return valueToIndex(inputContainer->at(index));
    }

    auto& indexToValue(std::size_t index) const{
        return outputContainer->at(indexToIndex(index));
    }

    auto& valueToValue(const K& inputValue) const{
        return outputContainer->at(valueToIndex(inputValue));
    }

    auto size(){
        return keyToIdxMap.size();
    }
};

/// @returns the inverse of a ConversionData instance.
template <typename hash_t, typename equals_t, typename conversion_data_t>
auto inverse(conversion_data_t& conv){
    std::unordered_map<typename conversion_data_t::V, std::size_t, hash_t, equals_t> keyToIdxMap;

    for (std::size_t idx = 0; idx < conv.size(); idx++){
        auto key = conv.inputContainer->at(idx);
        auto value = conv.valueToValue(key);
        keyToIdxMap.emplace(std::piecewise_construct, std::forward_as_tuple(value), std::forward_as_tuple(idx));
    }

    return ConversionData{
        conv.outputContainer,
        keyToIdxMap,
        conv.inputContainer
    };
}

/// @brief Creates a keyToIdxMap map each value in the input container to their index.
/// This represents a mapping from a value at index n in the input container to index n in the output container.
template <typename hash_t, typename equals_t, typename input_container_t>
auto create1To1(input_container_t& inputs){
    using K = typename input_container_t::value_type;
    std::unordered_map<K, std::size_t, hash_t, equals_t> keyToIdxMap;
    for (std::size_t idx = 0; idx < inputs.size(); idx++){
        keyToIdxMap.emplace(std::piecewise_construct, std::forward_as_tuple(inputs[idx]), std::forward_as_tuple(idx));
    }
    return keyToIdxMap;
}

/// @brief Creates an instance of the ConversionData class where the given indexMap parameter determines how the elements
/// depending on their index in the input container map to the elements of the output container depending on their index.
template <typename hash_t, typename equals_t, typename input_container_t, typename output_container_t, template <typename, typename> class index_map_t>
auto createFromIndexMap(input_container_t& inputs, output_container_t& outputs, index_map_t<std::size_t, std::size_t>& indexMap){
    using K = typename input_container_t::value_type;
    std::unordered_map<K, std::size_t, hash_t, equals_t> keyToIdxMap;
    for (std::size_t inputIdx = 0; inputIdx < inputs.size(); inputIdx++){
        auto outputIdx = indexMap.at(inputIdx);
        keyToIdxMap.emplace(std::piecewise_construct, std::forward_as_tuple(inputs[inputIdx]), std::forward_as_tuple(outputIdx));
    }
    return ConversionData{&inputs, std::move(keyToIdxMap), &outputs};
}

/// @brief Creates an instance of the ConversionData class where the mapping and output container is generated
/// by the provided function.
template <typename hash_t, typename equals_t, typename input_container_t, typename output_container_t, typename fn_t>
auto convert(input_container_t& inputs, fn_t fn, output_container_t& outputs){
    using K = typename input_container_t::value_type;
    std::transform(inputs.begin(), inputs.end(), std::back_inserter(outputs), fn);

    ConversionData conv{
        &inputs,
        create1To1<hash_t, equals_t>(inputs),
        &outputs
    };

    return conv;
}

}