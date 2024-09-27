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
#include <iterator>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace ActsExamples::Traccc::Common::Util {

/// @brief Class representing a mapping between two collections.
/// Functions for mapping from/to value and indexes are provided.
/// The class is useful as it allows you to obtain a reference to
/// the item mapped in the output collection.
template <typename It1, typename It2, typename Hash, typename Equals>
class ConversionData {
 public:
  using K = typename std::iterator_traits<It1>::value_type;
  using V = typename std::iterator_traits<It2>::value_type;

  ConversionData() = delete;

  const K& inputAt(std::size_t i) {
    It1 it = input_begin;
    std::advance(it, i);
    return *it;
  }

  std::size_t valueToIndex(const K& inputValue) const {
    return keyToIdxMap.at(inputValue);
  }

  std::size_t indexToIndex(std::size_t index) const {
    It1 i = input_begin;
    std::advance(i, index);
    return valueToIndex(*i);
  }

  std::conditional_t<std::is_lvalue_reference_v<V>, const V&, V> indexToValue(
      std::size_t index) const {
    It2 i = output_begin;
    std::advance(i, indexToIndex(index));
    return *i;
  }

  std::conditional_t<std::is_lvalue_reference_v<V>, const V&, V> valueToValue(
      const K& inputValue) const {
    It2 i = output_begin;
    std::advance(i, valueToIndex(inputValue));
    return *i;
  }

  std::size_t size() { return keyToIdxMap.size(); }

  ConversionData(It1 i_beg, It1 i_end, It2 o_beg, It2 o_end,
                 std::unordered_map<K, std::size_t, Hash, Equals>&& map)
      : input_begin(i_beg),
        input_end(i_end),
        output_begin(o_beg),
        output_end(o_end),
        keyToIdxMap(std::move(map)) {}

  template <typename _Hash, typename _Equals>
  ConversionData<It2, It1, _Hash, _Equals> invert(_Hash&& new_hash,
                                                  _Equals&& new_equals) {
    std::unordered_map<typename std::iterator_traits<It2>::value_type,
                       std::size_t, _Hash, _Equals>
        map(64, std::forward<_Hash>(new_hash),
            std::forward<_Equals>(new_equals));

    for (std::size_t idx = 0; idx < size(); idx++) {
      auto value = valueToValue(inputAt(idx));
      map.emplace(std::piecewise_construct, std::forward_as_tuple(value),
                  std::forward_as_tuple(idx));
    }

    return ConversionData<It2, It1, _Hash, _Equals>{
        output_begin, output_end, input_begin, input_end, std::move(map)};
  }

 private:
  It1 input_begin;
  It1 input_end;
  It2 output_begin;
  It2 output_end;
  std::unordered_map<K, std::size_t, Hash, Equals> keyToIdxMap;
};

/// @brief Creates an instance of the ConversionData class where the given indexMap parameter determines how the elements
/// depending on their index in the input container map to the elements of the
/// output container depending on their index.
template <std::forward_iterator It1, std::forward_iterator It2, typename Hash,
          typename Equals>
auto makeConversionFromIndexMap(
    It1 i_beg, It1 i_end, It2 o_beg, It2 o_end,
    const std::map<std::size_t, std::size_t>& indexMap, Hash&& hash = {},
    Equals&& equals = {}) {
  std::unordered_map<typename std::iterator_traits<It1>::value_type,
                     std::size_t, Hash, Equals>
      map(64, std::forward<Hash>(hash), std::forward<Equals>(equals));

  for (It1 i = i_beg; i != i_end; i++) {
    auto outputIdx = indexMap.at(std::distance(i_beg, i));
    map.emplace(std::piecewise_construct, std::forward_as_tuple(*i),
                std::forward_as_tuple(outputIdx));
  }
  return ConversionData<It1, It2, Hash, Equals>(i_beg, i_end, o_beg, o_end,
                                                std::move(map));
}

template <typename It1, typename It2, typename Hash, typename Equals>
ConversionData<It1, It2, Hash, Equals> makeConversionOneToOne(
    It1 i_beg, It1 i_end, It2 o_beg, It2 o_end, Hash&& hash = {},
    Equals&& equals = {}) {
  std::unordered_map<typename std::iterator_traits<It1>::value_type,
                     std::size_t, Hash, Equals>
      map(64, std::forward<Hash>(hash), std::forward<Equals>(equals));

  for (It1 i = i_beg; i != i_end; i++) {
    map.emplace(std::piecewise_construct, std::forward_as_tuple(*i),
                std::forward_as_tuple(std::distance(i_beg, i)));
  }

  return ConversionData<It1, It2, Hash, Equals>(i_beg, i_end, o_beg, o_end,
                                                std::move(map));
}

/// @brief Creates an instance of the ConversionData class where the mapping and output container is generated
/// by the provided function.
template <typename hash_t, typename equals_t, typename input_container_t,
          typename output_container_t, typename fn_t>
auto convert(input_container_t& inputs, fn_t fn, output_container_t& outputs) {
  std::ranges::transform(inputs, std::back_inserter(outputs), fn);

  return makeConversionOneToOne(inputs.cbegin(), inputs.cend(),
                                outputs.cbegin(), outputs.cend(), hash_t{},
                                equals_t{});
}

}  // namespace ActsExamples::Traccc::Common::Util
