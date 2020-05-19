// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/container/flat_map.hpp>
#include <cstddef>
#include <utility>

namespace FW {

/// Store elements that are identified by an index, e.g. in another container.
///
/// Each index can have zero or more associated elements. A typical case could
/// be to store all generating particles for a hit where the hit is identified
/// by its index in the hit container.
template <typename Value, typename Key = std::size_t>
using IndexMultimap = boost::container::flat_multimap<Key, Value>;

/// Invert the multimap, i.e. from a -> {b...} to b -> {a...}.
///
/// @note This assumes that the value in the initial multimap is itself a
///       sortable index-like object, as would be the case when mapping e.g.
///       hit ids to particle ids/ barcodes.
template <typename Value, typename Key>
inline IndexMultimap<Key, Value> invertIndexMultimap(
    const IndexMultimap<Value, Key>& multimap) {
  // switch key-value without enforcing the new ordering (linear copy)
  typename IndexMultimap<Key, Value>::sequence_type unordered;
  unordered.reserve(multimap.size());
  for (const auto& keyValue : multimap) {
    // value is now the key and the key is now the value
    unordered.emplace_back(keyValue.second, keyValue.first);
  }
  // adopting the unordered sequence will reestablish the correct order
  IndexMultimap<Key, Value> inverse;
  inverse.adopt_sequence(std::move(unordered));
  return inverse;
}

}  // namespace FW
