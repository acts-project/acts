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

namespace ActsExamples::Traccc::Common::Util {

/// @brief Helper function for matchMap to find the index of a match.
/// The index is return true the out parameter.
/// The function returns false if no match was found.
template <typename T, typename A, typename equal_fn_t>
auto findMatchIdx(std::size_t* out, const T& element,
                  const std::vector<std::size_t>& candidateIdxs,
                  const std::vector<T, A>& candidateVec,
                  const equal_fn_t& eqFn) {
  for (std::size_t i = 0; i < candidateIdxs.size(); i++) {
    auto idx = candidateIdxs[i];
    if (eqFn(element, candidateVec[idx])) {
      *out = i;
      return true;
    }
  }
  return false;
}

/// @brief Determines which pairs of indices in the two vectors contain items that are equal.
/// Each item in the vector 'from̈́' must be equal to at least 1 item in the
/// vector 'to'. Furthermore, if the flag 'bijection' is enabled each element in
/// ̈́'to' and 'from' can only be paired with exactly 1 element.
/// @note When 'bijection' is corresponds to the assumption that the two vectors contain exactly the same elements.
/// The ordering in the vectors may however be different (determining the
/// difference in the order is where this function can be helpful).
/// @param from a vector (the retutned map's domain).
/// @param to a vector (the returned map's codomain).
/// @param lshFn the hashing function used for locality sensitive hashing. Any items that could potential be equal should have the same hash value.
/// @param eqFn the function for determining if two items are equal.
/// @param bijection flag indicatiing whether the map should be a bijection.
/// @return a map: index ('from' vector) -> index ('to' vector).
template <typename T, typename A1, typename A2, typename hash_fn_t, typename equal_fn_t>
inline auto matchMap(const std::vector<T, A1>& from,
                     const std::vector<T, A2>& to,
                     const hash_fn_t& lshFn,
                     const equal_fn_t& eqFn,
                     const bool bijection = true) {
  // By composition of functions, we can combine the maps "index ('from' vector)
  // -> hash value" and "hash value -> index ('to' vector)" to obtain "index
  // ('from' vector) -> index ('to' vector)".

  if (bijection && from.size() != to.size()) {
    throw std::runtime_error(
        "Cannot create a bijection as domain and codomain do not have the same "
        "cardinality");
  }

  // We start by creating the map "hash value -> index ('to' vector)".

  // Since there can be collisions with the hash values
  // each hash value maps to a bucket of indices rather than a single index.

  std::unordered_map<std::size_t, std::vector<std::size_t>> map1;
  for (std::size_t toIdx = 0; toIdx < to.size(); toIdx++) {
    auto& toElement = to[toIdx];
    auto toElementHash = lshFn(toElement);
    map1[toElementHash].push_back(toIdx);
  }

  // We can build the map "index in 'from' -> index in 'to'" directly.
  std::map<std::size_t, std::size_t> res;
  for (std::size_t fromIdx = 0; fromIdx < from.size(); fromIdx++) {
    auto& fromElement = from[fromIdx];
    auto fromElementHash = lshFn(fromElement);
    // We now find the exact element to match fromElement with in the bucket.
    auto& candidateIdxs = map1[fromElementHash];
    std::size_t idx;
    if (candidateIdxs.empty() ||
        !findMatchIdx(&idx, fromElement, candidateIdxs, to, eqFn)) {
      throw std::runtime_error("Could not find a match for an element");
    }
    res[fromIdx] = candidateIdxs[idx];
    if (bijection) {
      candidateIdxs.erase(candidateIdxs.begin() + idx);
    }
  }
  return res;
}

}  // namespace ActsExamples::Traccc::Common::Util
