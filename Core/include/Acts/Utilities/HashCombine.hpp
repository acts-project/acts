// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <functional>

#include <boost/functional/hash.hpp>

namespace Acts {

/// Apply a bit-mixing function to break up patterns in integer hash values.
///
/// Many standard hash functions (including std::hash for integers) are the
/// identity, which leads to poor distribution when values are small or
/// sequential. This function applies avalanche mixing so that small changes
/// in the input produce large, uniformly-distributed changes in the output.
///
/// Uses the fmix64 finalizer from Austin Appleby's MurmurHash3.
/// @see https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp
///
/// @param x the value to mix
/// @return the mixed value
constexpr std::size_t hashMix(std::size_t x) {
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}

/// Combine a hash seed with one or more values.
///
/// Each value is passed through @c std::hash and then @c hashMix before being
/// folded into the seed via @c boost::hash_combine. This avoids the poor
/// distribution that @c boost::hash_combine alone exhibits when inputs are
/// small or patterned integers (where @c std::hash is typically the identity).
///
/// Multiple values can be combined in a single call:
/// @code
/// std::size_t seed = 0;
/// Acts::hashCombine(seed, a, b, c);
/// @endcode
///
/// @tparam T the type of the first value to combine
/// @tparam Rest the types of any additional values to combine
/// @param seed the hash seed to update in place
/// @param value the first value to fold into the seed
/// @param rest additional values to fold into the seed
template <typename T, typename... Rest>
constexpr void hashCombine(std::size_t& seed, const T& value, const Rest&... rest) {
  boost::hash_combine(seed, hashMix(std::hash<T>{}(value)));
  if constexpr (sizeof...(rest) > 0) {
    hashCombine(seed, rest...);
  }
}

}  // namespace Acts
