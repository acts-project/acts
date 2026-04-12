// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <type_traits>

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
constexpr std::uint64_t hashMix(std::uint64_t x) {
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}

/// Combine one or more values into a single hash.
///
/// Each value is passed through @c std::hash and @c hashMix, then folded
/// together using an asymmetric combine step (fold-left). This produces
/// well-distributed, order-dependent hashes even for small or sequential
/// integer inputs where @c std::hash is typically the identity.
///
/// @code
/// std::uint64_t h = Acts::hashMixAndCombine(a, b, c);
/// @endcode
///
/// @tparam T the type of the first value
/// @tparam Rest the types of any additional values
/// @param value the first value to hash
/// @param rest additional values to fold in
/// @return the combined hash
template <typename T, typename... Rest>
[[nodiscard]]
constexpr std::uint64_t hashMixAndCombine(const T& value, const Rest&... rest) {
  std::uint64_t seed = hashMix(std::hash<T>{}(value));
  if constexpr (sizeof...(rest) > 0) {
    // Fold-left combine adapted from boost::hash_combine's 32-bit formula,
    // extended to 64-bit with wider shifts. The XOR with shifted seed makes
    // the operation asymmetric, preserving argument order.
    // @see https://github.com/boostorg/container_hash/blob/e3cbbebc8a1f9833287c8eb52fb0484ba744646b/include/boost/container_hash/hash.hpp#L469-L473
    auto foldLeft = [&seed](const auto& v) {
      std::uint64_t h = hashMix(std::hash<std::decay_t<decltype(v)>>{}(v));
      seed ^= h + 0x9e3779b97f4a7c15ULL + (seed << 12) + (seed >> 4);
    };
    (foldLeft(rest), ...);
  }
  return seed;
}

}  // namespace Acts
