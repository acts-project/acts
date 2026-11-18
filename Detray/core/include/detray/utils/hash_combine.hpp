// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s)
#include <cstddef>
#include <functional>

namespace detray {

/// Mix the hash of @param value into the hash @param seed
///
/// @tparam T the type of the value (needs a @c std::hash specialization)
///
/// @returns the combined hash value
template <typename T>
constexpr std::size_t hash_combine(std::size_t seed, const T &value) {
  // Same mixing as boost::hash_combine (64-bit golden ratio constant)
  seed ^= std::hash<T>{}(value) + 0x9e3779b97f4a7c15ull + (seed << 6u) +
          (seed >> 2u);
  return seed;
}

}  // namespace detray
