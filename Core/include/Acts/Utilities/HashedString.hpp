// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <cstdint>

namespace Acts {
// Adapted from https://gist.github.com/Lee-R/3839813
namespace detail {
// FNV-1a 32bit hashing algorithm.
constexpr std::uint32_t fnv1a_32(char const* s, std::size_t count) {
  return ((count ? fnv1a_32(s, count - 1) : 2166136261u) ^ s[count]) *
         16777619u;
}
}  // namespace detail

constexpr std::uint32_t operator"" _hash(char const* s, std::size_t count) {
  return detail::fnv1a_32(s, count);
}

}  // namespace Acts
