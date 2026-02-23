// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <csignal>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <type_traits>
#include <utility>

namespace Acts {
/// @brief Type alias for hashed string representation
/// @details Represents a string as a compile-time hash value for efficient comparison
using HashedString = std::uint32_t;

// Adapted from https://gist.github.com/Lee-R/3839813
namespace detail {
// FNV-1a 32bit hashing algorithm.
constexpr HashedString fnv1a_32(char const* s, std::size_t count) {
  return count != 0u ? (fnv1a_32(s, count - 1) ^ s[count - 1]) * 16777619u
                     : 2166136261u;
}

constexpr HashedString fnv1a_32(std::string_view s) {
  return !s.empty() ? (fnv1a_32(s.substr(0, s.size() - 1)) ^ s[s.size() - 1]) *
                          16777619u
                    : 2166136261u;
}

// FNV-1a 64bit hashing algorithm.
constexpr std::uint64_t fnv1a_64(const char* data, std::size_t len) {
  constexpr std::uint64_t fnv_offset_basis = 0xcbf29ce484222325ULL;
  constexpr std::uint64_t fnv_prime = 0x100000001b3ULL;

  std::uint64_t hash = fnv_offset_basis;
  for (std::size_t i = 0; i < len; ++i) {
    hash ^= static_cast<std::uint64_t>(static_cast<unsigned char>(data[i]));
    hash *= fnv_prime;
  }
  return hash;
}

constexpr std::uint64_t fnv1a_64(std::string_view sv) {
  return fnv1a_64(sv.data(), sv.size());
}

constexpr int length(const char* str) {
  return *str != 0 ? 1 + length(str + 1) : 0;
}
}  // namespace detail

/// Compile-time hash of string literal
/// @param s String view to hash
/// @return Hashed string representation
consteval HashedString hashString(std::string_view s) {
  return detail::fnv1a_32(s);
}

/// Runtime hash of string
/// @param s String view to hash
/// @return Hashed string representation
constexpr HashedString hashStringDynamic(std::string_view s) {
  return detail::fnv1a_32(s);
}

namespace HashedStringLiteral {
constexpr HashedString operator""_hash(char const* s, std::size_t count) {
  return detail::fnv1a_32(s, count);
}

}  // namespace HashedStringLiteral

template <typename T>
std::uint64_t typeHash() {
  const static std::uint64_t value = detail::fnv1a_64(typeid(T).name());
  return value;
}

}  // namespace Acts
