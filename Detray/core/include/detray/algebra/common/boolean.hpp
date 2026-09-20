// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace detray {

namespace algebra::boolean {

/// Utilities for single booleans: default case
/// @{
constexpr bool any_of(bool b) {
  return b;
}
constexpr bool all_of(bool b) {
  return b;
}
constexpr bool none_of(bool b) {
  return !b;
}
constexpr std::size_t count(bool b) {
  return b ? 1u : 0u;
}

template <typename T, typename U>
  requires(std::is_scalar_v<T> && std::convertible_to<U, T>)
constexpr void set_if(T& value, const bool mask, const U new_value) {
  value = mask ? static_cast<T>(new_value) : value;
}

template <typename T>
  requires std::is_scalar_v<T>
constexpr void set_zero(T& value, const bool mask) {
  set_if(value, mask, static_cast<T>(0));
}

template <typename T>
  requires std::is_scalar_v<T>
constexpr void set_zero_inverted(T& value, const bool mask) {
  set_if(value, !mask, static_cast<T>(0));
}
///@}

}  // namespace algebra::boolean

namespace detail {

// Pull the boolean functions directly into the detray namespace
using namespace ::detray::algebra::boolean;

}  // namespace detail

}  // namespace detray
