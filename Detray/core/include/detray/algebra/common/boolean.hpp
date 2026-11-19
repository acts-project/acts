// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <concepts>
#include <cstddef>
#include <type_traits>

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
/// @param b Boolean to count.
/// @returns 1 if @p b is set, 0 otherwise.
DETRAY_HOST_DEVICE constexpr std::size_t count(bool b) {
  return b ? 1u : 0u;
}
/// @}

/// Masked assignments for single values
/// @{

/// Assign a new value if the mask is set
///
/// @tparam T Scalar type of the value to update.
/// @tparam U Type of the new value, must be convertible to @p T.
///
/// @param value Value to update.
/// @param mask Update @p value only if this is set.
/// @param new_value Value to assign.
template <typename T, typename U>
  requires(std::is_scalar_v<T> && std::convertible_to<U, T>)
DETRAY_HOST_DEVICE constexpr void set_if(T& value, const bool mask,
                                         const U new_value) {
  value = mask ? static_cast<T>(new_value) : value;
}

/// Set a value to zero if the mask is set
///
/// @tparam T Scalar type of the value to update.
///
/// @param value Value to update.
/// @param mask Set @p value to zero only if this is set.
template <typename T>
  requires std::is_scalar_v<T>
DETRAY_HOST_DEVICE constexpr void set_zero(T& value, const bool mask) {
  set_if(value, mask, static_cast<T>(0));
}

/// Set a value to zero if the mask is not set
///
/// @tparam T Scalar type of the value to update.
///
/// @param value Value to update.
/// @param mask Keep @p value only if this is set, otherwise set it to zero.
template <typename T>
  requires std::is_scalar_v<T>
DETRAY_HOST_DEVICE constexpr void set_zero_inverted(T& value, const bool mask) {
  set_if(value, !mask, static_cast<T>(0));
}
///@}

}  // namespace algebra::boolean

namespace detail {

// Pull the boolean functions directly into the detray namespace
using namespace ::detray::algebra::boolean;

}  // namespace detail

}  // namespace detray
