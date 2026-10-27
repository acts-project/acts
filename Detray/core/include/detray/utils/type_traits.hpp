// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"

// System include(s)
#include <limits>
#include <type_traits>

namespace detray::detail {

/// Extract the value type of e.g. a container.
/// @{
template <typename T>
struct get_value_type {
  /// The type itself when the argument is not a container or pointer.
  using type = T;
};

/// Extract the pointee type when the input is a raw pointer.
template <typename T>
struct get_value_type<T*> {
  /// The pointee type.
  using type = T;
};

/// Extract a container's `value_type` when it is not `void`.
template <typename container_t>
  requires(!std::same_as<typename std::remove_cvref_t<container_t>::value_type,
                         void>)
struct get_value_type<container_t> {
  /// The container's value type.
  using type = typename std::remove_cvref_t<container_t>::value_type;
};

template <typename T>
using get_value_t = typename get_value_type<T>::type;
/// @}

/// Extract the first type from a parameter pack.
/// @{
template <typename first_t = void, typename... other_t>
struct first_type {
  /// The first type in the pack, or `void` for an empty pack.
  using type = first_t;
};

template <typename... TYPES>
using first_t = typename first_type<TYPES...>::type;

/// Extract the first value from an index pack.
template <std::size_t first_t, std::size_t... other_t>
struct first_idx {
  /// The first index in the pack.
  static constexpr std::size_t value = first_t;
};

template <std::size_t... TYPES>
inline constexpr std::size_t first_idx_v = first_idx<TYPES...>::value;
/// @}

/// Get the position of a type in a template parameter pack
/// @{
template <typename T, typename... Ts>
struct get_type_pos {
  /// Unroll a parameter pack without using a tuple to find `T`.
  /// @tparam first_t current type being compared
  /// @tparam remaining_types types still to examine
  /// @return one-based position from the back, or the maximum `std::size_t` if no
  /// match is found
  template <typename first_t, typename... remaining_types>
  DETRAY_HOST_DEVICE static constexpr std::size_t type_pos_back() {
    if constexpr (!std::is_same_v<T, first_t>) {
      return type_pos_back<remaining_types...>();
    }
    if constexpr (std::is_same_v<T, first_t>) {
      return sizeof...(remaining_types) + 1;
    }
    return std::numeric_limits<std::size_t>::max();
  }

  /// Zero-based position of `T` in `Ts...`.
  static constexpr std::size_t value = sizeof...(Ts) - type_pos_back<Ts...>();
};

template <typename T, typename... Ts>
inline constexpr std::size_t get_type_pos_v = get_type_pos<T, Ts...>::value;
/// @}

}  // namespace detray::detail
