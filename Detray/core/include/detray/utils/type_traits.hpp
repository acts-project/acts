// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"

// System include(s)
#include <type_traits>

namespace detray::detail {

/// Extract the value type of e.g. a container.
/// @{
template <typename T>
struct get_value_type {
  using type = T;
};

template <typename T>
struct get_value_type<T*> {
  using type = T;
};

template <typename container_t>
  requires(!std::same_as<typename std::remove_cvref_t<container_t>::value_type,
                         void>)
struct get_value_type<container_t> {
  using type = typename std::remove_cvref_t<container_t>::value_type;
};

template <typename T>
using get_value_t = typename get_value_type<T>::type;
/// @}

/// Extract the first type from a parameter pack.
/// @{
template <typename first_t = void, typename... other_t>
struct first_type {
  using type = first_t;
};

template <typename... TYPES>
using first_t = typename first_type<TYPES...>::type;

/// Extract the first value from an index pack.
template <std::size_t first_t, std::size_t... other_t>
struct first_idx {
  static constexpr std::size_t value = first_t;
};

template <std::size_t... TYPES>
inline constexpr std::size_t first_idx_v = first_idx<TYPES...>::value;
/// @}

/// Get the position of a type in a template parameter pack
/// @{
template <typename T, typename... Ts>
struct get_type_pos {
  /// Unroll a parameter pack without using a tuple.
  ///
  /// @note Returns the position of the type counted from the back!
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

  static constexpr std::size_t value = sizeof...(Ts) - type_pos_back<Ts...>();
};

template <typename T, typename... Ts>
inline constexpr std::size_t get_type_pos_v = get_type_pos<T, Ts...>::value;
/// @}

}  // namespace detray::detail
