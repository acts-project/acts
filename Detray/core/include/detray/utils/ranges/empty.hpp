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
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"

namespace detray::ranges {

/// @brief Implements an empty view of a given type.
///
/// @see https://en.cppreference.com/w/cpp/ranges/empty_view
///
/// @tparam value_t type of the empty view.
template <typename value_t>
class empty_view : public detray::ranges::view_interface<empty_view<value_t>> {
 public:
  /// Iterator category: bidirectional
  using iterator_t = value_t*;

  /// @returns @c nullptr
  DETRAY_HOST_DEVICE
  static constexpr value_t* begin() noexcept { return nullptr; }

  /// @returns @c nullptr
  DETRAY_HOST_DEVICE
  static constexpr value_t* end() noexcept { return nullptr; }

  /// @returns @c nullptr
  DETRAY_HOST_DEVICE
  static constexpr value_t* data() noexcept { return nullptr; }

  /// @returns zero
  DETRAY_HOST_DEVICE
  static constexpr std::size_t size() noexcept { return 0; }

  /// @returns @c true
  DETRAY_HOST_DEVICE
  static constexpr bool empty() noexcept { return true; }

  /// No value contained in view
  /// @{
  DETRAY_HOST_DEVICE
  constexpr value_t front() const noexcept = delete;

  DETRAY_HOST_DEVICE
  constexpr value_t back() const noexcept = delete;

  DETRAY_HOST_DEVICE
  constexpr value_t operator[](const dindex) const noexcept = delete;
  /// @}
};

namespace views {

template <typename value_t>
using empty = empty_view<value_t>;

}  // namespace views

template <class T>
inline constexpr bool enable_borrowed_range<detray::ranges::empty_view<T>> =
    true;

}  // namespace detray::ranges
