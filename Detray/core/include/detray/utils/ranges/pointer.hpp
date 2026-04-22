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

/// @brief Struct that implements a pointer based view on a single element.
///
/// @tparam value_t type of the single element (outside of a container)
///
/// @note Does not take ownership of the value it operates on. Its lifetime
/// needs to be guaranteed throughout iteration or between iterations with the
/// same view instance.
/// @note Is not fit for lazy evaluation.
template <typename value_t>
class pointer_view
    : public detray::ranges::view_interface<pointer_view<value_t>> {
 public:
  using iterator_t = value_t*;

  /// Default constructor
  pointer_view() = default;

  /// Construct iterator from the single @param value - copy
  DETRAY_HOST_DEVICE constexpr explicit pointer_view(value_t& value)
      : m_value{&value} {}

  /// @return the pointer value
  DETRAY_HOST_DEVICE
  constexpr auto operator*() const -> const value_t& { return *m_value; }

  /// @return the pointer value
  DETRAY_HOST_DEVICE
  constexpr auto operator*() -> value_t& { return *m_value; }

  /// @returns value pointer.
  DETRAY_HOST_DEVICE
  constexpr auto begin() noexcept -> value_t* { return m_value; }

  /// @returns value pointer - const
  DETRAY_HOST_DEVICE
  constexpr auto begin() const noexcept -> const value_t* { return m_value; }

  /// @returns sentinel position.
  DETRAY_HOST_DEVICE
  constexpr auto end() noexcept -> value_t* { return m_value + 1; }

  /// @returns sentinel position - const
  DETRAY_HOST_DEVICE
  constexpr auto end() const noexcept -> const value_t* { return m_value + 1; }

  /// @returns a pointer to the beginning of the data
  DETRAY_HOST_DEVICE
  constexpr auto data() noexcept -> value_t* { return m_value; }

  /// @returns a pointer to the beginning of the data - const
  DETRAY_HOST_DEVICE
  constexpr auto data() const noexcept -> const value_t* { return m_value; }

  /// @returns the size of the pointer view, which is always 'one'.
  DETRAY_HOST_DEVICE
  static constexpr auto size() noexcept -> std::size_t { return 1; }

  /// @returns the value directly
  DETRAY_HOST_DEVICE
  constexpr auto front() noexcept -> value_t { return *m_value; }

  /// @returns the value directly
  DETRAY_HOST_DEVICE
  constexpr auto back() noexcept -> value_t { return *m_value; }

  /// @returns the value directly
  DETRAY_HOST_DEVICE constexpr auto operator[](const dindex) const -> value_t {
    return *m_value;
  }

 private:
  value_t* m_value{};
};

namespace views {

/// @brief interface type to construct a @c pointer_view with CTAD
template <typename value_t>
struct pointer : public detray::ranges::pointer_view<const value_t> {
  using base_type = detray::ranges::pointer_view<const value_t>;

  constexpr pointer() = default;

  template <typename deduced_value_t>
  DETRAY_HOST_DEVICE constexpr explicit pointer(deduced_value_t& value)
      : base_type(value) {}
};

// deduction guides

template <typename deduced_value_t>
DETRAY_HOST_DEVICE pointer(deduced_value_t&) -> pointer<deduced_value_t>;

}  // namespace views

}  // namespace detray::ranges
