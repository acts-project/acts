// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"

namespace detray::ranges {

/// @brief Struct that implements a view on a single element.
///
/// @see https://en.cppreference.com/w/cpp/ranges/single_view
///
/// @tparam value_t type of the single element (outside of a container)
///
/// @note Does not take ownership of the value it operates on. Its lifetime
/// needs to be guaranteed throughout iteration or between iterations with the
/// same view instance.
/// @note Is not fit for lazy evaluation.
template <typename value_t>
class single_view
    : public detray::ranges::view_interface<single_view<value_t>> {
 public:
  /// Pointer iterator over the single stored element.
  using iterator_t = value_t*;

  /// Default constructor
  constexpr single_view() = default;

  /// Construct the view by copying its single element.
  /// Store a copy of `value` as the view's only element.
  /// @param value value to copy into the view
  DETRAY_HOST_DEVICE constexpr explicit single_view(const value_t& value)
      : m_value{value} {}

  /// Construct the view by moving in its single element.
  /// Move `value` into the view's only element.
  /// @param value value to move into the view
  DETRAY_HOST_DEVICE constexpr explicit single_view(value_t&& value)
      : m_value{std::move(value)} {}

  /// Construct value in place from @p args
  /// @tparam Args argument types used to construct the value
  /// @param args arguments forwarded to the value constructor
  template <class... Args>
  DETRAY_HOST_DEVICE constexpr explicit single_view(std::in_place_t /*tag*/,
                                                    Args&&... args)
      : m_value{std::in_place, std::forward<Args>(args)...} {}

  /// Return the single stored value by value.
  /// @return the single value
  DETRAY_HOST_DEVICE
  constexpr auto value() const -> value_t { return m_value; }

  /// Return a const reference to the stored element.
  /// @return the single value
  DETRAY_HOST_DEVICE
  constexpr auto ref() const -> const value_t& { return m_value; }
  /// Return a mutable reference to the stored element.
  /// @returns A mutable reference to the element.
  DETRAY_HOST_DEVICE
  constexpr auto ref() -> value_t& { return m_value; }

  /// @returns value pointer.
  DETRAY_HOST_DEVICE
  constexpr auto begin() noexcept -> value_t* { return &m_value; }

  /// @returns value pointer - const
  DETRAY_HOST_DEVICE
  constexpr auto begin() const noexcept -> const value_t* { return &m_value; }

  /// @returns sentinel position.
  DETRAY_HOST_DEVICE
  constexpr auto end() noexcept -> value_t* { return &m_value + 1; }

  /// @returns sentinel position - const
  DETRAY_HOST_DEVICE
  constexpr auto end() const noexcept -> const value_t* { return &m_value + 1; }

  /// @returns a pointer to the beginning of the data
  DETRAY_HOST_DEVICE
  constexpr auto data() noexcept -> value_t* { return &m_value; }

  /// @returns a pointer to the beginning of the data - const
  DETRAY_HOST_DEVICE
  constexpr auto data() const noexcept -> const value_t* { return &m_value; }

  /// @returns the size of the single view, which is always 'one'.
  DETRAY_HOST_DEVICE
  static constexpr auto size() noexcept -> std::size_t { return 1; }

  /// @returns the value directly
  DETRAY_HOST_DEVICE
  constexpr auto front() noexcept -> value_t { return m_value; }

  /// @returns the value directly
  DETRAY_HOST_DEVICE
  constexpr auto back() noexcept -> value_t { return m_value; }

  /// @returns the value directly
  DETRAY_HOST_DEVICE constexpr auto operator[](const dindex /*unused*/) const
      -> value_t {
    return m_value;
  }

 private:
  value_t m_value{};
};

namespace views {

/// @brief interface type to construct a @c single_view with CTAD
template <typename value_t>
struct single : public detray::ranges::single_view<value_t> {
  /// Base implementation storing the single value.
  using base_type = detray::ranges::single_view<value_t>;

  constexpr single() = default;

  /// Construct the view by copying a value.
  /// @param value value to copy
  template <typename deduced_value_t>
    requires(!std::same_as<deduced_value_t, single>)
  DETRAY_HOST_DEVICE constexpr explicit single(const deduced_value_t& value)
      : base_type(value) {}

  /// Construct the view by moving a value.
  /// @param value value to move
  template <typename deduced_value_t>
    requires(!std::same_as<deduced_value_t, single>)
  DETRAY_HOST_DEVICE constexpr explicit single(deduced_value_t&& value)
      : base_type(std::forward<deduced_value_t>(value)) {}

  /// Construct the stored value in place.
  /// @tparam Args value constructor argument types
  /// @param args arguments forwarded to the value constructor
  template <class... Args>
  DETRAY_HOST_DEVICE constexpr explicit single(std::in_place_t /*unused*/,
                                               Args&&... args)
      : base_type(std::in_place /*tag*/, std::forward<Args>(args)... /*args*/) {
  }
};

// deduction guides

template <typename deduced_value_t>
DETRAY_HOST_DEVICE_DEDUCTION_GUIDE single(deduced_value_t)
    -> single<deduced_value_t>;

}  // namespace views

}  // namespace detray::ranges
