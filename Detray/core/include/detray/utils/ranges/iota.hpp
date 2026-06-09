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
#include "detray/definitions/indexing.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/ranges/detail/iterator_functions.hpp"
#include "detray/utils/ranges/ranges.hpp"

// System include(s)
#include <concepts>
#include <iterator>
#include <type_traits>

namespace detray::ranges {

/// @brief Range factory that produces a sequence of values.
///
/// @see https://en.cppreference.com/w/cpp/ranges/iota_view
///
/// @tparam incr_t the incrementable type that makes up the sequence
///
/// @note If given single value, does not do infinite iteration, but only jumps
///       to next value.
/// @note Is not fit for lazy evaluation.
template <std::incrementable incr_t>
class iota_view : public detray::ranges::view_interface<iota_view<incr_t>> {
 private:
  /// @brief Nested iterator to generate a range of values on demand.
  struct iterator {
    using difference_type = std::ptrdiff_t;
    using value_type = incr_t;
    using pointer = incr_t *;
    using reference = incr_t;
    using iterator_category =
        std::conditional_t<concepts::random_access_incrementable<incr_t> &&
                               std::totally_ordered<incr_t>,
                           detray::ranges::random_access_iterator_tag,
                           detray::ranges::bidirectional_iterator_tag>;

    constexpr iterator()
      requires std::default_initializable<incr_t>
    = default;

    /// Parametrized Constructor
    DETRAY_HOST_DEVICE constexpr explicit iterator(incr_t i) : m_i{i} {}

    /// Increment the index
    /// @{
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> iterator & {
      ++m_i;
      return *this;
    }

    DETRAY_HOST_DEVICE constexpr auto operator++(int) -> iterator {
      auto tmp(*this);
      ++(*this);
      return tmp;
    }
    DETRAY_HOST_DEVICE
    constexpr auto operator--() -> iterator & {
      --m_i;
      return *this;
    }

    DETRAY_HOST_DEVICE constexpr auto operator--(int) -> iterator {
      auto tmp(*this);
      --(*this);
      return tmp;
    }

    DETRAY_HOST_DEVICE constexpr iterator &operator+=(const difference_type n)
      requires concepts::random_access_incrementable<incr_t>
    {
      difference_type nv = static_cast<difference_type>(m_i) + n;
      // Check if the round trip between difference type and incrementor
      // type is valid.
      assert(static_cast<difference_type>(static_cast<incr_t>(nv)) == nv);
      m_i = static_cast<incr_t>(nv);
      return *this;
    }

    DETRAY_HOST_DEVICE friend constexpr iterator operator+(
        const iterator &i, const difference_type n)
      requires concepts::random_access_incrementable<incr_t>
    {
      return iterator(i.m_i + n);
    }

    DETRAY_HOST_DEVICE friend constexpr iterator operator+(
        const difference_type n, const iterator &i)
      requires concepts::random_access_incrementable<incr_t>
    {
      return iterator(i.m_i + n);
    }

    DETRAY_HOST_DEVICE constexpr iterator &operator-=(const difference_type n)
      requires concepts::random_access_incrementable<incr_t>
    {
      m_i -= n;
      return *this;
    }

    DETRAY_HOST_DEVICE friend constexpr iterator operator-(
        const iterator &i, const difference_type n)
      requires concepts::random_access_incrementable<incr_t>
    {
      return iterator(i.m_i - n);
    }

    DETRAY_HOST_DEVICE friend constexpr difference_type operator-(
        const iterator &i, const iterator &j)
      requires concepts::random_access_incrementable<incr_t>
    {
      return i.m_i - j.m_i;
    }
    /// @}

    /// @returns the current value in the sequence - copy
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const
        noexcept(std::is_nothrow_copy_constructible_v<incr_t>) -> incr_t {
      return m_i;
    }

    DETRAY_HOST_DEVICE
    constexpr incr_t operator[](const difference_type n) const
      requires concepts::random_access_incrementable<incr_t>
    {
      return static_cast<incr_t>(m_i + n);
    }

   private:
    DETRAY_HOST_DEVICE
    friend constexpr auto operator<=>(const iterator &lhs, const iterator &rhs)
      requires std::totally_ordered<incr_t>
    = default;

    /// Current value of sequence
    incr_t m_i{};
  };

  /// Start and end values of the sequence
  incr_t m_start;
  incr_t m_end;

 public:
  using iterator_t = iterator;

  /// Default constructor (only works if @c imrementable_t is default
  /// constructible)
  constexpr iota_view()
    requires std::default_initializable<incr_t>
  = default;

  /// Construct from an @param interval that defines start and end values.
  template <concepts::interval interval_t>
    requires(!std::same_as<interval_t, iota_view>)
  DETRAY_HOST_DEVICE constexpr explicit iota_view(interval_t &&interval)
      : m_start{detray::detail::get<0>(interval)},
        m_end{detray::detail::get<1>(interval)} {
    assert(m_start <= m_end);
  }

  /// Construct from a @param start start and @param end value.
  DETRAY_HOST_DEVICE constexpr iota_view(incr_t start, incr_t end)
      : m_start{start}, m_end{end} {
    assert(m_start <= m_end);
  }

  /// Construct from just a @param start value to represent a single value seq
  DETRAY_HOST_DEVICE
  constexpr explicit iota_view(incr_t start)
      : m_start{start}, m_end{start + 1} {
    assert(m_start <= m_end);
  }

  /// @return start position of range on container.
  DETRAY_HOST_DEVICE
  constexpr auto begin() const -> iterator_t { return iterator_t{m_start}; }

  /// @return sentinel of a sequence.
  DETRAY_HOST_DEVICE
  constexpr auto end() const -> iterator_t { return iterator_t{m_end}; }

  /// @returns the span of the sequence
  DETRAY_HOST_DEVICE
  constexpr auto size() const -> incr_t { return m_end - m_start; }
};

namespace views {

/// @brief interface type to construct a @c iota_view with CTAD
template <std::incrementable incr_t>
struct iota : public detray::ranges::iota_view<incr_t> {
  using base_type = detray::ranges::iota_view<incr_t>;

  constexpr iota()
    requires std::default_initializable<incr_t>
  = default;

  template <concepts::interval interval_t>
    requires(!std::same_as<interval_t, iota>)
  DETRAY_HOST_DEVICE constexpr explicit iota(interval_t &&interval)
      : base_type(std::forward<interval_t>(interval)) {}

  DETRAY_HOST_DEVICE constexpr iota(incr_t start, incr_t end)
      : base_type(start, end) {}

  DETRAY_HOST_DEVICE constexpr explicit iota(incr_t start) : base_type(start) {}
};

// deduction guides
template <concepts::interval interval_t>
DETRAY_HOST_DEVICE iota(interval_t &&interval)
    -> iota<std::remove_cvref_t<decltype(detray::detail::get<0>(interval))>>;

template <concepts::index I>
DETRAY_HOST_DEVICE iota(I start, I end) -> iota<I>;

template <concepts::index I>
DETRAY_HOST_DEVICE iota(I start) -> iota<I>;

}  // namespace views

}  // namespace detray::ranges
