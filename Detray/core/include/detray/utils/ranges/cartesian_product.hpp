// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"
#include "detray/utils/tuple.hpp"

// System include(s)
#include <tuple>  // < needed for structured binding of iterator values

namespace detray::ranges {

namespace detail {

template <std::input_iterator... T>
struct cartesian_product_iterator;

}

/// @brief Range adaptor that generates a cartesian product from the given input
/// ranges
///
/// @see https://en.cppreference.com/w/cpp/ranges/cartesian_product_view
/// @see
/// https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/include/std/ranges
template <detray::ranges::range... range_ts>
struct cartesian_product_view : public detray::ranges::view_interface<
                                    cartesian_product_view<range_ts...>> {
  using iterator_coll_t = dtuple<detray::ranges::iterator_t<range_ts>...>;
  using iterator_t = detray::ranges::detail::cartesian_product_iterator<
      detray::ranges::iterator_t<range_ts>...>;
  using value_type =
      dtuple<std::iter_reference_t<detray::ranges::iterator_t<range_ts>>...>;

  /// Default constructor
  constexpr cartesian_product_view() = default;

  /// Construct from a pack of @param ranges.
  DETRAY_HOST_DEVICE constexpr explicit cartesian_product_view(
      range_ts... ranges)
      : m_begins(detray::ranges::begin(ranges)...),
        m_ends(detray::ranges::end(ranges)...) {}

  /// @returns start position of range - const
  DETRAY_HOST_DEVICE
  constexpr auto begin() const -> iterator_t { return {m_begins, m_ends}; }

  /// @returns sentinel of the range - const
  DETRAY_HOST_DEVICE
  constexpr auto end() const -> iterator_t { return {m_ends, m_ends}; }

  /// @returns a pointer to the beginning of the data of the first underlying
  /// range - const
  DETRAY_HOST_DEVICE
  constexpr auto data() const { return &(*(detray::get<0>(m_begins()))); }

  /// @returns a pointer to the beginning of the data of the first underlying
  /// range - non-const
  DETRAY_HOST_DEVICE
  constexpr auto data() { return &(*(detray::get<0>(m_begins()))); }

  /// @returns product of the number elements of all ranges in the view
  template <int I = sizeof...(range_ts) - 1>
  DETRAY_HOST_DEVICE constexpr auto size(std::size_t s = 1u) const noexcept
      -> std::size_t {
    decltype(auto) begin = detray::get<I>(m_begins);
    decltype(auto) end = detray::get<I>(m_ends);

    s *= static_cast<std::size_t>(detray::ranges::distance(begin, end));

    if constexpr (I > 0) {
      return size<I - 1>(s);
    } else {
      return s;
    }
  }

 private:
  /// Start and end position of the subranges
  iterator_coll_t m_begins{};
  iterator_coll_t m_ends{};
};

namespace views {

/// @brief interface type to construct a @c cartesian_product_view with CTAD
template <detray::ranges::range... range_ts>
struct cartesian_product : public ranges::cartesian_product_view<range_ts...> {
  using base_type = ranges::cartesian_product_view<range_ts...>;

  constexpr cartesian_product() = default;

  template <detray::ranges::range... deduced_range_ts>
  DETRAY_HOST_DEVICE constexpr explicit cartesian_product(
      deduced_range_ts &&...ranges)
      : base_type(std::forward<deduced_range_ts>(ranges)...) {}
};

// deduction guides
template <detray::ranges::range... ranges_ts>
DETRAY_HOST_DEVICE cartesian_product(ranges_ts &&...ranges)
    -> cartesian_product<ranges_ts...>;

}  // namespace views

namespace detail {

/// @brief Iterator implementation for the cartesian product view
template <std::input_iterator... iterator_ts>
struct cartesian_product_iterator {
  using difference_type = std::ptrdiff_t;
  using value_type = std::tuple<std::iter_reference_t<iterator_ts>...>;
  using pointer = value_type *;
  using reference = value_type;
  // TODO: Adapt to the weakest iterator category in pack
  using iterator_category = detray::ranges::bidirectional_iterator_tag;

  /// Default constructor required by LegacyIterator trait
  constexpr cartesian_product_iterator() = default;

  /// Construct from a collection of @param begin and @param end positions
  DETRAY_HOST_DEVICE
  constexpr cartesian_product_iterator(dtuple<iterator_ts...> begins,
                                       dtuple<iterator_ts...> ends)
      : m_begins(begins), m_ends(ends), m_itrs(begins) {}

  /// @returns true if the last range iterators are equal.
  DETRAY_HOST_DEVICE constexpr bool operator==(
      const cartesian_product_iterator &rhs) const {
    return (detray::get<0>(m_itrs) == detray::get<0>(rhs.m_itrs));
  }

  /// Increment iterators.
  /// @{
  DETRAY_HOST_DEVICE constexpr auto operator++()
      -> cartesian_product_iterator & {
    unroll_increment();
    return *this;
  }

  DETRAY_HOST_DEVICE constexpr auto operator++(int)
      -> cartesian_product_iterator {
    auto tmp(*this);
    ++(*this);
    return tmp;
  }
  /// @}

  /// Decrement iterators.
  /// @{
  DETRAY_HOST_DEVICE constexpr auto operator--()
      -> cartesian_product_iterator & {
    unroll_decrement();
    return *this;
  }

  DETRAY_HOST_DEVICE constexpr auto operator--(int)
      -> cartesian_product_iterator {
    auto tmp(*this);
    ++(*this);
    return tmp;
  }
  /// @}

  /// @returns the structured binding of all current iterator values - const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) operator*() const {
    return unroll_values(
        std::make_integer_sequence<std::size_t, sizeof...(iterator_ts)>{});
  }

 private:
  /// Unroll the increment over all iterators and, for the inner iterators,
  /// reset to start a new iteration
  /// @see
  /// https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/include/std/ranges
  template <int I = sizeof...(iterator_ts) - 1>
  DETRAY_HOST_DEVICE constexpr void unroll_increment() {
    auto &itr = detray::get<I>(m_itrs);
    ++itr;
    if constexpr (I > 0) {
      if (itr == detray::get<I>(m_ends)) {
        itr = detray::get<I>(m_begins);
        unroll_increment<I - 1>();
      }
    }
  }

  /// Unroll the dencrement over all iterators and, for the inner iterators,
  /// reset to start a new iteration from the end
  /// @see
  /// https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/include/std/ranges
  template <int I = sizeof...(iterator_ts) - 1>
  DETRAY_HOST_DEVICE constexpr void unroll_decrement() {
    auto &itr = detray::get<I>(m_itrs);

    if constexpr (I > 0) {
      if (itr == detray::get<I>(m_begins)) {
        itr = detray::ranges::prev(detray::get<I>(m_ends));
        unroll_decrement<I - 1>();
      }
    }
    --itr;
  }

  /// Pack the return values
  /// @note uses @c std::tuple for structured binding
  template <std::size_t... I>
  DETRAY_HOST_DEVICE constexpr auto unroll_values(
      std::index_sequence<I...>) const {
    return std::tuple<std::iter_reference_t<iterator_ts>...>(
        *detray::get<I>(m_itrs)...);
  }

  /// Global range collection of begin and end iterators
  dtuple<iterator_ts...> m_begins;
  dtuple<iterator_ts...> m_ends;
  /// Current iterator states
  dtuple<iterator_ts...> m_itrs;
};

}  // namespace detail

}  // namespace detray::ranges
