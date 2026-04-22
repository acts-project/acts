// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/ranges/ranges.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::ranges {

namespace detail {

template <detray::ranges::range T>
  requires std::input_iterator<detray::detail::get_value_t<T>>
struct static_join_iterator;

}

/// @brief Range adaptor that joins different ranges of the same type (static)
///
/// @see https://en.cppreference.com/w/cpp/ranges/join_view
///
/// @tparam I the number of ranges in the join.
/// @tparam range_itr_t the iterator type of the ranges.
///
/// @note Static implementation: The number of ranges needs to be know at
/// compile time
/// @note Does not take ownership of the ranges it operates on. Their lifetime
/// needs to be guaranteed throughout iteration or between iterations with the
/// same join instance.
/// @note Is not fit for lazy evaluation.
/// @todo improve performance of e.g. @c operator+ and @c operator+=
template <std::size_t I, std::input_iterator range_itr_t>
struct static_join_view
    : public detray::ranges::view_interface<static_join_view<I, range_itr_t>> {
  using iterator_coll_t = darray<range_itr_t, I>;
  using iterator_t =
      detray::ranges::detail::static_join_iterator<iterator_coll_t>;
  using value_t = std::iter_value_t<iterator_t>;

  /// Default constructor
  constexpr static_join_view() = default;

  /// Construct from a pack of @param ranges.
  template <detray::ranges::range... ranges_t>
  DETRAY_HOST_DEVICE constexpr explicit static_join_view(ranges_t &&...ranges)
      : m_begins{detray::ranges::begin(std::forward<ranges_t>(ranges))...},
        m_ends{detray::ranges::end(std::forward<ranges_t>(ranges))...} {}

  /// Construct from a pack of @param ranges - const
  template <detray::ranges::range... ranges_t>
  DETRAY_HOST_DEVICE constexpr explicit static_join_view(
      const ranges_t &...ranges)
      : m_begins{detray::ranges::cbegin(ranges)...},
        m_ends{detray::ranges::cend(ranges)...} {}

  /// @return start position of range - const
  DETRAY_HOST_DEVICE
  constexpr auto begin() const -> iterator_t { return {m_begins, m_ends}; }

  /// @return sentinel of the range.
  DETRAY_HOST_DEVICE
  constexpr auto end() const -> iterator_t {
    // Build a joined itr from the last value in the iterator collection
    return {m_begins, m_ends, detray::detail::get<I - 1>(m_ends), I - 1};
  }

  /// @returns a pointer to the beginning of the data of the first underlying
  /// range - const
  DETRAY_HOST_DEVICE
  constexpr auto data() const -> const value_t * {
    return &(*(detray::detail::get<0>(m_begins())));
  }

  /// @returns a pointer to the beginning of the data of the first underlying
  /// range - non-const
  DETRAY_HOST_DEVICE
  constexpr auto data() -> value_t * {
    return &(*(detray::detail::get<0>(m_begins())));
  }

  /// @returns sum of the number elements of all ranges in the join
  DETRAY_HOST_DEVICE
  constexpr auto size() const noexcept -> std::size_t {
    std::size_t size{0};
    for (std::size_t i{0}; i < I; ++i) {
      const range_itr_t &begin = m_begins[i];
      const range_itr_t &end = m_ends[i];
      size += static_cast<std::size_t>(detray::ranges::distance(begin, end));
    }
    return size;
  }

  /// Start and end position of the subranges
  iterator_coll_t m_begins{};
  iterator_coll_t m_ends{};
};

namespace views {

/// @brief interface type to construct a @c static_join_view with CTAD
template <std::size_t I, std::input_iterator range_itr_t>
struct static_join : public ranges::static_join_view<I, range_itr_t> {
  using base_type = ranges::static_join_view<I, range_itr_t>;

  constexpr static_join() = default;

  template <detray::ranges::range... ranges_t>
  DETRAY_HOST_DEVICE constexpr explicit static_join(ranges_t &&...ranges)
      : base_type(std::forward<ranges_t>(ranges)...) {}
};

// deduction guides
template <detray::ranges::range... ranges_t>
DETRAY_HOST_DEVICE static_join(ranges_t &&...ranges)
    -> static_join<sizeof...(ranges_t),
                   typename detray::ranges::iterator_t<detray::detail::first_t<
                       std::remove_reference_t<ranges_t>...>>>;

}  // namespace views

namespace detail {

/// @brief Sequentially iterate through multiple ranges of the same type.
///
/// Once the sentinel of one range is reached, set the current iterator to the
/// next ranges 'begin' (or 'end' if decrementing)
///
/// @tparam iterator_coll_t type of iterator collection of ranges to be joined.
///         Can contain const iterators.
///
/// @note The iterator must not be typed on the current range index, so that
/// begin and sentinel type are the same.
template <detray::ranges::range iterator_coll_t>
  requires std::input_iterator<detray::detail::get_value_t<iterator_coll_t>>
struct static_join_iterator {
  using iterator_t = detray::detail::get_value_t<iterator_coll_t>;

  using difference_type = std::iter_difference_t<iterator_t>;
  using value_type = std::iter_value_t<iterator_t>;
  using pointer = typename std::iterator_traits<iterator_t>::pointer;
  using reference = std::iter_reference_t<iterator_t>;
  using iterator_category =
      typename std::iterator_traits<iterator_t>::iterator_category;

  /// Default constructor required by LegacyIterator trait
  constexpr static_join_iterator() = default;

  /// Construct from a collection of @param begin and @param  end positions
  DETRAY_HOST_DEVICE
  constexpr static_join_iterator(const iterator_coll_t &begins,
                                 const iterator_coll_t &ends)
      : m_begins(&begins), m_ends(&ends), m_iter{(*m_begins)[0]} {}

  /// Fully parametrized construction
  DETRAY_HOST_DEVICE
  constexpr static_join_iterator(const iterator_coll_t &begins,
                                 const iterator_coll_t &ends,
                                 iterator_t current, const std::size_t i)
      : m_begins(&begins), m_ends(&ends), m_iter{current}, m_idx{i} {}

  /// Increment current iterator and check for switch between ranges.
  /// @{
  DETRAY_HOST_DEVICE constexpr auto operator++() -> static_join_iterator & {
    ++m_iter;
    // Switch to next range in the collection
    if (constexpr std::size_t max_idx{
            sizeof(iterator_coll_t) / sizeof(iterator_t) - 1u};
        (m_iter == (*m_ends)[m_idx]) && (m_idx < max_idx)) {
      ++m_idx;
      m_iter = (*m_begins)[m_idx];
    }
    return *this;
  }

  DETRAY_HOST_DEVICE constexpr auto operator++(int) -> static_join_iterator {
    auto tmp(*this);
    ++(*this);
    return tmp;
  }
  /// @}

  /// Decrement current iterator and check for switch between ranges.
  /// @{
  DETRAY_HOST_DEVICE constexpr auto operator--() -> static_join_iterator &
    requires std::bidirectional_iterator<iterator_t>
  {
    if (m_iter != (*m_begins)[m_idx]) {
      // Normal case
      --m_iter;
    } else if (m_idx > 0u) {
      // Iterator has reached last valid position in this range during the
      // previous decrement. Now go to the end of the previous range
      --m_idx;
      m_iter = detray::ranges::prev((*m_ends)[m_idx]);
    }
    return *this;
  }
  /// @}

  DETRAY_HOST_DEVICE constexpr auto operator--(int) -> static_join_iterator
    requires std::bidirectional_iterator<iterator_t>
  {
    auto tmp(*this);
    ++(*this);
    return tmp;
  }

  /// @returns the single value that the iterator points to.
  DETRAY_HOST_DEVICE
  constexpr auto operator*() const -> decltype(auto) { return *m_iter; }

  /// @returns advance this iterator state by @param j.
  DETRAY_HOST_DEVICE constexpr auto operator+=(const difference_type j)
      -> static_join_iterator &
    requires std::random_access_iterator<iterator_t>
  {
    // walk through join to catch the switch between intermediate ranges
    if (difference_type i{j}; i >= difference_type{0}) {
      while (i--) {
        ++(*this);
      }
    } else {
      while (i++) {
        --(*this);
      }
    }
    return *this;
  }

  /// @returns advance this iterator state by @param j.
  DETRAY_HOST_DEVICE constexpr auto operator-=(const difference_type j)
      -> static_join_iterator &
    requires std::random_access_iterator<iterator_t>
  {
    m_iter += (-j);
    return *this;
  }

  /// @returns the value at a given position - const
  DETRAY_HOST_DEVICE constexpr decltype(auto) operator[](
      const difference_type i) const
    requires std::random_access_iterator<iterator_t>
  {
    difference_type offset{i - (m_iter - (*m_begins)[0])};
    return *(*this + offset);
  }

 private:
  /// @returns true if it points to the same value.
  DETRAY_HOST_DEVICE friend constexpr bool operator==(
      const static_join_iterator &lhs, const static_join_iterator &rhs) {
    return (lhs.m_iter == rhs.m_iter);
  }

  /// @returns decision of comparison operators
  DETRAY_HOST_DEVICE friend constexpr auto operator<=>(
      const static_join_iterator &lhs, const static_join_iterator &rhs)
    requires std::random_access_iterator<iterator_t>
  {
#if defined(__apple_build_version__)
    const auto l{lhs.m_iter};
    const auto r{rhs.m_iter};
    if (l < r || (l == r && l < r)) {
      return std::strong_ordering::less;
    }
    if (l > r || (l == r && l > r)) {
      return std::strong_ordering::greater;
    }
    return std::strong_ordering::equivalent;
#else
    return (lhs.m_iter <=> rhs.m_iter);
#endif
  }

  /// @returns an iterator advanced by @param j through the join.
  DETRAY_HOST_DEVICE friend constexpr auto operator+(
      const static_join_iterator &itr, const difference_type j)
      -> static_join_iterator
    requires std::random_access_iterator<iterator_t>
  {
    static_join_iterator<iterator_coll_t> tmp(itr);
    // walk through join to catch the switch between intermediate ranges
    if (difference_type i{j}; i >= 0) {
      while (i--) {
        ++tmp;
      }
    } else {
      while (i++) {
        --tmp;
      }
    }
    return tmp;
  }

  /// @returns an iterator advanced by @param j through the join.
  DETRAY_HOST_DEVICE friend constexpr auto operator+(
      const difference_type j, const static_join_iterator &itr)
      -> static_join_iterator
    requires std::random_access_iterator<iterator_t>
  {
    return itr + j;
  }

  /// @returns an iterator advanced by @param j through the join.
  DETRAY_HOST_DEVICE friend constexpr auto operator-(
      const static_join_iterator &itr, const difference_type j)
      -> static_join_iterator
    requires std::random_access_iterator<iterator_t>
  {
    return itr + (-j);
  }

  /// @returns the positional difference between two iterators
  DETRAY_HOST_DEVICE friend constexpr auto operator-(
      const static_join_iterator &lhs, const static_join_iterator &rhs)
      -> difference_type
    requires std::random_access_iterator<iterator_t>
  {
    const static_join_iterator l{lhs};
    const static_join_iterator r{rhs};

    if (l.m_idx == r.m_idx) {
      return l.m_iter - r.m_iter;
    }
    if (l.m_idx < r.m_idx) {
      // Negative distance
      difference_type diff{l.m_iter - (*l.m_ends)[l.m_idx]};
      for (std::size_t i{l.m_idx + 1u}; i < r.m_idx; ++i) {
        diff += (*l.m_begins)[i] - (*l.m_ends)[i];
      }
      diff += (*r.m_begins)[l.m_idx] - r.m_iter;
      return diff;
    } else {
      // Positive distance
      difference_type diff{l.m_iter - (*l.m_begins)[l.m_idx]};
      for (std::size_t i{l.m_idx - 1u}; i > r.m_idx; --i) {
        diff += (*l.m_ends)[i] - (*l.m_begins)[i];
      }
      diff += (*r.m_ends)[r.m_idx] - r.m_iter;
      return diff;
    }
  }

  /// Global range collection of begin and end iterators
  const iterator_coll_t *m_begins{nullptr};
  const iterator_coll_t *m_ends{nullptr};
  /// This is the actual iterator state that will be advanced during iteration
  iterator_t m_iter{};
  /// Index of the current range in the join
  std::size_t m_idx{0};
};

}  // namespace detail

}  // namespace detray::ranges
