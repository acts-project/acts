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
#include "detray/utils/ranges/empty.hpp"
#include "detray/utils/ranges/ranges.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::ranges {

namespace detail {

template <detray::ranges::range T>
struct join_iterator;

}

/// @brief Range adaptor that joins different ranges of the same type (static)
///
/// @see https://en.cppreference.com/w/cpp/ranges/join_view
///
/// @tparam range_t a range of the ranges that should be joined.
///
/// @note Static implementation: The number of ranges needs to be know at
/// compile time
/// @note Does not take ownership of the ranges it operates on. Their lifetime
/// needs to be guaranteed throughout iteration or between iterations with the
/// same join instance.
template <detray::ranges::range range_t>
struct join_view : public detray::ranges::view_interface<join_view<range_t>> {
  /// Iterator over the range of ranges
  using outer_iterator_t = detray::ranges::iterator_t<range_t>;
  // Type of a single range
  using outer_value_t =
      std::conditional_t<std::is_const_v<range_t>,
                         const detray::ranges::range_value_t<range_t>,
                         detray::ranges::range_value_t<range_t>>;
  // Iterator over a single range
  using inner_iterator_t = detray::ranges::iterator_t<outer_value_t>;

  using iterator_t = detray::ranges::detail::join_iterator<range_t>;
  using value_t = std::iter_value_t<iterator_t>;

  /// Default constructor
  constexpr join_view() = default;

  /// Construct from a range of @param ranges.
  template <detray::ranges::range R>
  DETRAY_HOST_DEVICE constexpr explicit join_view(R &&ranges)
      : m_begin{detray::ranges::begin(ranges)},
        m_end{detray::ranges::end(ranges)} {}

  /// @return start position of range - const
  DETRAY_HOST_DEVICE
  constexpr auto begin() const -> iterator_t { return {m_begin, m_end}; }

  /// @return sentinel of the range.
  DETRAY_HOST_DEVICE
  constexpr auto end() const -> iterator_t { return {m_end, m_end}; }

  /// @returns a pointer to the beginning of the data of the first underlying
  /// range - const
  DETRAY_HOST_DEVICE
  constexpr auto data() const -> const value_t * { return &(*(*m_begin)); }

  /// @returns a pointer to the beginning of the data of the first underlying
  /// range - non-const
  DETRAY_HOST_DEVICE
  constexpr auto data() -> value_t * { return &(*(*m_begin)); }

  /// @returns sum of the number elements of all ranges in the join
  DETRAY_HOST_DEVICE
  constexpr auto size() const noexcept -> std::size_t {
    std::size_t size{0u};
    for (outer_iterator_t itr = m_begin; itr != m_end; ++itr) {
      // subrange
      const auto begin = detray::ranges::begin(*itr);
      const auto end = detray::ranges::end(*itr);
      size += static_cast<std::size_t>(detray::ranges::distance(begin, end));
    }
    return size;
  }

  /// Start and end position of the subranges
  outer_iterator_t m_begin{};
  outer_iterator_t m_end{};
};

namespace views {

/// @brief interface type to construct a @c join_view with CTAD
template <detray::ranges::range range_t>
struct join : public ranges::join_view<range_t> {
  using base_type = detray::ranges::join_view<range_t>;

  constexpr join() = default;

  template <detray::ranges::range deduced_range_t>
  DETRAY_HOST_DEVICE constexpr explicit join(deduced_range_t &&ranges)
      : base_type(std::forward<deduced_range_t>(ranges)) {}

  /// Call operator for range composition
  template <detray::ranges::range deduced_range_t>
  DETRAY_HOST_DEVICE constexpr auto operator()(deduced_range_t &&ranges) {
    return detray::ranges::join_view<deduced_range_t>(
        std::forward<deduced_range_t>(ranges));
  }

  /// Copy assignment operator
  DETRAY_HOST_DEVICE
  join &operator=(const join &other) {
    base_type::operator=(other);
    return *this;
  }
};

// deduction guides
DETRAY_HOST_DEVICE join() -> join<detray::ranges::views::empty<dvector<int>>>;

template <detray::ranges::range R>
DETRAY_HOST_DEVICE join(R &&ranges) -> join<std::remove_reference_t<R>>;

}  // namespace views

namespace detail {

/// @brief Sequentially iterate through multiple ranges of the same type.
///
/// Once the sentinel of one range is reached, set the current iterator to the
/// next ranges 'begin' (or 'end' if decrementing)
///
/// @tparam range_t a range that contains the ranges to be joined.
template <detray::ranges::range range_t>
struct join_iterator {
  using outer_iterator_t =
      std::conditional_t<std::is_const_v<range_t>,
                         detray::ranges::const_iterator_t<range_t>,
                         detray::ranges::iterator_t<range_t>>;
  using outer_value_t = decltype(*std::declval<outer_iterator_t>());
  using inner_iterator_t = std::conditional_t<
      std::is_same_v<outer_iterator_t,
                     detray::ranges::const_iterator_t<range_t>>,
      detray::ranges::const_iterator_t<outer_value_t>,
      detray::ranges::iterator_t<outer_value_t>>;

  using iterator_t = inner_iterator_t;
  using difference_type = std::iter_difference_t<iterator_t>;
  using value_type = std::iter_value_t<iterator_t>;
  using pointer = typename std::iterator_traits<iterator_t>::pointer;
  using reference = std::iter_reference_t<iterator_t>;
  using iterator_category =
      typename std::iterator_traits<iterator_t>::iterator_category;

  /// Default constructor required by LegacyIterator trait
  constexpr join_iterator() = default;

  /// Construct from range of ranges ( @param begin and @param  end )
  DETRAY_HOST_DEVICE
  constexpr join_iterator(outer_iterator_t begin, outer_iterator_t end)
      : m_outer_begin(begin), m_outer_end(end), m_outer_itr(begin) {
    if (m_outer_itr != m_outer_end) {
      m_inner_itr = (*m_outer_itr).begin();
    } else {
      m_inner_itr = {};
    }
    next_inner();
  }

  /// Increment current iterator and check for switch between ranges.
  /// @{
  DETRAY_HOST_DEVICE auto operator++() -> join_iterator & {
    ++m_inner_itr;
    next_inner();

    return *this;
  }

  DETRAY_HOST_DEVICE constexpr auto operator++(int) -> join_iterator {
    auto tmp(*this);
    ++(*this);
    return tmp;
  }
  /// @}

  /// Decrement current iterator and check for switch between ranges.
  /// @{
  DETRAY_HOST_DEVICE constexpr auto operator--() -> join_iterator &
    requires std::bidirectional_iterator<inner_iterator_t> &&
             std::bidirectional_iterator<outer_iterator_t>
  {
    // If we are calling this at the end of the join iteration, go back into
    // the valid range
    if (m_outer_itr == m_outer_end) {
      --m_outer_itr;
    }

    previous_inner();
    --m_inner_itr;

    return *this;
  }

  DETRAY_HOST_DEVICE constexpr auto operator--(int) -> join_iterator
    requires std::bidirectional_iterator<inner_iterator_t> &&
             std::bidirectional_iterator<outer_iterator_t>
  {
    auto tmp(*this);
    --(*this);
    return tmp;
  }
  /// @}

  /// @returns the single value that the iterator points to - const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) operator*() const { return *m_inner_itr; }

  /// @returns advance this iterator state by @param j.
  DETRAY_HOST_DEVICE constexpr auto operator+=(const difference_type j)
      -> join_iterator &
    requires std::random_access_iterator<inner_iterator_t> &&
             std::random_access_iterator<outer_iterator_t>
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
      -> join_iterator &
    requires std::random_access_iterator<inner_iterator_t> &&
             std::random_access_iterator<outer_iterator_t>
  {
    m_inner_itr += (-j);
    return *this;
  }

  /// @returns the value at a given position - const
  DETRAY_HOST_DEVICE constexpr decltype(auto) operator[](
      const difference_type i) const
    requires std::random_access_iterator<inner_iterator_t> &&
             std::random_access_iterator<outer_iterator_t>
  {
    difference_type offset{i -
                           (m_inner_itr - detray::ranges::begin(*m_outer_itr))};
    return *(*this + offset);
  }

 private:
  /// @returns true if it points to the same value.
  DETRAY_HOST_DEVICE friend constexpr bool operator==(
      const join_iterator &lhs, const join_iterator &rhs) {
    return (lhs.m_inner_itr == rhs.m_inner_itr);
  }

  /// @returns false if it points to the same value (usually the global
  /// sentinel of the join).
  DETRAY_HOST_DEVICE friend constexpr bool operator!=(
      const join_iterator &lhs, const join_iterator &rhs) {
    return (lhs.m_outer_itr != rhs.m_outer_itr);
  }

  /// @returns relation according to comparison operators
  DETRAY_HOST_DEVICE friend constexpr auto operator<=>(const join_iterator &lhs,
                                                       const join_iterator &rhs)
    requires detray::ranges::random_access_iterator<inner_iterator_t> &&
             std::random_access_iterator<outer_iterator_t>
  {
#if defined(__apple_build_version__)
    const auto l_o_itr{lhs.m_outer_itr};
    const auto l_i_itr{lhs.m_inner_itr};
    const auto r_o_itr{rhs.m_outer_itr};
    const auto r_i_itr{rhs.m_inner_itr};
    if (l_o_itr == r_o_itr) {
      if (l_i_itr < r_i_itr || (l_i_itr == r_i_itr && l_i_itr < r_i_itr)) {
        return std::strong_ordering::less;
      }
      if (l_i_itr > r_i_itr || (l_i_itr == r_i_itr && l_i_itr > r_i_itr)) {
        return std::strong_ordering::greater;
      }
      return std::strong_ordering::equivalent;
    } else {
      if (l_o_itr < r_o_itr || (l_o_itr == r_o_itr && l_o_itr < r_o_itr)) {
        return std::strong_ordering::less;
      }
      if (l_o_itr > r_o_itr || (l_o_itr == r_o_itr && l_o_itr > r_o_itr)) {
        return std::strong_ordering::greater;
      }
      return std::strong_ordering::equivalent;
    }
#else
    if (lhs.m_outer_itr == rhs.m_outer_itr) {
      return lhs.m_inner_itr <=> rhs.m_inner_itr;
    } else {
      return lhs.m_outer_itr <=> rhs.m_outer_itr;
    }
#endif
  }

  /// @returns an iterator advanced by @param j through the join.
  DETRAY_HOST_DEVICE friend constexpr auto operator+(const join_iterator &itr,
                                                     const difference_type j)
      -> join_iterator
    requires std::random_access_iterator<inner_iterator_t> &&
             std::random_access_iterator<outer_iterator_t>
  {
    join_iterator<range_t> tmp(itr);
    // walk through join to catch the switch between intermediate ranges
    if (difference_type i{j}; i >= difference_type{0}) {
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
  DETRAY_HOST_DEVICE friend constexpr auto operator+(const difference_type j,
                                                     const join_iterator &itr)
      -> join_iterator
    requires std::random_access_iterator<inner_iterator_t> &&
             std::random_access_iterator<outer_iterator_t>
  {
    return itr + j;
  }

  /// @returns an iterator advanced by @param j through the join.
  DETRAY_HOST_DEVICE friend constexpr auto operator-(const join_iterator &itr,
                                                     const difference_type j)
      -> join_iterator
    requires std::random_access_iterator<inner_iterator_t> &&
             std::random_access_iterator<outer_iterator_t>
  {
    return itr + (-j);
  }

  /// @returns the positional difference between two iterators
  DETRAY_HOST_DEVICE friend constexpr auto operator-(const join_iterator &left,
                                                     const join_iterator &right)
      -> difference_type
    requires std::random_access_iterator<inner_iterator_t> &&
             std::random_access_iterator<outer_iterator_t>
  {
    const join_iterator lhs{left};
    const join_iterator rhs{right};
    outer_iterator_t tmp_outer_itr;

    if (tmp_outer_itr == rhs.m_outer_itr) {
      return lhs.m_inner_itr - rhs.m_inner_itr;
    }
    if ((tmp_outer_itr - rhs.m_outer_itr) < 0) {
      // Negative distance
      difference_type diff{lhs.m_inner_itr -
                           detray::ranges::end(*lhs.m_outer_itr)};
      for (tmp_outer_itr + 1; tmp_outer_itr != rhs.m_outer_itr;
           ++tmp_outer_itr) {
        diff += detray::ranges::end(*tmp_outer_itr) -
                detray::ranges::begin(*tmp_outer_itr);
      }
      diff += rhs.m_inner_itr - detray::ranges::end(*rhs.m_outer_itr);
      return diff;
    } else {
      // Positive distance
      difference_type diff{lhs.m_inner_itr -
                           detray::ranges::begin(*tmp_outer_itr)};
      for (tmp_outer_itr - 1; tmp_outer_itr != rhs.m_outer_itr;
           --tmp_outer_itr) {
        diff += detray::ranges::end(*tmp_outer_itr) -
                detray::ranges::begin(*tmp_outer_itr);
      }
      diff += rhs.m_inner_itr - detray::ranges::begin(*rhs.m_outer_itr);
      return diff;
    }
  }

  /// Find the first inner range that is not empty
  constexpr void next_inner() {
    if (m_outer_itr == m_outer_end) {
      return;
    }
    while (m_inner_itr == detray::ranges::end(*m_outer_itr)) {
      ++m_outer_itr;
      if (m_outer_itr != m_outer_end) {
        m_inner_itr = (*m_outer_itr).begin();
      } else {
        break;
      }
    }
  }

  /// Find the last inner range that is not empty
  constexpr void previous_inner() {
    // Get the start of the current inner range
    inner_iterator_t inner_begin = detray::ranges::begin(*m_outer_itr);

    // Iterator has reached last valid position in this range
    // during the previous decrement. Now go to the end of the
    // previous range
    while (m_inner_itr == inner_begin) {
      // No more inner ranges to try
      if (m_outer_itr == m_outer_begin) {
        m_inner_itr = detray::ranges::end(*m_outer_begin);
        return;
      }

      --m_outer_itr;

      inner_begin = detray::ranges::begin(*m_outer_itr);
      m_inner_itr = detray::ranges::end(*m_outer_itr);
    }
  }

  /// Global range collection begin and end (outer iterators)
  outer_iterator_t m_outer_begin{};
  outer_iterator_t m_outer_end{};
  /// Current range
  outer_iterator_t m_outer_itr{};
  /// Current iterators over the inner ranges
  inner_iterator_t m_inner_itr{};
};

}  // namespace detail

}  // namespace detray::ranges
