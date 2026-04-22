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
#include "detray/utils/ranges/empty.hpp"
#include "detray/utils/ranges/ranges.hpp"

// System include(s)
#include <type_traits>

namespace detray::ranges {

namespace detail {
/// @brief Implements a subrange by providing start and end iterators on
/// another range.
///
/// @see https://en.cppreference.com/w/cpp/ranges/subrange
///
/// @tparam range_t the iterable which to constrain to a subrange.
template <detray::ranges::range range_t>
class subrange_view
    : public detray::ranges::view_interface<subrange_view<range_t>> {
 public:
  using iterator_t = typename detray::ranges::iterator_t<range_t>;
  using const_iterator_t = typename detray::ranges::const_iterator_t<range_t>;
  using difference_t = typename detray::ranges::range_difference_t<range_t>;

  /// Default constructor
  constexpr subrange_view() = default;

  /// Construct from an @param start and @param end iterator pair.
  DETRAY_HOST_DEVICE constexpr subrange_view(iterator_t start, iterator_t end)
      : m_begin{start}, m_end{end} {}

  /// @return start position of range.
  DETRAY_HOST_DEVICE
  constexpr auto begin() -> iterator_t { return m_begin; }

  /// @return sentinel of the range.
  DETRAY_HOST_DEVICE
  constexpr auto end() -> iterator_t { return m_end; }

  /// @return start position of the range - const
  DETRAY_HOST_DEVICE
  constexpr auto begin() const -> const_iterator_t { return m_begin; }

  /// @return sentinel of the range.
  DETRAY_HOST_DEVICE
  constexpr auto end() const -> const_iterator_t { return m_end; }

  /// Equality operator
  ///
  /// @param rhs the subrange to compare with
  ///
  /// @returns whether the two subranges are equal
  DETRAY_HOST_DEVICE
  constexpr auto operator==(const subrange_view &rhs) const -> bool {
    return m_begin == rhs.m_begin && m_end == rhs.m_end;
  }

 private:
  /// Start and end position of the subrange
  iterator_t m_begin;
  iterator_t m_end;
};

}  // namespace detail

/// @brief interface type to construct a @c subrange_view with CTAD
template <detray::ranges::range range_t, typename index_range_t = dindex_range>
struct subrange : public detail::subrange_view<range_t> {
  using base_type = detray::ranges::detail::subrange_view<range_t>;

  using iterator_t = typename detray::ranges::iterator_t<range_t>;
  using difference_t = typename detray::ranges::range_difference_t<range_t>;

  /// Default constructor
  constexpr subrange() = default;

  /// Construct from a @param range and an index range @param pos.
  template <concepts::index index_t>
  DETRAY_HOST_DEVICE constexpr subrange(index_t start, index_t end)
      : m_pos{start, end} {}

  /// Construct from an @param start and @param end iterator pair.
  DETRAY_HOST_DEVICE constexpr subrange(iterator_t start, iterator_t end)
      : base_type{start, end} {}

  /// Construct from a @param range.
  template <detray::ranges::range deduced_range_t>
  DETRAY_HOST_DEVICE constexpr explicit subrange(deduced_range_t &&range)
      : base_type{detray::ranges::begin(range), detray::ranges::end(range)} {}

  /// Construct from a @param range and starting position @param pos. Used
  /// as an overload when only a single position is needed.
  template <detray::ranges::range deduced_range_t, concepts::index index_t>
    requires std::is_convertible_v<
        index_t, detray::ranges::range_difference_t<deduced_range_t>>
  DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range, index_t pos)
      : base_type{detray::ranges::next(detray::ranges::begin(range),
                                       static_cast<difference_t>(pos)),
                  detray::ranges::next(
                      detray::ranges::begin(range),
                      static_cast<difference_t>(pos) + difference_t{1})} {}

  /// Construct from a @param range and an index range @param pos.
  template <detray::ranges::range deduced_range_t,
            concepts::interval deduced_index_range_t>
  DETRAY_HOST_DEVICE constexpr subrange(deduced_range_t &&range,
                                        const deduced_index_range_t &pos)
      : base_type{get_itr<0>(range, pos), get_itr<1>(range, pos)} {}

  /// Call operator for range composition
  template <detray::ranges::range deduced_range_t>
  DETRAY_HOST_DEVICE constexpr auto operator()(deduced_range_t &&range) {
    return detail::subrange_view<deduced_range_t>{get_itr<0>(range, m_pos),
                                                  get_itr<1>(range, m_pos)};
  }

 private:
  /// @returns an iterator over the @param range at a position taken from
  /// the I-th element of @param pos
  template <std::size_t I, detray::ranges::range deduced_range_t,
            concepts::interval deduced_index_range_t>
  DETRAY_HOST_DEVICE constexpr auto get_itr(deduced_range_t &&range,
                                            const deduced_index_range_t &pos) {
    return detray::ranges::next(
        detray::ranges::begin(std::forward<deduced_range_t>(range)),
        static_cast<difference_t>(detray::detail::get<I>(pos)));
  }

  /// Index range that defines the subrange
  index_range_t m_pos{};
};

// deduction guides
DETRAY_HOST_DEVICE subrange()
    -> subrange<detray::ranges::views::empty<int>, dindex_range>;

template <concepts::index index_t>
DETRAY_HOST_DEVICE subrange(index_t start, index_t end)
    -> subrange<detray::ranges::views::empty<int>, darray<index_t, 2>>;

template <detray::ranges::range deduced_range_t>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range)
    -> subrange<std::remove_reference_t<deduced_range_t>, bool>;

template <detray::ranges::range deduced_range_t, concepts::index index_t>
  requires std::convertible_to<
      index_t, detray::ranges::range_difference_t<deduced_range_t>>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, index_t pos)
    -> subrange<std::remove_reference_t<deduced_range_t>, index_t>;

template <detray::ranges::range deduced_range_t,
          concepts::interval index_range_t>
DETRAY_HOST_DEVICE subrange(deduced_range_t &&range, index_range_t &&pos)
    -> subrange<std::remove_reference_t<deduced_range_t>,
                std::remove_cvref_t<index_range_t>>;

/// @see https://en.cppreference.com/w/cpp/ranges/borrowed_iterator_t
template <detray::ranges::range R, concepts::interval I>
using borrowed_subrange_t =
    std::conditional_t<borrowed_range<R>, detray::ranges::subrange<R, I>,
                       dangling>;

}  // namespace detray::ranges
