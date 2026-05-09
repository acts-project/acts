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
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/grid/detail/axis_bounds.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <utility>

namespace detray::axis::detail {

template <concepts::grid G, std::input_iterator I>
struct bin_iterator;

/// @returns the local bin indexer for the given @param search_window.
/// (cartesian product of the bin index ranges on the respective axes)
template <std::size_t... I>
DETRAY_HOST_DEVICE inline auto get_bin_indexer(
    const axis::multi_bin_range<sizeof...(I)> &search_window,
    std::index_sequence<I...> /*seq*/) {
  return detray::views::cartesian_product{
      detray::views::iota{detray::detail::get<I>(search_window)}...};
}

/// @brief Range adaptor that fetches grid bins according to a search window.
template <concepts::grid grid_t>
struct bin_view : public detray::ranges::view_interface<bin_view<grid_t>> {
  /// Cartesian product view over the local bin index sequences
  using bin_indexer_t = decltype(get_bin_indexer(
      std::declval<axis::multi_bin_range<grid_t::dim>>(),
      std::declval<std::make_index_sequence<grid_t::dim>>()));

  using iterator_t =
      bin_iterator<grid_t, detray::ranges::iterator_t<bin_indexer_t>>;
  using value_t = std::iter_value_t<iterator_t>;

  /// Default constructor
  constexpr bin_view() = default;

  /// Construct from a @param search_window of local bin index ranges and an
  /// underlying @param grid
  DETRAY_HOST_DEVICE constexpr explicit bin_view(
      const grid_t &grid, axis::multi_bin_range<grid_t::dim> &search_window)
      : m_grid{&grid},
        m_bin_indexer{get_bin_indexer(
            search_window,
            std::make_integer_sequence<std::size_t, grid_t::dim>{})} {}

  /// Copy constructor
  DETRAY_HOST_DEVICE
  constexpr bin_view(const bin_view &other)
      : m_grid{other.m_grid}, m_bin_indexer{other.m_bin_indexer} {}

  /// Default destructor
  ~bin_view() = default;

  /// Copy assignment operator
  DETRAY_HOST_DEVICE
  bin_view &operator=(const bin_view &other) {
    m_grid = other.m_grid;
    m_bin_indexer = other.m_bin_indexer;
    return *this;
  }

  /// @returns start position: first local bin index
  DETRAY_HOST_DEVICE
  constexpr auto begin() const -> iterator_t {
    return {m_grid, detray::ranges::begin(m_bin_indexer)};
  }

  /// @returns sentinel of the range: last local bin index
  DETRAY_HOST_DEVICE
  constexpr auto end() const -> iterator_t {
    return {m_grid, detray::ranges::end(m_bin_indexer)};
  }

  /// @returns number of all bins in the search area
  DETRAY_HOST_DEVICE
  constexpr auto size() const noexcept -> std::size_t {
    return m_bin_indexer.size();
  }

 private:
  /// The underlying grid that holds the bins
  const grid_t *m_grid{nullptr};
  /// How to index the bins in the search window (produces local indices)
  bin_indexer_t m_bin_indexer{};
};

/// @brief Iterate through the bin search area.
template <concepts::grid grid_t, std::input_iterator bin_indexer_t>
struct bin_iterator {
  using difference_type = std::ptrdiff_t;
  using value_type = typename grid_t::bin_type;
  using pointer = value_type *;
  using reference = value_type;
  using iterator_category = detray::ranges::bidirectional_iterator_tag;

  /// Default constructor required by LegacyIterator trait
  constexpr bin_iterator() = default;

  /// Construct from a bin indexing prescription @param bin_indexer and a
  /// @param grid
  DETRAY_HOST_DEVICE
  constexpr bin_iterator(const grid_t *grid, bin_indexer_t &&bin_indexer)
      : m_grid(grid), m_bin_indexer(std::move(bin_indexer)) {
    map_circular(*m_bin_indexer, m_lbin,
                 std::make_integer_sequence<std::size_t, grid_t::dim>{});
  }

  /// @returns true if it points to the same local bin.
  DETRAY_HOST_DEVICE constexpr bool operator==(const bin_iterator &rhs) const {
    return (m_bin_indexer == rhs.m_bin_indexer);
  }

  /// Increment to find next local bin index.
  DETRAY_HOST_DEVICE auto operator++() -> bin_iterator & {
    ++m_bin_indexer;

    // Get the correct local bin index
    map_circular(*m_bin_indexer, m_lbin,
                 std::make_integer_sequence<std::size_t, grid_t::dim>{});

    return *this;
  }

  /// Increment to find next local bin index (postfix)
  DETRAY_HOST_DEVICE constexpr bin_iterator operator++(int) {
    auto tmp(*this);
    ++(*this);
    return tmp;
  }

  /// Decrement to find previous local bin index.
  DETRAY_HOST_DEVICE auto operator--() -> bin_iterator & {
    --m_bin_indexer;

    // Get the correct local bin index
    map_circular(*m_bin_indexer, m_lbin,
                 std::make_integer_sequence<std::size_t, grid_t::dim>{});

    return *this;
  }

  /// Decrement to find previous local bin index (posfix)
  DETRAY_HOST_DEVICE constexpr bin_iterator operator--(int) {
    auto tmp(*this);
    --(*this);
    return tmp;
  }

  /// @returns the bin that corresponds to the current local bin index - const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) operator*() const {
    // Fetch the bin
    return m_grid->bin(m_lbin);
  }

 private:
  /// The iota range that is generated for circular axes does not map to their
  /// local bin index range, yet
  template <typename... Idx_t, std::size_t... I>
  DETRAY_HOST_DEVICE constexpr void map_circular(
      std::tuple<Idx_t...> index_tuple, typename grid_t::loc_bin_index &lbin,
      std::index_sequence<I...> /*seq*/) const {
    // Run the mapping for every axis in the grid
    (map_circular(m_grid->template get_axis<I>(), index_tuple, lbin), ...);
  }

  /// Map the local bin index for the phi axis to a periodic range and fill
  /// the local bin @param lbin
  template <typename axis_t, typename... Idx_t>
  DETRAY_HOST_DEVICE constexpr void map_circular(
      const axis_t &ax, std::tuple<Idx_t...> index_tuple,
      typename grid_t::loc_bin_index &lbin) const {
    constexpr auto loc_idx{
        static_cast<std::size_t>(axis_t::bounds_type::label)};

    if constexpr (axis_t::bounds_type::type == axis::bounds::e_circular) {
      lbin[loc_idx] = static_cast<dindex>(
          axis::circular<>{}.wrap(std::get<loc_idx>(index_tuple), ax.nbins()));
    } else {
      // All other axes start with a range that is already mapped
      lbin[loc_idx] = static_cast<dindex>(std::get<loc_idx>(index_tuple));
    }
  }

  /// Grid
  const grid_t *m_grid{nullptr};
  /// Bin indexing (cartesian product over local bin index ranges)
  bin_indexer_t m_bin_indexer{};
  /// Current local bin
  typename grid_t::loc_bin_index m_lbin{};
};

}  // namespace detray::axis::detail
