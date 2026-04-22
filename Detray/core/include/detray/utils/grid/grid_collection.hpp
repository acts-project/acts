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

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/grid/grid.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray {

/// @brief A collection of grids that can be moved to device as part of the
/// detector
///
/// @tparam grid_t The type of grid in this collection. Must be non-owning, so
///                that the grid collection can manage the underlying memory.
template <concepts::grid grid_t, typename = void>
class grid_collection {};

/// Specialization for @c detray::grid
template <typename axes_t, typename bin_t,
          template <std::size_t> class serializer_t>
  requires(!detray::grid_impl<axes_t, bin_t, serializer_t>::is_owning)
class grid_collection<detray::grid_impl<axes_t, bin_t, serializer_t>> {
  using grid_t = detray::grid_impl<axes_t, bin_t, serializer_t>;
  using const_grid_t =
      const detray::grid_impl<const axes_t, const bin_t, serializer_t>;
  using multi_axis_t = typename grid_t::axes_type;

 public:
  using value_type = grid_t;
  using size_type = dindex;

  /// Backend storage type for the grid
  using bin_container_type = typename grid_t::bin_container_type;
  /// Offsets into the bin edges container
  using edge_offset_container_type =
      typename multi_axis_t::edge_offset_container_type;
  /// Contains all bin edges for all axes
  using edges_container_type = typename multi_axis_t::edges_container_type;
  template <typename T>
  using vector_type = typename multi_axis_t::template vector_type<T>;

 private:
  /// @brief Iterator for the grid collection: Generates grids on the fly.
  struct iterator {
    using difference_type = std::ptrdiff_t;
    using value_type = grid_t;
    using pointer = grid_t *;
    using reference = grid_t &;
    using iterator_category = detray::ranges::bidirectional_iterator_tag;

    /// Parametrized Constructor from a grid collection and a start index
    DETRAY_HOST_DEVICE
    constexpr explicit iterator(const grid_collection &grid_coll, dindex i = 0u)
        : m_i{i}, m_grid_coll{grid_coll} {}

    constexpr iterator(const iterator &other) = default;
    constexpr iterator(iterator &&other) noexcept = default;
    constexpr iterator &operator=(const iterator &rhs) noexcept = default;

    /// @returns true if grid indices are the same
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const iterator &rhs) const -> bool {
      return (m_i == rhs.m_i);
    }

    /// Increment the grid index
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> iterator & {
      ++m_i;
      return *this;
    }

    /// Decrement the grid index
    DETRAY_HOST_DEVICE
    constexpr auto operator--() -> iterator & {
      --m_i;
      return *this;
    }

    /// @returns the grid instance for the current position in the
    /// collection
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const -> value_type { return m_grid_coll[m_i]; }

   private:
    /// Current value of sequence
    dindex m_i;
    /// The grid collection
    const grid_collection &m_grid_coll;
  };

 public:
  /// Vecmem based grid collection view type
  using view_type = dmulti_view<dvector_view<size_type>,
                                detail::get_view_t<bin_container_type>,
                                detail::get_view_t<edge_offset_container_type>,
                                detail::get_view_t<edges_container_type>>;

  /// Vecmem based grid collection view type
  using const_view_type =
      dmulti_view<dvector_view<const size_type>,
                  detail::get_view_t<const bin_container_type>,
                  detail::get_view_t<const edge_offset_container_type>,
                  detail::get_view_t<const edges_container_type>>;

  /// Vecmem based buffer type
  using buffer_type =
      dmulti_buffer<dvector_buffer<size_type>,
                    detail::get_buffer_t<bin_container_type>,
                    detail::get_buffer_t<edge_offset_container_type>,
                    detail::get_buffer_t<edges_container_type>>;

  /// Make grid collection default constructible: Empty
  grid_collection() = default;

  /// Create empty grid collection from specific vecmem memory resource
  DETRAY_HOST
  explicit grid_collection(vecmem::memory_resource *resource)
      : m_bin_offsets(resource),
        m_bins(resource),
        m_bin_edge_offsets(resource),
        m_bin_edges(resource) {}

  /// Create grid collection from existing data - move
  DETRAY_HOST_DEVICE
  grid_collection(vector_type<size_type> &&offs, bin_container_type &&bins,
                  edge_offset_container_type &&edge_offs,
                  edges_container_type &&edges)
      : m_bin_offsets(std::move(offs)),
        m_bins(std::move(bins)),
        m_bin_edge_offsets(std::move(edge_offs)),
        m_bin_edges(std::move(edges)) {}

  /// Device-side construction from a vecmem based view type
  template <concepts::device_view coll_view_t>
  DETRAY_HOST_DEVICE explicit grid_collection(coll_view_t &view)
      : m_bin_offsets(detail::get<0>(view.m_view)),
        m_bins(detail::get<1>(view.m_view)),
        m_bin_edge_offsets(detail::get<2>(view.m_view)),
        m_bin_edges(detail::get<3>(view.m_view)) {}

  /// Move constructor
  DETRAY_HOST_DEVICE grid_collection(grid_collection &&other) noexcept
      : m_bin_offsets(std::move(other.m_bin_offsets)),
        m_bins(std::move(other.m_bins)),
        m_bin_edge_offsets(std::move(other.m_bin_edge_offsets)),
        m_bin_edges(std::move(other.m_bin_edges)) {}

  /// Move assignment
  DETRAY_HOST_DEVICE grid_collection &operator=(
      grid_collection &&other) noexcept {
    if (this != &other) {
      m_bin_offsets = std::move(other.m_bin_offsets);
      m_bins = std::move(other.m_bins);
      m_bin_edge_offsets = std::move(other.m_bin_edge_offsets);
      m_bin_edges = std::move(other.m_bin_edges);
    }
    return *this;
  }

  /// @returns the number of grids in the collection - const
  DETRAY_HOST_DEVICE
  constexpr auto size() const noexcept -> dindex {
    return static_cast<dindex>(m_bin_offsets.size());
  }

  /// @returns an iterator that points to the first grid
  DETRAY_HOST_DEVICE
  constexpr auto begin() const { return iterator{*this, 0u}; }

  /// @returns an iterator that points to the coll. end
  DETRAY_HOST_DEVICE
  constexpr auto end() const { return iterator{*this, size()}; }

  /// @returns the number of grids in the collection - const
  DETRAY_HOST_DEVICE
  constexpr auto empty() const noexcept -> bool {
    return m_bin_offsets.empty();
  }

  /// @brief Resize the underlying containers
  /// @note Not defined! The amount of memory can differ for every grid
  DETRAY_HOST_DEVICE
  constexpr void resize(std::size_t) noexcept { /*Not defined*/ }

  /// @brief Reserve memory
  /// @note Not defined! The amount of memory can differ for every grid
  DETRAY_HOST_DEVICE
  constexpr void reserve(std::size_t) noexcept { /*Not defined*/ }

  /// Removes all data from the grid collection containers
  DETRAY_HOST_DEVICE
  constexpr void clear() noexcept {
    m_bin_offsets.clear();
    m_bins.clear();
    m_bin_edge_offsets.clear();
    m_bin_edges.clear();
  }

  /// Insert a number of grids
  /// @note Not defined! There is no grid iterator implementation
  template <typename... Args>
  DETRAY_HOST_DEVICE constexpr void insert(Args &&...) noexcept {
    /*Not defined*/
  }

  /// @returns the offsets for the grids in the bin storage - const
  DETRAY_HOST_DEVICE
  constexpr auto offsets() const -> const vector_type<size_type> & {
    return m_bin_offsets;
  }

  /// @returns the underlying bin content storage - const
  DETRAY_HOST_DEVICE
  constexpr auto bin_storage() const -> const bin_container_type & {
    return m_bins;
  }

  /// @returns the underlying axis boundary storage - const
  DETRAY_HOST_DEVICE
  constexpr auto axes_storage() const -> const edge_offset_container_type & {
    return m_bin_edge_offsets;
  }

  /// @returns the underlying bin edges storage - const
  DETRAY_HOST_DEVICE
  constexpr auto bin_edges_storage() const -> const edges_container_type & {
    return m_bin_edges;
  }

  /// Create grid from container pointers - const
  DETRAY_HOST_DEVICE
  auto operator[](const size_type i) const -> grid_t {
    const size_type axes_offset{grid_t::dim * i};
    return grid_t(&m_bins,
                  multi_axis_t(m_bin_edge_offsets, m_bin_edges, axes_offset),
                  m_bin_offsets[i]);
  }

  /// Create grid from container pointers with range check
  DETRAY_HOST_DEVICE
  auto at(const size_type i) const -> grid_t {
    const size_type axes_offset{grid_t::dim * i};
    return grid_t(&m_bins,
                  multi_axis_t(m_bin_edge_offsets, m_bin_edges, axes_offset),
                  m_bin_offsets.at(i));
  }

  /// @returns a vecmem view on the grid collection data - non-const
  DETRAY_HOST auto get_data() -> view_type {
    return view_type{detray::get_data(m_bin_offsets), detray::get_data(m_bins),
                     detray::get_data(m_bin_edge_offsets),
                     detray::get_data(m_bin_edges)};
  }

  /// @returns a vecmem view on the grid collection data - const
  DETRAY_HOST
  auto get_data() const -> const_view_type {
    return const_view_type{
        detray::get_data(m_bin_offsets), detray::get_data(m_bins),
        detray::get_data(m_bin_edge_offsets), detray::get_data(m_bin_edges)};
  }

  /// Add a new grid @param gr to the collection.
  /// @note this takes a data owning grid to transcribe the data from.
  template <typename other_grid_t>
    requires std::constructible_from<typename grid_t::template type<true>,
                                     other_grid_t>
  DETRAY_HOST constexpr auto push_back(const other_grid_t &gr) noexcept(false)
      -> void {
    // Current offset into the global bin storage for the new grid
    m_bin_offsets.push_back(static_cast<size_type>(m_bins.size()));

    // Add the bins of the new grid to the collection
    insert_bin_data(m_bins, gr.bins());

    // Add the bin edge offsets of the new grid to the collection
    // (how to lookup the axis bin edges)
    const auto &bin_edge_offsets = gr.axes().bin_edge_offsets();
    m_bin_edge_offsets.insert(m_bin_edge_offsets.end(),
                              bin_edge_offsets.begin(), bin_edge_offsets.end());

    // Current offset into the global bin edges storage
    auto bin_edges_offset{static_cast<dindex>(m_bin_edges.size())};

    // Update the bin edges index offset for the axes in the grid collection
    const auto start_idx{m_bin_edge_offsets.size() - grid_t::dim};
    for (std::size_t i = start_idx; i < m_bin_edge_offsets.size(); ++i) {
      auto &bin_entry_range = m_bin_edge_offsets.at(i);
      bin_entry_range.shift(bin_edges_offset);
    }

    // Add the bin edges of the new grid to the collection
    const auto &bin_edges = gr.axes().bin_edges();
    m_bin_edges.insert(m_bin_edges.end(), bin_edges.begin(), bin_edges.end());
  }

 private:
  /// Insert data into a vector of bins
  DETRAY_HOST void insert_bin_data(
      vector_type<typename grid_t::bin_type> &bin_data,
      const typename grid_t::template type<true>::bin_storage &grid_bins) {
    bin_data.insert(bin_data.end(), grid_bins.begin(), grid_bins.end());
  }

  /// Insert data into the backend containers of a grid with dynamic bin
  /// capacities
  template <typename container_t>
  DETRAY_HOST void insert_bin_data(
      detray::detail::dynamic_bin_container<bin_t, container_t> &bin_data,
      const grid_t::template type<true>::bin_storage &grid_bins) {
    bin_data.append(grid_bins);
  }

  /// Offsets for the respective grids into the bin storage
  vector_type<size_type> m_bin_offsets{};
  /// Contains the bin content for all grids
  bin_container_type m_bins{};
  /// Contains the offsets/no. bins for the bin edges of all axes of all grids
  edge_offset_container_type m_bin_edge_offsets{};
  /// Contains the bin edges for all grids
  edges_container_type m_bin_edges{};
};

}  // namespace detray
