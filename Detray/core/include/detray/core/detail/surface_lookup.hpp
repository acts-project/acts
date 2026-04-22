// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/identifier.hpp"
#include "detray/utils/logging.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <iostream>
#include <type_traits>

namespace detray {

/// General case: Brute force search for the corresponding sf-descriptor
struct default_searcher {
  template <typename source_link_contianer_t>
  auto operator()(const source_link_contianer_t &sf_container) {
    // Check that this searcher can be used on the passed surface container
    static_assert(
        std::is_same_v<decltype(sf_container[0].source), std::uint64_t>,
        "Source link searcher not compatible with detector");

    // Cannot assume any sorting
    for (const auto &sf : sf_container) {
      if (sf.source == m_source) {
        return sf;
      }
    }

    return typename source_link_contianer_t::value_type{};
  }

  /// The query source link
  std::uint64_t m_source;
};

/// Couple the surface descriptor to a source link
template <typename sf_desc_t>
struct source_link : sf_desc_t {
  source_link() = default;

  source_link(sf_desc_t sf_desc, std::uint64_t src)
      : sf_desc_t{sf_desc}, source{src} {}

  std::uint64_t source{detail::invalid_value<std::uint64_t>()};
};

/// @brief Wraps a vector-like container that holds the surface descriptors of a
/// detector and makes them searchable by index and source link.
///
/// @tparam sf_desc_t The surface descriptor type
/// @tparam container_t The type of container to use for the descriptor
/// collection.
template <typename sf_desc_t,
          template <typename...> class container_t = dvector>
class surface_lookup {
 public:
  /// Underlying container type that can handle vecmem views
  using base_type = container_t<source_link<sf_desc_t>>;
  using size_type = typename base_type::size_type;
  using value_type = typename base_type::value_type;
  using iterator = typename base_type::iterator;
  using const_iterator = typename base_type::const_iterator;

  /// Vecmem view types
  using view_type = detail::get_view_t<container_t<source_link<sf_desc_t>>>;
  using const_view_type =
      detail::get_view_t<const container_t<source_link<sf_desc_t>>>;
  using buffer_type = detail::get_buffer_t<container_t<source_link<sf_desc_t>>>;

  /// Empty container
  constexpr surface_lookup() = default;

  /// Construct with a specific memory resource @param resource
  /// (host-side only)
  template <typename allocator_t = vecmem::memory_resource>
    requires(!concepts::device_view<allocator_t>)
  DETRAY_HOST explicit surface_lookup(allocator_t &resource)
      : m_container(&resource) {}

  /// Copy Construct with a specific memory resource @param resource
  /// (host-side only)
  template <typename allocator_t = vecmem::memory_resource,
            typename C = container_t<source_link<sf_desc_t>>>
    requires std::is_same_v<C, std::vector<source_link<sf_desc_t>>>
  DETRAY_HOST explicit surface_lookup(allocator_t &resource,
                                      const source_link<sf_desc_t> &arg)
      : m_container(&resource, arg) {}

  /// Construct from the container @param view . Mainly used device-side.
  template <concepts::device_view container_view_t>
  DETRAY_HOST_DEVICE explicit surface_lookup(container_view_t &view)
      : m_container(view) {}

  /// @returns the size of the underlying container
  DETRAY_HOST_DEVICE
  constexpr auto size() const noexcept -> dindex {
    return static_cast<dindex>(m_container.size());
  }

  /// @returns true if the underlying container is empty
  DETRAY_HOST_DEVICE
  constexpr auto empty() const noexcept -> bool { return m_container.empty(); }

  /// Reserve memory of size @param n for a given geometry context
  DETRAY_HOST void reserve(std::size_t n) { m_container.reserve(n); }

  /// Resize the underlying container to @param n for a given geometry context
  DETRAY_HOST void resize(std::size_t n) { m_container.resize(n); }

  /// Removes and destructs all elements in the container.
  DETRAY_HOST void clear() { m_container.clear(); }

  /// @returns the collections iterator at the start position.
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) begin() { return m_container.begin(); }

  /// @returns the collections iterator sentinel.
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) end() { return m_container.end(); }

  /// @returns the collections iterator at the start position - const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) begin() const { return m_container.begin(); }

  /// @returns the collections iterator sentinel - const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) end() const { return m_container.end(); }

  /// @returns the reverse iterator at the start position - const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) rbegin() const { return m_container.rbegin(); }

  /// @returns the reverse iterator sentinel - const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) rend() const { return m_container.rend(); }

  /// Elementwise access. Needs @c operator[] for storage type - non-const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) operator[](const std::size_t i) {
    assert(i < m_container.size());
    return m_container[i];
  }

  /// Elementwise access. Needs @c operator[] for storage type - const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) operator[](const std::size_t i) const {
    assert(i < m_container.size());
    return m_container[i];
  }

  /// @returns context based access to an element (also range checked)
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) at(const dindex i) noexcept {
    assert(i < m_container.size());
    return m_container.at(i);
  }

  /// @returns context based access to an element (also range checked) - const
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) at(const dindex i) const noexcept {
    assert(i < m_container.size());
    return m_container.at(i);
  }

  /// @returns the surface descriptor according to the global surface index
  /// @param sf_index
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) search(dindex sf_index) const {
    assert(sf_index < m_container.size());
    return m_container[sf_index];
  }

  /// @returns the surface descriptor according to the surface identifier
  /// @param geo_id
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) search(geometry::identifier geo_id) const {
    return search(geo_id.index());
  }

  /// @returns the surface descriptor according to the searcher passed as
  /// @param source_searcher
  template <typename searcher_t = default_searcher>
  DETRAY_HOST_DEVICE constexpr decltype(auto) search(
      searcher_t &&source_searcher) const {
    return source_searcher(m_container);
  }

  /// Add a new element to the collection
  ///
  /// @param sf_desc the surface descriptor
  /// @param src the source index
  DETRAY_HOST constexpr auto push_back(sf_desc_t sf_desc,
                                       std::uint64_t src) noexcept(false)
      -> void {
    m_container.push_back({sf_desc, src});
  }

  /// Add a new element to the collection - copy
  ///
  /// @param sf_link the detray source link
  DETRAY_HOST constexpr auto push_back(source_link<sf_desc_t> sf_link) noexcept(
      false) -> void {
    m_container.push_back(sf_link);
  }

  /// Insert a surface descriptor @param sf_desc and its source index
  /// @param src into the container
  DETRAY_HOST void insert(
      sf_desc_t sf_desc,
      std::uint64_t src =
          detail::invalid_value<std::uint64_t>()) noexcept(false) {
    insert({sf_desc, src});
  }

  /// Insert a source link @param sf_link at the position of its surface
  /// index.
  DETRAY_HOST void insert(source_link<sf_desc_t> sf_link) noexcept(false) {
    if (detail::is_invalid_value(sf_link.index())) {
      DETRAY_ERROR_HOST("Invalid surface descriptor: " << sf_link);
    }
    if (m_container.size() <= sf_link.index()) {
      m_container.resize(sf_link.index() + 1u);
    }
    m_container.at(sf_link.index()) = sf_link;
  }

  /// @return the view on the underlying container - non-const
  DETRAY_HOST auto get_data() -> view_type {
    return detray::get_data(m_container);
  }

  /// @return the view on the underlying container - const
  DETRAY_HOST auto get_data() const -> const_view_type {
    return detray::get_data(m_container);
  }

 private:
  /// The underlying container implementation
  base_type m_container;
};

}  // namespace detray
