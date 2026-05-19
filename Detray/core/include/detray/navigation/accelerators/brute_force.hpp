// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Detray include(s).
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/algorithms.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/navigation/accelerators/search_window.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <type_traits>

namespace detray {

/// @brief A collection of brute force surface finders, callable by index.
///
/// This class fulfills all criteria to be used in the detector @c multi_store .
///
/// @tparam value_t the entry type in the collection (e.g. surface descriptors).
/// @tparam container_t the types of underlying containers to be used.
template <class value_t, typename container_t = host_container_types>
class brute_force_collection {
 public:
  template <typename T>
  using vector_type = typename container_t::template vector_type<T>;
  using size_type = dindex;

  /// A nested surface finder that returns all surfaces in a range (brute
  /// force). This type will be returned when the surface collection is
  /// queried for the surfaces of a particular volume.
  struct brute_forcer
      : public detray::ranges::subrange<const vector_type<value_t>> {
    using base = detray::ranges::subrange<const vector_type<value_t>>;
    using value_type = value_t;
    using query_type = bool;

    /// Default constructor
    brute_forcer() = default;

    /// Constructor from @param surface_range - move
    DETRAY_HOST_DEVICE constexpr brute_forcer(
        const vector_type<value_t>& surfaces, const dindex_range& range)
        : base(surfaces, range) {}

    /// @returns the complete surface range of the search volume
    DETRAY_HOST_DEVICE constexpr auto search(query_type /*p*/) const {
      return *this;
    }

    /// @returns the complete surface range of the search volume
    template <concepts::arithmetic window_size_t>
    DETRAY_HOST_DEVICE constexpr auto search(
        query_type /*p*/, search_window<window_size_t, 2> /*win_size*/) const {
      return *this;
    }

    /// @returns the complete surface range of the search volume
    template <typename detector_t, typename track_t,
              concepts::arithmetic window_size_t>
    DETRAY_HOST_DEVICE constexpr auto search(
        const detector_t& /*det*/,
        const typename detector_t::volume_type& /*volume*/,
        const track_t& /*track*/,
        const search_window<window_size_t, 2>& /*win_size*/,
        const typename detector_t::geometry_context& /*ctx*/) const {
      DETRAY_DEBUG_HOST("Brute force search...");
      return *this;
    }

    /// @returns the surface at a given index @param i - const
    DETRAY_HOST_DEVICE constexpr value_t at(const dindex i) const {
      return (*this)[i];
    }
    /// @returns the surface at a given index @param i - non-const
    DETRAY_HOST_DEVICE constexpr value_t& at(const dindex i) {
      return (*this)[i];
    }

    /// @returns an iterator over all surfaces in the data structure
    DETRAY_HOST_DEVICE constexpr auto all() const { return *this; }

    /// @returns an iterator over all surfaces in the data structure
    DETRAY_HOST_DEVICE auto all() { return *this; }

    /// @returns a string stream that prints the brute forcer details
    DETRAY_HOST
    friend std::ostream& operator<<(std::ostream& os, const brute_forcer& bf) {
      for (const auto& entry : bf.all()) {
        os << entry << std::endl;
      }

      return os;
    }
  };

  using value_type = brute_forcer;

  using view_type = dmulti_view<dvector_view<size_type>, dvector_view<value_t>>;
  using const_view_type =
      dmulti_view<dvector_view<const size_type>, dvector_view<const value_t>>;
  using buffer_type =
      dmulti_buffer<dvector_buffer<size_type>, dvector_buffer<value_t>>;

  /// Default constructor
  constexpr brute_force_collection() {
    // Start of first subrange
    m_offsets.push_back(0u);
  };

  /// Constructor from memory resource
  DETRAY_HOST
  explicit constexpr brute_force_collection(vecmem::memory_resource* resource)
      : m_offsets(resource), m_surfaces(resource) {
    // Start of first subrange
    m_offsets.push_back(0u);
  }

  /// Constructor from memory resource
  DETRAY_HOST
  explicit constexpr brute_force_collection(vecmem::memory_resource& resource)
      : brute_force_collection(&resource) {}

  /// Device-side construction from a vecmem based view type
  template <concepts::device_view coll_view_t>
  DETRAY_HOST_DEVICE explicit brute_force_collection(coll_view_t& view)
      : m_offsets(detail::get<0>(view.m_view)),
        m_surfaces(detail::get<1>(view.m_view)) {}

  /// @returns access to the volume offsets - const
  DETRAY_HOST const auto& offsets() const { return m_offsets; }

  /// @returns access to the volume offsets
  DETRAY_HOST auto& offsets() { return m_offsets; }

  /// @returns number of surface collections (at least on per volume) - const
  DETRAY_HOST_DEVICE
  constexpr auto size() const noexcept -> size_type {
    // The start index of the first range is always present
    return static_cast<dindex>(m_offsets.size()) - 1u;
  }

  /// @note outside of navigation, the number of elements is unknown
  DETRAY_HOST_DEVICE
  constexpr auto empty() const noexcept -> bool {
    return size() == size_type{0};
  }

  /// @return access to the surface container - const.
  DETRAY_HOST_DEVICE
  auto all() const -> const vector_type<value_t>& { return m_surfaces; }

  /// @return access to the surface container - non-const.
  DETRAY_HOST_DEVICE
  auto all() -> vector_type<value_t>& { return m_surfaces; }

  /// Create brute force surface finder from surface container - const
  DETRAY_HOST_DEVICE
  auto operator[](const size_type i) const -> value_type {
    return {m_surfaces, dindex_range{m_offsets[i], m_offsets[i + 1u]}};
  }

  /// Create brute force surface finder from surface container - const
  DETRAY_HOST_DEVICE
  auto operator[](const size_type i) -> value_type {
    return {m_surfaces, dindex_range{m_offsets[i], m_offsets[i + 1u]}};
  }

  /// Add a new surface collection
  template <detray::ranges::range sf_container_t>
    requires std::is_same_v<typename sf_container_t::value_type, value_t>
  DETRAY_HOST auto push_back(const sf_container_t& surfaces) noexcept(false)
      -> void {
    m_surfaces.reserve(m_surfaces.size() + surfaces.size());
    m_surfaces.insert(m_surfaces.end(), surfaces.begin(), surfaces.end());
    // End of this range is the start of the next range
    m_offsets.push_back(static_cast<dindex>(m_surfaces.size()));
  }

  /// Remove surface from collection
  DETRAY_HOST auto erase(typename vector_type<value_t>::iterator pos) noexcept(
      false) {
    // Remove one element
    auto next = m_surfaces.erase(pos);

    // Update the upper bound of the range and all following ranges
    const auto idx{static_cast<std::size_t>(pos - m_surfaces.begin())};
    auto offset = detray::upper_bound(m_offsets.begin(), m_offsets.end(), idx);
    for (auto itr = offset; itr != m_offsets.end(); ++itr) {
      --(*itr);
    }

    return next;
  }

  /// @return the view on the brute force finders - non-const
  DETRAY_HOST
  constexpr auto get_data() noexcept -> view_type {
    return view_type{detray::get_data(m_offsets), detray::get_data(m_surfaces)};
  }

  /// @return the view on the brute force finders - const
  DETRAY_HOST
  constexpr auto get_data() const noexcept -> const_view_type {
    return const_view_type{detray::get_data(m_offsets),
                           detray::get_data(m_surfaces)};
  }

 private:
  /// Offsets for the respective volumes into the surface storage
  vector_type<size_type> m_offsets{};
  /// The storage for all surface handles
  vector_type<value_t> m_surfaces{};
};

}  // namespace detray
