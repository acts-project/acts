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
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/grid/axis.hpp"
#include "detray/utils/grid/detail/axis_helpers.hpp"
#include "detray/utils/grid/detail/bin_storage.hpp"
#include "detray/utils/grid/populators.hpp"
#include "detray/utils/grid/serializers.hpp"
#include "detray/utils/ranges.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray {

/// @brief An N-dimensional grid for object storage.
///
/// @tparam axes_t the types of grid axes
/// @tparam bin_t type of bin in the (global) bin storage.
/// @tparam serializer_t how to serialize axis-local bin indices into global bin
///                      indices in the grid backend storage and vice versa.
template <typename axes_t, typename bin_t,
          template <std::size_t> class serializer_t = simple_serializer>
class grid_impl {
 public:
  /// Grid dimension
  static constexpr unsigned int dim = axes_t::dim;

  /// Single value in a bin entry
  using bin_type = bin_t;
  using value_type = typename bin_type::entry_type;

  template <std::size_t DIM>
  using serializer_type = serializer_t<DIM>;

  /// The type of the multi-axis is tied to the type of the grid: a non-
  /// owning grid holds a non-owning multi-axis member.
  using axes_type = axes_t;
  using glob_bin_index = dindex;
  using loc_bin_index = typename axes_type::loc_bin_index;
  using local_frame_type = typename axes_type::local_frame_type;
  using algebra_type = typename axes_type::algebra_type;
  using scalar_type = dscalar<algebra_type>;
  using point_type = typename axes_type::point_type;

  static constexpr bool is_owning{axes_type::is_owning};

  /// Backend storage type for the grid
  using bin_storage =
      detray::detail::bin_storage<is_owning, bin_type,
                                  typename axes_type::container_types>;
  using bin_container_type = typename bin_storage::bin_container_type;

  /// Vecmem based grid view type
  using view_type = dmulti_view<typename bin_storage::view_type,
                                typename axes_type::view_type>;
  /// Vecmem based grid view type - const
  using const_view_type = dmulti_view<typename bin_storage::const_view_type,
                                      typename axes_type::const_view_type>;
  /// Vecmem based buffer type
  using buffer_type = dmulti_buffer<typename bin_storage::buffer_type,
                                    typename axes_type::buffer_type>;

  /// Find the corresponding (non-)owning grid type
  template <bool owning>
  using type =
      grid_impl<typename axes_t::template type<owning>, bin_t, serializer_t>;

  /// Make grid default constructible: Empty grid with empty axis
  grid_impl() = default;

  /// Create empty grid with empty axes from specific vecmem memory resource
  DETRAY_HOST
  explicit grid_impl(vecmem::memory_resource &resource)
      : m_bins(resource), m_axes(resource) {}

  /// Create grid with well defined @param axes and @param bins_data - move
  DETRAY_HOST_DEVICE
  grid_impl(bin_container_type &&bin_data, axes_type &&axes)
      : m_bins(std::move(bin_data)), m_axes(std::move(axes)) {}

  /// Create grid from container pointers - non-owning (both grid and axes)
  DETRAY_HOST_DEVICE
  grid_impl(const bin_container_type *bin_data_ptr, const axes_type &axes,
            const dindex offset = 0)
      : m_bins(*bin_data_ptr, offset, axes.nbins()), m_axes(axes) {}

  /// Create grid from container pointers - non-owning (both grid and axes)
  DETRAY_HOST_DEVICE
  grid_impl(bin_container_type *bin_data_ptr, axes_type &axes,
            const dindex offset = 0u)
      : m_bins(*bin_data_ptr, offset, axes.nbins()), m_axes(axes) {}

  /// Create grid from container pointers - non-owning (both grid and axes)
  // TODO: const correctness
  DETRAY_HOST_DEVICE
  grid_impl(const bin_container_type *bin_data_ptr, axes_type &&axes,
            const dindex offset = 0)
      : m_bins(*(const_cast<bin_container_type *>(bin_data_ptr)), offset,
               axes.nbins()),
        m_axes(std::move(axes)) {}

  /// Create grid from container pointers - non-owning (both grid and axes)
  DETRAY_HOST_DEVICE
  grid_impl(bin_container_type *bin_data_ptr, axes_type &&axes,
            const dindex offset = 0u)
      : m_bins(*bin_data_ptr, offset, axes.nbins()), m_axes(std::move(axes)) {}

  /// Device-side construction from a vecmem based view type
  template <concepts::device_view grid_view_t>
  DETRAY_HOST_DEVICE explicit grid_impl(grid_view_t &view)
      : m_bins(detray::detail::get<0>(view.m_view)),
        m_axes(detray::detail::get<1>(view.m_view)) {}

  /// @returns the multi-axis used by the grid - const
  DETRAY_HOST_DEVICE
  auto axes() const -> const axes_type & { return m_axes; }

  /// @returns the grid local coordinate system
  DETRAY_HOST_DEVICE
  static constexpr auto get_local_frame() -> local_frame_type { return {}; }

  /// Transform a point in global cartesian coordinates to bound coordinates
  ///
  /// @param trf the placement transform of the grid (e.g. from a volume or
  ///            a surface).
  /// @param p   the point in global coordinates
  /// @param d   direction of a track at position p
  ///
  /// @returns a point in the coordinate system that is spanned by the grid's
  /// axes.
  template <concepts::transform3D transform3_t, concepts::point3D point3_t,
            concepts::vector3D vector3_t>
  DETRAY_HOST_DEVICE point_type project(const transform3_t &trf,
                                        const point3_t &p,
                                        const vector3_t &d) const {
    return grid_impl::get_local_frame().global_to_local(trf, p, d);
  }

  /// @returns an axis object, corresponding to the index.
  template <std::size_t index>
  DETRAY_HOST_DEVICE constexpr auto get_axis() const {
    return m_axes.template get_axis<index>();
  }

  /// @returns an axis object, corresponding to the label.
  template <axis::label L>
  DETRAY_HOST_DEVICE constexpr auto get_axis() const {
    return m_axes.template get_axis<L>();
  }

  /// @returns an axis object of the given type.
  template <typename axis_t>
  DETRAY_HOST_DEVICE constexpr axis_t get_axis() const {
    return m_axes.template get_axis<axis_t>();
  }

  /// @returns the total number of bins in the grid
  DETRAY_HOST_DEVICE constexpr auto nbins() const -> dindex {
    return m_axes.nbins();
  }

  /// @returns the total number of values in the grid
  /// @note this has to query every bin for the number of elements
  DETRAY_HOST_DEVICE constexpr auto size() const -> dindex {
    return static_cast<dindex>(all().size());
  }

  /// @returns an instance of the grid serializer
  static constexpr auto serializer() -> serializer_t<dim> { return {}; }

  /// @returns a local multi-bin index from a global bin index @param gid
  constexpr auto deserialize(const glob_bin_index gid) const -> loc_bin_index {
    return serializer()(axes(), gid);
  }

  /// @returns a global bin index from a local bin index @param mbin
  constexpr auto serialize(const loc_bin_index &mbin) const -> glob_bin_index {
    return serializer()(axes(), mbin);
  }

  /// @returns the full range of bins - const
  DETRAY_HOST_DEVICE
  auto bins() const -> const bin_storage & { return m_bins; }

  /// @returns the full range of bins
  DETRAY_HOST_DEVICE
  auto bins() -> bin_storage & { return m_bins; }

  /// @returns the iterable view of the bin content
  /// @{
  /// @param gbin the global bin index - const
  DETRAY_HOST_DEVICE
  decltype(auto) bin(const glob_bin_index gbin) const { return m_bins[gbin]; }

  /// @param gbin the global bin index
  DETRAY_HOST_DEVICE
  decltype(auto) bin(const glob_bin_index gbin) { return m_bins[gbin]; }

  /// @param mbin the multi-index of bins over all axes - const
  DETRAY_HOST_DEVICE
  decltype(auto) bin(const loc_bin_index &mbin) const {
    return bin(serialize(mbin));
  }

  /// @param mbin the multi-index of bins over all axes
  DETRAY_HOST_DEVICE
  decltype(auto) bin(const loc_bin_index &mbin) { return bin(serialize(mbin)); }

  /// Find the value of a single bin - const
  ///
  /// @param p is point in the local (bound) frame
  ///
  /// @return the iterable view of the bin content
  DETRAY_HOST_DEVICE decltype(auto) bin(const point_type &p) const {
    return bin(m_axes.bins(p));
  }

  /// Find the value of a single bin
  ///
  /// @param p is point in the local (bound) frame
  ///
  /// @return the iterable view of the bin content
  DETRAY_HOST_DEVICE decltype(auto) bin(const point_type &p) {
    return bin(m_axes.bins(p));
  }

  /// @param indices the single indices corresponding to a multi_bin
  template <concepts::index... I>
    requires(sizeof...(I) == dim)
  DETRAY_HOST_DEVICE decltype(auto) bin(I... indices) const {
    return bin(loc_bin_index{indices...});
  }
  /// @}

  /// Access a single entry in a bin from the global bin index, as well as
  /// the index of the entry in the bin.
  ///
  /// @param idx the index of a specific grid entry
  ///
  /// @returns a single bin entry.
  /// @{
  DETRAY_HOST_DEVICE
  decltype(auto) at(const glob_bin_index gbin, const dindex idx) {
    return bin(gbin)[idx];
  }
  DETRAY_HOST_DEVICE
  decltype(auto) at(const glob_bin_index gbin, const dindex idx) const {
    return bin(gbin)[idx];
  }
  DETRAY_HOST_DEVICE
  decltype(auto) at(const loc_bin_index &mbin, const dindex idx) {
    return bin(mbin)[idx];
  }
  DETRAY_HOST_DEVICE
  decltype(auto) at(const loc_bin_index &mbin, const dindex idx) const {
    return bin(mbin)[idx];
  }
  /// @}

  /// @returns a view over the flattened bin content by joining the bin ranges
  DETRAY_HOST_DEVICE auto all() { return detray::views::join(bins()); }

  /// @returns a view over the flattened bin content by joining the bin ranges
  DETRAY_HOST_DEVICE auto all() const { return detray::views::join(bins()); }

  /// Poupulate a bin with a single one of its corresponding values @param v
  /// @{
  /// @param mbin the multi bin index to be populated
  template <typename populator_t, typename V = value_type>
  DETRAY_HOST_DEVICE void populate(const loc_bin_index &mbin, V &&v) {
    populator_t{}(bin(mbin), std::forward<V>(v));
  }

  /// @param gbin the global bin index to be populated
  template <typename populator_t, typename V = value_type>
  DETRAY_HOST_DEVICE void populate(const glob_bin_index gbin, V &&v) {
    populator_t{}(bin(gbin), std::forward<V>(v));
  }

  /// @param p the point in local coordinates that defines the bin to be
  ///          populated
  template <typename populator_t, typename V = value_type>
  DETRAY_HOST_DEVICE void populate(const point_type &p, V &&v) {
    populator_t{}(bin(m_axes.bins(p)), std::forward<V>(v));
  }
  /// @}

  /// @returns view of a grid, including the grids multi_axis. Also valid if
  /// the value type of the grid is cv qualified (then value_t propagates
  /// quialifiers) - non-const
  template <bool owning = is_owning>
    requires owning
  DETRAY_HOST auto get_data() -> view_type {
    return view_type{detray::get_data(m_bins), detray::get_data(m_axes)};
  }

  /// @returns view of a grid, including the grids multi_axis. Also valid if
  /// the value type of the grid is cv qualified (then value_t propagates
  /// quialifiers) - const
  template <bool owning = is_owning>
    requires owning
  DETRAY_HOST auto get_data() const -> const_view_type {
    return const_view_type{detray::get_data(m_bins), detray::get_data(m_axes)};
  }

  /// Equality comparison
  ///
  /// @param rhs the right-hand side of the comparison
  ///
  /// @note grids could have different bin storage ranges, but could still be
  /// identical, hence compare the actual grid content
  ///
  /// @returns whether the two grids are equal
  DETRAY_HOST_DEVICE constexpr auto operator==(const grid_impl &rhs) const
      -> bool {
    // Check axes: they need to be identical
    if (m_axes != rhs.m_axes) {
      return false;
    }
    // Loop over global bin index and compare the two
    for (glob_bin_index i = 0; i < nbins(); ++i) {
      if (bin(i) != rhs.bin(i)) {
        return false;
      }
    }
    return true;
  }

  /// @returns a string stream that prints the grid details
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os, const grid_impl &gr) {
    os << gr.axes() << std::endl;

    for (glob_bin_index gbin = 0; gbin < gr.nbins(); ++gbin) {
      os << "bin " << gbin << ": ";
      for (const auto &bin_content : gr.bin(gbin)) {
        os << bin_content << ", ";
      }
      os << std::endl;
    }

    return os;
  }

 private:
  /// Struct that contains the grid's data state
  bin_storage m_bins{};
  /// The axes of the grid
  axes_type m_axes{};
};

/// Type alias for easier construction
template <concepts::algebra algebra_t, typename axes_t, typename bin_t,
          template <std::size_t> class serializer_t = simple_serializer,
          typename containers = host_container_types, bool ownership = true>
using grid =
    grid_impl<coordinate_axes<axes_t, algebra_t, ownership, containers>, bin_t,
              simple_serializer>;

}  // namespace detray
