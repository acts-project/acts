// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/navigation/accelerators/concepts.hpp"
#include "detray/navigation/accelerators/search_window.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/grid/detail/axis_helpers.hpp"
#include "detray/utils/grid/detail/bin_view.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/grid/grid_bins.hpp"
#include "detray/utils/grid/grid_collection.hpp"

namespace detray {

/// @brief An N-dimensional spatial grid for geometry object searches.
template <concepts::grid grid_t>
class spatial_grid_impl : public grid_t {
  using base_grid = grid_t;
  using frame_t = typename grid_t::local_frame_type;
  using algebra_t = typename grid_t::algebra_type;

 public:
  using value_type = typename base_grid::value_type;
  using query_type = typename base_grid::point_type;

  /// Find the corresponding (non-)owning grid type
  template <bool is_owning>
  using type = spatial_grid_impl<typename grid_t::template type<is_owning>>;

  using mask_type =
      detray::mask<concentric_cylinder2D, algebra_t, std::uint8_t>;

  /// Use all of the grid constructors
  using base_grid::base_grid;

  /// Construct from existing grid - move
  ///
  /// @param mask_values additional mask boundary values that cannot be
  /// inferred from the axis spans, e.g. the concentric cylinder radius
  template <concepts::scalar... Args>
  DETRAY_HOST_DEVICE explicit constexpr spatial_grid_impl(base_grid &&gr,
                                                          Args &&...mask_values)
      : base_grid(std::move(gr)), m_mask{get_mask_from_axes(mask_values...)} {}

  /// Construct from existing base grid - copy
  template <concepts::scalar... Args>
  DETRAY_HOST_DEVICE explicit constexpr spatial_grid_impl(base_grid &gr,
                                                          Args &&...mask_values)
      : base_grid(gr), m_mask{get_mask_from_axes(mask_values...)} {}

  /// Construct from existing base grid and mask - copy
  DETRAY_HOST_DEVICE constexpr spatial_grid_impl(base_grid &&gr,
                                                 const mask_type &mask)
      : base_grid(std::move(gr)), m_mask{mask} {}

  /// @return the mask of the (virtual) reference surface
  DETRAY_HOST_DEVICE constexpr const mask_type &mask() const { return m_mask; }

  /// Set a new mask @param m
  /// @TODO: Check against result of @c get_mask_from_axes()
  DETRAY_HOST_DEVICE constexpr void mask(const mask_type &m) { m_mask = m; }

  /// Find the value of a single bin - const
  ///
  /// @param p is point in the local (bound) frame
  ///
  /// @return the iterable view of the bin content
  DETRAY_HOST_DEVICE decltype(auto) search(const query_type &p) const {
    DETRAY_DEBUG_HOST(" -> lookup pos: " << p);
    return this->bin(p);
  }

  /// Find the value of a single bin
  ///
  /// @param p is point in the local (bound) frame
  ///
  /// @return the iterable view of the bin content
  DETRAY_HOST_DEVICE decltype(auto) search(const query_type &p) {
    DETRAY_DEBUG_HOST(" -> lookup pos: " << p);
    return this->bin(p);
  }

  /// Straight line intersection with the grid
  ///
  /// @returns the point of intersection in grid local coordinates
  DETRAY_HOST_DEVICE query_type
  local_point_of_intersection(const detray::detail::ray<algebra_t> &tangential,
                              const dtransform3D<algebra_t> &trf) const {
    // Intersect the (virtual) reference surface of the grid to find
    // the correct bin
    using intersector_t =
        ray_intersector_impl<frame_t, algebra_t, intersection::contains_pos>;

    constexpr intersector_t intersector{};
    typename intersector_t::result_type result{};
    constexpr auto overstep_tol{
        -std::numeric_limits<dscalar<algebra_t>>::max()};

    // The cylinder intersector requires the radius from the mask
    if constexpr (concepts::cylindrical<frame_t>) {
      result = intersector.point_of_intersection(tangential, trf, m_mask,
                                                 overstep_tol);
    } else {
      result = intersector.point_of_intersection(tangential, trf, overstep_tol);
    }

    // Retrieve the closest intersection point and check its valid
    typename intersector_t::point_type intr_point;
    bool is_inside{true};
    if constexpr (intersector_t::n_solutions == 1) {
      // Intersection with the grid was found
      if (detray::detail::all_of(result.is_valid())) [[likely]] {
        intr_point = result.point;
      } else {
        is_inside = false;
      }
    } else {
      // Closest intersection with the grid was found
      if (detray::detail::all_of(result[0].is_valid())) [[likely]] {
        intr_point = result[0].point;
      } else {
        is_inside = false;
      }
    }
    // Use ray origin as query position
    if (!is_inside) [[unlikely]] {
      // In case the intersection results come in local coordinates
      if constexpr (std::same_as<decltype(intr_point), query_type>) {
        intr_point =
            frame_t::global_to_local(trf, tangential.pos(), tangential.dir());
      } else {
        intr_point = tangential.pos();
      }
    }

    // Most intersectors return global positions -> project to grid axes
    if constexpr (std::same_as<frame_t, concentric_cylindrical2D<algebra_t>>) {
      // The concentric cylinder intersector only returns the z-pos
      return {vector::phi(tangential.pos(result.path)), intr_point[1]};
    } else if constexpr (std::same_as<decltype(intr_point), query_type>) {
      return intr_point;
    } else {
      // Project the intersection point into the local grid frame
      return this->project(trf, intr_point, tangential.dir());
    }
  }

  /// @brief Return a neighborhood of values from the grid
  ///
  /// The lookup is done with a search window around the bin
  ///
  /// @param p is point in the local frame
  /// @param win_size size of the binned/scalar search window
  ///
  /// @return the sequence of values
  template <concepts::arithmetic window_size_t>
  DETRAY_HOST_DEVICE auto search(
      const query_type &p,
      const search_window<window_size_t, 2> &win_size) const {
    DETRAY_DEBUG_HOST(" -> lookup pos: " << p);
    DETRAY_DEBUG_HOST(" -> search window: [" << win_size[0] << ", "
                                             << win_size[1] << "]");

    // Return iterable over bins in the search window
    auto bin_ranges = this->axes().bin_ranges(p, win_size);
    auto search_area = axis::detail::bin_view(*this, bin_ranges);

    // Join the respective bins to a single iteration
    return detray::views::join(std::move(search_area));
  }

  /// Interface for the navigator
  template <typename detector_t, typename track_t,
            concepts::arithmetic window_size_t>
  DETRAY_HOST_DEVICE auto search(
      const detector_t &det, const typename detector_t::volume_type &volume,
      const track_t &track, const search_window<window_size_t, 2> &win_size,
      const typename detector_t::geometry_context &ctx) const {
    // Placement of the grid (same as volume)
    const auto &trf = det.transform_store().at(volume.transform(), ctx);

    query_type loc_pos{};
    // For 2-dimensional grids, project the track position along the track
    // direction to find the optimal grid bin
    if constexpr (base_grid::dim == 2) {
      if constexpr (concepts::cylindrical<frame_t>) {
        DETRAY_DEBUG_HOST("2D spatial grid: " << m_mask);
      } else {
        DETRAY_DEBUG_HOST("2D spatial grid: " << DETRAY_TYPENAME(frame_t));
      }

      // Tangential to the current track parameters
      const detray::detail::ray<algebra_t> tangential{track};

      // Intersect the (virtual) reference surface of the grid to find
      // the correct bin
      loc_pos = local_point_of_intersection(tangential, trf);
    } else {
      loc_pos = this->project(trf, track.pos(), track.dir());
    }

    // Grid lookup
    return search(loc_pos, win_size);
  }

 private:
  /// @returns a mask that has boundaries which match the grid axis spans
  template <typename... Args>
  DETRAY_HOST_DEVICE mask_type get_mask_from_axes(Args &&...mask_values) {
    constexpr auto inv_vol_link{detray::detail::invalid_value<std::uint8_t>()};

    /// @param mask_values contains the cylinder radius
    if constexpr (concepts::cylindrical<frame_t>) {
      static_assert(sizeof...(Args) == 1,
                    "Spatial cylinder grid: Only the cylinder radius "
                    "needs to be provided externally");

      constexpr auto inf{std::numeric_limits<dscalar<algebra_t>>::max()};
      return mask_type{inv_vol_link, mask_values..., -inf, inf};
    } else {
      /// @TODO: Implement axes to mask conversion for the other shapes
      return {};
    }
  }

  /// Struct that contains the grid's data state
  mask_type m_mask{};
};

template <concepts::algebra algebra_t, typename axes_t, typename bin_t,
          template <std::size_t> class serializer_t = simple_serializer,
          typename containers = host_container_types, bool ownership = true>
using spatial_grid = spatial_grid_impl<
    grid_impl<coordinate_axes<axes_t, algebra_t, ownership, containers>, bin_t,
              simple_serializer>>;

namespace concepts {

template <class SG>
concept spatial_grid =
    grid<SG> && (surface_accelerator<SG> || volume_accelerator<SG>) &&
    requires(const SG &sgr) {
      typename SG::mask_type;

      { sgr.mask() } -> std::same_as<const typename SG::mask_type &>;
    };

}  // namespace concepts

/// Accelerator collection specialization for @c detray::spatial_grid_impl
///
/// @note holds an additional container for mask values
template <concepts::grid grid_t>
  requires(!spatial_grid_impl<grid_t>::is_owning)
class grid_collection<spatial_grid_impl<grid_t>>
    : public grid_collection<grid_t> {
  // Use a normal grid collection for the grid related data
  using base_collection = grid_collection<grid_t>;

  using frame_t = typename grid_t::local_frame_type;
  using mask_t = typename spatial_grid_impl<grid_t>::mask_type;
  using scalar_t = dscalar<typename grid_t::algebra_type>;

  template <typename T>
  using vector_t = typename base_collection::template vector_type<T>;

 public:
  using size_type = typename base_collection::size_type;
  using value_type = spatial_grid_impl<grid_t>;

  /// Vecmem based grid collection view type
  using view_type =
      dmulti_view<detail::get_view_t<base_collection>, dvector_view<mask_t>>;

  /// Vecmem based grid collection view type
  using const_view_type = dmulti_view<detail::get_view_t<const base_collection>,
                                      dvector_view<const mask_t>>;

  /// Vecmem based buffer type
  using buffer_type = dmulti_buffer<detail::get_buffer_t<base_collection>,
                                    dvector_buffer<mask_t>>;

  /// Make spatial grid collection default constructible: Empty
  grid_collection() = default;

  /// Create empty spatial grid collection from specific vecmem memory
  /// resource
  DETRAY_HOST
  explicit grid_collection(vecmem::memory_resource *resource)
      : base_collection(resource), m_mask_values(resource) {}

  /// Create spatial grid collection from existing data - move
  DETRAY_HOST_DEVICE
  grid_collection(base_collection &&grid_coll, vector_t<scalar_t> &&mask_values)
      : base_collection(std::move(grid_coll)),
        m_mask_values(std::move(mask_values)) {}

  /// Device-side construction from a vecmem based view type
  template <concepts::device_view coll_view_t>
  DETRAY_HOST_DEVICE explicit grid_collection(coll_view_t &view)
      : base_collection(detail::get<0>(view.m_view)),
        m_mask_values(detail::get<1>(view.m_view)) {}

  /// Move constructor
  /// @note the base class constructor only moves the base class members
  DETRAY_HOST_DEVICE grid_collection(grid_collection &&other) noexcept
      : base_collection(std::move(other)),
        m_mask_values(std::move(other.m_mask_values)) {}

  /// Move assignment
  DETRAY_HOST_DEVICE grid_collection &operator=(
      grid_collection &&other) noexcept {
    if (this != &other) {
      m_mask_values = std::move(other.m_mask_values);
      base_collection::operator=(std::move(other));
    }
    return *this;
  }

  /// Create spatial grid acceleration structure from underlying grid data
  DETRAY_HOST_DEVICE
  constexpr auto operator[](const size_type i) const
      -> spatial_grid_impl<grid_t> {
    if constexpr (concepts::cylindrical<frame_t>) {
      assert(static_cast<dindex>(m_mask_values.size()) == this->size());
      return spatial_grid_impl<grid_t>(base_collection::operator[](i),
                                       m_mask_values[i]);
    } else {
      return spatial_grid_impl<grid_t>(base_collection::operator[](i));
    }
  }

  /// Create spatial grid acceleration structure from underlying grid data
  DETRAY_HOST_DEVICE
  constexpr auto at(const size_type i) const -> spatial_grid_impl<grid_t> {
    if constexpr (concepts::cylindrical<frame_t>) {
      assert(static_cast<dindex>(m_mask_values.size()) == this->size());
      return spatial_grid_impl<grid_t>(base_collection::at(i),
                                       m_mask_values[i]);
    } else {
      return spatial_grid_impl<grid_t>(base_collection::at(i));
    }
  }

  /// @returns a vecmem view on the spatial grid collection data - non-const
  DETRAY_HOST auto get_data() -> view_type {
    return view_type{base_collection::get_data(),
                     detray::get_data(m_mask_values)};
  }

  /// @returns a vecmem view on the spatial grid collection data - const
  DETRAY_HOST
  auto get_data() const -> const_view_type {
    return const_view_type{base_collection::get_data(),
                           detray::get_data(m_mask_values)};
  }

  /// Add a new grid @param gr to the collection.
  /// @note this takes a data owning grid to transcribe the data from.
  template <concepts::spatial_grid other_grid_t>
    requires std::constructible_from<
        typename spatial_grid_impl<grid_t>::template type<true>, other_grid_t>
  DETRAY_HOST constexpr auto push_back(const other_grid_t &gr) noexcept(false)
      -> void {
    // Copy over the base grid data
    base_collection::push_back(gr);

    // Add the additional mask values
    m_mask_values.push_back(gr.mask());
  }

 private:
  /// Container that holds the additional mask parameters for every grid
  vector_t<mask_t> m_mask_values{};
};

}  // namespace detray
